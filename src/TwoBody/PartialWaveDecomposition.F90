module PartialWaveDecomposition
  use, intrinsic :: iso_c_binding
  use MyLibrary, only: pi
  implicit none
  public :: PWD
  public :: identity_helicity_rep, sigma_dot_sigma, sigma1q_sigma2q, &
      & sigma_plus_q, sigma_plus_QQ, sigma_minus_q, sigma_minus_QQ, sigma_cross_q, sigma_cross_QQ, &
      & helicity_expectation_value_1, helicity_expectation_value_x, &
      & helicity_expectation_value_y, helicity_expectation_value_z, helicity_expectation_value_sigma, &
      & get_zmesh, get_y1
  private

  type :: PWD
    integer, private :: NMesh = 100 ! angular mesh number
    integer, private :: Jmax = 8    ! max of angular momentum
    integer, private :: rank = 0    ! operator rank
    character(:), allocatable, private :: Regulator
    integer, private :: RegulatorPower = 2  ! Regulator power n : exp( - (q^2/ Lam^2)^n )
    real(8), private :: RegulatorLambda = 500.d0 ! Cutoff in unit of MeV

    integer, private :: num_helicities = 4 ! ++, +-, -+, --
    integer :: lams_to_idx(-1:1,-1:1)
    integer :: lam1(4) = [ 1, 1,-1,-1]
    integer :: lam2(4) = [ 1,-1, 1,-1]

    real(8), allocatable :: ZMesh(:), WMesh(:) ! z=cos(t), w=weight
    real(8), allocatable :: dj(:,:,:,:) ! z, m1, m2, j
  contains
    procedure :: GetJmax
    procedure :: GetNMesh
    procedure :: GetNumHeli
    procedure :: GetHelicities
    procedure :: GetIndexHeli
    procedure :: GetRank
    procedure :: init_pwd
    procedure :: do_pwd
    procedure :: show_pwd_options
    procedure :: fin_pwd
  end type PWD

  type :: ovlap_lams
    real(8), allocatable :: v(:,:)
  contains
    procedure :: init_ovlap_lams
    procedure :: release_ovlap_lams
  end type ovlap_lams

  type :: ovlap_j
    type(ovlap_lams), allocatable :: v(:)
    integer :: jmin, jmax
  contains
    procedure :: init_ovlap_j
    procedure :: release_ovlap_j
  end type ovlap_j

  type :: ovlap_ls
    type(ovlap_j), allocatable :: v(:,:)
    integer :: lmin, lmax, smin, smax
  contains
    procedure :: init_ovlap_ls
    procedure :: release_ovlap_ls
    procedure :: get_ovlap_lsj
  end type ovlap_ls

  type(ovlap_ls) :: ovlaps
  real(8), allocatable :: sigma_1_store(:,:,:)
  real(8), allocatable :: sigma_x_store(:,:,:)
  real(8), allocatable :: sigma_y_store(:,:,:)
  real(8), allocatable :: sigma_z_store(:,:,:)
  real(8), allocatable :: sigma_m_store(:,:,:,:)
  real(8), allocatable :: ZMesh(:)
  real(8), allocatable :: y1store(:,:) ! i, m
contains

  function GetJmax(this) result(r)
    class(PWD), intent(in) :: this
    integer :: r
    r = this%Jmax
  end function GetJmax

  function GetRank(this) result(r)
    class(PWD), intent(in) :: this
    integer :: r
    r = this%rank
  end function GetRank

  function GetNMesh(this) result(r)
    class(PWD), intent(in) :: this
    integer :: r
    r = this%NMesh
  end function GetNMesh

  function GetNumHeli(this) result(r)
    class(PWD), intent(in) :: this
    integer :: r
    r = this%num_helicities
  end function GetNumHeli

  function GetHelicities(this, idx) result(r)
    class(PWD), intent(in) :: this
    integer, intent(in) :: idx
    integer :: r(2)
    r = [this%lam1(idx), this%lam2(idx)]
  end function GetHelicities

  function GetIndexHeli(this, lam1, lam2) result(r)
    class(PWD), intent(in) :: this
    integer, intent(in) :: lam1, lam2
    integer :: r
    r = this%lams_to_idx(lam1, lam2)
  end function GetIndexHeli

  subroutine show_pwd_options(this)
    class(PWD), intent(in) :: this
    write(*,"(a)") "# Option for PWD"
    write(*,"(a,i3,a,i4)") "# Jmax: ", this%GetJmax(), &
        &", NMesh for partial wave decomposition: ", this%GetNMesh()
    if(this%Regulator=="None") then
      write(*,"(2a)") "# Regulator: ", trim(this%Regulator)
    else
      write(*,"(3a,f8.2,a,i3)") "# Regulator: ", trim(this%Regulator), &
          & ", Lambda: ", this%RegulatorLambda, " MeV, Power: ", this%RegulatorPower
    end if
  end subroutine show_pwd_options

  subroutine init_pwd(this, NMesh, Jmax, rank, Regulator, RegulatorPower, RegulatorLambda)
    use MyLibrary, only: gauss_legendre, spherical_harmonics0
    class(PWD), intent(inout) :: this
    integer, intent(in), optional :: NMesh, Jmax, rank, RegulatorPower
    real(8), intent(in), optional :: RegulatorLambda
    character(*), intent(in), optional :: Regulator
    integer :: i_heli, i_mesh, lam1, lam2
    real(8) :: z

    this%Regulator="None"
    if(present(NMesh)) this%NMesh = NMesh
    if(present(Jmax)) this%Jmax = Jmax
    if(present(rank)) this%rank = rank
    if(present(Regulator)) this%Regulator = Regulator
    if(present(RegulatorPower)) this%RegulatorPower = RegulatorPower
    if(present(RegulatorLambda)) this%RegulatorLambda = RegulatorLambda

    allocate(this%dj(this%NMesh, min(-1,-this%Jmax):max(1,this%Jmax), -1:1, 0:this%Jmax) )
    this%dj(:,:,:,:) = 0.d0
    call gauss_legendre(-1.d0, 1.d0, this%ZMesh, this%WMesh, this%NMesh)
    call set_wigner_function(this)
    call ovlaps%init_ovlap_ls(this%GetJmax())
    this%lams_to_idx(:,:) = -1
    do i_heli = 1, this%GetNumHeli()
      this%lams_to_idx(this%lam1(i_heli), this%lam2(i_heli)) = i_heli
    end do
    allocate(ZMesh(this%NMesh))
    allocate(y1store(this%NMesh, -1:1))
    allocate(sigma_1_store(this%NMesh, -1:1, -1:1))
    allocate(sigma_x_store(this%NMesh, -1:1, -1:1))
    allocate(sigma_y_store(this%NMesh, -1:1, -1:1))
    allocate(sigma_z_store(this%NMesh, -1:1, -1:1))
    allocate(sigma_m_store(-1:1, this%NMesh, -1:1, -1:1))

    ! Note that i is not taken into account in sigma_y
    ZMesh(:) = this%ZMesh(:)
    do i_heli = 1, this%GetNUmHeli()
      lam1 = this%lam1(i_heli)
      lam2 = this%lam2(i_heli)
      do i_mesh = 1, this%NMesh
        z = this%ZMesh(i_mesh)
        sigma_1_store(i_mesh,lam1,lam2) = 0.5d0*(dble(abs(lam1+lam2))*sqrt(0.5d0*(1.d0+z)) + dble(lam1-lam2)*sqrt(0.5d0*(1.d0-z)))
        sigma_x_store(i_mesh,lam1,lam2) = 0.5d0*(dble(abs(lam1-lam2))*sqrt(0.5d0*(1.d0+z)) + dble(lam1+lam2)*sqrt(0.5d0*(1.d0-z)))
        sigma_y_store(i_mesh,lam1,lam2) = 0.5d0*(dble(lam2-lam1)*sqrt(0.5d0*(1.d0+z)) + dble(abs(lam1+lam2))*sqrt(0.5d0*(1.d0-z)))
        sigma_z_store(i_mesh,lam1,lam2) = 0.5d0*(dble(lam1+lam2)*sqrt(0.5d0*(1.d0+z)) - dble(abs(lam1-lam2))*sqrt(0.5d0*(1.d0-z)))
        sigma_m_store(-1,i_mesh,lam1,lam2) = (sigma_x_store(i_mesh,lam1,lam2) + sigma_y_store(i_mesh,lam1,lam2))/sqrt(2.d0)
        sigma_m_store( 0,i_mesh,lam1,lam2) =  sigma_z_store(i_mesh,lam1,lam2)
        sigma_m_store( 1,i_mesh,lam1,lam2) =-(sigma_x_store(i_mesh,lam1,lam2) - sigma_y_store(i_mesh,lam1,lam2))/sqrt(2.d0)
      end do
    end do

    do i_heli = -1, 1
      do i_mesh = 1, this%NMesh
        z = this%ZMesh(i_mesh)
        y1store(i_mesh,i_heli) = spherical_harmonics0(1,i_heli,z)
      end do
    end do

  end subroutine init_pwd

  subroutine fin_pwd(this)
    class(PWD), intent(inout) :: this
    this%NMesh = 100
    this%Jmax = 8
    this%rank = 0
    this%Regulator = "None"
    this%RegulatorPower = 2
    this%RegulatorLambda = 500.d0
    call ovlaps%release_ovlap_ls()
    deallocate(sigma_1_store)
    deallocate(sigma_x_store)
    deallocate(sigma_y_store)
    deallocate(sigma_z_store)
    deallocate(sigma_m_store)
    deallocate(y1store)
    deallocate(ZMesh)
    deallocate(this%dj)
    deallocate(this%ZMesh)
    deallocate(this%WMesh)
  end subroutine fin_pwd

  function do_pwd(this, op, pbra, lbra, sbra, jbra, pket, lket, sket, jket) result(r)
    class(PWD), intent(in) :: this
    real(8), intent(inout) :: op(this%NMesh,-this%Rank:this%Rank,this%num_helicities,this%num_helicities) ! z, mu, lam, lam
    real(8), intent(in) :: pbra, pket ! relative momenta in unit of MeV
    integer, intent(in) :: lbra, sbra, jbra, lket, sket, jket
    integer :: ibra, iket
    integer :: lams_bra(2), lams_ket(2)
    real(8) :: r
    complex(8), parameter :: i_unit = (0.d0, 1.d0)

    call multiply_regulator(this, op, pbra, pket)
    r = 0.d0
    do ibra = 1, this%GetNumHeli()
      lams_bra = this%GetHelicities(ibra)
      do iket = 1, this%GetNumHeli()
        lams_ket = this%GetHelicities(iket)
        !r = r + projection_to_j(this, op, ibra, iket, jbra, jket) * &
        !    & overlap_helicity_lsj(lbra,sbra,lams_bra(1),lams_bra(2),jbra) * &
        !    & overlap_helicity_lsj(lket,sket,lams_ket(1),lams_ket(2),jket)
        r = r + projection_to_j(this, op, ibra, iket, jbra, jket) * &
            & ovlaps%get_ovlap_lsj(lbra,sbra,lams_bra(1),lams_bra(2),jbra) * &
            & ovlaps%get_ovlap_lsj(lket,sket,lams_ket(1),lams_ket(2),jket)
      end do
    end do
    if(mod(abs(lbra-lket),2) == 0) then
      r = r * dble(i_unit ** (lbra-lket)) ! symmetric
    else
      r = r * dble(i_unit ** (lbra-lket-1)) ! skew, in unit of i
    end if
  end function do_pwd

  function projection_to_j(this, op, ibra, iket, jbra, jket) result(r)
    use MyLibrary, only: dcg
    type(PWD), intent(in) :: this
    real(8), intent(in) :: op(this%NMesh,-this%Rank:this%Rank,this%num_helicities,this%num_helicities) ! z, mu, lam, lam
    integer, intent(in) :: ibra, iket, jbra, jket
    real(8) :: r
    real(8) :: s
    integer :: lams_bra(2), lams_ket(2), lam_bra, lam_ket
    integer :: mu, i
    lams_bra = this%GetHelicities(ibra)
    lams_ket = this%GetHelicities(iket)
    lam_bra = (lams_bra(1) - lams_bra(2))/2
    lam_ket = (lams_ket(1) - lams_ket(2))/2

    r = 0.d0
    if(abs(lam_ket) > jket) return
    if(abs(lam_bra) > jbra) return
    do mu = -this%GetRank(), this%GetRank()
      if(abs(mu+lam_ket) > jbra) cycle
      s = 0.d0
      do i = 1, this%GetNMesh()
        s = s + this%WMesh(i) * this%dj(i,mu+lam_ket, lam_bra, jbra) * op(i,mu,ibra,iket)
      end do
      r = r + dcg(2*jket, 2*lam_ket, 2*this%GetRank(), 2*mu, 2*jbra, 2*(mu+lam_ket)) * s
    end do
    r = r * 2.d0 * pi * sqrt(dble(2*jket+1))
  end function projection_to_j

  function overlap_helicity_lsj(l, s, lam1, lam2, j) result(r)
    use MyLibrary, only: dcg
    integer, intent(in) :: l, s, lam1, lam2, j
    real(8) :: r
    integer :: lam
    r = 0.d0
    lam = lam1-lam2
    if(2*s < abs(lam)) return
    if(2*j < abs(lam)) return
    r = sqrt( dble(2*l+1) / dble(2*j+1) ) * &
        & dcg(2*l, 0, 2*s, lam, 2*j, lam) * dcg(1, lam1, 1, -lam2, 2*s, lam)
  end function overlap_helicity_lsj

  subroutine init_ovlap_ls(this, jmax)
    class(ovlap_ls), intent(inout) :: this
    integer, intent(in) :: jmax
    integer :: l, s, j
    this%lmin = 0
    this%lmax = jmax+1
    this%smin = 0
    this%smax = 1
    allocate(this%v(this%lmin:this%lmax, this%smin:this%smax))
    do l = this%lmin, this%lmax
      do s = this%smin, this%smax
        call this%v(l,s)%init_ovlap_j(l,s)
      end do
    end do
  end subroutine init_ovlap_ls

  subroutine release_ovlap_ls(this)
    class(ovlap_ls), intent(inout) :: this
    integer :: l, s
    do l = this%lmin, this%lmax
      do s = this%smin, this%smax
        call this%v(l,s)%release_ovlap_j()
      end do
    end do
    deallocate(this%v)
  end subroutine release_ovlap_ls

  subroutine init_ovlap_j(this, l, s)
    class(ovlap_j), intent(inout) :: this
    integer, intent(in) :: l, s
    integer :: j
    this%jmin = abs(l-s)
    this%jmax = l+s
    allocate(this%v(this%jmin:this%jmax))
    do j = this%jmin, this%jmax
      call this%v(j)%init_ovlap_lams(l,s,j)
    end do
  end subroutine init_ovlap_j

  subroutine release_ovlap_j(this)
    class(ovlap_j), intent(inout) :: this
    integer :: j
    do j = this%jmin, this%jmax
      call this%v(j)%release_ovlap_lams()
    end do
    deallocate(this%v)
  end subroutine release_ovlap_j

  subroutine init_ovlap_lams(this, l, s, j)
    class(ovlap_lams), intent(inout) :: this
    integer, intent(in) :: l, s, j
    integer :: lam1, lam2
    allocate(this%v(-1:1, -1:1))
    this%v(:,:) = 0.d0
    do lam1 = -1, 1, 2
      do lam2 = -1, 1, 2
        this%v(lam1,lam2) = overlap_helicity_lsj(l, s, lam1, lam2, j)
      end do
    end do
  end subroutine init_ovlap_lams

  subroutine release_ovlap_lams(this)
    class(ovlap_lams), intent(inout) :: this
    deallocate(this%v)
  end subroutine release_ovlap_lams

  function get_ovlap_lsj(this, l, s, lam1, lam2, j) result(r)
    class(ovlap_ls), intent(in) :: this
    integer, intent(in) :: l, s, j, lam1, lam2
    real(8) :: r
    r = this%v(l,s)%v(j)%v(lam1,lam2)
  end function get_ovlap_lsj

  function get_y1(idx_mesh,m) result(r)
    integer, intent(in) :: idx_mesh, m
    real(8) :: r
    r = y1store(idx_mesh,m)
  end function get_y1

  function get_zmesh(idx_mesh) result(r)
    integer, intent(in) :: idx_mesh
    real(8) :: r
    r = ZMesh(idx_mesh)
  end function get_zmesh
  !
  ! helicity functions
  !
  function identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh) result(r)
    ! 1
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8) :: r
    r = sigma_1_store(idx_mesh,lam1,lam3) * sigma_1_store(idx_mesh,-lam2,-lam4)
  end function identity_helicity_rep

  function sigma_dot_sigma(lam1, lam2, lam3, lam4, idx_mesh) result(r)
    ! s1 . s2
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8) :: r
    r =   sigma_x_store(idx_mesh,lam1,lam3) * sigma_x_store(idx_mesh,-lam2,-lam4) - &
        & sigma_y_store(idx_mesh,lam1,lam3) * sigma_y_store(idx_mesh,-lam2,-lam4) + &
        & sigma_z_store(idx_mesh,lam1,lam3) * sigma_z_store(idx_mesh,-lam2,-lam4)
  end function sigma_dot_sigma

  function sigma_plus_q(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! (s1 + s2) . q
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout * dble(lam1-lam2) - pin * dble(lam3-lam4)) * &
        & identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh)
  end function sigma_plus_q

  function sigma_plus_qq(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! (s1 + s2) . Q  with Q = pout + pin
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout * dble(lam1-lam2) + pin * dble(lam3-lam4)) * &
        & identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh)
  end function sigma_plus_qq

  function sigma_minus_q(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! (s1 - s2) . q
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout * dble(lam1+lam2) - pin * dble(lam3+lam4)) * &
        & identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh)
  end function sigma_minus_q

  function sigma_minus_qq(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! (s1 - s2) . Q  with Q = pout + pin
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout * dble(lam1+lam2) + pin * dble(lam3+lam4)) * &
        & identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh)
  end function sigma_minus_qq

  function sigma_cross_q(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! i(s1 x s2) . q
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout*ZMesh(idx_mesh) - pin) * ( &
        & helicity_expectation_value_x( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_y(-lam2, -lam4, idx_mesh) - &
        & helicity_expectation_value_y( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_x(-lam2, -lam4, idx_mesh) ) + &
        & pout * sqrt(1.d0-ZMesh(idx_mesh)**2) * (&
        & helicity_expectation_value_y( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_z(-lam2, -lam4, idx_mesh) - &
        & helicity_expectation_value_z( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_y(-lam2, -lam4, idx_mesh) )
    r = r * (-1.d0)
  end function sigma_cross_q

  function sigma_cross_qq(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! i(s1 x s2) . Q with Q = pout + pin
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = (pout*ZMesh(idx_mesh) + pin) * ( &
        & helicity_expectation_value_x( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_y(-lam2, -lam4, idx_mesh) - &
        & helicity_expectation_value_y( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_x(-lam2, -lam4, idx_mesh) ) + &
        & pout * sqrt(1.d0-ZMesh(idx_mesh)**2) * (&
        & helicity_expectation_value_y( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_z(-lam2, -lam4, idx_mesh) - &
        & helicity_expectation_value_z( lam1,  lam3, idx_mesh) * &
        & helicity_expectation_value_y(-lam2, -lam4, idx_mesh) )
    r = r * (-1.d0)
  end function sigma_cross_qq

  function sigma1q_sigma2q(lam1,lam2,lam3,lam4,pout,pin,idx_mesh) result(r)
    ! (s1 . q) (s2 . q)
    integer, intent(in) :: lam1, lam2, lam3, lam4, idx_mesh
    real(8), intent(in) :: pout, pin
    real(8) :: r
    r = -(pout*dble(lam1) - pin*dble(lam3)) * (pout*dble(lam2) - pin*dble(lam4)) * &
        & identity_helicity_rep(lam1, lam2, lam3, lam4, idx_mesh)
  end function sigma1q_sigma2q

  function helicity_expectation_value_1(lam1, lam2, idx_mesh) result(r)
    integer, intent(in) :: lam1, lam2, idx_mesh
    real(8) :: r
    r = sigma_1_store(idx_mesh,lam1,lam2)
  end function helicity_expectation_value_1

  function helicity_expectation_value_x(lam1, lam2, idx_mesh) result(r)
    integer, intent(in) :: lam1, lam2, idx_mesh
    real(8) :: r
    r = sigma_x_store(idx_mesh,lam1,lam2)
  end function helicity_expectation_value_x

  function helicity_expectation_value_y(lam1, lam2, idx_mesh) result(r)
    ! Note: i is not taken into account
    integer, intent(in) :: lam1, lam2, idx_mesh
    real(8) :: r
    r = sigma_y_store(idx_mesh,lam1,lam2)
  end function helicity_expectation_value_y

  function helicity_expectation_value_z(lam1, lam2, idx_mesh) result(r)
    integer, intent(in) :: lam1, lam2, idx_mesh
    real(8) :: r
    r = sigma_z_store(idx_mesh,lam1,lam2)
  end function helicity_expectation_value_z

  function helicity_expectation_value_sigma(lam1, lam2, idx_mesh, m) result(r)
    integer, intent(in) :: lam1, lam2, idx_mesh, m
    real(8) :: r
    r = sigma_m_store(m,idx_mesh,lam1,lam2)
  end function helicity_expectation_value_sigma

  subroutine multiply_regulator(this, op, pbra, pket)
    type(PWD), intent(in) :: this
    real(8), intent(inout) :: op(this%NMesh,-this%Rank:this%Rank,this%num_helicities,this%num_helicities) ! z, mu, lam, lam
    real(8), intent(in) :: pbra, pket
    integer :: i
    real(8) :: q2
    select case(this%Regulator)
    case("Local", "local")
      do i = 1, this%GetNMesh()
        q2 = pbra**2 + pket**2 - 2.d0 * pbra * pket * this%ZMesh(i)
        op(i,:,:,:) = op(i,:,:,:) * exp(-(q2 / this%RegulatorLambda**2)**this%RegulatorPower )
      end do
    case("NonLocal", "non-local", "Nonlocal")
      op(:,:,:,:) = op(:,:,:,:) * &
          & exp( - (pbra**2 /this%RegulatorLambda**2)**this%RegulatorPower ) * &
          & exp( - (pket**2 /this%RegulatorLambda**2)**this%RegulatorPower )
    case("None", "none")
    case default
      write(*,*) "Unknown regulator: ", trim(this%Regulator)
      stop
    end select
  end subroutine multiply_regulator

  subroutine set_wigner_function(this)
    type(PWD), intent(inout) :: this
    integer :: J, m1, m2, ith
    real(8) :: z
    do J = 0, this%Jmax
      do m1 = -J, J
        do m2 = max(-1, -J), min(1, J)
          do ith = 1, this%NMesh
            z = acos(this%ZMesh(ith)) ! cos theta
            this%dj(ith, m1, m2, J) = d_func(J, m1, m2, z)
          end do
        end do
      end do
    end do
  end subroutine set_wigner_function

  function ln_binom(n, k) result(r)
    use MyLibrary, only: ln_gamma
    integer, intent(in) :: n, k
    real(8) :: r
    r = ln_gamma(dble(n+1)) - ln_gamma(dble(k+1)) - ln_gamma(dble(n-k+1))
  end function ln_binom

  function Jacobi_poly(alpha, beta, n, x) result(r)
    use MyLibrary, only: ln_gamma
    integer, intent(in) :: n
    real(8), intent(in) :: alpha, beta, x
    real(8) :: r, fact
    integer :: m
    r = 0.d0
    do m = 0, n
      fact = ln_binom(n,m) + ln_gamma(alpha+beta+dble(n+m+1)) - ln_gamma(alpha+dble(m+1))
      r = r + exp(fact) * ((x-1.d0)*0.5d0)**m
    end do
    fact = ln_gamma(alpha+dble(n+1)) - ln_gamma(dble(n+1)) - ln_gamma(alpha+beta+dble(n+1))
    r = r * exp(fact)
  end function Jacobi_poly

  function d_func(j, m1, m2, theta) result(r)
    integer, intent(in) :: j, m1, m2
    real(8), intent(in) :: theta
    real(8) :: r, ln_fac, fac
    integer :: k, a, b, lam
    k = min(j+m2, j-m2, j+m1, j-m1)
    if(k==j+m2) then
      a = m1 - m2
      lam=m1 - m2
    else if(k==j-m2) then
      a = m2 - m1
      lam=0
    else if(k==j+m1) then
      a = m2 - m1
      lam=0
    else if(k==j-m1) then
      a = m1 - m2
      lam=m1 - m2
    end if
    b = 2*j - 2*k - a
    ln_fac = 0.5d0*ln_binom(2*j-k, k+a) - 0.5d0*ln_binom(k+b, b)
    fac = (-1.d0)**lam * exp(ln_fac) * sin(0.50d0*theta)**a * cos(0.5d0*theta)**b
    r = fac * jacobi_poly(dble(a), dble(b), k, cos(theta))
  end function d_func
end module PartialWaveDecomposition

