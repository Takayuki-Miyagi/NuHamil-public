module NNNForceLocal
  !$ use omp_lib
  use StoreCouplings
  use LinAlgLib
  implicit none

  type :: PWChan
    integer :: l12, s12, j12, t12, l3, j3
    integer, allocatable :: n12(:), n3(:), ns2i(:,:)
    integer :: Nho
  contains
    procedure :: InitPWChan
    procedure :: FinPWChan
    generic :: init => InitPWChan
    generic :: fin => FinPWChan
  end type PWChan

  type :: Coordinates
    integer :: i1, i2
    real(8) :: x1, x2, w1, w2
  end type Coordinates

  type, private :: ModelSpace
    ! Coordinate
    type(Coordinates), allocatable :: xis(:)
    integer :: Nxi
    integer :: Nxis
    ! HO radial quantum number
    integer :: Nmax
    real(8) :: hw
    ! partial wave channels
    type(PWChan), allocatable :: alpha(:)
    integer, private :: Nalpha
    integer, allocatable :: LSJTlj2idx(:,:,:,:,:,:)
    integer :: jtot, ptot, ttot
  contains
    procedure :: InitModelSpace
    procedure :: FinModelSpace
    procedure :: PrintModelSpace
    generic :: init => InitModelSpace
    generic :: Fin => FinModelSpace
    generic :: prt => PrintModelSpace
  end type ModelSpace

  type, extends(DMat) :: vmat
    logical :: zero=.true.
  end type vmat

  type :: NNNIntLocal
    type(vmat), allocatable :: MatCh(:,:)
    type(ModelSpace) :: ms
    real(8) :: lec = 1.d0
    real(8) :: lambda = 500.d0
  contains
    procedure :: InitNNNIntLocal
    procedure :: FinNNNIntLocal
    procedure :: ReleaseNNNIntLocal
    procedure :: SetNNNIntLocal
    procedure :: GetMEFromIndex
    procedure :: GetMEFromQuantumNumbers
    generic :: init => InitNNNIntLocal
    generic :: fin => FinNNNIntLocal
    generic :: release => ReleaseNNNIntLocal
    generic :: set => SetNNNIntLocal
    generic :: get => GetMEFromIndex, GetMEFromQuantumNumbers
  end type NNNIntLocal

  type :: r_dependence
    real(8), allocatable :: v(:)
  end type r_dependence
  type(r_dependence) :: z0
  type(r_dependence), allocatable :: zx(:), fk(:), fkx(:,:)

  type(CGzeroStore) :: CGs
  type(SixJsStore) :: lsj12, lsj3, jjx, kkxy
  type(NineJsStore) :: ls12, ls3, spin12

  integer :: NMesh_r = 80
  real(8) :: rmin = 0.d0, rmax = 15.d0
  real(8), allocatable :: r_mesh(:), rw_mesh(:)
  real(8), allocatable :: radial_ho_wf(:,:,:)

  integer :: NMesh_p = 100
  real(8) :: pmin = 0.d0, pmax = 8.d0
  real(8), allocatable :: p_mesh(:), pw_mesh(:)

  integer :: NMesh_cos = 80
  real(8), allocatable :: cos_mesh(:), cosw_mesh(:)
contains

  subroutine FinNNNIntLocal(this)
    class(NNNIntLocal), intent(inout) :: this
    call release_arrays(this%ms)
    call this%ms%fin()
  end subroutine FinNNNIntLocal

  subroutine ReleaseNNNIntLocal(this)
    class(NNNIntLocal), intent(inout) :: this
    integer :: chbra, chket
    do chbra = 1, this%ms%Nalpha
      do chket = 1, this%ms%Nalpha
        call this%MatCh(chbra,chket)%fin()
      end do
    end do
    deallocate(this%MatCh)
  end subroutine ReleaseNNNIntLocal

  subroutine InitNNNIntLocal(this, Jtot, Ptot, Ttot, Nmax, hw, NMesh_in, xmin_in, xmax_in, lambda_in, &
        & regulator_power_in)
    class(NNNIntLocal), intent(inout) :: this
    integer, intent(in) :: Jtot, Ptot, Ttot, Nmax
    real(8), intent(in) :: hw
    integer, intent(in), optional :: NMesh_in, regulator_power_in
    real(8), intent(in), optional :: lambda_in, xmin_in, xmax_in
    integer :: regulator_power=2
    real(8) :: time

    if(present(lambda_in)) this%lambda = lambda_in
    if(present(regulator_power_in)) regulator_power = regulator_power_in
    call this%ms%init(Jtot, Ptot, Ttot, Nmax, hw, NMesh_in, xmin_in, xmax_in)
    time = omp_get_wtime()
    call precalculations(Nmax,Jtot,hw,this%ms%xis,this%lambda,regulator_power)
#ifndef MPI
    write(*,"(a,f12.6,a)") "# Local 3NF, precalculations: ", omp_get_wtime()-time, " sec"
#endif
  end subroutine InitNNNIntLocal

  function GetMEFromIndex(this, chbra, chket, bra, ket) result(me)
    class(NNNIntLocal), intent(in) :: this
    integer, intent(in) :: chbra, chket, bra, ket
    real(8) :: me
    me = 0.d0
    if(this%MatCH(chbra,chket)%zero) return
    me = this%MatCh(chbra,chket)%m(bra,ket)
  end function GetMEFromIndex

  function GetMEFromQuantumNumbers(this, n12, l12, s12, j12, t12, n3, l3, j3, &
        & n45, l45, s45, j45, t45, n6, l6, j6) result(me)
    class(NNNIntLocal), intent(in) :: this
    integer, intent(in) :: n12, l12, s12, j12, t12, n3, l3, j3
    integer, intent(in) :: n45, l45, s45, j45, t45, n6, l6, j6
    real(8) :: me
    integer :: chbra, chket, bra, ket

    me = 0.d0
    chbra = this%ms%LSJTlj2idx(l12,s12,j12,t12,l3,(j3+1)/2)
    chket = this%ms%LSJTlj2idx(l45,s45,j45,t45,l6,(j6+1)/2)
    if(this%MatCH(chbra,chket)%zero) return
    if(chbra * chket == 0) then
      write(*,*) "Warning in GetMEFromQuantumNumbers channel, something might be wrong"
      return
    end if
    bra = this%ms%alpha(chbra)%ns2i(n12,n3)
    ket = this%ms%alpha(chket)%ns2i(n45,n6)
    if(bra * ket == 0) then
      write(*,*) "Warning in GetMEFromQuantumNumbers radial, something might be wrong"
      return
    end if
    me = this%get(chbra,chket,bra,ket)
  end function GetMEFromQuantumNumbers

  subroutine SetNNNIntLocal(this, term_name, lec)
    class(NNNIntLocal), intent(inout) :: this
    character(*), intent(in) :: term_name
    real(8), intent(in), optional :: lec
    integer :: chbra, chket
    integer :: Jtot, Ttot
    real(8) :: time
    this%lec = 1.d0
    if(present(lec)) this%lec = lec
    time = omp_get_wtime()
    Jtot = this%ms%Jtot
    Ttot = this%ms%Ttot
    allocate(this%MatCh(this%ms%Nalpha, this%ms%Nalpha))
    do chbra = 1, this%ms%Nalpha
      do chket = 1, chbra
        select case(term_name)
        case("C1","c1")
          call set_two_pion_exchange_c1(this, chbra, chket)
        case("C3","c3")
          call set_two_pion_exchange_c3(this, chbra, chket)
        case("C4","c4")
          call set_two_pion_exchange_c4(this, chbra, chket)
        case("CD","cD","Cd","cd")
          call set_one_pion_exchange_cd(this, chbra, chket)
        case("CE","cE","Ce","ce")
          call set_contact_ce(this, chbra, chket)
        end select
      end do
    end do
#ifndef MPI
    write(*,"(3a,f12.6,a)") "# Local 3NF, set ", trim(term_name), ": ", omp_get_wtime()-time, " sec"
#endif
  end subroutine SetNNNIntLocal

  subroutine set_two_pion_exchange_c1(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (42) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: gamma_function, sjs, f_pi, hc, m_pi, g_A => g_A_GT
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    integer :: v, r, ymin, ymax, y, k3, k4, x
    real(8) :: spin_factor, isospin_factor, prefact, xi1, xi2, phase
    real(8) :: fvr, fy, fk3k4, fx
    real(8) :: vvr, vy
    type(DVec) :: vr, vk, vx, vi
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif
    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)
    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12; l3  = chbra%l3; j3  = chbra%j3
    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12; l6  = chket%l3; j6  = chket%j3
    Jtot = ms%Jtot; Ttot = ms%Ttot

    isospin_factor = tau2_dot_tau3(t12,t45,Ttot)
    if(abs(isospin_factor) < 1.d-8) return
    spin_factor = sqrt( dble(2*s12+1) * dble(2*s45+1) ) * sjs(2*s12, 2*s45, 2, 1, 1, 1)
    if(abs(spin_factor) < 1.d-8) return
    prefact = this%lec * m_pi**2 * g_A**2 * (hc**4) * 1.d-3 / (f_pi**4 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*l45+1) * dble(2*l6+1) )

    phase = (-1.d0)**( (Jtot-j3)/2 + s12 + j45)
    call vr%zeros(ms%Nxis)
    call vk%zeros(ms%Nxis)
    call vx%zeros(ms%Nxis)
    call vi%zeros(ms%Nxis)

#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif

    vk%v = 0.d0
    do k3 = 0, 1
      k4 = 1 - k3
      fk3k4 = sqrt( dble(2*k3+1) * dble(2*k4+1) ) * &
          & sqrt( gamma_function(4.d0) / &
          & ( gamma_function(dble(2*k3+2)) * gamma_function(dble(2*k4+2))) )

      vx%v = 0.d0
      do x = 0, ms%Nmax+2

        vvr = 0.d0
        do r = max(abs(l3-l6), abs(x-k4)), min(l3+l6, x+k4)
          if(mod(x+k4+r,2)==1) cycle
          do v = abs(l12-l45), l12+l45
            if(mod(v+l12+l45,2)==1) cycle
            if(mod(r+l3+l6,2)==1) cycle
            fvr = sqrt( dble(2*v+1) * dble(2*r+1) ) * &
                & CGs%get(v,l45,l12) * CGs%get(r,l6,l3) * CGs%get(x,k4,r)
            if(abs(fvr) < 1.d-16) cycle

            ymin = max( abs(v-1), abs(j12-j45), abs(r-1), abs(j3-j6)/2, abs(x-k3) )
            ymax = min(    (v+1),    (j12+j45),    (r+1),    (j3+j6)/2,    (x+k3) )
            if(ymin > ymax) cycle
            vy = 0.d0
            do y = ymin, ymax
              if(mod(x+k3+y,2)==1) cycle
              fy = (-1.d0)**y * sqrt( dble(2*y+1) ) * CGs%get(y,1,v) * &
                  & ls12%get(2*v, 2, 2*y, 2*l12, 2*s12, 2*j12, 2*l45, 2*s45, 2*j45) * &
                  & ls3%get( 2*r, 2, 2*y, 2*l3,      1,   j3,  2*l6,      1,   j6 ) * &
                  & jjx%get( 2*y, 2*j12, 2*j45, Jtot, j6, j3) * &
                  & CGs%get(x,k3,y) * &
                  & kkxy%get(2*k4, 2, 2*k3, 2*y, 2*x, 2*r)

              vy = vy + fy
            end do
            vvr = vvr + vy * fvr
          end do
        end do
        fx = dble(2*x+1) * vvr
        if(abs(fx) < 1.d-16) cycle
        vi%v = 0.d0
        do i = 1, ms%Nxis
          i1 = ms%xis(i)%i1
          i2 = ms%xis(i)%i2
          xi1 = ms%xis(i)%x1 / dsqrt(2.d0)
          xi2 = ms%xis(i)%x2 * dsqrt(1.5d0)
          vi%v(i) = (xi1**k3) * (xi2**k4) * fk(1)%v(i1) * fkx(1,x)%v(i)
        end do
        vx = vx + vi * fx
      end do
      vk = vk + vx * fk3k4
    end do
    vr = vk * (-36.d0 * prefact * spin_factor * isospin_factor * phase)
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      call vk%fin()
      call vx%fin()
      call vi%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
    call vk%fin()
    call vx%fin()
    call vi%fin()
  end subroutine set_two_pion_exchange_c1

  subroutine set_two_pion_exchange_c3(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (47) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: gamma_function, sjs, f_pi, hc, g_A => g_A_GT
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    integer :: k1, k2, k3, k4, v, r, y, zmin, zmax, z, x
    real(8) :: spin_factor, isospin_factor, prefact, xi1, xi2, phase
    real(8) :: fk1k2, fvr, fy, fz, fk3k4, fx, vz, vy, vvr
    type(DVec) :: vr, vk, vk3k4, vx, vi
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif
    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)
    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12; l3  = chbra%l3; j3  = chbra%j3
    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12; l6  = chket%l3; j6  = chket%j3
    Jtot = ms%Jtot; Ttot = ms%Ttot

    isospin_factor = tau2_dot_tau3(t12,t45,Ttot)
    if(abs(isospin_factor) < 1.d-8) return
    spin_factor = sqrt( dble(2*s12+1) * dble(2*s45+1) ) * sjs(2*s12, 2*s45, 2, 1, 1, 1)
    if(abs(spin_factor) < 1.d-8) return
    prefact = this%lec * g_A**2 * (hc**6) * 1.d-3 / (f_pi**4 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*l45+1) * dble(2*l6+1) )
    phase = (-1.d0)**( (Jtot-j3)/2 + s12 + j45)
    call vr%zeros(ms%Nxis)
    call vk%zeros(ms%Nxis)
    call vk3k4%zeros(ms%Nxis)
    call vx%zeros(ms%Nxis)
    call vi%zeros(ms%Nxis)
#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif
    vk%v = 0.d0
    do k1 = 0, 2, 2
      do k2 = 0, 2, 2
        fk1k2 = ((-1.d0)**((k1+k2)/2) ) * sqrt( dble(2*k1+1) * dble(2*k2+1) ) * &
            & CGs%get(1,1,k1) * CGs%get(1,1,k2)
        if(abs(fk1k2) < 1.d-16) cycle

        vk3k4%v = 0.d0
        do k3 = 0, k2
          k4 = k2 - k3
          fk3k4 = sqrt( dble(2*k3+1) * dble(2*k4+1) ) * &
              & sqrt( gamma_function( dble(2*k2+2) ) / &
              & ( gamma_function( dble(2*k3+2) ) * gamma_function( dble(2*k4+2) ) ) )

          vx%v = 0.d0
          do x = 0, ms%Nmax+2

            vvr = 0.d0
            do r = max(abs(l3-l6), abs(x-k4)), min(l3+l6, x+k4)
              if(mod(x+k4+r,2)==1) cycle
              if(mod(r+l6+l3,2)==1) cycle
              do v = abs(l12-l45), l12+l45
                if(mod(v+l45+l12,2)==1) cycle
                fvr =  ((-1.d0)**(v+r)) * sqrt( dble(2*v+1) * dble(2*r+1) ) * &
                  & CGs%get(v,l45,l12) * CGs%get(r,l6,l3) * CGs%get(x,k4,r)
                if(abs(fvr) < 1.d-16) cycle

                vy = 0.d0
                do y = max(abs(v-k1), abs(r-k2), abs(x-k3)), min(v+k1,r+k2,x+k3)
                  if(mod(x+k3+y,2)==1) cycle
                  fy = sqrt(dble(2*y+1)) * CGs%get(y,k1,v) * CGs%get(x,k3,y) * &
                    & kkxy%get(2*k3, 2*k4, 2*k2, 2*r, 2*y, 2*x)
                  if(abs(fy) < 1.d-16) cycle

                  zmin = max( abs(v-1), abs(j12-j45), abs(r-1), abs(j3-j6)/2, abs(y-1) )
                  zmax = min(    (v+1),    (j12+j45),    (r+1),    (j3+j6)/2,    (y+1) )
                  if(zmin > zmax) cycle
                  vz = 0.d0
                  do z = zmin, zmax
                    fz = ((-1.d0)**z) * dble(2*z+1) * &
                        & ls12%get(2*v, 2, 2*z, 2*l12, 2*s12, 2*j12, 2*l45, 2*s45, 2*j45) * &
                        & ls3%get( 2*r, 2, 2*z, 2*l3,      1,   j3,  2*l6,      1,   j6 ) * &
                        & jjx%get( 2*z, 2*j12, 2*j45, Jtot, j6, j3) * &
                        & kkxy%get(2, 2, 2*k1, 2*v, 2*y, 2*z) * &
                        & kkxy%get(2, 2, 2*k2, 2*r, 2*y, 2*z)
                    vz = vz + fz
                  end do
                  vy = vy + vz * fy
                end do
                vvr = vvr + vy * fvr
              end do
            end do
            fx = dble(2*x+1) * vvr
            if(abs(fx) < 1.d-16) cycle
            vi%v = 0.d0
            do i = 1, ms%Nxis
              i1 = ms%xis(i)%i1
              i2 = ms%xis(i)%i2
              xi1 = ms%xis(i)%x1 / dsqrt(2.d0)
              xi2 = ms%xis(i)%x2 * dsqrt(1.5d0)
              vi%v(i) = (xi1**k3) * (xi2**k4) * fk(k1)%v(i1) * fkx(k2,x)%v(i)
            end do
            vx = vx + vi * fx
          end do
          vk3k4 = vk3k4 + vx * fk3k4
        end do
        vk = vk + vk3k4 * fk1k2
      end do
    end do
    vr = vk * (18.d0 * prefact * spin_factor * isospin_factor * phase)

#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      call vk%fin()
      call vk3k4%fin()
      call vx%fin()
      call vi%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
    call vk%fin()
    call vk3k4%fin()
    call vx%fin()
    call vi%fin()
  end subroutine set_two_pion_exchange_c3

  subroutine set_two_pion_exchange_c4(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (49) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: gamma_function, sjs, f_pi, hc, g_A => g_A_GT
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    integer :: k1, k2, k3, k4, k5, k5min, k5max, v, r, y, zmin, zmax, z, x
    real(8) :: isospin_factor, prefact, xi1, xi2, phase
    real(8) :: fk1k2, fvr, fy, fz, fk3k4, fk5, fx
    real(8) :: vk5, vz, vy, vvr
    type(DVec) :: vr, vk, vk3k4, vx, vi
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif
    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)
    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12; l3  = chbra%l3; j3  = chbra%j3
    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12; l6  = chket%l3; j6  = chket%j3
    Jtot = ms%Jtot; Ttot = ms%Ttot

    isospin_factor = tau_cross(t12,t45,Ttot)
    if(abs(isospin_factor) < 1.d-8) return
    prefact = this%lec * g_A**2 * (hc**6) * 1.d-3 / (f_pi**4 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*s12+1) * dble(2*s45+1) * dble(2*l45+1) * dble(2*l6+1) )
    phase = (-1.d0)**( (Jtot+j3)/2 + j45)
    call vr%zeros(ms%Nxis)
    call vk%zeros(ms%Nxis)
    call vk3k4%zeros(ms%Nxis)
    call vx%zeros(ms%Nxis)
    call vi%zeros(ms%Nxis)

#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif
    do k1 = 0, 2, 2
      do k2 = 0, 2, 2
        fk1k2 = ((-1.d0)**((k1+k2)/2) ) * sqrt( dble(2*k1+1) * dble(2*k2+1) ) * &
            & CGs%get(1,1,k1) * CGs%get(1,1,k2)
        if(abs(fk1k2) < 1.d-16) cycle

        vk3k4%v = 0.d0
        do k3 = 0, k2
          k4 = k2 - k3
          fk3k4 = sqrt( dble(2*k3+1) * dble(2*k4+1) ) * &
              & sqrt( gamma_function( dble(2*k2+2) ) / &
              & ( gamma_function( dble(2*k3+2) ) * gamma_function( dble(2*k4+2) ) ) )

          vx%v = 0.d0
          do x = 0, ms%Nmax+2

            vvr = 0.d0
            do r = max(abs(l3-l6), abs(x-k4)), min(l3+l6,x+k4)
              if(mod(x+k4+r,2)==1) cycle
              if(mod(r+l6+l3,2)==1) cycle
              do v = abs(l12-l45), l12+l45
                if(mod(v+l45+l12,2)==1) cycle
                fvr = ((-1.d0)**(v+r)) * sqrt( dble(2*v+1) * dble(2*r+1) ) * &
                    & CGs%get(v,l45,l12) * CGs%get(r,l6,l3) * CGs%get(x,k4,r)
                if(abs(fvr) < 1.d-16) cycle

                vy = 0.d0
                do y = max(abs(v-k1), abs(r-k2), abs(x-k3)), min(v+k1,r+k2,x+k3)
                  if(mod(x+k3+y,2)==1) cycle
                  fy = sqrt(dble(2*y+1)) * CGs%get(y,k1,v) * CGs%get(x,k3,y) * kkxy%get(2*k3, 2*k4, 2*k2, 2*r, 2*y, 2*x)
                  if(abs(fy) < 1.d-16) cycle

                  zmin = max( abs(j12-j45), abs(r-1), abs(j3-j6)/2, abs(y-1) )
                  zmax = min(    (j12+j45),    (r+1),    (j3+j6)/2,    (y+1) )
                  if(zmin > zmax) cycle
                  vz = 0.d0
                  do z = zmin, zmax
                    fz =  ((-1.d0)**z) * dble(2*z+1) * &
                        & ls3%get( 2*r, 2, 2*z, 2*l3,      1,   j3,  2*l6,      1,   j6 ) * &
                        & jjx%get( 2*z, 2*j12, 2*j45, Jtot, j6, j3) * &
                        & kkxy%get(2, 2, 2*k2, 2*r, 2*y, 2*z)
                    if(abs(fz) < 1.d-16) cycle

                    k5min = max( abs(s12-s45), abs(v-z), abs(k1-1) )
                    k5max = min(    (s12+s45),    (v+z),    (k1+1), 2 )
                    if(k5min > k5max) cycle

                    vk5 = 0.d0
                    do k5 = k5min, k5max
                      fk5 = dble(2*k5+1) * &
                          & ls12%get(2*v, 2*k5, 2*z, 2*l12, 2*s12, 2*j12, 2*l45, 2*s45, 2*j45) * &
                          & spin12%get( 1, 1, 2*s12, 1, 1, 2*s45, 2, 2, 2*k5) * &
                          & kkxy%get(2*k5, 2, 2*k1, 2*y, 2*v, 2*z) * &
                          & kkxy%get(2, 2, 2, 2, 2*k5, 2*k1)
                      vk5 = vk5 + fk5
                    end do
                    vz = vz + vk5 * fz
                  end do
                  vy = vy + vz * fy
                end do
                vvr = vvr + vy * fvr
              end do
            end do
            fx = dble(2*x+1) * vvr
            if(abs(fx) < 1.d-16) cycle
            vi%v = 0.d0
            do i = 1, ms%Nxis
              i1 = ms%xis(i)%i1
              i2 = ms%xis(i)%i2
              xi1 = ms%xis(i)%x1 / dsqrt(2.d0)
              xi2 = ms%xis(i)%x2 * dsqrt(1.5d0)
              vi%v(i) = (xi1**k3) * (xi2**k4) * fk(k1)%v(i1) * fkx(k2,x)%v(i)
            end do
            vx = vx + vi * fx
          end do
          vk3k4 = vk3k4 + vx * fk3k4
        end do
        vk = vk + vk3k4 * fk1k2
      end do
    end do
    vr = vk * (-324.d0 * prefact * isospin_factor * phase)

#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      call vk%fin()
      call vk3k4%fin()
      call vx%fin()
      call vi%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
    call vk%fin()
    call vk3k4%fin()
    call vx%fin()
    call vi%fin()
  end subroutine set_two_pion_exchange_c4

  subroutine set_one_pion_exchange_cd(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (38) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: lambda_chi, f_pi, hc, sjs, g_A => g_A_GT
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    real(8) :: spin_factor, isospin_factor, prefact, phase
    real(8) :: fkk, fv, fx, fz
    real(8) :: vv, vz
    integer :: xmin, xmax, x, v, k, zmin, zmax, z
    type(DVec) :: vr, vk, vx, vi
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif
    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)
    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12; l3  = chbra%l3; j3  = chbra%j3
    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12; l6  = chket%l3; j6  = chket%j3
    Jtot = ms%Jtot; Ttot = ms%Ttot

    isospin_factor = tau2_dot_tau3(t12, t45, Ttot)
    if(abs(isospin_factor) < 1.d-8) return
    spin_factor = sqrt( dble(2*s12+1) * dble(2*s45+1) ) * sjs(2*s12, 2*s45, 2, 1, 1, 1)
    if(abs(spin_factor) < 1.d-8) return
    prefact = this%lec * 9.d0 * g_a / (f_pi**4 * lambda_chi) * hc**6 * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*l45+1) * dble(2*l6+1) )
    phase = (-1.d0)**( (Jtot-j3)/2 + s12 + j45)
    call vr%zeros(ms%Nxis)
    call vk%zeros(ms%Nxis)
    call vx%zeros(ms%Nxis)
    call vi%zeros(ms%Nxis)
#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif

    vk%v = 0.d0
    do k = 0, 2, 2
      fkk = sqrt(dble(2*k+1)) * (-1.d0)**(k/2) * CGs%get(1,1,k)
      if(abs(fkk) < 1.d-16) cycle

        xmin = max( abs(v-k), abs(l3-l6))
        xmax = min(    (v+k),    (l3+l6))
        vx%v = 0.d0
        do x = abs(l3-l6), l3+l6
          if(mod(x+l3+l6,2)==1) cycle

          vv = 0.d0
          do v = max(abs(l12-l45), abs(x-k)), min(l12+l45,x+k)
            if(mod(x+k+v,2)==1) cycle
            if(mod(l12+l45+v,2) == 1) cycle
            fv = (-1.d0)**v * sqrt(dble(2*v+1)) * CGs%get(v, l45, l12) * CGs%get(x, k, v)
            if(abs(fv) < 1.d-16) cycle

            zmin = max( abs(j12-j45), abs(v-1), abs(x-1), abs(j3-j6)/2 )
            zmax = min(    (j12+j45),    (v+1),    (x+1),    (j3+j6)/2 )
            if(zmin > zmax) cycle
            vz = 0.d0
            do z = zmin, zmax
              fz = dble(2*z+1) * &
                  & ls12%get(2*v, 2, 2*z, 2*l12, 2*s12, 2*j12, 2*l45, 2*s45, 2*j45) * &
                  & ls3%get( 2*x, 2, 2*z, 2*l3 ,     1,   j3,  2*l6,      1,   j6) * &
                  & jjx%get( 2*z, 2*j12, 2*j45, Jtot, j6, j3) * &
                  & kkxy%get(2, 2, 2*K, 2*v, 2*x, 2*z)
              vz = vz + fz
            end do
            vv = vv + vz * fv
          end do
          fx = dble(2*x+1) * CGs%get(x, l6, l3) * vv
          if(abs(fx) < 1.d-16) cycle
          vi%v = 0.d0
          do i = 1, ms%Nxis
            i1 = ms%xis(i)%i1
            i2 = ms%xis(i)%i2
            vi%v(i) = fk(k)%v(i1) * zx(x)%v(i)
          end do
          vx = vx + vi * fx
        end do
      vk = vk + vx * fkk
    end do
    vr = vk * (-1.d0 * isospin_factor * spin_factor * phase * prefact)
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      call vk%fin()
      call vx%fin()
      call vi%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
    call vk%fin()
    call vx%fin()
    call vi%fin()
  end subroutine set_one_pion_exchange_cd

  subroutine set_one_pion_exchange_cd_sigma1(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (39) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: lambda_chi, f_pi, hc, g_A => g_A_GT
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    real(8) :: isospin_factor, prefact, vtot, vv, vx
    integer :: xmin, xmax, x, v, k
    type(DVec) :: vr
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif
    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)
    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12; l3  = chbra%l3; j3  = chbra%j3
    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12; l6  = chket%l3; j6  = chket%j3
    Jtot = ms%Jtot; Ttot = ms%Ttot

    if(t12 /= t45) return
    isospin_factor = tau1_dot_tau2(t12)
    if(abs(isospin_factor) < 1.d-8) return
    prefact = this%lec / (f_pi**4 * lambda_chi) * hc**6
    call vr%zeros(ms%Nxis)

#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif
    !$omp parallel
    !$omp do private(i, i1, i2, vtot, k, vv, v, xmin, xmax, vx, x)
    do i = 1, ms%Nxis
      i1 = ms%xis(i)%i1
      i2 = ms%xis(i)%i2
      vtot = 0.d0

      do k = abs(s12-s45), abs(s12+s45)
        if(k == 1) cycle
        vv = 0.d0
        do v = abs(l12-l45), l12+l45
          if(mod(l12+l45+v,2) == 1) cycle
          xmin = max( abs(v-k), abs(l3-l6), abs(j12-j45), abs(j3-j6)/2)
          xmax = min(    (v+k),    (l3+l6),    (j12+j45),    (j3+j6)/2)
          if(xmin > xmax) cycle
          vx = 0.d0
          do x = xmin, xmax
            if(mod(x+k+v,2)==1) cycle
            if(mod(x+l3+l6,2)==1) cycle

            vx = vx + dble(2*x+1) * &
                & jjx%get(  2*x, 2*j12, 2*j45, Jtot, j6, j3) * &
                & lsj3%get( 2*x, 2*l6,  2*l3,  1, j3, j6) * &
                & ls12%get( 2*v, 2*k, 2*x, 2*l12, 2*s12, 2*j12, 2*l45, 2*s45, 2*j45) * &
                & CGs%get(x, k, v) * CGs%get(x, l6, l3) * &
                & fk(k)%v(i1) * zx(x)%v(i)

          end do
          vv = vv + sqrt(dble(2*v+1)) * CGs%get(v, l45, l12) * vx
        end do
        vtot = vtot + sqrt(dble(2*k+1)) * (-1.d0)**(k/2) * CGs%get(1,1,k) * &
            & spin12%get(1, 1, 2, 1, 1, 2, 2*s12, 2*s45, 2*k) * vv
      end do
      vr%v(i) = - 9.d0 * vtot * isospin_factor * &
          & (-1.d0)**( s12+s45+j45+(Jtot+j3+j6+1)/2+l3 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*s12+1) * dble(2*s45+1) * dble(2*l45+1) * dble(2*l6+1) ) * &
          & prefact * g_a
    end do
    !$omp end do
    !$omp end parallel
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
  end subroutine set_one_pion_exchange_cd_sigma1

  subroutine set_contact_ce_old(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (15) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: lambda_chi, f_pi, hc
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    real(8) :: isospin_factor, v, prefact
    integer :: xmin, xmax, x
    type(DVec) :: vr
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif

    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)

    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12
    l3  = chbra%l3; j3  = chbra%j3

    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12
    l6  = chket%l3; j6  = chket%j3

    Jtot = ms%Jtot
    Ttot = ms%Ttot

    if(s12 /= s45) return
    isospin_factor = tau2_dot_tau3(t12, t45, Ttot)
    if(abs(isospin_factor) < 1.d-8) return

    prefact = this%lec / (f_pi**4 * lambda_chi) * hc**6
    call vr%zeros(size(ms%xis))
    xmin = max( abs(l12-l45), abs(j12-j45), abs(j3-j6)/2, abs(l3-l6) )
    xmax = min(    (l12+l45),    (j12+j45),    (j3+j6)/2,    (l3+l6) )
    if(xmin > xmax) return
#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif
    !$omp parallel
    !$omp do private(i, i1, i2, v, x)
    do i = 1, size(ms%xis)
      i1 = ms%xis(i)%i1
      i2 = ms%xis(i)%i2
      v = 0.d0
      do x = xmin, xmax
        if(mod(l12+l45+x,2) == 1) cycle
        if(mod(l3+l6+x,2) == 1) cycle
        v = v + (-1.d0)**x * dble(2*x+1) * &
            & lsj12%get(2*x, 2*l45, 2*l12, 2*s12, 2*j12, 2*j45) * &
            & jjx%get(  2*x, 2*j12, 2*j45, Jtot, j6, j3) * &
            & lsj3%get( 2*x, 2*l6,  2*l3,  1, j3, j6) * &
            & CGs%get(l45, x, l12) * CGs%get(l6, x, l3) * &
            & z0%v(i1) * zx(x)%v(i)

      end do
      vr%v(i) = 6.d0 * v * (-1.d0)**( (Jtot-1+j6-j3)/2 +l12+l3+s12 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*l45+1) * dble(2*l6+1) ) * prefact * isospin_factor
    end do
    !$omp end do
    !$omp end parallel
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
  end subroutine set_contact_ce_old

  subroutine set_contact_ce(this, ch_bra, ch_ket)
    !
    ! For details, see Eq. (15) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: lambda_chi, f_pi, hc
    class(NNNIntLocal), intent(inout), target :: this
    integer, intent(in) :: ch_bra, ch_ket
    type(ModelSpace), pointer :: ms
    type(PWChan), pointer :: chbra, chket
    integer :: Jtot, Ttot
    integer :: l12, s12, j12, t12, l3, j3
    integer :: l45, s45, j45, t45, l6, j6
    integer :: i, i1, i2
    real(8) :: isospin_factor, v, prefact
    integer :: xmin, xmax, x
    type(DVec) :: vr
#ifdef DebugNNNForceLocal
    real(8) :: time
#endif

    ms => this%ms
    chbra => ms%alpha(ch_bra)
    chket => ms%alpha(ch_ket)

    l12 = chbra%l12; s12 = chbra%s12; j12 = chbra%j12; t12 = chbra%t12
    l3  = chbra%l3; j3  = chbra%j3

    l45 = chket%l12; s45 = chket%s12; j45 = chket%j12; t45 = chket%t12
    l6  = chket%l3; j6  = chket%j3

    Jtot = ms%Jtot
    Ttot = ms%Ttot

    if(s12 /= s45) return
    isospin_factor = tau2_dot_tau3(t12, t45, Ttot)
    if(abs(isospin_factor) < 1.d-8) return

    prefact = this%lec / (f_pi**4 * lambda_chi) * hc**6
    call vr%zeros(size(ms%xis))
    xmin = max( abs(l12-l45), abs(j12-j45), abs(j3-j6)/2, abs(l3-l6) )
    xmax = min(    (l12+l45),    (j12+j45),    (j3+j6)/2,    (l3+l6) )
    if(xmin > xmax) return
#ifdef DebugNNNForceLocal
    time = omp_get_wtime()
#endif
    !$omp parallel
    !$omp do private(i, i1, i2, v, x)
    do i = 1, size(ms%xis)
      i1 = ms%xis(i)%i1
      i2 = ms%xis(i)%i2
      v = 0.d0
      do x = xmin, xmax
        if(mod(l12+l45+x,2) == 1) cycle
        if(mod(l3+l6+x,2) == 1) cycle
        v = v + (-1.d0)**x * dble(2*x+1) * &
            & lsj12%get(2*x, 2*l45, 2*l12, 2*s12, 2*j12, 2*j45) * &
            & jjx%get(  2*x, 2*j12, 2*j45, Jtot, j6, j3) * &
            & lsj3%get( 2*x, 2*l6,  2*l3,  1, j3, j6) * &
            & CGs%get(l45, x, l12) * CGs%get(l6, x, l3) * &
            & z0%v(i1) * zx(x)%v(i)
      end do
      vr%v(i) = 6.d0 * v * (-1.d0)**( (Jtot-1+j6-j3)/2 +l12+l3+s12 ) * &
          & sqrt( dble(2*j12+1) * dble(2*j45+1) * dble(j3+1) * dble(j6+1) * &
          & dble(2*l45+1) * dble(2*l6+1) ) * prefact * isospin_factor
    end do
    !$omp end do
    !$omp end parallel
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Set coordinate space 3NF: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
#endif
    if(maxval(abs(vr%v)) < 1.d-16) then
      call vr%fin()
      return
    end if
    this%MatCh(ch_bra, ch_ket)%DMat = transform_xis_to_ho(ms%Nxis,ms%xis,chbra,chket,vr)
    this%MatCh(ch_ket, ch_bra)%DMat = this%MatCh(ch_bra, ch_ket)%DMat%T()
    this%MatCh(ch_bra, ch_ket)%zero = .false.
    this%MatCh(ch_ket, ch_bra)%zero = .false.
#ifdef DebugNNNForceLocal
    write(*,"(a,f12.6,a)") "Transform Coordinate space -> HO space: ", omp_get_wtime() - time, " sec"
#endif
    call vr%fin()
  end subroutine set_contact_ce

  function transform_xis_to_ho_element(nxis, xis, chbra, chket, v) result(mat)
    integer, intent(in) :: nxis
    type(Coordinates), intent(in) :: xis(:)
    type(PWChan), intent(in) :: chbra, chket
    type(DVec), intent(in) :: v
    type(DMat) :: mat
    integer :: bra, ket, i, i1, i2
    integer :: n12, l12, n3, l3, n45, l45, n6, l6
    real(8) :: s
    call mat%zeros(chbra%Nho, chket%Nho)
    do bra = 1, chbra%Nho
      n12 = chbra%n12(bra)
      n3  = chbra%n3( bra)
      l12 = chbra%l12
      l3  = chbra%l3
      do ket = 1, chket%Nho
        n45 = chket%n12(ket)
        n6  = chket%n3( ket)
        l45 = chket%l12
        l6  = chket%l3
        s = 0.d0
        do i = 1, nxis
          i1 = xis(i)%i1
          i2 = xis(i)%i2
          s = s + rw_mesh(i1) * rw_mesh(i2) * &
              & radial_ho_wf(i1,n12,l12) * radial_ho_wf(i2,n3,l3) * &
              & radial_ho_wf(i1,n45,l45) * radial_ho_wf(i2,n6,l6) * &
              & v%v(i)
        end do
        mat%m(bra,ket) = s
      end do
    end do
  end function transform_xis_to_ho_element

  function transform_xis_to_ho(nxis,xis,chbra,chket,v) result(mat)
    integer, intent(in) :: nxis
    type(Coordinates), intent(in) :: xis(:)
    type(PWChan), intent(in) :: chbra, chket
    type(DVec), intent(in) :: v
    type(DMat) :: mat, m, ovlp_bra, ovlp_ket
    real(8) :: c
    integer :: i, i1, i2

    ovlp_bra = get_overlap_xis_ho(nxis, xis, chbra)
    ovlp_ket = get_overlap_xis_ho(nxis, xis, chket)
    m = ovlp_bra%t()
    !$omp parallel
    !$omp do private(i, i1, i2, c)
    do i = 1, v%n_size
      i1 = xis(i)%i1
      i2 = xis(i)%i2
      c = v%v(i) * rw_mesh(i1) * rw_mesh(i2)
      m%m(:,i) = m%m(:,i) * c
    end do
    !$omp end do
    !$omp end parallel
    mat = m * ovlp_ket
    call ovlp_bra%fin()
    call ovlp_ket%fin()
    call m%fin()
  end function transform_xis_to_ho

  function transform_xis_to_ho_old(nxis,xis,chbra,chket,v) result(mat)
    integer, intent(in) :: nxis
    type(Coordinates), intent(in) :: xis(:)
    type(PWChan), intent(in) :: chbra, chket
    type(DVec), intent(in) :: v
    type(DMat) :: mat, mat_r, ovlp_bra, ovlp_ket
    integer :: i, i1, i2

    ovlp_bra = get_overlap_xis_ho(nxis, xis, chbra)
    ovlp_ket = get_overlap_xis_ho(nxis, xis, chket)
    call mat_r%zeros(v%n_size, v%n_size)
    do i = 1, v%n_size
      i1 = xis(i)%i1
      i2 = xis(i)%i2
      mat_r%m(i,i) = v%v(i) * rw_mesh(i1) * rw_mesh(i2)
    end do
    mat = ovlp_bra%t() * mat_r * ovlp_ket
    call ovlp_bra%fin()
    call ovlp_ket%fin()
    call mat_r%fin()
  end function transform_xis_to_ho_old

  function get_overlap_xis_ho(nxis,xis,ch) result(mat)
    integer, intent(in) :: nxis
    type(Coordinates), intent(in) :: xis(:)
    type(PWChan), intent(in) :: ch
    type(DMat) :: mat, eye
    integer :: iHO, n12, n3, i, i1, i2

    call mat%ini(nxis, ch%Nho)
    !$omp parallel
    !$omp do private(iHO, n12, n3, i, i1, i2)
    do iHO = 1, ch%Nho
      n12 = ch%n12(iHO)
      n3  = ch%n3( iHO)

      do i = 1, nxis
        i1 = xis(i)%i1
        i2 = xis(i)%i2
        mat%m(i,iHO) = radial_ho_wf(i1,n12,ch%l12) * radial_ho_wf(i2,n3,ch%l3)
      end do
    end do
    !$omp end do
    !$omp end parallel

    return
    ! norm check
    call mat%ini(nxis, ch%Nho)
    call eye%zeros(nxis, nxis)
    !$omp parallel
    !$omp do private(iHO, n12, n3, i, i1, i2)
    do iHO = 1, ch%Nho
      n12 = ch%n12(iHO)
      n3  = ch%n3( iHO)

      do i = 1, nxis
        i1 = xis(i)%i1
        i2 = xis(i)%i2
        mat%m(i,iHO) = radial_ho_wf(i1,n12,ch%l12) * radial_ho_wf(i2,n3,ch%l3)
      end do
    end do
    !$omp end do
    !$omp end parallel
    do i = 1, nxis
      i1 = xis(i)%i1
      i2 = xis(i)%i2
      eye%m(i,i) = rw_mesh(i1) * rw_mesh(i2)
    end do
    eye = mat%t() * eye * mat
    call eye%prt("eye")
  end function get_overlap_xis_ho

  function tau1_dot_tau2(t) result(r)
    use MyLibrary, only: sjs
    integer, intent(in) :: t
    real(8) :: r
    r = sjs(1,1,2*t,1,1,2) * (-1.d0)**(t-1)
  end function tau1_dot_tau2

  function tau2_dot_tau3(t12, t45, Ttot) result(r)
    use MyLibrary, only: sjs
    integer, intent(in) :: t12, t45, Ttot
    real(8) :: r
    r = sqrt( dble(2*t12+1) * dble(2*t45+1) ) * &
        & sjs(2*t12, 2*t45, 2, 1, 1, 1) * &
        & sjs(2*t12, 2*t45, 2, 1, 1, Ttot) * &
        & (-1.d0) ** (t12+t45+(Ttot+1)/2)
  end function tau2_dot_tau3

  function tau_cross(t12, t45, Ttot) result(r)
    use MyLibrary, only: sjs, snj
    integer, intent(in) :: t12, t45, Ttot
    real(8) :: r
    r = sqrt( dble(2*t12+1) * dble(2*t45+1) ) * &
        & sjs(2*t12, 2*t45, 2, 1, 1, Ttot) * &
        & snj(1, 1, 2*t45, 1, 1, 2*t12, 2, 2, 2) * &
        & (-1.d0) ** (t45+(Ttot+1)/2)
  end function tau_cross

  subroutine FinPWChan(this)
    class(PWChan), intent(inout) :: this
    deallocate(this%n12)
    deallocate(this%n3)
    deallocate(this%ns2i)
  end subroutine FinPWChan

  subroutine InitPWChan(this, Nmax, l12, s12, j12, t12, l3, j3)
    class(PWChan), intent(inout) :: this
    integer, intent(in) :: Nmax, l12, s12, j12, t12, l3, j3
    integer :: cnt, N, n12, n3
    this%l12 = l12
    this%s12 = s12
    this%j12 = j12
    this%t12 = t12
    this%l3  = l3
    this%j3  = j3

    cnt = 0
    do N = 0, (Nmax-l12-l3)/2
      do n12 = 0, N
        n3 = N - n12
        cnt = cnt + 1
      end do
    end do

    this%Nho = cnt
    allocate(this%n12(this%Nho))
    allocate(this%n3(this%Nho))
    allocate(this%ns2i(0:Nmax/2, 0:Nmax/2))
    this%ns2i(:,:) = 0

    cnt = 0
    do N = 0, (Nmax-l12-l3)/2
      do n12 = 0, N
        n3 = N - n12
        cnt = cnt + 1
        this%n12(cnt) = n12
        this%n3( cnt) = n3
        this%ns2i(n12,n3) = cnt
      end do
    end do
  end subroutine InitPWChan

  subroutine FinModelSpace(this)
    class(ModelSpace), intent(inout) :: this
    integer :: ch
    do ch = 1, this%Nalpha
      call this%alpha(ch)%fin()
    end do
    deallocate(this%xis)
    deallocate(this%LSJTlj2idx)
  end subroutine FinModelSpace

  subroutine InitModelSpace(this, Jtot, Ptot, Ttot, Nmax, hw, Nmesh_in, xmin_in, xmax_in)
    use MyLibrary, only: triag, gauss_legendre
    class(ModelSpace), intent(inout) :: this
    integer, intent(in) :: Jtot, Ptot, Ttot, Nmax
    real(8), intent(in) :: hw
    integer, intent(in), optional :: NMesh_in
    real(8), intent(in), optional :: xmin_in, xmax_in
    integer :: j12, s12, l12, t12, l3, j3, cnt
    integer :: i1, i2, i

    this%Jtot = Jtot
    this%Ptot = Ptot
    this%Ttot = Ttot
    this%Nmax = Nmax
    this%hw = hw
    if(present(NMesh_in)) NMesh_r = NMesh_in
    if(present(xmin_in)) rmin = xmin_in
    if(present(xmax_in)) rmax = xmax_in
    allocate(this%LSJTlj2idx(0:Nmax,0:1,0:Nmax+1,0:1,0:Nmax,Nmax+1))
    this%LSJTlj2idx(:,:,:,:,:,:) = 0
    cnt = 0
    do j12 = 0, Nmax+1
      do s12 = 0, 1
        if(j12==Nmax+1 .and. s12==0) cycle
        do l12 = abs(j12-s12), (j12+s12)
          do t12 = 0, 1
            if((-1)**(l12+s12+t12) == 1) cycle
            if( triag(2*t12, 1, Ttot)) cycle
            do l3 = abs( abs(Jtot-2*j12)-1 )/2, (Jtot+2*j12+1)/2
              if(l12 + l3 > Nmax) cycle
              if((-1)**(l12+l3) /= Ptot) cycle
              do j3 = abs(2*l3-1), 2*l3+1, 2
                if( triag(j3, 2*j12, Jtot)) cycle
                cnt = cnt + 1
              end do
            end do
          end do
        end do
      end do
    end do
    this%Nalpha = cnt
    if(this%Nalpha < 1) return
    allocate(this%alpha(this%Nalpha))

    cnt = 0
    do j12 = 0, Nmax+1
      do s12 = 0, 1
        if(j12==Nmax+1 .and. s12==0) cycle
        do l12 = abs(j12-s12), (j12+s12)
          do t12 = 0, 1
            if((-1)**(l12+s12+t12) == 1) cycle
            if( triag(2*t12, 1, Ttot)) cycle
            do l3 = abs( abs(Jtot-2*j12)-1 )/2, (Jtot+2*j12+1)/2
              if(l12 + l3 > Nmax) cycle
              if((-1)**(l12+l3) /= Ptot) cycle
              do j3 = abs(2*l3-1), 2*l3+1, 2
                if( triag(j3, 2*j12, Jtot)) cycle
                cnt = cnt + 1
                this%LSJTlj2idx(l12,s12,j12,t12,l3,(j3+1)/2) = cnt
                call this%alpha(cnt)%init(Nmax, l12, s12, j12, t12, l3, j3)
              end do
            end do
          end do
        end do
      end do
    end do

    this%Nxi = NMesh_r
    this%Nxis = NMesh_r**2
    allocate(this%xis(this%Nxis))
    allocate(r_mesh(NMesh_r))
    allocate(rw_mesh(NMesh_r))
    call gauss_legendre(rmin, rmax, r_mesh, rw_mesh, NMesh_r)
    do i = 1, this%Nxis
      i1 = (i-1)/NMesh_r + 1
      i2 = mod(i-1,NMesh_r)+1
      this%xis(i)%i1 = i1
      this%xis(i)%i2 = i2
      this%xis(i)%x1 = r_mesh(i1)
      this%xis(i)%x2 = r_mesh(i2)
      this%xis(i)%w1 = rw_mesh(i1)
      this%xis(i)%w2 = rw_mesh(i2)
    end do
    deallocate(r_mesh)
    deallocate(rw_mesh)
  end subroutine InitModelSpace

  subroutine PrintModelSpace(ms)
    class(ModelSpace), intent(in) :: ms
    integer :: idx

    write(*,"(a,i2,a,i2,a,i2,a)") "Three-Body Model Space (Mom Space): J=", &
      & ms%Jtot, "/2, P=", ms%Ptot, ", ", ms%Ttot, "/2"

    write(*,"(a)") "angular-momentum channels: "
    write(*,"(7x,a,3x,a,3x,a,3x,a,3x,a,4x,a,2x,a,3x,a,3x,a,3x,a)") &
        &"idx", "l12", "s12", "j12", "t12", "l3", "2*j3", "2*J", "Par", "2*T"
    do idx = 1, ms%Nalpha
      write(*,"(4x,10i6)") idx, ms%alpha(idx)%l12, ms%alpha(idx)%s12, &
          & ms%alpha(idx)%j12, ms%alpha(idx)%t12, ms%alpha(idx)%l3, &
          & ms%alpha(idx)%j3, ms%Jtot, ms%Ptot, ms%Ttot
    end do
  end subroutine PrintModelSpace

  subroutine precalculations(Nmax,Jtot,hw,xis,lambda,power)
    use MyLibrary, only: triag, gauss_legendre, m_red_pn, hc, ho_radial_wf_norm, pi
    integer, intent(in) :: Nmax, Jtot, power
    real(8), intent(in) :: hw, lambda
    type(Coordinates), intent(in) :: xis(:)
    integer :: n, l, i
    real(8) :: nu
    ! ( j1 j2 | j12 )
    ! (  0  0 |   0 )
    call CGs%init(0,2*Nmax+2,0,2*Nmax+2)

    ! {   X l12' l12 }
    ! { S12 j12  j12'}
    call lsj12%init(0,2*Nmax+4,.false., 0,2*Nmax,.false., 0,2,.false., j12dmax_in=2*Nmax, j23dmax_in=2*Nmax+2)

    ! {   X l3' l3 }
    ! { 1/2 j3  j3'}
    call lsj3%init(0,2*Nmax+4,.false., 0,2*Nmax,.false., 1,1,.true., j12dmax_in=2*Nmax, j23dmax_in=2*Nmax+1)

    ! {   X j12 j12' }
    ! {Jtot j3'   j3 }
    call jjx%init(0,2*Nmax+4,.false., 0,2*Nmax+2,.false., Jtot,Jtot,.true., j12dmax_in=2*Nmax+2, j23dmax_in=2*Nmax+1)

    ! {  K1  K2  X }
    ! {   Y   R  Z }
    call kkxy%init(0,4,.false., 0,4,.false., 0,2*Nmax+4,.false.)

    ! {   X    K    Z }
    ! { l12  S12  j12 }
    ! { l12' S12' j12'}
    call ls12%init(0,2*Nmax+4,.false., 0,4,.false., 0,2*Nmax,.false., 0,2,.false., j13dmax_in=2*Nmax, jdmax_in=2*Nmax+2)

    ! {   X    K    Z }
    ! {  l3  1/2   j3 }
    ! {  l3' 1/2   j3'}
    call ls3%init(0,2*Nmax+4,.false., 0,4,.false., 0,2*Nmax,.false., 1,1,.true., j13dmax_in=2*Nmax, jdmax_in=2*Nmax+1)

    ! {  1/2 1/2   S1 }
    ! {  1/2 1/2   S2 }
    ! {   S3  S4    S }
    call spin12%init(1,1,.true., 1,1,.true., 1,1,.true., 1,1,.true.)

    call gauss_legendre(rmin, rmax, r_mesh, rw_mesh, NMesh_r)
    allocate(radial_ho_wf(NMesh_r, 0:Nmax/2, 0:Nmax))
    !nu = 2.d0 * amp * amn * hw / ( (amp + amn) * hc**2 )
    nu = 2.d0 * m_red_pn * hw / hc**2
    radial_ho_wf(:,:,:) = 0.d0
    do n = 0, Nmax/2
      do l = 0, Nmax
        do i = 1, NMesh_r
          radial_ho_wf(i,n,l) = ho_radial_wf_norm(n,l,nu,r_mesh(i))
        end do
      end do
    end do

    call gauss_legendre(pmin, pmax, p_mesh, pw_mesh, NMesh_p)
    call gauss_legendre(-1.d0, 1.d0, cos_mesh, cosw_mesh, NMesh_cos)
    call init_z0_function(lambda, power)
    call init_zx_function(lambda, power, xis, Nmax+2)
    call init_fk_function(lambda, power)
    call init_fkx_function(lambda, power, xis, Nmax+2)

  end subroutine precalculations

  subroutine release_arrays(ms)
    type(ModelSpace), intent(in) :: ms
    integer :: L, k
    call CGs%fin()
    call lsj12%fin()
    call lsj3%fin()
    call jjx%fin()
    call kkxy%fin()
    call ls12%fin()
    call ls3%fin()
    call spin12%fin()
    deallocate(r_mesh)
    deallocate(rw_mesh)
    deallocate(radial_ho_wf)
    deallocate(z0%v)
    do k = 0, 2
      do L = 0, ms%Nmax+2
        deallocate(fkx(k,L)%v)
      end do
      deallocate(fk(k)%v)
    end do
    do L = 0, ms%Nmax+2
      deallocate(zx(L)%v)
    end do
    deallocate(fkx)
    deallocate(fk)
    deallocate(zx)


  end subroutine release_arrays

  subroutine init_z0_function(lambda,power)
    !
    ! For details, see Eq. (13) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: spherical_bessel, pi, hc
    real(8), intent(in) :: lambda
    integer, intent(in) :: power
    integer :: i, j
    real(8) :: s, r, p, w
    allocate(z0%v(NMesh_r))
    z0%v(:) = 0.d0
    !$omp parallel
    !$omp do private(i, r, s, j, p, w)
    do i = 1, NMesh_r
      r = r_mesh(i) * sqrt(2.d0)
      s = 0.d0
      do j = 1, NMesh_p
        p = p_mesh(j)
        w = pw_mesh(j)
        s = s + w * p**2 * spherical_bessel(0, p*r) * local_regulator(p*hc, lambda, power)
      end do
      z0%v(i) = s / (2.d0 * pi**2)
    end do
    !$omp end do
    !$omp end parallel
  end subroutine init_z0_function

  subroutine init_zx_function(lambda,power,xis, L)
    !
    ! For details, see Eq. (16) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: spherical_bessel, pi, hc
    real(8), intent(in) :: lambda
    integer, intent(in) :: power, L
    type(Coordinates), intent(in) :: xis(:)
    integer :: i, j, x
    real(8) :: s, r1, r2, p, w, a

    allocate(zx(0:L))
    do x = 0, L
      allocate(zx(x)%v(size(xis)))
      zx(x)%v(:) = 0.d0
      a = exp(-200.d0 / dble(x) * log(10.d0) + dble(2*x+1)/dble(x) * log(dble(2*x+1)) - &
          & dble(2*x+1)/dble(x) - log(dble(x)) + 1 - log(2.d0) )
      !$omp parallel
      !$omp do private(i, r1, r2, s, j, p, w)
      do i = 1, size(xis)
        r1 = xis(i)%x1 / sqrt(2.d0)
        r2 = xis(i)%x2 * sqrt(1.5d0)

        s = 0.d0
        do j = 1, NMesh_p
          p = p_mesh(j)
          w = pw_mesh(j)
          if( x > 42 .and. (r1*p < 1.d-4 .or. r2*p < 1.d-4) ) cycle
          s = s + w * p**2 * spherical_bessel(x, r1*p) * spherical_bessel(x, r2*p) * &
              & local_regulator(p*hc, lambda, power)
        end do
        zx(x)%v(i) = s / (2.d0 * pi**2)
      end do
      !$omp end do
      !$omp end parallel

    end do
  end subroutine init_zx_function

  subroutine init_fk_function(lambda,power)
    !
    ! For details, see Eq. (39) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: spherical_bessel, pi, hc, m_pi
    real(8), intent(in) :: lambda
    integer, intent(in) :: power
    integer :: i, j, k, kk
    real(8) :: s, r, p, w, mpi
    mpi = m_pi / hc
    allocate(fk(0:2))
    kk = 0
    do k = 0, 2
      if(k==0 .or. k==2) kk = 4
      if(k==1) kk = 3
      allocate(fk(k)%v(NMesh_r))
      fk(k)%v(:) = 0.d0
      !$omp parallel
      !$omp do private(i, r, s, j, p, w)
      do i = 1, NMesh_r
        r = r_mesh(i) * sqrt(2.d0)
        s = 0.d0
        do j = 1, NMesh_p
          p = p_mesh(j)
          w = pw_mesh(j)
          s = s + w * p**kk * spherical_bessel(k, p*r) * &
              & local_regulator(p*hc, lambda, power) / (p**2 + mpi**2)
        end do
        fk(k)%v(i) = s / (2.d0 * pi**2)
      end do
      !$omp end do
      !$omp end parallel
    end do
  end subroutine init_fk_function

  subroutine init_fkx_function(lambda,power,xis, L)
    !
    ! For details, see Eqs. (34), (36), and (44) in P. Navrátil, Few-Body Syst. 41, 117 (2007).
    !
    use MyLibrary, only: spherical_bessel, pi, hc, m_pi, legendre_polynomial
    real(8), intent(in) :: lambda
    integer, intent(in) :: power, L
    type(Coordinates), intent(in) :: xis(:)
    integer :: i, j, x
    real(8) :: s, r1, r2, r, p, w, costh, wcosth, mpi

    mpi = m_pi / hc
    allocate(fkx(0:2,0:L))

    ! k = 0
    do x = 0, L
      allocate(fkx(0,x)%v( size(xis) ) )
      fkx(0,x)%v(:) = 0.d0

      !$omp parallel
      !$omp do private(i, r1, r2, s, j, p, w)
      do i = 1, size(xis)
        r1 = xis(i)%x1 / sqrt(2.d0)
        r2 = xis(i)%x2 * sqrt(1.5d0)

        s = 0.d0
        do j = 1, NMesh_p
          p = p_mesh(j)
          w = pw_mesh(j)
          s = s + w * p**4 * spherical_bessel(x, p*r1) * &
              & spherical_bessel(x, p*r2) * &
              & local_regulator(p*hc, lambda, power) / (p**2 + mpi**2)
        end do
        fkx(0,x)%v(i) = s / (2.d0 * pi**2)
      end do
      !$omp end do
      !$omp end parallel
    end do

    ! k = 1
    do x = 0, L
      allocate(fkx(1,x)%v( size(xis) ) )
      fkx(1,x)%v(:) = 0.d0

      !$omp parallel
      !$omp do private(i, r1, r2, s, costh, wcosth, r)
      do i = 1, size(xis)
        r1 = xis(i)%x1 / sqrt(2.d0)
        r2 = xis(i)%x2 * sqrt(1.5d0)

        s = 0.d0
        do j = 1, NMesh_cos
          costh = cos_Mesh(j)
          wcosth = cosw_Mesh(j)
          r = sqrt(r1**2 + r2**2 - 2.d0*r1*r2*costh)
          s = s + wcosth * legendre_polynomial(x,costh) * f1_func(r, lambda, power) / r
        end do
        fkx(1,x)%v(i) = s * 0.5d0
      end do
      !$omp end do
      !$omp end parallel
    end do

    ! k = 2
    do x = 0, L
      allocate(fkx(2,x)%v( size(xis) ) )
      fkx(2,x)%v(:) = 0.d0

      !$omp parallel
      !$omp do private(i, r1, r2, s, costh, wcosth, r)
      do i = 1, size(xis)
        r1 = xis(i)%x1 / sqrt(2.d0)
        r2 = xis(i)%x2 * sqrt(1.5d0)

        s = 0.d0
        do j = 1, NMesh_cos
          costh = cos_Mesh(j)
          wcosth = cosw_Mesh(j)
          r = sqrt(r1**2 + r2**2 - 2.d0*r1*r2*costh)
          s = s + wcosth * legendre_polynomial(x,costh) * f2_func(r, lambda, power) / r**2
        end do
        fkx(2,x)%v(i) = s * 0.5d0
      end do
      !$omp end do
      !$omp end parallel
    end do

  end subroutine init_fkx_function

  function local_regulator(q, lambda, n) result(f)
    real(8), intent(in) :: q, lambda
    integer, intent(in) :: n
    real(8) :: f, x
    x = ( q / lambda )**n
    f = exp( - x**2 )
  end function local_regulator

  function f1_func(r, lambda, n) result(f)
    use MyLibrary, only: spherical_bessel, pi, hc, m_pi
    real(8), intent(in) :: r, lambda
    integer, intent(in) :: n
    real(8) :: f
    real(8) :: p, w, mpi
    integer :: i

    mpi = m_pi / hc
    f = 0.d0
    do i = 1, NMesh_p
      p = p_mesh(i)
      w = pw_mesh(i)
      f = f + w * p**3 * spherical_bessel(1, r*p) * &
          & local_regulator(p*hc, lambda, n) / (p**2 + mpi**2)
    end do
    f = f / (2.d0 * pi**2)
  end function f1_func

  function f2_func(r, lambda, n) result(f)
    use MyLibrary, only: spherical_bessel, pi, hc, m_pi
    real(8), intent(in) :: r, lambda
    integer, intent(in) :: n
    real(8) :: f
    real(8) :: p, w, mpi
    integer :: i

    mpi = m_pi / hc
    f = 0.d0
    do i = 1, NMesh_p
      p = p_mesh(i)
      w = pw_mesh(i)
      f = f + w * p**4 * spherical_bessel(2, r*p) * &
          & local_regulator(p*hc, lambda, n) / (p**2 + mpi**2)
    end do
    f = f / (2.d0 * pi**2)
  end function f2_func
end module NNNForceLocal
