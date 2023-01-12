module NNForceIsospin
  use omp_lib
  use Profiler, only: timer
  use ClassSys
  use LinAlgLib
  use MPIfunction, only: myrank
  use TwoBodyRelativeSpace
  implicit none

  public :: NNForceHOIsospin

  private :: FinNNForceHOIsospin
  private :: InitNNForceHOIsospin
  private :: SetNNForceHOIsospin
  private :: check_nnforce_iso_ho_space
  private :: renorm_ho_space_isospin
  private :: HOSRGChannel
  private :: LeeSuzukiChannel

  type :: NNForceHOIsospin
    type(DMat), allocatable :: MatCh(:)
    type(TwoBodyRelSpaceSpinIsoHOBasis), pointer :: ms
    logical :: is_init = .false.
  contains
    procedure :: InitNNForceHOIsospin
    procedure :: FinNNForceHOIsospin
    procedure :: SetNNForceHOIsospin

    generic :: init => InitNNForceHOIsospin
    generic :: fin => FinNNForceHOIsospin
  end type NNForceHOIsospin
contains
  subroutine FinNNForceHOIsospin(this)
    class(NNForceHOIsospin), intent(inout) :: this
    integer :: ich
    do ich = 1, this%ms%NChan
      call this%MatCh(ich)%fin()
    end do
    deallocate(this%MatCh)
    this%ms => null()
    this%is_init = .false.
  end subroutine FinNNForceHOIsospin

  subroutine InitNNForceHOIsospin(this, ms)
    class(NNForceHOIsospin), intent(inout) :: this
    type(TwoBodyRelSpaceSpinIsoHOBasis), target :: ms
    type(TwoBodyRelChanSpinIsoHOBasis), pointer :: tbc
    integer :: ch, n

    this%ms => ms
    allocate(this%MatCh(ms%NChan))
    do ch = 1, ms%NChan
      tbc => ms%GetChannel(ch)
      n = tbc%GetNumberStates()
      call this%MatCh(ch)%zeros(n,n)
    end do
    this%is_init = .true.
  end subroutine InitNNForceHOIsospin

  subroutine SetNNForceHOIsospin(this, U, params, T_total, Tz_total, verbose)
    !
    ! Note: T_total and Tz_total are double
    !
    use NNForce, only: NNForceHO
    use NuHamilInput, only: InputParameters
    class(NNForceHOIsospin), intent(inout) :: this
    type(NNForceHOIsospin), intent(inout) :: U
    type(InputParameters), intent(in) :: params
    integer, intent(in), optional :: T_total, Tz_total
    logical, intent(in), optional :: verbose
    type(TwoBodyRelSpaceSpinIsoHOBasis), pointer :: two
    type(TwoBodyRelChanSpinIsoHOBasis), pointer :: tbc
    integer :: na, tt, tz, np, nn
    real(8) :: cpp, cnn, cpn
    type(NNForceHO) :: vpn, Upn
    type(TwoBodyRelSpaceSpinHOBasis) :: twopn
    type(TwoBodyRelChanSpinHOBasis), pointer :: tbc_pn, tbc_pp, tbc_nn
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    type(InputParameters) :: input
    integer :: ich, ichpp, ichnn, ichpn
    integer :: j, p, s, t, n
    integer :: bra, ket
    integer :: n1, l1, n2, l2
    integer :: bnpn, knpn, bnpp, knpp, bnnn, knnn

    if(.not. this%is_init) then
      write(*,*) "Initialization error in ", __FILE__, " at ", __LINE__
      return
    end if

    if(.not. U%is_init) then
      write(*,*) "Initialize NNForceHOIsospin before calling SetNNForceHOIsospin"
      write(*,*) "Initialization error in ", __FILE__, " at ", __LINE__
      return
    end if

    cpp = 1.d0 / 3.d0
    cnn = 1.d0 / 3.d0
    cpn = 1.d0 / 3.d0

    if(present(T_total) .and. .not. present(Tz_total)) then
      write(*,*) "Argument warning in ", __FILE__, " at ", __LINE__
    end if

    if(.not. present(T_total) .and. present(Tz_total)) then
      write(*,*) "Argument warning in ", __FILE__, " at ", __LINE__
    end if
    if(present(T_total) .and. present(Tz_total)) then
      NA = params%particle_rank
      tt = T_total
      np = (NA - Tz_total) / 2
      nN = (NA + Tz_total) / 2
      tz = nn - np
      cpp = dble((na - tz) * (na - tz - 2)) / dble( 3 * na * (na - 2) + tt * (tt + 2) )
      cnn = dble((na + tz) * (na + tz - 2)) / dble( 3 * na * (na - 2) + tt * (tt + 2) )
      cpn = dble(na * (na - 2) + tt * (tt + 2) - 2 * tz ** 2) / dble( 3 * na * (na - 2) + tt * (tt + 2) )
    end if
    if(cpp< -1.d-8 .or. cnn < -1.d-8 .or. cpn < -1.d-8) then
      write(*,*) cpp, cnn, cpn
      write(*,*) "Error in ", __FILE__, " at ", __LINE__
      stop
    end if

    !write(*,"(a)") "# ratio NN int of T=1 chann: "
    !write(*,"(a,f8.6,a,f8.6,a,f8.6)") "# pp: ", cpp, ", nn: ", cnn, ", pn: ", cpn

    two => this%ms
    call twopn%init(two%hw, two%Nmax, two%Jmax)

    call vpn%init(twopn)
    call Upn%init(twopn)
    input = params
    if(input%renorm_space2%val == 'ho') input%renorm = "bare"
    call vpn%SetNNForceHO(Upn, input)
    call Upn%fin()


    do ich = 1, two%GetNumberChannels()
      tbc => two%GetChannel(ich)
      n = tbc%GetNumberStates()
      call U%MatCh(ich)%eye(n)
    end do

    do ich = 1, two%NChan
      tbc => two%GetChannel(ich)
      j = tbc%GetJ()
      p = tbc%GetParity()
      s = tbc%GetS()
      t = tbc%GetT()

      if(t == 0) then
        ichpn = twopn%GetIndex(j,p,s,0)
        if( ichpn == 0 ) cycle
        tbc_pn => twopn%GetChannel(ichpn)
        !$omp parallel
        !$omp do private(bra, ho_bra, n1, l1, bnpn, ket, ho_ket, n2, l2, knpn)
        do bra = 1, tbc%GetNumberStates()
          ho_bra => tbc%getp(bra)
          n1 = ho_bra%GetN()
          l1 = ho_bra%GetL()
          bnpn = tbc_pn%GetIndex(n1,l1)
          do ket = 1, bra
            ho_ket => tbc%getp(ket)
            n2 = ho_ket%GetN()
            l2 = ho_ket%GetL()
            knpn = tbc_pn%GetIndex(n2,l2)

            this%MatCh(ich)%m(bra,ket) = vpn%MatCh(ichpn)%m(bnpn,knpn)
            this%MatCh(ich)%m(ket,bra) = vpn%MatCh(ichpn)%m(bnpn,knpn)

          end do
        end do
        !$omp end do
        !$omp end parallel
      elseif(t == 1) then
        ichpp = twopn%GetIndex(j,p,s,-1)
        ichpn = twopn%GetIndex(j,p,s, 0)
        ichnn = twopn%GetIndex(j,p,s, 1)
        tbc_pp => twopn%GetChannel(ichpp)
        tbc_pn => twopn%GetChannel(ichpn)
        tbc_nn => twopn%GetChannel(ichnn)

        !$omp parallel
        !$omp do private(bra, ho_bra, n1, l1, bnpp, bnpn, bnnn, ket, &
        !$omp &  ho_ket, n2, l2, knpp, knpn, knnn)
        do bra = 1, two%jpst(ich)%GetNumberStates()
          ho_bra => tbc%getp(bra)
          n1 = ho_bra%GetN()
          l1 = ho_bra%GetL()
          bnpp = tbc_pp%GetIndex(n1,l1)
          bnpn = tbc_pn%GetIndex(n1,l1)
          bnnn = tbc_nn%GetIndex(n1,l1)
          do ket = 1, bra
            ho_ket => tbc%getp(ket)
            n2 = ho_ket%GetN()
            l2 = ho_ket%GetL()
            knpp = tbc_pp%GetIndex(n2,l2)
            knpn = tbc_pn%GetIndex(n2,l2)
            knnn = tbc_nn%GetIndex(n2,l2)

            this%MatCh(ich)%m(bra,ket) = (vpn%MatCh(ichpp)%m(bnpp,knpp) * cpp + &
                & vpn%MatCh(ichpn)%m(bnpn,knpn) * cpn + &
                & vpn%MatCh(ichnn)%m(bnnn,knnn) * cnn)
            this%MatCh(ich)%m(ket,bra) = this%MatCh(ich)%m(bra,ket)

          end do
        end do
        !$omp end do
        !$omp end parallel
      end if
    end do
#ifdef NNForceDebug
    call check_nnforce_iso_ho_space(this, U, two%GetFrequency())
#endif
    if(params%renorm%val /= 'bare' .and. params%renorm_space2%val == 'ho') then
      call renorm_ho_space_isospin(this, U, params%renorm, params%lambda, &
          & params%particle_rank, params%N_ls, params%srg_generator, params%Nmax_srg_edge, verbose)
#ifdef NNForceDebug
      call check_nnforce_iso_ho_space(this, U, two%GetFrequency())
#endif
    end if
  end subroutine SetNNForceHOIsospin

  subroutine check_nnforce_iso_ho_space(vnn, U, hw, verbose)
    use MyLibrary, only: hc, m_red_pn
    use TwoBodyRelativeSpace, only: TwoBodyRelSpaceSpinIsoHOBasis
    type(NNForceHOIsospin), intent(in) :: vnn, U
    type(TwoBodyRelSpaceSpinIsoHOBasis), pointer :: two
    type(TwoBodyRelChanSpinIsoHOBasis), pointer :: tbc
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    real(8), intent(in) :: hw
    logical, intent(in), optional :: verbose
    integer :: n2max, j, p, t, s, ich, n, m
    integer :: bra, ket, n1, l1, n2, l2
    real(8) :: kine, r, f
    type(DMat) :: h, tkin, rm
    type(DVec) :: v
    type(EigenSolSymD) :: sol
    two => vnn%ms
    N2max = two%GetNmax()
    j = 1
    p = 1
    t = 0
    f = hc ** 2 / (m_red_pn * hw)
    do s = 0, 1
      ich = two%GetIndex(j, p, s, t)
      if(ich < 1) cycle
      tbc => two%GetChannel(ich)
      n = two%jpst(ich)%GetNumberStates()
      call H%zeros(n,n); call tkin%zeros(n,n)
      call rm%zeros(n,n)
      do bra = 1, n
        ho_bra => tbc%getp(bra)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        if(2 * n1 + l1 > n2max) cycle
        do ket = 1, n
          ho_ket => tbc%getp(ket)
          n2 = ho_ket%GetN()
          l2 = ho_ket%GetL()
          if(2 * n2 + l2 > n2max) cycle
          kine = 0.d0
          if(l1 /= l2) cycle
          if(abs(n1 - n2) > 1) cycle
          if(n1 == n2) kine = dble(2 * n1 + l1) + 1.5d0
          if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
          if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          tkin%m(bra, ket) = kine * hw * 0.5d0

          r = 0.d0
          if(n1 == n2) r = dble(2 * n1 + l1) + 1.5d0
          if(n1 == n2-1) r = -dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
          if(n1 == n2+1) r = -dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          rm%m(bra, ket) = r * f
        end do
      end do
      H = tkin + vnn%MatCh(ich)
      rm = U%MatCh(ich)%T() * rm * U%MatCh(ich)
      call tkin%fin()
      call sol%init(h)
      call sol%DiagSym(h)
      if(myrank == 0) then
        m = min(5, n)
        call v%zeros(m)
        v%v(:) = sol%eig%v(:m)
        if(present(verbose)) write(*,'(a)') 'd Eigen values in CheckNNForceHOSpace (Isospin form)'
        if(present(verbose)) write(*,'(a,i3,a,i3,a,i3,a,i4)') 'J = ', j, ',  P = ', p, ',  T = ', t, ",  Nmax = ", tbc%GetNmax()
        call v%prt()
        call v%fin()
        if(present(verbose)) write(*, '(a, f12.6)') 'rm = ', &
            & sqrt(dot_product(sol%vec%m(:,1), matmul(rm%m, sol%vec%m(:,1))) / 4.d0)
      end if
      call sol%fin()
      call H%fin()
    end do
  end subroutine check_nnforce_iso_ho_space

  subroutine renorm_ho_space_isospin(vnn, Trs, renorm, lambda, A, N_ls, generator, Nmax_srg_edge, verbose)
    use MyLibrary, only: hc, m_red_pn
    class(NNForceHOIsospin), intent(inout) :: vnn
    type(NNForceHOIsospin), intent(inout) :: Trs
    type(str), intent(in) :: renorm
    real(8), intent(in) :: lambda
    integer, intent(in) :: A, N_ls, Nmax_srg_edge
    type(str), intent(in) :: generator
    logical, intent(in), optional :: verbose
    type(TwoBodyRelSpaceSpinIsoHOBasis), pointer :: two
    integer :: ich, Nmax
    real(8) :: lc, alpha, hw

    if(myrank == 0) then
      if(present(verbose)) write(*,'(a)') '## NN (isospin formalism) Renormalization in HO space ...'
    end if
    two => vnn%ms
    Nmax = two%Nmax
    hw = two%hw
    if(renorm%val == 'srg') then
      lc = lambda * hc / dsqrt(2.d0 * m_red_pn)
      alpha = (1.d0 / lc) ** (4.d0)
      do ich = 1, two%NChan
        call HOSRGChannel(vnn%MatCh(ich),Trs%MatCh(ich),vnn%ms%jpst(ich),alpha,hw,generator, Nmax_srg_edge)
      end do
    elseif(renorm%val == 'ls') then
      do ich = 1, two%NChan
        call LeeSuzukiChannel(vnn%MatCh(ich),Trs%MatCh(ich),vnn%ms%jpst(ich),hw, N_ls, A)
      end do
    end if
  end subroutine renorm_ho_space_isospin

  subroutine HOSRGChannel(v, u, ms, alpha, hw, generator_type, Nmax_srg_edge)
    use OperatorDefinitions, only: CalcMERel, OperatorDef
    use Renormalization, only: SRGSolver
    type(DMat), intent(inout) :: v, u
    type(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: ms
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    real(8), intent(inout) :: alpha, hw
    type(str), intent(in) :: generator_type
    integer, intent(in) :: Nmax_srg_edge
    integer :: n, n_zero
    integer :: bra, ket, n1, l1, n2, l2
    real(8) :: kine
    type(DMat) :: h, t, generator
    type(SRGSolver) :: sol
    integer, allocatable :: rhs_zero_term(:)
    type(OperatorDef) :: Gen

    n = ms%GetNumberStates()
    if(n < 1) return
    call t%zeros(n,n)
    do bra = 1, n
      ho_bra => ms%getp(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        ho_ket => ms%getp(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        kine = 0.d0
        if(l1 /= l2) cycle
        if(abs(n1 - n2) > 1) cycle
        if(n1 == n2) kine = dble(2 * n1 + l1) + 1.5d0
        if(n1 == n2-1) kine = dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
        if(n1 == n2+1) kine = dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
        t%m(bra, ket) = kine * hw * 0.5d0
      end do
    end do

    ! Generator
    call generator%zeros(n,n)
    call Gen%initOpDef(generator_type, .false.)
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        ho_ket => ms%GetP(ket)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        generator%m(bra,ket) = CalcMERel(Gen, [n1,l1,ms%GetS(),ms%GetJ(),ms%GetT()], [n2,l2,ms%GetS(),ms%GetJ(),ms%GetT()])
      end do
    end do

    ! Neumann boundary for SRG flow
    n_zero = 0
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if( 2*n1+l1 < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
    end do
    allocate(rhs_zero_term(n_zero))
    n_zero = 0
    do bra = 1, n
      ho_bra => ms%GetP(bra)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if( 2*n1+l1 < ms%GetNmax()-Nmax_srg_edge ) cycle
      n_zero = n_zero + 1
      rhs_zero_term(n_zero) = bra
    end do

    h = t + v
    !call sol%init(h, 'HUflow')
    call sol%init(h, 'Hflow')
    if( n_zero > 0 ) call sol%SRGFlow(h, generator, alpha, rhs_zero_term)
    if( n_zero == 0) call sol%SRGFlow(h, generator, alpha)
    v = sol%h - t
    u = sol%U

    call sol%fin()
    call generator%fin()
    call t%fin()
    call h%fin()
  end subroutine HOSRGChannel

  subroutine LeeSuzukiChannel(vin, u, ms, hw, N_ls, A)
    use Renormalization, only: LeeSuzukiSolver
    type(DMat), intent(inout) :: vin, u
    type(TwoBodyRelChanSpinIsoHOBasis), intent(in) :: ms
    type(HarmonicOscillator), pointer :: ho_bra, ho_ket
    real(8), intent(in) :: hw
    integer, intent(in) :: N_ls, A
    integer :: n
    integer :: n1, l1, bra
    integer :: n2, l2, ket
    integer :: np, nq
    integer :: NN, cnt, cntb, cntk
    integer, allocatable :: reord(:), back(:)
    real(8) :: rpot
    type(DMat) :: h, hd, v, ho
    type(LeeSuzukiSolver) :: sol

    n = ms%GetNumberStates()
    if(n < 1) return
    allocate(reord(n), back(n))
    cnt = 0
    do NN = 0, ms%GetNmax()
      do bra = 1, n
        ho_bra => ms%getp(bra)
        n1 = ho_bra%GetN()
        l1 = ho_bra%GetL()
        if(NN /= 2 * n1 + l1) cycle
        cnt = cnt + 1
        reord(bra) = cnt
        back(cnt) = bra
      end do
    end do
    np = 0; nq = 0
    do bra = 1, n
      cnt = reord(bra)
      ho_bra => ms%getp(cnt)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      if(2 * n1 + l1 <= N_ls) np = np + 1
      if(2 * n1 + l1 >  N_ls) nq = nq + 1
    end do

    call hd%zeros(n,n); call h%zeros(n,n); call v%zeros(n,n)
    call ho%zeros(n,n)

    do bra = 1, n
      cntb = reord(bra)
      ho_bra => ms%getp(cntb)
      n1 = ho_bra%GetN()
      l1 = ho_bra%GetL()
      do ket = 1, n
        cntk = reord(ket)
        ho_ket => ms%getp(cntk)
        n2 = ho_ket%GetN()
        l2 = ho_ket%GetL()
        rpot = 0.d0
        if(l1 /= l2) then
          if(abs(n1 - n2) > 1) then
            if(n1 == n2) rpot = dble(2 * n1 + l1) + 1.5d0
            if(n1 == n2-1) rpot = -dsqrt(dble(n1 + 1) * (dble(n1 + l1) + 1.5d0))
            if(n1 == n2+1) rpot = -dsqrt(dble(n2 + 1) * (dble(n2 + l2) + 1.5d0))
          end if
        end if
        ho%m(bra, ket) = - rpot * hw / dble(A)
        v%m(bra, ket) = vin%m(cntb, cntk)
      end do
      hd%m(bra, bra) = dble(2 * n1 + l1) * hw
    end do

    h = hd + v + ho
    call sol%init(h)
    call sol%LeeSuzuki(0, np, nq, h)
    v = sol%h - hd - ho
    do bra = 1, n
      cntb = back(bra)
      do ket = 1, n
        cntk = back(ket)
        vin%m(bra,ket) = v%m(cntb, cntk)
        U%m(bra,ket) = sol%expS%m(cntb, cntk)
      end do
    end do
    call sol%fin()
    call h%fin(); call hd%fin(); call v%fin()
    deallocate(reord, back)
  end subroutine LeeSuzukiChannel
end module NNForceIsospin

!program test
!  use Profiler, only: timer
!  use TwoBodyRelativeSpace
!  use NNForceIsospin
!  use MyLibrary
!  !type(TwoBodyRelativeSpaceSpinMeshBasis), target :: two
!  type(TwoBodyRelSpaceSpinIsoHOBasis) :: two
!  type(NNForceHOIsospin) :: vho, u
!
!  call timer%init()
!
!  call two%init(20.d0,20,1)
!  call vho%init(two)
!  call u%init(two)
!  call vho%SetParams("N3LO_EMN500","srg","ho",2.d0,.True.,200,8)
!  call vho%SetNNForceHOIsospin(u, 500, 0.d0, 25.d0, .True., 2, 20, 0, 0, &
!      & 1.d0/3.d0, 1.d0/3.d0, 1.d0/3.d0)
!
!  call u%fin()
!  call vho%fin()
!  call two%fin()
!
!  call timer%fin()
!end program test
