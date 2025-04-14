module OneBodyLabOps
  use omp_lib
  use ClassSys
  use MyLibrary, only: g_A, m_pi, g_pi, f_pi, m_nucleon
  use Profiler, only: timer
  use LinAlgLib
  use SingleParticleState
  use OperatorDefinitions
  implicit none

  private
  public :: OneBodyLabOp

  type, extends(OperatorDef) :: OneBodyLabOp
    type(DMat) :: mat
    type(Orbits), pointer :: sps
    real(8) :: hw
  contains
    procedure :: InitOneBodyLabOp
    procedure :: InitOneBodyLabOpFromString
    procedure :: FinOneBodyLabOp
    procedure :: SetOneBodyLabOp
    procedure :: GetFileName
    procedure :: WriteOperator
    procedure :: ReducedToNormal
    procedure :: NormalToReduced
    !procedure :: ReadOperator
    generic :: init => InitOneBodyLabOp, InitOneBodyLabOpFromString
    generic :: fin => FinOneBodyLabOp
    generic :: set => SetOneBodyLabOp
  end type OneBodyLabOp

  integer, parameter :: tz_phase = -1 ! particle physic convention -> nuclear physics convention
  real(8), allocatable :: rnl_bra(:,:,:), rnl_ket(:,:,:,:), pbra(:), pket(:,:), wbra(:), wket(:,:)
  integer :: NMesh = 120
  real(8), parameter :: kmin = 0.d0, kmax = 10.d0
  !logical, parameter :: use_parametric_form_factor = .true.
  real(8), parameter :: Lambda_A = 1040.d0 ! Phys. Rev. D 86, 103511 (2012)
  real(8), parameter :: d18 = g_A * (1.d0 - g_pi * f_pi / (m_nucleon * g_A)) / m_pi**2 ! This should be small, because of Goldberger-Treimann relation
  real(8), parameter :: d22 = 2.d0 * g_A / Lambda_A**2
  !real(8), parameter :: d18 = 0.d0
  !real(8), parameter :: d22 = 0.d0
  logical, parameter :: use_parametric_form_factor = .false.
  complex(8), parameter :: i_unit = (0.d0, 1.d0)
contains
  subroutine FinOneBodyLabOp(this)
    class(OneBodyLabOp), intent(inout) :: this
    call this%mat%fin()
    this%sps => null()
  end subroutine FinOneBodyLabOp

  subroutine InitOneBodyLabOpFromString(this, sps, hw, opname)
    class(OneBodyLabOp), intent(inout) :: this
    type(Orbits), intent(in), target :: sps
    real(8), intent(in) :: hw
    type(str), intent(in) :: opname
    call this%InitOpDef(opname, .true.)
    call this%init(sps, this%GetOpJ(), this%GetOpP(), this%GetOpZ(), hw)
    if(this%GetOpJ()/=0 .or. this%GetOpP()/=1 .or. this%GetOpZ()/=0) then
      call this%SetReduced(.false.)
    end if
  end subroutine InitOneBodyLabOpFromString

  subroutine InitOneBodyLabOp(this, sps, jr, pr, tzr, hw)
    class(OneBodyLabOp), intent(inout) :: this
    type(Orbits), intent(in), target :: sps
    integer, intent(in) :: jr, pr, tzr
    real(8), intent(in) :: hw
    this%sps => sps
    this%hw = hw
    call this%InitOpDef(.true., jr, pr, tzr)
    call this%mat%zeros(sps%norbs, sps%norbs)
    if(this%GetOpJ()/=0 .or. this%GetOpP()/=1 .or. this%GetOpZ()/=0) then
      call this%SetReduced(.false.)
    end if
  end subroutine InitOneBodyLabOp

  subroutine SetOneBodyLabOp(this, isospin, NMeshMultipole)
    class(OneBodyLabOp), intent(inout) :: this
    integer, intent(in), optional :: isospin
    integer, intent(in), optional :: NMeshMultipole
    type(Orbits), pointer :: sps
    integer :: a, b
    type(SingleParticleOrbit), pointer :: oa, ob
    type(sys) :: s
    type(str) :: opname

    if(present(NMeshMultipole)) NMesh = NMeshMultipole
    if( s%find(this%GetOpName(), s%str('L5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('M5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tel5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tmag5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('M_1B')) .or. &
      & s%find(this%GetOpName(), s%str('L_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tel_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tmag_1B')) ) then
      call this%SetReduced(.true.)
      call set_multipole(this, isospin=isospin)
      return
    end if

    sps => this%sps
    do a = 1, sps%norbs
      oa => sps%GetOrbit(a)
      do b =1, sps%norbs
        ob => sps%GetOrbit(b)
        this%mat%m(a,b) = CalcMEOneBody( this%OperatorDef, &
            & [oa%n,oa%l,oa%j,oa%z], [ob%n,ob%l,ob%j,ob%z], isospin)
      end do
    end do
  end subroutine SetOneBodyLabOp

  function GetFileName(this, filename) result(f)
    class(OneBodyLabOp), intent(in) :: this
    type(str), intent(in) :: filename
    type(sys) :: s
    type(str) :: f, OpName

    OpName = this%GetOpName()
    if(filename%val /= "default") then
      f = OpName + "_" + filename
      return
    end if
    f = OpName + '_OBME'
    !f = f + '_hw' + s%str(this%hw) + '_emax' + s%str(this%sps%emax) + ".me2j.gz"
    f = f + '_hw' + s%str(this%hw) + '_emax' + s%str(this%sps%emax) + ".snt"
  end function GetFileName

  subroutine WriteOperator(this, f)
    class(OneBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti

    ti = omp_get_wtime()
    if(this%GetOpJ()==0 .and. this%GetOpP()==1 .and. this%GetOpZ()==0 .and. this%reduced_me()) then
      call this%ReducedToNormal()
    end if
    call write_operator(this, f)
    call timer%Add(s%str('Write to file'), omp_get_wtime() - ti)
  end subroutine WriteOperator

  subroutine write_operator(this, f)
    class(OneBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(sys) :: s

    if(s%find(f,s%str(".snt"))) then
      call write_operator_kshell_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".snt"))) then
      call write_operator_snt(this, f)
      return
    end if

    if(s%find(f,s%str(".me2j.gz"))) then
      call write_operator_gzip(this, f)
      return
    end if

    write(*,"(a)") "Unknown file format, writing file assuming readable snt format"
    call write_operator_snt(this, f)
  end subroutine write_operator

  subroutine write_operator_snt(this, f)
    class(OneBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(Orbits), pointer :: sps
    integer :: a, b, cnt
    integer :: wunit=22
    character(255) :: header
    type(str) :: OpName
    sps => this%sps
    open(wunit, file=f%val, status='replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(wunit,"(a)") trim(header)
    write(wunit,"(4i5)") sps%norbs/2, sps%norbs/2, 0, 0
    do a = 1, sps%norbs
      write(wunit, '(5i5)') a, sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j, sps%orb(a)%z
    end do
    write(wunit,'(a)') '####### one-body term'
    cnt = 0
    do a = 1, sps%norbs
      do b = 1, a
        if(abs(this%mat%m(a,b)) < 1.d-50) cycle
        cnt = cnt + 1
      end do
    end do
    write(wunit,'(i5, i3)') cnt, 0
    write(wunit,'(a)') '### a, b, <a| op |b> or < a|| op || b >'
    do a = 1, sps%norbs
      do b = 1, a
        if(abs(this%mat%m(a,b)) < 1.d-50) cycle
        write(wunit, '(2i5, es18.8)') a, b, this%mat%m(a,b)
      end do
    end do
    write(wunit, '(i15,i3)') 0, 0
    write(wunit, '(a)') '##### two-body term #####'
    close(wunit)
  end subroutine write_operator_snt

  subroutine write_operator_kshell_snt(this, f)
    class(OneBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(Orbits), pointer :: sps
    type(Orbits) :: sps_k
    type(SingleParticleOrbit), pointer :: oa, ob
    integer :: a_k, b_k, a, b, cnt
    integer :: wunit=22
    character(255) :: header
    type(str) :: OpName

    sps => this%sps
    call sps_k%init(sps%emax, mode="kshell")
    open(wunit, file = f%val, status = 'replace')
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    write(wunit,"(a)") trim(header)
    write(wunit,'(4i5)') sps_k%norbs/2, sps_k%norbs/2, 0, 0
    do a = 1, sps_k%norbs
      write(wunit, '(5i5)') a, sps_k%orb(a)%n, sps_k%orb(a)%l, sps_k%orb(a)%j, sps_k%orb(a)%z
    end do
    cnt = 0
    do a_k = 1, sps_k%norbs
!      do b_k = 1, a_k
      do b_k = 1, sps_k%norbs
        oa => sps_k%GetOrbit(a_k)
        ob => sps_k%GetOrbit(b_k)
        a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
        b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
        if(abs(this%mat%m(a,b)) < 1.d-50) cycle
        cnt = cnt + 1
      end do
    end do
    write(wunit,'(i5, i3)') cnt, 0
    write(wunit,'(a)') '### a, b, <a| op |b> or < a|| op || b >'
    do a_k = 1, sps_k%norbs
!      do b_k = 1, a_k
      do b_k = 1, sps_k%norbs
        oa => sps_k%GetOrbit(a_k)
        ob => sps_k%GetOrbit(b_k)
        a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
        b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
        if(abs(this%mat%m(a,b)) < 1.d-50) cycle
        write(wunit, '(2i5, es18.8)') a_k, b_k, this%mat%m(a,b)
      end do
    end do
    write(wunit, '(i15,i3)') 0, 0
    write(wunit, '(a)') '##### two-body term #####'
    close(wunit)
    call sps_k%fin()
  end subroutine write_operator_kshell_snt

  subroutine write_operator_gzip(this, f)
    use, intrinsic :: iso_c_binding
    use MyLibrary, only: gzip_open, gzip_writeline, gzip_close, triag
    class(OneBodyLabOp), intent(in) :: this
    type(str), intent(in) :: f
    type(Orbits), pointer :: sps
    type(OrbitsIsospin) :: remap
    integer :: a, b
    type(SingleParticleOrbitIsospin), pointer :: oa, ob
    character(255) :: header
    integer :: ap, bp
    integer :: an, bn
    character(kind=c_char, len=100) :: bufferm
    character(kind=c_char, len=80) :: buffer1
    type(c_ptr) :: fp, err
    type(str) :: OpName
    sps => this%sps
    call remap%init(sps%emax, sps%lmax)
    OpName = this%GetOpName()
    header = "# " // OpName%val // " operator generated by NuHamil (Tokyo code), "
#ifdef VERSION
    header = trim(header) // trim(VERSION)
#endif
    fp = gzip_open(f%val,"wt")
    err = gzip_writeline(fp, trim(header), len_trim(header))

    write(bufferm,"(6i4)") this%GetOpJ(), this%GetOpP(), this%GetOpZ(), sps%emax, -1, sps%lmax
    err = gzip_writeline(fp, trim(bufferm), len_trim(bufferm))

    write(bufferm,'(f16.8)') 0.d0
    err = gzip_writeline(fp, trim(bufferm), len_trim(bufferm))

    do a = 1, remap%norbs
      oa => remap%GetOrbit(a)
      ap = remap%iso2pn(sps,a,-1)
      an = remap%iso2pn(sps,a, 1)
      do b = 1, remap%norbs
        ob => remap%GetOrbit(b)
        bp = remap%iso2pn(sps,b,-1)
        bn = remap%iso2pn(sps,b, 1)
        if((-1)**(oa%l+ob%l) * this%GetOpP() /= 1) cycle
        if( triag(oa%j,ob%j,2*this%GetOpJ() ) ) cycle
        write(buffer1,"(4es18.8)") this%mat%m(ap,bp), &
            & this%mat%m(an,bn), this%mat%m(an,bp), this%mat%m(ap,bn)
        err = gzip_writeline(fp, trim(buffer1), len_trim(buffer1))
      end do
    end do
    err = gzip_close(fp)
    call remap%fin()
  end subroutine write_operator_gzip

  subroutine ReducedToNormal(this)
    class(OneBodyLabOp), intent(inout) :: this
    integer :: i, j
    type(SingleParticleOrbit), pointer :: oi

    if(this%GetOpJ()/=0 .or. this%GetOpP()/=1 .or. this%GetOpZ()/=0) return
    if(.not. this%reduced_me()) return
    write(*,*) "OneBody: reduced -> normal"
    do i = 1, this%sps%norbs
      oi => this%sps%GetOrbit(i)
      do j = 1, this%sps%norbs
        this%mat%m(i,j) = this%mat%m(i,j) / sqrt(dble(oi%j+1))
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine ReducedToNormal

  subroutine NormalToReduced(this)
    class(OneBodyLabOp), intent(inout) :: this
    integer :: i, j
    type(SingleParticleOrbit), pointer :: oi

    if(this%GetOpJ()/=0 .or. this%GetOpP()/=1 .or. this%GetOpZ()/=0) return
    if(this%reduced_me()) return
    write(*,*) "OneBody: normal -> reduced"
    do i = 1, this%sps%norbs
      oi => this%sps%GetOrbit(i)
      do j = 1, this%sps%norbs
        this%mat%m(i,j) = this%mat%m(i,j) * sqrt(dble(oi%j+1))
      end do
    end do
    call this%SetReduced(.true.)
  end subroutine NormalToReduced

  subroutine set_multipole(op, isospin)
    use MyLibrary, only: triag, hc, m_nucleon, gauss_legendre, ho_radial_wf_norm
    abstract interface
      subroutine func(res, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
        real(8), intent(inout) :: res(:,:)
        real(8), intent(in) :: Q
        integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
        integer, optional :: isospin
      end subroutine func
    end interface

    type(OneBodyLabOp), intent(inout) :: op
    integer, intent(in), optional :: isospin
    type(Orbits), pointer :: sps
    real(8) :: hw, par, phase
    integer :: i, j, n, l
    type(sys) :: s
    type(str) :: opname
    real(8), allocatable :: p_tmp(:), w_tmp(:), mom_mat(:,:)
    procedure(func), pointer :: f_ptr
    integer :: lbra, jbra, zbra, lket, jket, zket
    integer :: lambda_phase
    integer :: a, b
    type(SingleParticleOrbit), pointer :: oa, ob


    hw = op%hw
    sps => op%sps

    par = hc**2 / (m_nucleon * hw)
    call gauss_legendre(kmin, kmax, pbra, wbra, NMesh)
    allocate(pket(NMesh,NMesh), wket(NMesh,NMesh))
    allocate(rnl_bra(NMesh, 0:sps%emax/2, 0:sps%emax))
    allocate(rnl_ket(NMesh, NMesh, 0:sps%emax/2, 0:sps%emax))
    do i = 1, NMesh
      call gauss_legendre(abs(pbra(i)-op%GetQ()/hc), pbra(i)+op%GetQ()/hc, p_tmp, w_tmp, NMesh)
      pket(i,:) = p_tmp(:)
      wket(i,:) = w_tmp(:)
      do n = 0, sps%emax/2
        do l = 0, sps%emax
        if(2*n + l > sps%emax) cycle
          rnl_bra(i,n,l) = ho_radial_wf_norm(n, l, par, pbra(i)) * (-1.d0)**n
          do j = 1, NMesh
            rnl_ket(j,i,n,l) = ho_radial_wf_norm(n, l, par, pket(i,j)) * (-1.d0)**n
          end do
        end do
      end do
    end do

    opname = op%GetOpName()
    lambda_phase = op%GetOpJ()
    f_ptr => null()
    if(s%find(op%GetOpName(), s%str('UTM5_1B')    )) f_ptr => M5_UT_momentum_mat
    if(s%find(op%GetOpName(), s%str('M-1M5_1B')   )) f_ptr => M5_m1_momentum_mat
    if(s%find(op%GetOpName(), s%str('L5_1B')      )) f_ptr => L5_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tel5_1B')    )) f_ptr => Tel5_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tmag5_1B')   )) f_ptr => Tmag5_momentum_mat
    if(s%find(op%GetOpName(), s%str('ISL5_1B')    )) f_ptr => L5_is_momentum_mat
    if(s%find(op%GetOpName(), s%str('ISTel5_1B')  )) f_ptr => Tel5_is_momentum_mat
    if(s%find(op%GetOpName(), s%str('ISTmag5_1B') )) f_ptr => Tmag5_is_momentum_mat
    if(s%find(op%GetOpName(), s%str('M_1B')       )) f_ptr => M_momentum_mat
    if(s%find(op%GetOpName(), s%str('StaticM_1B') )) f_ptr => M_static_momentum_mat
    if(s%find(op%GetOpName(), s%str('M-1M_1B')    )) f_ptr => M_m1_momentum_mat
    if(s%find(op%GetOpName(), s%str('M-2M_1B')    )) f_ptr => M_m2_momentum_mat
    if(s%find(op%GetOpName(), s%str('SOM_1B')     )) f_ptr => M_so_momentum_mat
    if(s%find(op%GetOpName(), s%str('DFM_1B')     )) f_ptr => M_df_momentum_mat
    if(s%find(op%GetOpName(), s%str('L_1B')       )) f_ptr => L_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tel_1B')     )) f_ptr => Tel_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tmag_1B')    )) f_ptr => Tmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvL_1B')   )) f_ptr => ConvL_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvTmag_1B'))) f_ptr => ConvTmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('SpinTmag_1B'))) f_ptr => SpinTmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvTel_1B') )) f_ptr => ConvTel_momentum_mat
    if(s%find(op%GetOpName(), s%str('SpinTel_1B') )) f_ptr => SpinTel_momentum_mat

    if(s%find(op%GetOpName(), s%str('M_1B')    )) lambda_phase = op%GetOpJ()
    if(s%find(op%GetOpName(), s%str('Tel_1B')  )) lambda_phase = op%GetOpJ()+1
    if(s%find(op%GetOpName(), s%str('Tmag_1B') )) lambda_phase = op%GetOpJ()+1
    if(s%find(op%GetOpName(), s%str('L_1B')    )) lambda_phase = op%GetOpJ()+1
    if(s%find(op%GetOpName(), s%str('M5_1B')   )) lambda_phase = op%GetOpJ()
    if(s%find(op%GetOpName(), s%str('Tel5_1B') )) lambda_phase = op%GetOpJ()+1
    if(s%find(op%GetOpName(), s%str('Tmag5_1B'))) lambda_phase = op%GetOpJ()+1
    if(s%find(op%GetOpName(), s%str('L5_1B')   )) lambda_phase = op%GetOpJ()+1

    if(.not. associated(f_ptr)) then
      write(*,*) "Unknown operator name"
    end if

    allocate(mom_mat(NMesh,NMesh))
    do zbra = -1, 1, 2
      do zket = -1, 1, 2
        do lbra = 0, sps%emax
          do lket = 0, sps%emax
            do jbra = abs(2*lbra-1), 2*lbra+1, 2
              do jket = abs(2*lket-1), 2*lket+1, 2
                if(triag(jbra, jket, 2*op%GetOpJ())) cycle
                if((-1)**(lbra+lket) /= op%GetOpP()) cycle
                if(abs(zbra-zket)/2 /= op%GetOpZ()) cycle

                call f_ptr(mom_mat, op%GetQ(), op%GetOpJ(), lbra, jbra, zbra, lket, jket, zket, isospin)
                if(mod(abs(lbra-lket-lambda_phase),2)==0) phase = dble(i_unit ** (lbra - lket - lambda_phase))
                if(mod(abs(lbra-lket-lambda_phase),2)==1) phase = dble(i_unit ** (lbra - lket - lambda_phase - 1)) 
                mom_mat(:,:) = mom_mat(:,:) * phase
                call transform_to_ho(op, mom_mat, lbra, jbra, zbra, lket, jket, zket)

              end do
            end do
          end do
        end do
      end do
    end do
    deallocate(mom_mat)
    deallocate(pbra, wbra, pket, wket, rnl_bra, rnl_ket)
  end subroutine set_multipole

  subroutine transform_to_ho(this, mom_mat, lbra, jbra, zbra, lket, jket, zket)
    use MyLibrary, only: hc
    type(OneBodyLabOp), intent(inout) :: this
    real(8), intent(in) :: mom_mat(:,:)
    integer, intent(in) :: lbra, jbra, zbra, lket, jket, zket
    type(Orbits), pointer :: sps
    integer :: ibra, iket, nbra, nket, idx_bra, idx_ket
    real(8) :: me

    sps => this%sps
    do nbra = 0, sps%emax/2
      do nket = 0, sps%emax/2
        if(2*nbra + lbra > sps%emax) cycle
        if(2*nket + lket > sps%emax) cycle

        me = 0.d0
        do ibra = 1, NMesh
          do iket = 1, NMesh
            me = me + wbra(ibra) * wket(ibra,iket) * pbra(ibra) * pket(ibra,iket) * &
              & rnl_bra(ibra, nbra, lbra) * rnl_ket(iket,ibra,nket,lket) * &
              & mom_mat(ibra, iket) * hc**3
          end do
        end do
        idx_bra = sps%nljz2idx(nbra, lbra, jbra, zbra)
        idx_ket = sps%nljz2idx(nket, lket, jket, zket)
        if(idx_bra==0) write(*,*) nbra, lbra, jbra, zbra
        if(idx_ket==0) write(*,*) nket, lket, jket, zket
        this%mat%m(idx_bra,idx_ket) = me
      end do
    end do
  end subroutine transform_to_ho

  !
  ! mom-space formalism
  !
  subroutine M_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8), allocatable :: tmp1(:,:), tmp2(:,:), tmp3(:,:)
    allocate(tmp1(NMesh,NMesh), tmp2(NMesh,NMesh), tmp3(NMesh, NMesh))
    call M_static_momentum_mat(tmp1, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    call M_so_momentum_mat(tmp2, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    call M_df_momentum_mat(tmp3, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !call M_m1_momentum_mat(tmp2, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !call M_m2_momentum_mat(tmp3, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    mat = tmp1 + tmp2 + tmp3
    deallocate(tmp1, tmp2, tmp3)
  end subroutine M_momentum_mat

  subroutine L_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    call ConvL_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
  end subroutine L_momentum_mat

  subroutine Tel_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8), allocatable :: tmp(:,:)
    allocate(tmp(NMesh,NMesh))
    call ConvTel_momentum_mat(tmp, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    mat(:,:) = tmp(:,:)
    call SpinTel_momentum_mat(tmp, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    mat(:,:) = mat(:,:) + tmp(:,:)
    deallocate(tmp)
  end subroutine Tel_momentum_mat

  subroutine Tmag_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8), allocatable :: tmp(:,:)
    allocate(tmp(NMesh,NMesh))
    call ConvTmag_momentum_mat(tmp, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    mat(:,:) = tmp(:,:)
    call SpinTmag_momentum_mat(tmp, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    mat(:,:) = mat(:,:) + tmp(:,:)
    deallocate(tmp)
  end subroutine Tmag_momentum_mat

  subroutine M5_UT_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    fact = g_A / (8.d0 * pi) / (Q**2 + m_pi**2) * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    fact = make_beta_hermitian(fact, 1, tz_phase*(zbra-zket)/2) 
    mat = func_axial_charge_ut(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M5_UT_momentum_mat

  subroutine M5_m1_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    fact = (-1.d0) * g_A / (8.d0 * pi) * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    fact = make_beta_hermitian(fact, 1, tz_phase*(zbra-zket)/2) 
    mat = func_axial_charge_m1(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M5_m1_momentum_mat

  subroutine L5_is_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! J = sigma/2
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_1
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==1) return
    end if
    fact = 1.d0 / (8.d0 * pi) * tau_1(zbra,zket)
    if(abs(fact) < 1.d-8) return
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * (-1.d0) * sqrt(dble(rankJ+1)/dble(2*rankJ+1))
    if(rankJ>0) mat = mat + func_axial_vector(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ)/dble(2*rankJ+1))
    mat = mat * fact
  end subroutine L5_is_momentum_mat

  subroutine Tel5_is_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! J = sigma/2
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_1
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==1) return
    end if
    fact = 1.d0 / (8.d0 * pi) * tau_1(zbra,zket)
    if(abs(fact) < 1.d-8) return
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * sqrt(dble(rankJ)/dble(2*rankJ+1))
    if(rankJ>0) mat = mat + func_axial_vector(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ+1)/dble(2*rankJ+1))
    mat = mat * fact
  end subroutine Tel5_is_momentum_mat

  subroutine Tmag5_is_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! J = sigma/2
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_1
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==1) return
    end if
    fact = 1.d0 / (8.d0 * pi) * tau_1(zbra,zket)
    if(abs(fact) < 1.d-8) return
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ, rankJ, Q) * fact * (-1.d0)
  end subroutine Tmag5_is_momentum_mat

  subroutine L5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! LO: J = (-gA/2) sigma tau
    ! For DM, you have to multiply (-1/gA)
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    fact = 1.d0 / (4.d0 * pi) * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    fact = fact * ( -g_A + g_A * Q**2 / (Q**2 + m_pi**2) - 2.d0 * d18 * m_pi**2 * Q**2 / (Q**2 + m_pi**2)) * 0.5d0
    fact = make_beta_hermitian(fact, 1, tz_phase*(zbra-zket)/2) 
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * (-1.d0) * sqrt(dble(rankJ+1)/dble(2*rankJ+1))
    if(rankJ>0) mat = mat + func_axial_vector(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ)/dble(2*rankJ+1))
    mat = mat * fact
  end subroutine L5_momentum_mat

  subroutine Tel5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! LO: J = (-gA/2) sigma tau
    ! For DM, you have to multiply (-1/gA)
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    fact = 1.d0 / (4.d0 * pi) * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    fact = fact * (-g_A+d22*Q**2) * 0.5d0
    fact = make_beta_hermitian(fact, 1, tz_phase*(zbra-zket)/2) 
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * sqrt(dble(rankJ)/dble(2*rankJ+1))
    if(rankJ>0) mat = mat + func_axial_vector(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ+1)/dble(2*rankJ+1))
    mat = mat * fact
  end subroutine Tel5_momentum_mat

  subroutine Tmag5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! LO: J = (-gA/2) sigma tau
    ! For DM, you have to multiply (-1/gA)
    !
    use MyLibrary, only: hc, m_pi, g_A, pi, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    fact = 1.d0 / (4.d0 * pi) * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    fact = fact * (-g_A+d22*Q**2) * 0.5d0
    fact = make_beta_hermitian(fact, 1, tz_phase*(zbra-zket)/2) 
    mat = func_axial_vector(lbra, jbra, lket, jket, rankJ, rankJ, Q) * fact * (-1.d0)
  end subroutine Tmag5_momentum_mat

  subroutine M_static_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE, fact, tmp
    integer :: ibra, iket
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    fact = GE / (4.d0 * pi)
    mat = func_charge_static(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M_static_momentum_mat

  subroutine M_m1_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! spin-orbit contribution
    !
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, m_nucleon, tau_1, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GM, fact, ang_fact, tmp
    integer :: ibra, iket, K
    mat(:,:) = 0.d0
    GM = 0.d0
    if(.not. present(isospin)) then
      GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GM = GM + tmp * GM_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GM = tmp * GM_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GM) < 1.d-8) return
    fact = GM / (8.d0 * pi * m_nucleon**2)
    mat = func_charge_m1(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M_m1_momentum_mat

  subroutine M_m2_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, m_nucleon, tau_1, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE, tmp
    real(8) :: fact1, fact2
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    fact1 = (-1.d0) * GE * Q**2 / (32.d0 * pi * m_nucleon**2)
    fact2 = (-1.d0) * GE / (16.d0 * pi * m_nucleon**2)
    fact2 = 0.d0
    mat = func_charge_static(lbra, jbra, lket, jket, rankJ, Q) * fact1 + &
        & func_charge_m1(lbra, jbra, lket, jket, rankJ, Q) * fact2
  end subroutine M_m2_momentum_mat

  subroutine M_so_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    !
    ! spin-orbit contribution
    !
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, m_nucleon, tau_1, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GM, GE, fact, ang_fact, tmp
    integer :: ibra, iket, K
    mat(:,:) = 0.d0
    GE = 0.d0
    GM = 0.d0
    if(.not. present(isospin)) then
      GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GM = GM + tmp * GM_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GM = tmp * GM_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if

    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if

    fact = (GM-0.5d0*GE) / (8.d0 * pi * m_nucleon**2)
    mat = func_charge_m1(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M_so_momentum_mat

  subroutine M_df_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, m_nucleon, tau_1, tau_m, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE
    real(8) :: fact, tmp
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    fact = (-1.d0) * GE * Q**2 / (32.d0 * pi * m_nucleon**2)
    mat = func_charge_static(lbra, jbra, lket, jket, rankJ, Q) * fact
  end subroutine M_df_momentum_mat


  subroutine ConvL_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE, fact, tmp
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1) * 0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    mat = func_vector_conv(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * sqrt(dble(rankJ+1) / dble(2*rankJ+1)) * (-1.d0)
    if(rankJ > 0) mat = mat + func_vector_conv(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ) / dble(2*rankJ+1))
    mat(:,:) = mat(:,:) * GE / (4.d0 * pi * (m_proton+m_neutron))
  end subroutine ConvL_momentum_mat

  subroutine ConvTel_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE, fact, tmp
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1) * 0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    mat = func_vector_conv(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * sqrt(dble(rankJ) / dble(2*rankJ+1))
    if(rankJ > 0) mat = mat + func_vector_conv(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ+1) / dble(2*rankJ+1))
    mat(:,:) = mat(:,:) * GE / (4.d0 * pi * (m_proton+m_neutron))
  end subroutine ConvTel_momentum_mat

  subroutine ConvTmag_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GE, fact, tmp
    mat(:,:) = 0.d0
    GE = 0.d0
    if(.not. present(isospin)) then
      GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GE = GE + tmp * GE_form_factor_iso(Q,1) * 0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GE = tau_1(zbra,zket) * GE_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GE = tmp * GE_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GE) < 1.d-8) return
    mat = func_vector_conv(lbra, jbra, lket, jket, rankJ, rankJ, Q)
    mat(:,:) = mat(:,:) * GE / (4.d0 * pi * (m_proton+m_neutron)) * (-1.d0)
  end subroutine ConvTmag_momentum_mat

  subroutine SpinTel_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GM, fact, tmp
    integer :: k, ibra, iket

    mat(:,:) = 0.d0
    GM = 0.d0
    if(.not. present(isospin)) then
      GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GM = GM + tmp * GM_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GM = tmp * GM_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GM) < 1.d-8) return
    mat = func_vector_spin(lbra, jbra, lket, jket, rankJ+1, rankJ, Q) * sqrt(dble(rankJ) / dble(2*rankJ+1))
    if(rankJ > 0) mat = mat + func_vector_spin(lbra, jbra, lket, jket, rankJ-1, rankJ, Q) * sqrt(dble(rankJ+1) / dble(2*rankJ+1))
    mat(:,:) = mat(:,:) * GM / (4.d0 * pi * (m_proton+m_neutron)) 
  end subroutine SpinTel_momentum_mat

  subroutine SpinTmag_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin)
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv, make_beta_hermitian
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: GM, fact, tmp
    integer :: k, ibra, iket

    mat(:,:) = 0.d0
    GM = 0.d0
    if(.not. present(isospin)) then
      GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0)*0.5d0
      tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
      GM = GM + tmp * GM_form_factor_iso(Q,1)*0.5d0
    end if
    if(present(isospin)) then
      if(isospin==0) then
        GM = tau_1(zbra,zket) * GM_form_factor_iso(Q,0) * 0.5d0
      else if(isospin==1) then
        tmp = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
        tmp = make_beta_hermitian(tmp, 1, tz_phase*(zbra-zket)/2) 
        GM = tmp * GM_form_factor_iso(Q,1) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(GM) < 1.d-8) return
    mat = func_vector_spin(lbra, jbra, lket, jket, rankJ, rankJ, Q)
    mat(:,:) = mat(:,:) * GM / (4.d0 * pi * (m_proton+m_neutron)) * (-1.d0)
  end subroutine SpinTmag_momentum_mat

  function func_charge_static(lbra, jbra, lket, jket, lambda, Q) result(res)
    use MyLibrary, only: hc, sjs
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, lambda
    real(8), allocatable :: res(:,:)
    integer :: ibra, iket

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do ibra = 1, NMesh
      do iket = 1, NMesh
        res(ibra,iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, lambda)
      end do
    end do
    res(:,:) = res(:,:) * (-1.d0)**(lbra+lambda+(1+jket)/2) * &
        & sqrt(dble(jbra+1)*dble(jket+1)) * sjs(jbra, 2*lbra, 1, 2*lket, jket, 2*lambda)
  end function func_charge_static

  function func_charge_m1_(lbra, jbra, lket, jket, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk, fks, flam
    integer :: K, ibra, iket, lam1, lam2, k1, k2

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do K = abs(lbra-lket), lbra+lket
      fk = (-1.d0)**(lambda+K) * sqrt(dble(2*K+1)) * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda)
      if(abs(fk) < 1.d-8) cycle
      do lam1 = 0, lambda
        lam2 = lambda - lam1
        flam = (-1.d0)**lam2 * sqrt(dble(2*lam1+1)*dble(2*lam2+1)) * &
            & sqrt(4.d0*pi*gamma_function(dble(2*lambda+2)) / &
            & (gamma_function(dble(2*lam1+2)) * gamma_function(dble(2*lam2+2))))
        do k1 = abs(lam1-1), lam1+1
          do k2 = max(abs(lam2-1), abs(K-k1)), min(lam2+1, K+k1)
            fks = snj(2*lam1, 2*lam2, 2*lambda, 2, 2, 2, 2*k1, 2*k2, 2*K) * &
                & dcg(2*lam1, 0, 2, 0, 2*k1, 0) * dcg(2*lam2, 0, 2, 0, 2*k2, 0)
            if(abs(fks) < 1.d-8) cycle

            do ibra = 1, NMesh
              do iket = 1, NMesh
                res(ibra,iket) = res(ibra,iket) + (fk * flam * fks) * &
                    & general_integral_(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, k1, k2, K, Q) * &
                    & (pbra(ibra)*hc)**(lam1+1) * (pket(ibra,iket)*hc)**(lam2+1) / Q**lambda
              end do
            end do
          end do
        end do
      end do
    end do
    res(:,:) = res(:,:) * 6.d0 * sqrt(dble(jbra+1)*dble(jket+1)*dble(2*lambda+1))
  end function func_charge_m1_

  function func_charge_m1(lbra, jbra, lket, jket, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk, fks, flam, pb, pk
    integer :: K, ibra, iket, lam1, lam2, k1, k2

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do K = abs(lbra-lket), lbra+lket
      fk = sqrt(dble(2*K+1)) * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda)
      if(abs(fk) < 1.d-8) cycle
      do ibra = 1, NMesh
        pb = pbra(ibra)*hc
        do iket = 1, NMesh
          pk = pket(ibra,iket)*hc
          res(ibra,iket) = res(ibra,iket) + &
            & pb * pk * general_integral(pb, lbra, pk, lket, 1, 1, 1, lambda, K, Q) * fk
        end do
      end do
    end do
    res(:,:) = res(:,:) * (2.d0/sqrt(3.d0)) * sqrt(dble(jbra+1)*dble(jket+1)) * (4.d0 * pi)
  end function func_charge_m1

  function func_vector_spin(lbra, jbra, lket, jket, kappa, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, m_proton, m_neutron
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, kappa, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk
    integer :: K, ibra, iket
    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do K = max(abs(kappa-1), abs(lambda-1)), min(kappa+1, lambda+1)
      fk = sjs(2*kappa, 2, 2*K, 2, 2*lambda, 2) * dcg(2*kappa, 0, 2, 0, 2*K, 0) * &
          & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda)
      if(abs(fk) < 1.d-8) cycle
      do ibra = 1, NMesh
        do iket = 1, NMesh
          res(ibra,iket) = res(ibra,iket) + fk * angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, K) 
        end do
      end do
    end do
    res(:,:) = res(:,:) * 6.d0 * sqrt(dble(jbra+1)*dble(jket+1)*dble(2*lambda+1)*dble(2*kappa+1)) * Q
    res(:,:) = res(:,:) * (-1.d0)**(1+lambda+kappa)
  end function func_vector_spin

  ! it seems this one is not numerically stable for high lambda
  function func_vector_conv_(lbra, jbra, lket, jket, kappa, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, m_proton, m_neutron, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, kappa, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk, fact
    integer :: k1, k2, K, ibra, iket

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    fact = (-1.d0)**(lbra+(1+jket)/2+lambda) * sqrt(dble(jbra+1)*dble(jket+1)*dble(2*kappa+1)) * &
        & sjs(jbra, 2*lbra, 1, 2*lket, jket, 2*lambda)
    do k1 = 0, kappa
      k2 = kappa - k1
      do K = max(abs(1-k1), abs(k2-lambda)), min(1+k1, k2+lambda)
        fk = (-1.d0)**k2 * sqrt(dble(2*k1+1)) * &
            & sqrt(4.d0*pi*gamma_function(dble(2*kappa+2)) / (gamma_function(dble(2*k1+2)) * gamma_function(dble(2*k2+2)))) * &
            & sjs(2, 2*k1, 2*K, 2*k2, 2*lambda, 2*kappa) * dcg(2, 0, 2*k1, 0, 2*K, 0)
        if(abs(fk) < 1.d-8) cycle
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fk * fact * &
                & general_integral_(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, K, k2, lambda, Q) * &
                & (pbra(ibra)*hc)**(k1+1) * (pket(ibra,iket)*hc)**k2 / Q**kappa 
          end do
        end do
      end do
    end do

    fact = (-1.d0)**(lbra+(3+jket)/2) * sqrt(dble(jbra+1)*dble(jket+1)*dble(2*kappa+1)) * &
        & sjs(jbra, 2*lbra, 1, 2*lket, jket, 2*lambda)
    do k1 = 0, kappa
      k2 = kappa - k1
      do K = max(abs(1-k2), abs(k1-lambda)), min(1+k2, k1+lambda)
        fk = (-1.d0)**k1 * sqrt(dble(2*k2+1)) * &
            & sqrt(4.d0*pi*gamma_function(dble(2*kappa+2)) / (gamma_function(dble(2*k1+2)) * gamma_function(dble(2*k2+2)))) * &
            & sjs(2, 2*k2, 2*K, 2*k1, 2*lambda, 2*kappa) * dcg(2, 0, 2*k2, 0, 2*K, 0)
        if(abs(fk) < 1.d-8) cycle
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fk * fact * &
                & general_integral_(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, k1, K, lambda, Q) * &
                & (pbra(ibra)*hc)**k1 * (pket(ibra,iket)*hc)**(k2+1) / Q**kappa
          end do
        end do
      end do
    end do
  end function func_vector_conv_

  function func_vector_conv(lbra, jbra, lket, jket, kappa, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, m_proton, m_neutron, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, kappa, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk, fact, pb, pk
    integer :: k1, k2, K, ibra, iket

    allocate(res(NMesh,NMesh))
    fact = (-1.d0)**(lbra+(1+jket)/2+lambda) * sqrt(dble(jbra+1)*dble(jket+1)) * &
        & sjs(jbra, 2*lbra, 1, 2*lket, jket, 2*lambda) 
    fact = fact * 4.d0 * pi / sqrt(3.d0) * (-1.d0)**(1+kappa+lambda)
    res(:,:) = 0.d0
    do ibra = 1, NMesh
      pb = pbra(ibra)*hc
      do iket = 1, NMesh
        pk = pket(ibra,iket)*hc
        res(ibra,iket) = pb * general_integral(pb, lbra, pk, lket, 1, 0, 1, kappa, lambda, Q) + &
          & pk * general_integral(pb, lbra, pk, lket, 0, 1, 1, kappa, lambda, Q)
      end do
    end do
    res(:,:) = res(:,:) * fact
  end function func_vector_conv

  function func_axial_charge_ut(lbra, jbra, lket, jket, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, m_proton, m_neutron, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: fk, fact
    integer :: k1, k2, K, ibra, iket

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do K = max(abs(lbra-lket), abs(lambda-1)), min(lbra+lket, lambda+1)
      fk = (-1.d0)**K * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda) * &
          & dcg(2*lambda, 0, 2, 0, 2*K, 0)
      if(abs(fk) < 1.d-8) cycle
      do ibra = 1, NMesh
        do iket = 1, NMesh
          res(ibra,iket) = res(ibra,iket) + fk * angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, K)
        end do
      end do
    end do
    res(:,:) = res(:,:) * (-1.d0)**lambda * sqrt(6.d0*dble(jbra+1)*dble(jket+1)*dble(2*lambda+1)) * Q
  end function func_axial_charge_ut

  function func_axial_charge_m1(lbra, jbra, lket, jket, lambda, Q) result(res)
    use MyLibrary, only: hc, dcg, sjs, snj, pi, m_proton, m_neutron, gamma_function
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, lambda
    real(8), allocatable :: res(:,:)
    real(8) :: flam, fk, fx
    integer :: ibra, iket, lam1, lam2, K, X

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do lam1 = 0, lambda
      lam2 = lambda - lam1
      flam = sqrt(dble(2*lam1+1)) * &
          & sqrt(4.d0*pi*gamma_function(dble(2*lambda+2)) / &
          & (gamma_function(dble(2*lam1+2)) * gamma_function(dble(2*lam2+2))))
      do K = max(abs(lbra-lket), abs(lambda-1)), min(lbra+lket, lambda+1)
        fk = (-1.d0)**K * sqrt(dble(2*K+1)) * &
            & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda)
        if(abs(fk) < 1.d-8) cycle
        do X = max(abs(lam1-1), abs(lam2-K)), min(lam1+1, lam2+K)
          fx = sjs(2*lam1, 2*lam2, 2*lambda, 2*K, 2, 2*X) * dcg(2*lam1, 0, 2, 0, 2*X, 0)
          if(abs(fx) < 1.d-8) cycle

          do ibra = 1, NMesh
            do iket = 1, NMesh
              res(ibra,iket) = res(ibra,iket) + (fk * flam * fx) * &
                  & general_integral_(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, X, lam2, K, Q) * &
                  & (pbra(ibra)*hc)**(lam1+1) * (pket(ibra,iket)*hc)**lam2 / Q**lambda
            end do
          end do

        end do
      end do
    end do

    do lam1 = 0, lambda
      lam2 = lambda - lam1
      flam = (-1.d0)**lam2 * sqrt(dble(2*lam2+1)) * &
          & sqrt(4.d0*pi*gamma_function(dble(2*lambda+2)) / &
          & (gamma_function(dble(2*lam1+2)) * gamma_function(dble(2*lam2+2))))
      do K = max(abs(lbra-lket), abs(lambda-1)), min(lbra+lket, lambda+1)
        fk = (-1.d0)**K * sqrt(dble(2*K+1)) * &
            & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*K, 2, 2*lambda)
        if(abs(fk) < 1.d-8) cycle
        do X = max(abs(lam1-K), abs(lam2-1)), min(lam1+K, lam2+1)
          fx = sjs(2*lam1, 2*lam2, 2*lambda, 2, 2*K, 2*X) * dcg(2*lam2, 0, 2, 0, 2*X, 0)
          if(abs(fx) < 1.d-8) cycle

          do ibra = 1, NMesh
            do iket = 1, NMesh
              res(ibra,iket) = res(ibra,iket) + (fk * flam * fx) * &
                  & general_integral_(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, lam1, X, K, Q) * &
                  & (pbra(ibra)*hc)**(lam1+1) * (pket(ibra,iket)*hc)**lam2 / Q**lambda
            end do
          end do

        end do
      end do
    end do
    res(:,:) = res(:,:) * (-1.d0) * 0.5d0 * sqrt(6.d0*dble(jbra+1)*dble(jket+1)*dble(2*lambda+1))
  end function func_axial_charge_m1

  function func_axial_vector(lbra, jbra, lket, jket, kappa, lambda, Q) result(res)
    use MyLibrary, only: hc, snj
    real(8), intent(in) :: Q
    integer, intent(in) :: lbra, jbra, lket, jket, kappa, lambda
    real(8), allocatable :: res(:,:)
    integer :: ibra, iket

    allocate(res(NMesh,NMesh))
    res(:,:) = 0.d0
    do ibra = 1, NMesh
      do iket = 1, NMesh
        res(ibra,iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, kappa)
      end do
    end do
    res(:,:) = res(:,:) * sqrt(6.d0*dble(jbra+1)*dble(jket+1)*dble(2*lambda+1)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*kappa, 2, 2*lambda)
  end function func_axial_vector

  function red_mat_Ys(pbra, lbra, jbra, pket, lket, jket, k, s, rank, Q) result(res)
    use MyLibrary, only: snj
    real(8), intent(in) :: pbra, pket, Q
    integer, intent(in) :: lbra, jbra, lket, jket, k, s, rank
    real(8) :: res
    res = sqrt(dble(6*(jbra+1)*(jket+1)*(2*rank+1))) * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*k, 2*s, 2*rank) * &
        & angular_part(pbra, pket, Q, lbra, lket, k)
  end function red_mat_Ys

  function angular_part(pbra, pket, Q, lbra, lket, L) result(res)
    use MyLibrary, only: pi, triag
    real(8), intent(in) :: pbra, pket, Q
    integer, intent(in) :: lbra, lket, L
    real(8) :: res
    res = 0.d0
    if(triag(lbra, lket, L)) return
    res = sqrt(4.d0 * pi) * (2.d0*pi) * (-1.d0)**(lbra+L) * Yfunc(pbra, lbra, pket, lket, Q, L) / (pbra * pket * Q) 
  end function angular_part

  function general_integral_(pbra, lbra, pket, lket, k1, k2, k, Q) result(res)
    use MyLibrary, only: pi, triag, dcg, sjs
    real(8), intent(in) :: pbra, pket, Q
    integer, intent(in) :: lbra, lket, k1, k2, k
    real(8) :: res
    integer :: lam
    res = 0.d0
    if(triag(lbra,lket,k)) return
    do lam = max(abs(lbra-k1), abs(lket-k2)), min(lbra+k1, lket+k2)
      res = res + dcg(2*lbra, 0, 2*k1, 0, 2*lam, 0) * dcg(2*lket, 0, 2*k2, 0, 2*lam, 0) * &
          & sjs(2*lbra, 2*k1, 2*lam, 2*k2, 2*lket, 2*k) * &
          & Yfunc(pbra, lam, pket, lam, Q, 0) / sqrt(dble(2*lam+1))
    end do
    res = res * (-1.d0)**(lket+k) * &
        & sqrt(dble(2*lbra+1)*dble(2*lket+1)*dble(2*k1+1)*dble(2*k2+1)*dble(2*k+1)) * &
        & 2.d0 * pi / (pbra * pket * Q)
  end function general_integral_

  function general_integral(pbra, lbra, pket, lket, k1, k2, k, lambda, kappa, Q) result(res)
    use MyLibrary, only: pi, triag, dcg, snj
    real(8), intent(in) :: pbra, pket, Q
    integer, intent(in) :: lbra, lket, k1, k2, k, lambda, kappa
    real(8) :: res, fact
    integer :: alpha1, alpha2
    res = 0.d0
    do alpha1 = abs(lbra-k1), lbra+k1
      do alpha2 = abs(lket-k2), lket+k2
        fact = snj(2*k1, 2*alpha1, 2*lbra, 2*k2, 2*alpha2, 2*lket, 2*K, 2*lambda, 2*kappa) * &
          & dcg(2*lbra, 0, 2*k1, 0, 2*alpha1, 0) * dcg(2*lket, 0, 2*k2, 0, 2*alpha2, 0) 
        if(abs(fact) < 1.d-16) cycle
        !res = res + fact * angular_part(pbra, pket, Q, alpha1, alpha2, lambda) * (-1.d0)**(lket+alpha1+kappa) 
        res = res + fact * angular_part(pbra, pket, Q, alpha1, alpha2, lambda) * (-1.d0)**(lbra+alpha1)
      end do
    end do
    res = res * sqrt(dble(2*lbra+1) * dble(2*lket+1) * dble(2*k1+1) * dble(2*k2+1) * dble(2*k+1) * dble(2*kappa+1)) / (4.d0 * pi)
  end function general_integral

  function Yfunc(p1, l1, p2, l2, q, L) result(res)
    use MyLibrary, only: dcg, spherical_harmonics0
    integer, intent(in) :: l1, l2, L
    real(8), intent(in) :: p1, p2, q
    integer :: m
    real(8) :: cost1, cost2
    real(8) :: res

    res = 0.d0
    if(abs(p1-p2)>q) return
    if(    p1+p2 <q) return
    cost1 = (p1**2 - p2**2 + q**2) / (2.d0*p1*q)
    cost2 = (p1**2 - p2**2 - q**2) / (2.d0*p2*q)

    do m = -min(l1, l2), min(l1, l2)
      res = res + dcg(2*l1, -2*m, 2*l2, 2*m, 2*L, 0) * spherical_harmonics0(l1,-m,cost1) * spherical_harmonics0(l2,m,cost2)
    end do
  end function Yfunc

  function GE_form_factor(Q, tz) result(res)
    !
    ! Z. Ye, J. Arrington, R. J. Hill, and G. Lee, Phys. Lett. B 777, 8 (2018).
    !
    real(8), intent(in) :: Q
    integer, intent(in) :: tz
    real(8) :: z, res
    integer :: i
    real(8), parameter :: t0 = -7.0d5, t_cut = 4.d0*(139.57)**2
    real(8), parameter :: ap(13) = [0.239163298067,  -1.10985857441,  1.44438081306,  &
        & 0.479569465603,  -2.28689474187,  1.12663298498,  1.25061984354,  &
        & -3.63102047159,  4.08221702379,  0.504097346499,  -5.08512046051,  &
        & 3.96774254395,  -0.981529071103]
    real(8), parameter :: an(13) = [0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, &
        & 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461, &
        &-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0]

    if(use_parametric_form_factor) then
      z = (sqrt(t_cut+Q**2) - sqrt(t_cut-t0)) / ((sqrt(t_cut+Q**2) + sqrt(t_cut-t0)))
      res = 0.d0
      do i = 1, 13
        if(tz==-1) res = res + ap(i) * z**(i-1)
        if(tz== 1) res = res + an(i) * z**(i-1)
      end do
      return
    end if

    ! LO
    if(tz==-1) res = 1.d0
    if(tz== 1) res = 0.d0
  end function GE_form_factor

  function GE_form_factor_iso(Q, t) result(res)
    real(8), intent(in) :: Q
    integer, intent(in) :: t
    real(8) :: res
    if(t==0) res = GE_form_factor(Q,-1) + GE_form_factor(Q,1)
    if(t==1) res = GE_form_factor(Q,-1) - GE_form_factor(Q,1)
  end function GE_form_factor_iso

  function GM_form_factor(Q, tz) result(res)
    !
    ! Z. Ye, J. Arrington, R. J. Hill, and G. Lee, Phys. Lett. B 777, 8 (2018).
    !
    use MyLibrary, only: gs, gv
    real(8), intent(in) :: Q
    integer, intent(in) :: tz
    real(8) :: z, res
    integer :: i
    real(8), parameter :: t0 = -7.0d5, t_cut = 4.d0*(139.57)**2
    real(8), parameter :: ap(13) = [0.264142994136, -1.09530612212, &
        & 1.21855378178, 0.661136493537, -1.40567892503, -1.35641843888, &
        & 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096, &
        & 8.70740306757, -5.70699994375, 1.28081437589]
    real(8), parameter :: an(13) = [0.257758326959,-1.079540642058, &
        & 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, &
        & 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,&
        & -0.712072389946, 0.0, 0.0]

    if(use_parametric_form_factor) then
      z = (sqrt(t_cut+Q**2) - sqrt(t_cut-t0)) / ((sqrt(t_cut+Q**2) + sqrt(t_cut-t0)))
      res = 0.d0
      do i = 1, 13
        if(tz==-1) res = res + ap(i) * z**(i-1)
        if(tz== 1) res = res + an(i) * z**(i-1)
      end do
      if(tz==-1) res = (gs+gv)*0.5d0*res
      if(tz== 1) res = (gs-gv)*0.5d0*res
      return
    end if

    ! LO
    if(tz==-1) res = (gs + gv) * 0.5d0
    if(tz== 1) res = (gs - gv) * 0.5d0
  end function GM_form_factor

  function GM_form_factor_iso(Q, t) result(res)
    real(8), intent(in) :: Q
    integer, intent(in) :: t
    real(8) :: res
    if(t==0) res = GM_form_factor(Q,-1) + GM_form_factor(Q,1)
    if(t==1) res = GM_form_factor(Q,-1) - GM_form_factor(Q,1)
  end function GM_form_factor_iso

end module OneBodyLabOps
