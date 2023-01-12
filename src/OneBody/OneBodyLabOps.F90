module OneBodyLabOps
  use omp_lib
  use ClassSys
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
    logical :: reduced = .true.
  contains
    procedure :: InitOneBodyLabOp
    procedure :: InitOneBodyLabOpFromString
    procedure :: FinOneBodyLabOp
    procedure :: SetOneBodyLabOp
    procedure :: GetFileName
    procedure :: WriteOperator
    !procedure :: ReadOperator
    generic :: init => InitOneBodyLabOp, InitOneBodyLabOpFromString
    generic :: fin => FinOneBodyLabOp
    generic :: set => SetOneBodyLabOp
  end type OneBodyLabOp

  integer, parameter :: tz_phase = -1 ! particle physic convention -> nuclear physics convention
  real(8), allocatable :: rnl_bra(:,:,:), rnl_ket(:,:,:,:), pbra(:), pket(:,:), wbra(:), wket(:,:)
  integer, parameter :: NMesh = 40
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
  end subroutine InitOneBodyLabOp

  subroutine SetOneBodyLabOp(this, isospin, e_charge)
    class(OneBodyLabOp), intent(inout) :: this
    integer, intent(in), optional :: isospin
    real(8), intent(in), optional :: e_charge
    type(Orbits), pointer :: sps
    integer :: a, b
    type(SingleParticleOrbit), pointer :: oa, ob
    type(sys) :: s

    if( s%find(this%GetOpName(), s%str('L5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('M5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('M_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tel5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tmag5_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tel_1B')) .or. &
      & s%find(this%GetOpName(), s%str('Tmag_1B')) ) then
      !call set_1b_multipole(this, isospin=isospin)
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
            & [oa%n,oa%l,oa%j,oa%z], [ob%n,ob%l,ob%j,ob%z], isospin, e_charge)
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
    f = f + '_hw' + s%str(this%hw) + '_emax' + s%str(this%sps%emax) + ".me2j.gz"
    !f = f + '_hw' + s%str(this%hw) + '_emax' + s%str(this%sps%emax) + ".op.snt"
  end function GetFileName

  subroutine WriteOperator(this, f)
    class(OneBodyLabOp), intent(inout) :: this
    type(str), intent(in) :: f
    type(sys) :: s
    real(8) :: ti

    ti = omp_get_wtime()
    if(this%GetOpJ()==0 .and. this%GetOpP()==1 .and. this%GetOpZ()==0 .and. this%reduced_me()) then
      call reduced_to_non_reduced(this)
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
        if(abs(this%mat%m(a,b)) < 1.d-8) cycle
        cnt = cnt + 1
      end do
    end do
    write(wunit,'(i5, i3)') cnt, 0
    write(wunit,'(a)') '### a, b, <a| op |b> or < a|| op || b >'
    do a = 1, sps%norbs
      do b = 1, a
        if(abs(this%mat%m(a,b)) < 1.d-8) cycle
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
      do b_k = 1, a_k
        oa => sps_k%GetOrbit(a_k)
        ob => sps_k%GetOrbit(b_k)
        a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
        b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
        if(abs(this%mat%m(a,b)) < 1.d-8) cycle
        cnt = cnt + 1
      end do
    end do
    write(wunit,'(i5, i3)') cnt, 0
    write(wunit,'(a)') '### a, b, <a| op |b> or < a|| op || b >'
    do a_k = 1, sps_k%norbs
      do b_k = 1, a_k
        oa => sps_k%GetOrbit(a_k)
        ob => sps_k%GetOrbit(b_k)
        a = sps%nljz2idx(oa%n, oa%l, oa%j, oa%z)
        b = sps%nljz2idx(ob%n, ob%l, ob%j, ob%z)
        if(abs(this%mat%m(a,b)) < 1.d-8) cycle
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

  subroutine reduced_to_non_reduced(this)
    type(OneBodyLabOp), intent(inout) :: this
    integer :: i, j
    type(SingleParticleOrbit), pointer :: oi

    if(this%GetOpJ()/=0 .or. this%GetOpP()/=1 .or. this%GetOpZ()/=0) return
    write(*,*) "OneBody: reduced -> non-reduced"
    do i = 1, this%sps%norbs
      oi => this%sps%GetOrbit(i)
      do j = 1, this%sps%norbs
        this%mat%m(i,j) = this%mat%m(i,j) / sqrt(dble(oi%j+1))
      end do
    end do
    call this%SetReduced(.false.)
  end subroutine reduced_to_non_reduced


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
    real(8) :: hw, par
    integer :: i, j, n, l
    type(sys) :: s
    real(8), allocatable :: p_tmp(:), w_tmp(:), mom_mat(:,:)
    procedure(func), pointer :: f_ptr
    integer :: lbra, jbra, zbra, lket, jket, zket

    hw = op%hw
    sps => op%sps

    par = hc**2 / (m_nucleon * hw)
    call gauss_legendre(0.d0, 10.d0, pbra, wbra, NMesh)
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
    
    if(s%find(op%GetOpName(), s%str('M5_1B'))) f_ptr => M5_momentum_mat
    if(s%find(op%GetOpName(), s%str('L5_1B'))) f_ptr => L5_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tel5_1B'))) f_ptr => Tel5_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tmag5_1B'))) f_ptr => Tmag5_momentum_mat
    if(s%find(op%GetOpName(), s%str('M_1B'))) f_ptr => M_momentum_mat
    if(s%find(op%GetOpName(), s%str('L_1B'))) f_ptr => L_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tel_1B'))) f_ptr => Tel_momentum_mat
    if(s%find(op%GetOpName(), s%str('Tmag_1B'))) f_ptr => Tmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvL_1B'))) f_ptr => ConvTmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvTmag_1B'))) f_ptr => ConvTmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('SpinTmag_1B'))) f_ptr => SpinTmag_momentum_mat
    if(s%find(op%GetOpName(), s%str('ConvTel_1B'))) f_ptr => ConvTel_momentum_mat
    if(s%find(op%GetOpName(), s%str('SpinTel_1B'))) f_ptr => SpinTel_momentum_mat
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
  subroutine M5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
  end subroutine M5_momentum_mat

  subroutine L5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, m_pi, g_A, triag, snj, pi, tau_m
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    integer :: ibra, iket
    real(8) :: tfact, fact

    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if

    tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
    if(.not. triag(lbra, lket, rankJ+1)) then
      fact = -sqrt(dble(6*(jbra+1)*(jket+1)*(rankJ+1))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ+1) * fact
        end do
      end do
    end if

    if((.not. triag(lbra, lket, rankJ-1)) .and. rankJ>0) then
      fact = sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = mat(ibra,iket) + angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ-1) * fact
        end do
      end do
    end if
    fact = (-0.5d0) * g_A * (1.d0 - Q**2 / (Q**2 + m_pi**2)) * tfact * (-1.d0)**((lket-lbra)/2)*(-1.d0)**(rankJ/2) / (4.d0 * pi)
    mat(:,:) = mat(:,:) * fact
  end subroutine L5_momentum_mat

  subroutine Tel5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, g_A, triag, snj, pi, tau_m
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    integer :: ibra, iket
    real(8) :: tfact, fact

    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if

    tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
    if(.not. triag(lbra, lket, rankJ+1)) then
      fact = sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ+1) * fact
        end do
      end do
    end if

    if((.not. triag(lbra, lket, rankJ-1)) .and. rankJ>0) then
      fact = sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = mat(ibra,iket) + angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ-1) * fact
        end do
      end do
    end if
    fact = (-0.5d0) * g_A * tfact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0 * pi)
    mat(:,:) = mat(:,:) * fact
  end subroutine Tel5_momentum_mat

  subroutine Tmag5_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, g_A, triag, snj, pi, tau_m
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    integer :: ibra, iket
    real(8) :: tfact, fact
    mat(:,:) = 0.d0
    if(present(isospin)) then
      if(isospin==0) return
    end if
    if(triag(lbra, lket, rankJ)) return
    tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
    fact = sqrt(dble(6*(jbra+1)*(jket+1)*(2*rankJ+1))) * &
      & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*rankJ, 2, 2*rankJ) 
    do ibra = 1, NMesh
      do iket = 1, NMesh
        mat(ibra, iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ) * fact
      end do
    end do
    fact = (-0.5d0) * g_A * tfact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0 * pi)
    mat(:,:) = mat(:,:) * fact
  end subroutine Tmag5_momentum_mat

  subroutine M_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: tfact, fact
    integer :: ibra, iket
    mat(:,:) = 0.d0
    tfact = (tau_1(zbra,zket) + tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)) * 0.5d0
    if(present(isospin)) then
      if(isospin==0) then
        tfact = tau_1(zbra,zket) * 0.5d0
      else if(isospin==1) then
        tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    fact = tfact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0 * pi) 
    fact = fact * sqrt(dble((jbra+1)*(jket+1))) * sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ) * (-1.d0)**((1+jket)/2+lbra+rankJ)
    do ibra = 1, NMesh
      do iket = 1, NMesh
        mat(ibra, iket) = angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, rankJ) * fact
      end do
    end do
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

  subroutine ConvL_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
  end subroutine ConvL_momentum_mat

  subroutine ConvTel_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: tfact, fact
    real(8), allocatable :: tmp1(:,:), tmp2(:,:)
    mat(:,:) = 0.d0
    tfact = tau_1(zbra,zket) + tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(present(isospin)) then
      if(isospin==0) then
        tfact = tau_1(zbra,zket)
      else if(isospin==1) then
        tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    allocate(tmp1(NMesh,NMesh), tmp2(NMesh,NMesh))
    call alpha_m1(tmp1)
    call alpha_p1(tmp2)
    fact = 1.d0 / (sqrt(dble(2*rankJ+1)) * 2.d0 * (m_proton + m_neutron))
    fact = fact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0*pi) * tfact
    mat(:,:) = (tmp1(:,:) + tmp2(:,:)) * fact
    deallocate(tmp1, tmp2)
  contains
    function fact_mom(l1, l2, l) result(res)
      use MyLibrary, only: gamma_function
      integer, intent(in) :: l1, l2, l
      real(8) :: res
      real(8) :: a1, a2, a
      a1 = dble(2*l1+2); a2 = dble(2*l2+2); a  = dble(2*l+2)
      res = sqrt(4*pi*gamma_function(a) / (gamma_function(a1)*gamma_function(a2))) * (-1.d0)**l2
    end function fact_mom

    subroutine alpha_m1(res) 
      real(8), intent(inout) :: res(:,:)
      integer :: lam1, lam2, lam
      real(8), allocatable :: tmp1(:,:), tmp2(:,:)
      res(:,:) = 0.d0
      if(rankJ-1<0) return
      lam = rankJ-1
      allocate(tmp1(NMesh,NMesh), tmp2(NMesh,NMesh))
      do lam1 = 0, lam
        lam2 = lam-lam1
        call beta_bra(tmp1, lam1, lam2, lam)
        call beta_ket(tmp2, lam1, lam2, lam)
        res(:,:) = res(:,:) + fact_mom(lam1, lam2, lam) * (tmp1 + tmp2)
      end do
      res(:,:) = res(:,:) * sqrt(dble( (rankJ+1)*(2*lam+1) ))
      deallocate(tmp1, tmp2)
    end subroutine alpha_m1

    subroutine alpha_p1(res) 
      real(8), intent(inout) :: res(:,:)
      integer :: lam1, lam2, lam
      real(8), allocatable :: tmp1(:,:), tmp2(:,:)
      res(:,:) = 0.d0
      if(rankJ-1<0) return
      lam = rankJ+1
      allocate(tmp1(NMesh,NMesh), tmp2(NMesh,NMesh))
      do lam1 = 0, lam
        lam2 = lam-lam1
        call beta_bra(tmp1, lam1, lam2, lam)
        call beta_ket(tmp2, lam1, lam2, lam)
        res(:,:) = res(:,:) + fact_mom(lam1, lam2, lam) * (tmp1 + tmp2)
      end do
      res(:,:) = res(:,:) * sqrt(dble( rankJ*(2*lam+1) ))
      deallocate(tmp1, tmp2)
    end subroutine alpha_p1

    subroutine beta_bra(res, lam1, lam2, lam) 
      integer, intent(in) :: lam1, lam2, lam
      real(8), intent(inout) :: res(:,:)
      integer :: k, ibra, iket
      real(8) :: fact
      res = 0.d0
      do k = abs(lam1-1), lam1+1, 2
        fact =  sqrt(dble(2*lam1+1)) * &
          & sjs(2, 2*lam1, 2*K, 2*lam2, 2*rankJ, 2*lam) * dcg(2, 0, 2*lam1, 0, 2*k, 0) * &
          & sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ)
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fact * &
              & (pbra(ibra)*hc)**(lam1+1) * (pket(ibra,iket)*hc)**lam2 / Q**lam * &
              & general_integral(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, k, lam2, rankJ, Q)
          end do
        end do
      end do
    end subroutine beta_bra

    subroutine beta_ket(res, lam1, lam2, lam) 
      integer, intent(in) :: lam1, lam2, lam
      real(8), intent(inout) :: res(:,:)
      integer :: k, ibra, iket
      real(8) :: fact
      res = 0.d0
      do k = abs(lam1-1), lam1+1, 2
        fact =  sqrt(dble(2*lam1+1)) * &
          & sjs(2, 2*lam2, 2*K, 2*lam1, 2*rankJ, 2*lam) * dcg(2, 0, 2*lam2, 0, 2*k, 0) * &
          & sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ)
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fact * &
              & (pbra(ibra)*hc)**lam1 * (pket(ibra,iket)*hc)**(lam2+1) / Q**lam * &
              & general_integral(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, lam1, k, rankJ, Q)
          end do
        end do
      end do
    end subroutine beta_ket
  end subroutine ConvTel_momentum_mat

  subroutine ConvTmag_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: tfact, fact
    mat(:,:) = 0.d0
    tfact = tau_1(zbra,zket) + tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(present(isospin)) then
      if(isospin==0) then
        tfact = tau_1(zbra,zket)
      else if(isospin==1) then
        tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    call alpha_0(mat)
    fact = sqrt(dble(2*rankJ+1)) / (2.d0 * (m_proton + m_neutron))
    fact = fact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0*pi) * tfact
    mat(:,:) = mat(:,:) * fact
  contains
    function fact_mom(l1, l2, l) result(res)
      use MyLibrary, only: gamma_function
      integer, intent(in) :: l1, l2, l
      real(8) :: res
      real(8) :: a1, a2, a
      a1 = dble(2*l1+2); a2 = dble(2*l2+2); a  = dble(2*l +2)
      res = sqrt(4*pi*gamma_function(a) / (gamma_function(a1)*gamma_function(a2))) * (-1.d0)**l2
    end function fact_mom

    subroutine alpha_0(res) 
      integer :: lam1, lam2
      real(8), intent(inout) :: res(:,:)
      real(8), allocatable :: tmp1(:,:), tmp2(:,:)
      res = 0.d0
      if(rankJ<0) return
      allocate(tmp1(NMesh,NMesh), tmp2(NMesh,NMesh))
      do lam1 = 0, rankJ
        lam2 = rankJ-lam1
        call beta_bra(tmp1, lam1, lam2)
        call beta_ket(tmp2, lam1, lam2)
        res(:,:) = res(:,:) + fact_mom(lam1, lam2, rankJ) * (tmp1(:,:) + tmp2(:,:))
      end do
      deallocate(tmp1, tmp2)
    end subroutine alpha_0

    subroutine beta_bra(res, lam1, lam2) 
      integer, intent(in) :: lam1, lam2
      real(8), intent(inout) :: res(:,:)
      integer :: k, ibra, iket
      real(8) :: fact
      res(:,:) = 0.d0
      do k = abs(lam1-1), lam1+1, 2
        fact =  sqrt(dble(2*lam1+1)) * &
          & sjs(2, 2*lam1, 2*K, 2*lam2, 2*rankJ, 2*rankJ) * dcg(2, 0, 2*lam1, 0, 2*k, 0) * &
          & sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ)
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fact * &
              & (pbra(ibra)*hc)**(lam1+1) * (pket(ibra,iket)*hc)**lam2 / Q**rankJ * &
              & general_integral(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, k, lam2, rankJ, Q)
          end do
        end do
      end do
    end subroutine beta_bra

    subroutine beta_ket(res, lam1, lam2) 
      integer, intent(in) :: lam1, lam2
      real(8), intent(inout) :: res(:,:)
      integer :: k, ibra, iket
      real(8) :: fact
      res(:,:) = 0.d0
      do k = abs(lam2-1), lam2+1, 2
        fact =  sqrt(dble(2*lam2+1)) * &
          & sjs(2, 2*lam2, 2*K, 2*lam1, 2*rankJ, 2*rankJ) * dcg(2, 0, 2*lam2, 0, 2*k, 0) * &
          & sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ) * (-1.d0)**(lam1+k-rankJ)
        do ibra = 1, NMesh
          do iket = 1, NMesh
            res(ibra,iket) = res(ibra,iket) + fact * &
              & (pbra(ibra)*hc)**lam1 * (pket(ibra,iket)*hc)**(lam2+1) / Q**rankJ * &
              & general_integral(pbra(ibra)*hc, lbra, pket(ibra,iket)*hc, lket, lam1, k, rankJ, Q)
          end do
        end do
      end do
    end subroutine beta_ket
  end subroutine ConvTmag_momentum_mat

  subroutine SpinTel_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: tfact, fact
    integer :: k, ibra, iket

    mat(:,:) = 0.d0
    tfact = (gs * tau_1(zbra,zket) + gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)) * 0.5d0
    if(present(isospin)) then
      if(isospin==0) then
        tfact = gs * tau_1(zbra,zket) * 0.5d0
      else if(isospin==1) then
        tfact = gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return

    if(rankJ-1 >= 0) then
      do k = abs(rankJ-2), rankJ, 2
        fact = sqrt(dble(6 * (rankJ+1) * (2*rankJ-1)) / dble(2*rankJ+1)) * &
          & sjs(2*(rankJ-1), 2, 2*k, 2, 2*rankJ, 2) * dcg(2*(rankJ-1), 0, 2, 0, 2*k, 0) * &
          & sqrt(dble(6*(jbra+1)*(jket+1)*(2*rankJ+1))) * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*k, 2, 2*rankJ) 
        do ibra = 1, NMesh
          do iket = 1, NMesh
            mat(ibra, iket) = mat(ibra, iket) + &
              & angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, k) * fact
          end do
        end do
      end do
    end if

    do k = rankJ, rankJ+2, 2
      fact = sqrt(dble(6 * rankJ * (2*rankJ+31)) / dble(2*rankJ+1)) * &
        & sjs(2*(rankJ+1), 2, 2*k, 2, 2*rankJ, 2) * dcg(2*(rankJ+1), 0, 2, 0, 2*k, 0) * &
        & sqrt(dble(6*(jbra+1)*(jket+1)*(2*rankJ+1))) * snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*k, 2, 2*rankJ) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = mat(ibra, iket) + &
            & angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, k) * fact
        end do
      end do
    end do
    fact = tfact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0*pi) * (-1.d0) / (m_proton+m_neutron) * Q
    mat(:,:) = mat(:,:) * fact
  end subroutine SpinTel_momentum_mat

  subroutine SpinTmag_momentum_mat(mat, Q, rankJ, lbra, jbra, zbra, lket, jket, zket, isospin) 
    use MyLibrary, only: hc, triag, dcg, sjs, snj, pi, tau_1, tau_m, m_proton, m_neutron, gs, gv
    real(8), intent(inout) :: mat(:,:)
    real(8), intent(in) :: Q
    integer, intent(in) :: rankJ, lbra, jbra, zbra, lket, jket, zket
    integer, optional :: isospin
    real(8) :: tfact, fact
    integer :: k, ibra, iket

    mat(:,:) = 0.d0
    tfact = (gs * tau_1(zbra,zket) + gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)) * 0.5d0
    if(present(isospin)) then
      if(isospin==0) then
        tfact = gs * tau_1(zbra,zket) * 0.5d0
      else if(isospin==1) then
        tfact = gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return

    do k = abs(rankJ-1), rankJ+1, 2
      fact = snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*k, 2, 2*rankJ) * &
        & sjs(2*rankJ, 2, 2*k, 2, 2*rankJ, 2) * dcg(2*rankJ, 0, 2, 0, 2*k, 0) 
      do ibra = 1, NMesh
        do iket = 1, NMesh
          mat(ibra, iket) = mat(ibra, iket) + &
            & angular_part(pbra(ibra)*hc, pket(ibra,iket)*hc, Q, lbra, lket, k) * fact
        end do
      end do
    end do
    fact = dble(6 * (2*rankJ+1)) * sqrt(dble((jbra+1)*(jket+1)))
    fact = fact * tfact * (-1.d0)**((lket-lbra)/2) * (-1.d0)**(rankJ/2) / (4.d0*pi) * (-1.d0) / (m_proton+m_neutron) * Q
    mat(:,:) = mat(:,:) * fact
  end subroutine SpinTmag_momentum_mat










  subroutine set_1b_multipole(op, isospin)
    use MyLibrary, only: triag, hc, m_nucleon, gauss_legendre, ho_radial_wf_norm
    type(OneBodyLabOp), intent(inout) :: op
    integer, intent(in), optional :: isospin
    type(Orbits), pointer :: sps
    real(8) :: hw
    integer :: i, j, n, l
    type(SingleParticleOrbit), pointer :: oi, oj
    type(sys) :: s
    real(8) :: me, par
    real(8), allocatable :: p_tmp(:), w_tmp(:)

    hw = op%hw
    sps => op%sps

    par = hc**2 / (m_nucleon * hw)
    call gauss_legendre(0.d0, 10.d0, pbra, wbra, NMesh)
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

    !$omp parallel
    !$omp do private(i,j,oi,oj,me) schedule(dynamic)
    do i = 1, sps%norbs
      do j = 1, i
        oi => sps%GetOrbit(i)
        oj => sps%GetOrbit(j)
        if(triag(oi%j, oj%j, 2*op%GetOpJ())) cycle
        if((-1)**(oi%l+oj%l) /= op%GetOpP()) cycle
        if(abs(oi%z-oj%z)/2 /= op%GetOpZ()) cycle

        if(s%find(op%GetOpName(), s%str('L5_1Br'))) then
          me = me_L5_rspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('L5_1B'))) then
          me = me_L5_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('Tel5_1Br'))) then
          me = me_Tel5_rspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('Tel5_1B'))) then
          me = me_Tel5_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('Tmag5_1Br'))) then
          me = me_Tmag5_rspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('Tmag5_1B'))) then
          me = me_Tmag5_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z)
        else if(s%find(op%GetOpName(), s%str('Tel_1B'))) then
          me = me_Tel_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z, isospin)
        else if(s%find(op%GetOpName(), s%str('SpinTmag_1B'))) then
          me = me_Tmag_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z, &
            & isospin, c_conv=0.d0, c_spin=1.d0)
        else if(s%find(op%GetOpName(), s%str('ConvTmag_1B'))) then
          me = me_Tmag_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z, &
            & isospin, c_conv=1.d0, c_spin=0.d0)
        else if(s%find(op%GetOpName(), s%str('Tmag_1B'))) then
          me = me_Tmag_pspace(op%GetQ(), op%GetOpJ(), hw, oi%n, oi%l, oi%j, oi%z, oj%n, oj%l, oj%j, oj%z, isospin)
        end if
        me = 0.d0
        op%mat%m(i,j) = me
        op%mat%m(j,i) = me * (-1.d0)**((oi%j-oj%j)/2)
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(pbra, wbra, pket, wket, rnl_bra, rnl_ket)
  end subroutine set_1b_multipole

  !
  ! r-space formalism; for details, see P. Klos, J. Menéndez, D. Gazit, and A. Schwenk, Phys. Rev. D 88, 083516 (2013).
  !
  function me_L5_rspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, r, par, fact_iso
    real(8), allocatable :: rmesh(:), rwmesh(:)
    integer :: NMesh, i

    NMesh=300
    call gauss_legendre(0.d0, 10.d0, rmesh, rwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    par = 1.d0 / par
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do i = 1, NMesh
      r = rmesh(i)
      res = res + rwmesh(i) * &
          & ho_radial_wf_norm(nbra, lbra, par, r) * &
          & ho_radial_wf_norm(nket, lket, par, r) * &
          & L5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket)
    end do
    res = res * fact_iso
    deallocate(rmesh, rwmesh)
  end function me_L5_rspace

  function me_Tel5_rspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, r, par, fact_iso
    real(8), allocatable :: rmesh(:), rwmesh(:)
    integer :: NMesh, i

    NMesh=300
    call gauss_legendre(0.d0, 10.d0, rmesh, rwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    par = 1.d0 / par
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do i = 1, NMesh
      r = rmesh(i)
      res = res + rwmesh(i) * &
          & ho_radial_wf_norm(nbra, lbra, par, r) * &
          & ho_radial_wf_norm(nket, lket, par, r) * &
          & Tel5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket)
    end do
    res = res * fact_iso
    deallocate(rmesh, rwmesh)
  end function me_Tel5_rspace

  function me_Tmag5_rspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, r, par, fact_iso
    real(8), allocatable :: rmesh(:), rwmesh(:)
    integer :: NMesh, i

    NMesh=300
    call gauss_legendre(0.d0, 10.d0, rmesh, rwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    par = 1.d0 / par
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do i = 1, NMesh
      r = rmesh(i)
      res = res + rwmesh(i) * &
          & ho_radial_wf_norm(nbra, lbra, par, r) * &
          & ho_radial_wf_norm(nket, lket, par, r) * &
          & Tmag5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket)
    end do
    res = res * fact_iso
    deallocate(rmesh, rwmesh)
  end function me_Tmag5_rspace

  function L5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket) result(res)
    use MyLibrary, only: snj, tjs, spherical_bessel, pi, g_A, m_pi, hc
    real(8), intent(in) :: Q, r
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res, me
    me = sqrt(dble((jbra+1)*(jket+1)*(rankJ+1)*(2*lbra+1)*(2*lket+1)*(2*rankJ+3))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) * &
        & (-1.d0)**lbra * sqrt(6.d0 / (4.d0*pi)) * tjs(2*lbra, 2*(rankJ+1), 2*lket, 0, 0, 0) * &
        & spherical_bessel(rankJ+1, Q*r/hc)
    if(rankJ>0) then
      me = me + sqrt(dble((jbra+1)*(jket+1)*rankJ*(2*lbra+1)*(2*lket+1)*(2*rankJ-1))) * &
          & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) * &
          & (-1.d0)**lbra * sqrt(6.d0 / (4.d0*pi)) * tjs(2*lbra, 2*(rankJ-1), 2*lket, 0, 0, 0) * &
          & spherical_bessel(rankJ-1, Q*r/hc)
    end if
    res = -me * 0.5d0 * g_A * (1.d0 - Q**2 / (Q**2 + m_pi**2))
  end function L5_lo_rspace

  function Tel5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket) result(res)
    use MyLibrary, only: snj, tjs, spherical_bessel, pi, g_A, hc
    real(8), intent(in) :: Q, r
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res, me
    me = -sqrt(dble((jbra+1)*(jket+1)*rankJ*(2*lbra+1)*(2*lket+1)*(2*rankJ+3))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) * &
        & (-1.d0)**lbra * sqrt(6.d0 / (4.d0*pi)) * tjs(2*lbra, 2*(rankJ+1), 2*lket, 0, 0, 0) * &
        & spherical_bessel(rankJ+1, Q*r/hc)
    if(rankJ>0) then
      me = me + sqrt(dble((jbra+1)*(jket+1)*(rankJ+1)*(2*lbra+1)*(2*lket+1)*(2*rankJ-1))) * &
          & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) * &
          & (-1.d0)**lbra * sqrt(6.d0 / (4.d0*pi)) * tjs(2*lbra, 2*(rankJ-1), 2*lket, 0, 0, 0) * &
          & spherical_bessel(rankJ-1, Q*r/hc)
    end if
    res = -me * 0.5d0 * g_A
  end function Tel5_lo_rspace

  function Tmag5_lo_rspace(Q, rankJ, r, lbra, jbra, lket, jket) result(res)
    use MyLibrary, only: snj, tjs, spherical_bessel, pi, g_A, hc, triag
    real(8), intent(in) :: Q, r
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res, me
    res = 0.d0
    me = sqrt(dble((jbra+1)*(jket+1)*(2*lbra+1)*(2*lket+1))) * dble(2*rankJ+1) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*rankJ, 2, 2*rankJ) * &
        & (-1.d0)**lbra * sqrt(6.d0 / (4.d0*pi)) * tjs(2*lbra, 2*rankJ, 2*lket, 0, 0, 0) * &
        & spherical_bessel(rankJ, Q*r/hc)
    res = -me * 0.5d0 * g_A
  end function Tmag5_lo_rspace

  !
  ! mom-space formalism
  !
  function me_L5_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, pbra, pket, par, fact_iso
    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
    integer :: NMesh, ibra, iket

    NMesh=40
    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do ibra = 1, NMesh
      pbra = pmesh(ibra)
      call gauss_legendre(abs(pbra-Q/hc), pbra+Q/hc, pp, ww, NMesh)
      do iket = 1, NMesh
        pket = pp(iket)
        res = res + pwmesh(ibra) * ww(iket) * pbra * pket * &
            & ho_radial_wf_norm(nbra, lbra, par, pbra)*(-1.d0)**nbra * &
            & ho_radial_wf_norm(nket, lket, par, pket)*(-1.d0)**nket * &
            & L5_lo_pspace(Q, rankJ, pbra*hc, lbra, jbra, pket*hc, lket, jket) * hc**3
      end do
    end do
    res = res * fact_iso
    deallocate(pmesh, pwmesh)
  end function me_L5_pspace

  function me_Tel5_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, pbra, pket, par, fact_iso
    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
    integer :: NMesh, ibra, iket

    NMesh=40
    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do ibra = 1, NMesh
      pbra = pmesh(ibra)
      call gauss_legendre(abs(pbra-Q/hc), pbra+Q/hc, pp, ww, NMesh)
      do iket = 1, NMesh
        pket = pp(iket)
        if(abs(pbra-pket)*hc>Q) cycle
        if(   (pbra+pket)*hc<Q) cycle
        res = res + pwmesh(ibra) * ww(iket) * pbra * pket * &
            & ho_radial_wf_norm(nbra, lbra, par, pbra)*(-1.d0)**nbra * &
            & ho_radial_wf_norm(nket, lket, par, pket)*(-1.d0)**nket * &
            & Tel5_lo_pspace(Q, rankJ, pbra*hc, lbra, jbra, pket*hc, lket, jket) * hc**3
      end do
    end do
    res = res * fact_iso
    deallocate(pmesh, pwmesh)
  end function me_Tel5_pspace

  function me_Tmag5_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm, tau_m
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    real(8) :: res, pbra, pket, par, fact_iso
    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
    integer :: NMesh, ibra, iket

    NMesh=40
    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    res = 0.d0
    fact_iso = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(abs(fact_iso)<1.d-8) return
    do ibra = 1, NMesh
      pbra = pmesh(ibra)
      call gauss_legendre(abs(pbra-Q/hc), pbra+Q/hc, pp, ww, NMesh)
      do iket = 1, NMesh
        pket = pp(iket)
        if(abs(pbra-pket)*hc>Q) cycle
        if(   (pbra+pket)*hc<Q) cycle
        res = res + pwmesh(ibra) * ww(iket) * pbra * pket * &
            & ho_radial_wf_norm(nbra, lbra, par, pbra)*(-1.d0)**nbra * &
            & ho_radial_wf_norm(nket, lket, par, pket)*(-1.d0)**nket * &
            & Tmag5_lo_pspace(Q, rankJ, pbra*hc, lbra, jbra, pket*hc, lket, jket) * hc**3
      end do
    end do
    res = res * fact_iso
    deallocate(pmesh, pwmesh)
  end function me_Tmag5_pspace

  function L5_lo_pspace(Q, rankJ, pbra, lbra, jbra, pket, lket, jket) result(res)
    use MyLibrary, only: dcg, snj, pi, g_A, m_pi, hc, triag
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res, me
    res = 0.d0
    me = 0.d0
    if(.not. triag(lbra, lket, rankJ+1)) &
        & me =-sqrt(dble(6*(jbra+1)*(jket+1)*(rankJ+1))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) * &
        & angular_part(pbra, pket, Q, lbra, lket, rankJ+1)
    if((.not. triag(lbra, lket, rankJ-1)) .and. rankJ>0) &
        & me = me + sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) * &
        & angular_part(pbra, pket, Q, lbra, lket, rankJ-1)
    res = -me * 0.5d0 * g_A * (1.d0 - Q**2 / (Q**2 + m_pi**2))
    !res = res * (-1.d0)**(rankJ/2) * (-1.d0)**((lbra-lket)/2) / (4.d0*pi)
    res = res * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi)
  end function L5_lo_pspace

  function Tel5_lo_pspace(Q, rankJ, pbra, lbra, jbra, pket, lket, jket) result(res)
    use MyLibrary, only: dcg, snj, pi, g_A, hc, triag
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res, me

    res = 0.d0
    me = 0.d0
    if(.not. triag(lbra, lket, rankJ+1)) &
        & me = sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) * &
        & angular_part(pbra, pket, Q, lbra, lket, rankJ+1)
    if((.not. triag(lbra, lket, rankJ-1)) .and. rankJ>0) &
        & me = me + sqrt(dble(6*(jbra+1)*(jket+1)*(rankJ+1))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) * &
        & angular_part(pbra, pket, Q, lbra, lket, rankJ-1)
    res = me * (-0.5d0) * g_A
    !res = res * (-1.d0)**(rankJ/2) * (-1.d0)**((lbra-lket)/2) / (4.d0*pi)
    res = res * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi)
  end function Tel5_lo_pspace

  function Tmag5_lo_pspace(Q, rankJ, pbra, lbra, jbra, pket, lket, jket) result(res)
    use MyLibrary, only: dcg, snj, pi, g_A, triag
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
    real(8) :: res

    res = 0.d0
    if(triag(lbra, lket, rankJ)) return
    res = sqrt(dble(6*(jbra+1)*(jket+1)*(2*rankJ+1))) * &
        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*rankJ, 2, 2*rankJ) * &
        & angular_part(pbra, pket, Q, lbra, lket, rankJ) * (-0.5d0) * g_A
    !res = res * (-1.d0)**(rankJ/2) * (-1.d0)**((lbra-lket)/2) / (4.d0*pi)
    res = res * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi)
  end function Tmag5_lo_pspace

!  function me_L_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket) result(res)
!    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm
!    real(8), intent(in) :: Q, hw
!    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
!    real(8) :: res, pbra, pket, par, fact_iso
!    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
!    integer :: NMesh, ibra, iket, cnt
!
!    NMesh=40
!    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
!    par = hc**2 / (m_nucleon * hw)
!    res = 0.d0
!    if(zbra /= zket) return
!    do ibra = 1, NMesh
!      pbra = pmesh(ibra)
!      call gauss_legendre(abs(pbra-Q/hc), pbra+Q/hc, pp, ww, NMesh)
!      do iket = 1, NMesh
!        pket = pp(iket)
!        res = res + pwmesh(ibra) * ww(iket) * pbra * pket * &
!            & ho_radial_wf_norm(nbra, lbra, par, pbra)*(-1.d0)**nbra * &
!            & ho_radial_wf_norm(nket, lket, par, pket)*(-1.d0)**nket * &
!            & L_lo_pspace(Q, rankJ, pbra*hc, lbra, jbra, pket*hc, lket, jket) * hc**3
!      end do
!    end do
!    fact_iso = -1.d0
!    if(zbra==1) fact_iso = 1.d0
!    res = res * fact_iso
!    deallocate(pmesh, pwmesh)
!  end function me_L_pspace

  function me_Tel_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket, isospin) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    integer, intent(in), optional :: isospin
    real(8) :: res, pbra, pket, par
    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
    integer :: NMesh, ibra, iket


    NMesh=40
    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    res = 0.d0
    do ibra = 1, NMesh
      pbra = pmesh(ibra)
      call gauss_legendre(abs(pbra-Q/hc), pbra+Q/hc, pp, ww, NMesh)
      do iket = 1, NMesh
        pket = pp(iket)
        res = res + pwmesh(ibra) * ww(iket) * pbra * pket * &
            & ho_radial_wf_norm(nbra, lbra, par, pbra)*(-1.d0)**nbra * &
            & ho_radial_wf_norm(nket, lket, par, pket)*(-1.d0)**nket * &
            & Tel_lo_pspace(Q, rankJ, pbra*hc, lbra, jbra, zbra, pket*hc, lket, jket, zket, isospin) * hc**3
      end do
    end do
    deallocate(pmesh, pwmesh)
  end function me_Tel_pspace

  function me_Tmag_pspace(Q, rankJ, hw, nbra, lbra, jbra, zbra, nket, lket, jket, zket, isospin, c_conv, c_spin) result(res)
    use MyLibrary, only: gauss_legendre, hc, m_nucleon, ho_radial_wf_norm
    real(8), intent(in) :: Q, hw
    integer, intent(in) :: rankJ, nbra, lbra, jbra, zbra, nket, lket, jket, zket
    integer, intent(in), optional :: isospin
    real(8), intent(in), optional :: c_conv, c_spin
    real(8) :: res, p_bra, p_ket, par
    real(8), allocatable :: pmesh(:), pwmesh(:), pp(:), ww(:)
    integer :: ibra, iket, NMesh

    res = 0.d0
    NMesh=40
    call gauss_legendre(0.d0, 10.d0, pmesh, pwmesh, NMesh)
    par = hc**2 / (m_nucleon * hw)
    do ibra = 1, NMesh
      p_bra = pmesh(ibra)
      call gauss_legendre(abs(p_bra-Q/hc), p_bra+Q/hc, pp, ww, NMesh)
      do iket = 1, NMesh
        p_ket = pp(iket)
        res = res + pwmesh(ibra) * ww(iket) * p_bra * p_ket * &
            & ho_radial_wf_norm(nbra, lbra, par, p_bra)*(-1.d0)**nbra * &
            & ho_radial_wf_norm(nket, lket, par, p_ket)*(-1.d0)**nket * &
            & Tmag_lo_pspace(Q, rankJ, p_bra*hc, lbra, jbra, zbra, p_ket*hc, lket, jket, zket, isospin, &
            & c_conv=c_conv, c_spin=c_spin) * hc**3
      end do
      deallocate(pp,ww)
    end do
    deallocate(pmesh, pwmesh)
  end function me_Tmag_pspace

!  function L_lo_pspace(Q, rankJ, pbra, lbra, jbra, pket, lket, jket) result(res)
!    use MyLibrary, only: dcg, snj, pi, g_A, m_pi, hc, triag
!    real(8), intent(in) :: Q, pbra, pket
!    integer, intent(in) :: rankJ, lbra, jbra, lket, jket
!    integer :: k
!    real(8) :: res, me
!    res = 0.d0
!    me = 0.d0
!    if(.not. triag(lbra, lket, rankJ+1)) &
!        & me =-sqrt(dble(6*(jbra+1)*(jket+1)*(rankJ+1))) * &
!        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ+1), 2, 2*rankJ) * &
!        & angular_part(pbra, pket, Q, lbra, lket, rankJ+1)
!    if((.not. triag(lbra, lket, rankJ-1)) .and. rankJ>0) &
!        & me = me + sqrt(dble(6*(jbra+1)*(jket+1)*rankJ)) * &
!        & snj(2*lbra, 1, jbra, 2*lket, 1, jket, 2*(rankJ-1), 2, 2*rankJ) * &
!        & angular_part(pbra, pket, Q, lbra, lket, rankJ-1)
!    res = -me * 0.5d0 * g_A * (1.d0 - Q**2 / (Q**2 + m_pi**2))
!    res = res * (-1.d0)**(rankJ/2) * (-1.d0)**((lbra-lket)/2) / (4.d0*pi)
!  end function L_lo_pspace

  function Tel_lo_pspace(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) result(res)
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    real(8) :: res
    res = Tel_lo_pspace_spin(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) + &
        & Tel_lo_pspace_conv(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin)
  end function Tel_lo_pspace

  function Tmag_lo_pspace(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin, &
      & c_conv, c_spin) result(res)
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    real(8), intent(in), optional :: c_conv, c_spin
    real(8) :: res, c1, c2
    c1 = 1.d0; c2 = 1.d0
    if(present(c_conv)) c1 = c_conv
    if(present(c_spin)) c2 = c_spin
    res = Tmag_lo_pspace_conv(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin)*c1 + &
      &   Tmag_lo_pspace_spin(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin)*c2 
  end function Tmag_lo_pspace

  function Tel_lo_pspace_spin(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) result(res)
    use MyLibrary, only: dcg, snj, pi, hc, triag, m_proton, m_neutron, sjs, dcg, gv, gs, tau_1, tau_m
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    integer :: k
    real(8) :: res, me, tfact

    res = 0.d0
    me = 0.d0
    tfact = (gs * tau_1(zbra,zket) + gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)) * 0.5d0
    if(present(isospin)) then
      if(isospin==0) then
        tfact = gs * tau_1(zbra,zket) * 0.5d0
      else if(isospin==1) then
        tfact = gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    if(rankJ-1 >= 0) then
      do k = abs(rankJ-2), rankJ, 2
        me = me + sqrt(dble(6 * (rankJ+1) * (2*rankJ-1)) / dble(2*rankJ+1)) * &
            & sjs(2*(rankJ-1), 2, 2*k, 2, 2*rankJ, 2) * dcg(2*(rankJ-1), 0, 2, 0, 2*k, 0) * &
            & red_mat_Ys(pbra, lbra, jbra, pket, lket, jket, k, 1, rankJ, Q)
      end do
    end if

    do k = rankJ, rankJ+2, 2
      me = me + sqrt(dble(6 * rankJ * (2*rankJ+3)) / dble(2*rankJ+1)) * &
          & sjs(2*(rankJ+1), 2, 2*k, 2, 2*rankJ, 2) * dcg(2*(rankJ+1), 0, 2, 0, 2*k, 0) * &
          & red_mat_Ys(pbra, lbra, jbra, pket, lket, jket, k, 1, rankJ, Q)
    end do
    res = me / (m_proton + m_neutron) * Q
    res = res * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi) * tfact * (-1.d0)
  end function Tel_lo_pspace_spin

  function Tmag_lo_pspace_spin(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) result(res)
    use MyLibrary, only: dcg, snj, pi, triag, m_proton, m_neutron, sjs, dcg, gv, gs, tau_1, tau_m
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    integer :: k
    real(8) :: res, me, tfact

    res = 0.d0
    me = 0.d0
    tfact = (gs * tau_1(zbra,zket) + gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)) * 0.5d0
    if(present(isospin)) then
      if(isospin==0) then
        tfact = gs * tau_1(zbra,zket) * 0.5d0
      else if(isospin==1) then
        tfact = gv * tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase) * 0.5d0
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    do k = abs(rankJ-1), rankJ+1, 2
      me = me + sqrt(dble(6 * (2*rankJ+1))) * &
          & sjs(2*rankJ, 2, 2*k, 2, 2*rankJ, 2) * dcg(2*rankJ, 0, 2, 0, 2*k, 0) * &
          & red_mat_Ys(pbra, lbra, jbra, pket, lket, jket, k, 1, rankJ, Q)
    end do
    res = me / (m_proton + m_neutron) * Q
    res = res * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi) * tfact * (-1.d0)
    !if(lbra==0 .and. lket==0 .and. jbra==1 .and. jket==1) write(*,*) res * pbra * pket
  end function Tmag_lo_pspace_spin

  function Tel_lo_pspace_conv(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) result(me)
    use MyLibrary, only: dcg, snj, pi, hc, triag, m_proton, m_neutron, sjs, dcg, tau_1, tau_m
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    real(8) :: me, tfact

    me = 0.d0
    tfact = tau_1(zbra,zket) + tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(present(isospin)) then
      if(isospin==0) then
        tfact = tau_1(zbra,zket)
      else if(isospin==1) then
        tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    me = (alpha_m1() + alpha_p1()) / sqrt(dble(2*rankJ+1))
    me = me / (2.d0 * (m_proton + m_neutron))
    me = me * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi) * tfact
  contains

    function fact_mom(l1, l2, l) result(res)
      use MyLibrary, only: gamma_function
      integer, intent(in) :: l1, l2, l
      real(8) :: res
      real(8) :: a1, a2, a
      a1 = dble(2*l1+2); a2 = dble(2*l2+2); a  = dble(2*l+2)
      res = sqrt(4*pi*gamma_function(a) / (gamma_function(a1)*gamma_function(a2))) * (-1.d0)**l2
    end function fact_mom

    function alpha_m1() result(res)
      integer :: lam1, lam2, lam
      real(8) :: res
      res = 0.d0
      if(rankJ-1<0) return
      lam = rankJ-1
      do lam1 = 0, lam
        lam2 = lam-lam1
        res = res + fact_mom(lam1, lam2, lam) * (beta_bra(lam1, lam2, lam) + beta_ket(lam1, lam2, lam))
      end do
      res = res * sqrt(dble( (rankJ+1)*(2*lam+1) ))
    end function alpha_m1

    function alpha_p1() result(res)
      integer :: lam1, lam2, lam
      real(8) :: res
      res = 0.d0
      lam = rankJ+1
      do lam1 = 0, lam
        lam2 = lam-lam1
        res = res + fact_mom(lam1, lam2, lam) * (beta_bra(lam1, lam2, lam) + beta_ket(lam1, lam2, lam))
      end do
      res = res * sqrt(dble( rankJ*(2*lam+1) ))
    end function alpha_p1

    function beta_bra(lam1, lam2, lam) result(res)
      integer, intent(in) :: lam1, lam2, lam
      real(8) :: res
      integer :: k
      res = 0.d0
      do k = abs(lam1-1), lam1+1, 2
        res = res + sjs(2, 2*lam1, 2*K, 2*lam2, 2*rankJ, 2*lam) * dcg(2, 0, 2*lam1, 0, 2*k, 0) * integral(k,lam2)
      end do
      res = res * sqrt(dble(2*lam1+1)) * pbra**(lam1+1) * pket**lam2 / Q**lam
    end function beta_bra

    function beta_ket(lam1, lam2, lam) result(res)
      integer, intent(in) :: lam1, lam2, lam
      real(8) :: res
      integer :: k
      res = 0.d0
      do k = abs(lam2-1), lam2+1, 2
        res = res + sjs(2, 2*lam2, 2*K, 2*lam1, 2*rankJ, 2*lam) * dcg(2, 0, 2*lam2, 0, 2*k, 0) * integral(lam1,k)
      end do
      res = res * sqrt(dble(2*lam2+1)) * pbra**lam1 * pket**(lam2+1) / Q**lam
    end function beta_ket

    function integral(k1, k2) result(res)
      integer, intent(in) :: k1, k2
      real(8) :: res
      res = sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ) * &
          & general_integral(pbra, lbra, pket, lket, k1, k2, rankJ, Q)
    end function integral
  end function Tel_lo_pspace_conv

  function Tmag_lo_pspace_conv(Q, rankJ, pbra, lbra, jbra, zbra, pket, lket, jket, zket, isospin) result(me)
    use MyLibrary, only: dcg, snj, pi, triag, m_proton, m_neutron, sjs, dcg, tau_1, tau_m
    real(8), intent(in) :: Q, pbra, pket
    integer, intent(in) :: rankJ, lbra, jbra, lket, jket, zbra, zket
    integer, intent(in), optional :: isospin
    real(8) :: me, tfact

    me = 0.d0
    tfact = tau_1(zbra,zket) + tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
    if(present(isospin)) then
      if(isospin==0) then
        tfact = tau_1(zbra,zket)
      else if(isospin==1) then
        tfact = tau_m(zbra,zket,tz_phase*(zbra-zket)/2,phase=tz_phase)
      else
        write(*,*) "Error:", __FILE__, __LINE__
      end if
    end if
    if(abs(tfact) < 1.d-8) return
    me = alpha_0() * sqrt(dble(2*rankJ+1))
    me = me / (2.d0 * (m_proton + m_neutron))
    me = me * (-1.d0)**((lket-lbra+rankJ)/2) / (4.d0*pi) * tfact
  contains
    function fact_mom(l1, l2, l) result(res)
      use MyLibrary, only: gamma_function
      integer, intent(in) :: l1, l2, l
      real(8) :: res
      real(8) :: a1, a2, a
      a1 = dble(2*l1+2); a2 = dble(2*l2+2); a  = dble(2*l +2)
      res = sqrt(4*pi*gamma_function(a) / (gamma_function(a1)*gamma_function(a2))) * (-1.d0)**l2
    end function fact_mom

    function alpha_0() result(res)
      integer :: lam1, lam2
      real(8) :: res
      res = 0.d0
      if(rankJ<0) return
      do lam1 = 0, rankJ
        lam2 = rankJ-lam1
        res = res + fact_mom(lam1, lam2, rankJ) * &
            (beta_bra(lam1, lam2) + beta_ket(lam1, lam2))
      end do
    end function alpha_0

    function beta_bra(lam1, lam2) result(res)
      integer, intent(in) :: lam1, lam2
      real(8) :: res
      integer :: k
      res = 0.d0
      do k = abs(lam1-1), lam1+1, 2
        res = res + sjs(2, 2*lam1, 2*K, 2*lam2, 2*rankJ, 2*rankJ) * dcg(2, 0, 2*lam1, 0, 2*k, 0) * &
            & integral(k,lam2)
      end do
      res = res * sqrt(dble(2*lam1+1)) * pbra**(lam1+1) * pket**lam2 / Q**rankJ
    end function beta_bra

    function beta_ket(lam1, lam2) result(res)
      integer, intent(in) :: lam1, lam2
      real(8) :: res
      integer :: k
      res = 0.d0
      do k = abs(lam2-1), lam2+1, 2
        res = res + sjs(2, 2*lam2, 2*K, 2*lam1, 2*rankJ, 2*rankJ) * dcg(2, 0, 2*lam2, 0, 2*k, 0) * &
            & integral(lam1,k) * (-1.d0)**(lam1+K-rankJ)
      end do
      res = res * sqrt(dble(2*lam2+1)) * pbra**lam1 * pket**(lam2+1) / Q**rankJ
    end function beta_ket

    function integral(k1, k2) result(res)
      integer, intent(in) :: k1, k2
      real(8) :: res
      res = sqrt(dble(jbra+1)*dble(jket+1)) * (-1.d0)**((1+jket)/2+lbra+rankJ) * &
          & sjs(2*lbra, jbra, 1, jket, 2*lket, 2*rankJ) * &
          & general_integral(pbra, lbra, pket, lket, k1, k2, rankJ, Q)
    end function integral
  end function Tmag_lo_pspace_conv

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
    if(triag(lbra,lket,L)) return
    res = sqrt(4.d0 * pi) * (2.d0*pi) * (-1.d0)**(lbra+L) * Yfunc(pbra, lbra, pket, lket, Q, L) / &
        & (pbra * pket * Q)
  end function angular_part

  function general_integral(pbra, lbra, pket, lket, k1, k2, k, Q) result(res)
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
        & sqrt(dble( (2*lbra+1)*(2*lket+1)*(2*k1+1)*(2*k2+1)*(2*k+1) )) * &
        & 2.d0 * pi / (pbra * pket * Q)
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
      res = res + dcg(2*l1, -2*m, 2*l2, 2*m, 2*L, 0) * &
          & spherical_harmonics0(l1,-m,cost1) * spherical_harmonics0(l2,m,cost2)
    end do
  end function Yfunc
end module OneBodyLabOps
