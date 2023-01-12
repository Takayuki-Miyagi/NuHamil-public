module ClassSys
  implicit none

  public :: sys
  public :: str
  public :: iList
  public :: cList
  public :: dList
  public :: imap
  public :: cmap
  public :: dmap
  public :: assignment(=)
  public :: operator(+)
  private

  type :: sys
  contains
    procedure :: mkdir
    procedure :: sys_i2str
    procedure :: sys_d2str
    procedure :: sys_f2str
    procedure :: sys_char2str
    procedure :: sys_str2i
    procedure :: sys_str2d
    procedure :: sys_str2f
    procedure :: find
    procedure :: split
    procedure :: isfile_stop
    procedure :: isfile_func
    procedure :: sort
    procedure :: sorted_index
    generic :: str => sys_i2str, sys_d2str, sys_f2str, sys_char2str
    generic :: dbl => sys_str2d
    generic :: flt => sys_str2f
    generic :: itg => sys_str2i
    generic :: isfile => isfile_stop, isfile_func
  end type sys

  type :: str
    character(:), allocatable :: val
  end type str

  type :: iList
    integer, allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_ilist
    procedure :: append => append_ilist
    procedure :: prt => print_ilist
  end type iList

  type :: cList
    type(str), allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_clist
    procedure :: append => append_clist
    procedure :: prt => print_clist
  end type cList

  type :: dList
    real(8), allocatable :: i(:)
  contains
    procedure :: GetLen => GetLength_dlist
    procedure :: append => append_dlist
    procedure :: prt => print_dlist
  end type dList

  type :: imap
    integer, allocatable :: val(:)
    type(str), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_imap
    procedure :: append => append_imap
    procedure :: prt => print_imap
    procedure :: Get => Get_imap
    procedure :: Search => Search_imap
  end type imap

  type :: cmap
    type(str), allocatable :: val(:)
    type(str), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_cmap
    procedure :: append => append_cmap
    procedure :: prt => print_cmap
    procedure :: Get => Get_cmap
    procedure :: Search => Search_cmap
  end type cmap

  type :: dmap
    real(8), allocatable :: val(:)
    type(str), allocatable :: key(:)
  contains
    procedure :: GetLen => GetLength_dmap
    procedure :: append => append_dmap
    procedure :: prt => print_dmap
    procedure :: Get => Get_dmap
    procedure :: Search => Search_dmap
  end type dmap

  interface assignment(=)
    module procedure :: char2str
    module procedure :: str2char
    module procedure :: str2str
    module procedure :: i2str
    module procedure :: d2str
    module procedure :: f2str
    module procedure :: str2i
    module procedure :: str2d
    module procedure :: str2f
  end interface assignment(=)

  interface operator(+)
    module procedure :: combine_str_str
    module procedure :: combine_str_char
    module procedure :: combine_char_str
  end interface operator(+)
contains

  subroutine mkdir(this, dir)
    class(sys), intent(in) :: this
    character(*), intent(in) :: dir
    character(512) :: comm
    integer :: is, ie, idxdir
    integer :: lnblnk

    idxdir = lnblnk(dir)
    is = 1; ie = is + len('if [ ! -d ')
    write(comm(is:ie), '(a)') 'if [ ! -d '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1
    ie = is + len(' ]; then (echo "Creating directory, ')
    write(comm(is:ie), '(a)') ' ]; then (echo "Creating directory, '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + len('"; mkdir -p ') - 1
    write(comm(is:ie), '(a)') '"; mkdir -p '
    is = ie + 1; ie = is + idxdir - 1
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + 4
    write(comm(is:ie), '(a)') '); fi'
    !write(*,'(a)') comm(1:ie)
    call system(comm(1:ie))
  end subroutine mkdir

  function sys_char2str(this, cha) result(string)
    class(sys), intent(in) :: this
    character(*), intent(in) :: cha
    type(str) :: string
    string = cha
  end function sys_char2str

  function sys_i2str(this, i) result(string)
    class(sys), intent(in) :: this
    integer, intent(in) :: i
    type(str) :: string
    string = i
  end function sys_i2str

  function sys_d2str(this, d) result(string)
    class(sys), intent(in) :: this
    real(8), intent(in) :: d
    type(str) :: string
    string = d
  end function sys_d2str

  function sys_f2str(this, f) result(string)
    class(sys), intent(in) :: this
    real(4), intent(in) :: f
    type(str) :: string
    string = f
  end function sys_f2str

  function sys_str2i(this, string) result(i)
    class(sys), intent(in) :: this
    type(str), intent(in) :: string
    integer :: i
    i = string
  end function sys_str2i

  function sys_str2d(this, string) result(d)
    class(sys), intent(in) :: this
    type(str), intent(in) :: string
    real(8) :: d
    d = string
  end function sys_str2d

  function sys_str2f(this, string) result(f)
    class(sys), intent(in) :: this
    type(str), intent(in) :: string
    real(4) :: f
    f = string
  end function sys_str2f

  subroutine i2str(string, i)
    integer, intent(in) :: i
    type(str), intent(out) :: string
    character(range(i)+3) :: tmp
    write(tmp,*) i
    string = trim(adjustl(tmp))
  end subroutine i2str

  subroutine str2i(i, string)
    integer, intent(out) :: i
    type(str), intent(in) :: string
    character(:), allocatable :: tmp
    tmp = string
    read(tmp,*) i
  end subroutine str2i

  subroutine str2f(f, string)
    real(4), intent(out) :: f
    type(str), intent(in) :: string
    character(:), allocatable :: tmp
    tmp = string
    read(tmp,*) f
  end subroutine str2f

  subroutine str2d(d, string)
    real(8), intent(out) :: d
    type(str), intent(in) :: string
    character(:), allocatable :: tmp
    tmp = string
    read(tmp,*) d
  end subroutine str2d

  subroutine d2str(string, r)
    real(8), intent(in) :: r
    type(str), intent(out) :: string
    character(100) :: rst
    integer :: l
    l = int(log10(abs(r) + 1.d-8))
    if(l >= 1) then
      write(rst, *) int(r)
    else if(l < 1 .and. l >= 0) then
      write(rst, '(f10.2)') r
    else if(l < 0 .and. l >= -1) then
      write(rst, '(f10.3)') r
    else if(l < -1 .and. l >= -2) then
      write(rst, '(f10.4)') r
    else if(l < -2 .and. l >= -3) then
      write(rst, '(f10.5)') r
    else if(l < -3 .and. l >= -4) then
      write(rst, '(f10.6)') r
    else if(l < -4 .and. l >= -5) then
      write(rst, '(f10.7)') r
    else if(l < -5 .and. l >= -6) then
      write(rst, '(f10.8)') r
    else if(l < -6 .and. l >= -7) then
      write(rst, '(f10.9)') r
    else
      rst = "0"
    end if
    string= trim(adjustl(rst))
  end subroutine d2str

  subroutine f2str(string, r)
    real(4), intent(in) :: r
    type(str), intent(out) :: string
    character(100) :: rst
    integer :: l
    l = int(log10(abs(r) + 1.d-8))
    if(l >= 1) then
      write(rst, *) int(r)
    else if(l < 1 .and. l >= 0) then
      write(rst, '(f10.2)') r
    else if(l < 0 .and. l >= -1) then
      write(rst, '(f10.3)') r
    else if(l < -1 .and. l >= -2) then
      write(rst, '(f10.4)') r
    else if(l < -2 .and. l >= -3) then
      write(rst, '(f10.5)') r
    else if(l < -3 .and. l >= -4) then
      write(rst, '(f10.6)') r
    else if(l < -4 .and. l >= -5) then
      write(rst, '(f10.7)') r
    else if(l < -5 .and. l >= -6) then
      write(rst, '(f10.8)') r
    else if(l < -6 .and. l >= -7) then
      write(rst, '(f10.9)') r
    else
      rst = "0"
    end if
    string= trim(adjustl(rst))
  end subroutine f2str

  subroutine char2str(string, cha)
    character(*), intent(in) :: cha
    type(str), intent(out) :: string
    string%val = trim(cha)
  end subroutine char2str

  subroutine str2char(cha, string)
    type(str), intent(in) :: string
    character(:), allocatable, intent(out) :: cha
    cha = string%val
  end subroutine str2char

  subroutine str2str(str1, str2)
    type(str), intent(in) :: str2
    type(str), intent(out) :: str1
    str1%val = str2%val
  end subroutine str2str

  function combine_str_str(str1, str2) result(str3)
    type(str), intent(in) :: str1, str2
    type(str) :: str3
    str3%val = trim(str1%val) // trim(str2%val)
  end function combine_str_str

  function combine_str_char(str1, cha2) result(str3)
    type(str), intent(in) :: str1
    character(*), intent(in) :: cha2
    type(str) :: str2, str3
    str2 = cha2
    str3 = str1 + str2
  end function combine_str_char

  function combine_char_str(cha1, str2) result(str3)
    character(*), intent(in) :: cha1
    type(str), intent(in) :: str2
    type(str) :: str1, str3
    str1 = cha1
    str3 = str1 + str2
  end function combine_char_str

  logical function find(this, str1, key) result(r)
    class(sys), intent(in) :: this
    type(str), intent(in) :: str1, key
    integer :: i
    i = index(str1%val, key%val)
    if(i == 0) r = .false.
    if(i /= 0) r = .true.
  end function find

  subroutine split(this, string, key, splitted)
    class(sys), intent(in) :: this
    type(str), intent(in) :: string, key
    type(str), allocatable, intent(out) :: splitted(:)
    integer :: len_key, n, i, n_elm
    integer, allocatable :: spos(:), epos(:)
    len_key = len(key%val)
    n_elm = 0
    n = 0
    do i = 1, len(string%val)-len_key+1
      if(string%val(i:i+len_key-1) == key%val) n = n + 1
    end do
    n_elm = n + 1
    allocate(splitted(n_elm))
    allocate(spos(n_elm))
    allocate(epos(n_elm))
    spos(1) = 1
    epos(n_elm) = len(string%val)
    n = 0
    do i = 1, len(string%val)-len_key+1
      if(string%val(i:i+len_key-1) == key%val) then
        n = n + 1
        spos(n+1) = i+len_key
        epos(n) = i-1
      end if
    end do

    do i = 1, n_elm-1
      splitted(i) = trim( adjustl (string%val(spos(i):epos(i)) ))
    end do
    splitted(n_elm) = trim(adjustl(string%val(spos(n_elm):)))
    deallocate(spos, epos)
  end subroutine split

  function isfile_stop(this, f, msg) result(ex)
    class(sys), intent(in) :: this
    type(str), intent(in) :: f, msg
    logical :: ex

    inquire(file=f%val, exist=ex)
    if(.not. ex) then
      write(*,'(3a)') trim(f%val), ' is not found! In ', trim(msg%val)
      stop
    end if

  end function isfile_stop

  function isfile_func(this, f) result(ex)
    class(sys), intent(in) :: this
    type(str), intent(in) :: f
    logical :: ex
    inquire(file=f%val, exist=ex)
  end function isfile_func

  function GetLength_ilist(this) result(n)
    class(iList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_iList

  subroutine append_iList(this, i)
    class(iList), intent(inout) :: this
    integer, intent(in) :: i
    integer, allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_iList

  subroutine print_iList(this, iunit)
    class(iList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,*) this%i
    end if
  end subroutine print_iList

  function GetLength_clist(this) result(n)
    class(cList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_cList

  subroutine append_cList(this, i)
    class(cList), intent(inout) :: this
    type(str), intent(in) :: i
    type(str), allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_cList

  subroutine print_cList(this, iunit)
    class(cList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      do i = 1, n
        write(unt,'(2a)',advance='no') trim(this%i(i)%val), ' '
      end do
      write(unt,*)
    end if
  end subroutine print_cList

  function GetLength_dlist(this) result(n)
    class(dList), intent(in) :: this
    integer :: n
    if(allocated(this%i)) then
      n = size(this%i)
    else
      n = 0
    end if
  end function GetLength_dList

  subroutine append_dList(this, i)
    class(dList), intent(inout) :: this
    real(8), intent(in) :: i
    real(8), allocatable :: list(:)
    integer :: n
    if(allocated(this%i)) then
      n = this%GetLen()
      allocate(list(n+1))
      list(1:n) = this%i(:)
      list(n+1) = i
      deallocate(this%i)
      allocate(this%i(n+1))
      this%i = list
      deallocate(list)
    else
      allocate(this%i(1))
      this%i(1) = i
    end if
  end subroutine append_dList

  subroutine print_dList(this, iunit)
    class(dList), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,*) this%i
    end if
  end subroutine print_dList

  function GetLength_imap(this) result(n)
    class(imap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_imap

  subroutine append_imap(this, key, val)
    class(imap), intent(inout) :: this
    type(str), intent(in) :: key
    integer, intent(in) :: val
    integer :: n, pos
    type(str), allocatable :: keys(:)
    integer, allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key(1) = key
      this%val(1) = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_imap

  function Get_imap(this, key) result(r)
    class(imap), intent(in) :: this
    type(str), intent(in) :: key
    integer :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key%val)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_imap

  subroutine print_imap(this, iunit)
    class(imap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(2a,i5,a)',advance='no') trim(this%key(i)%val), ':', this%val(i), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_imap

  function Search_imap(this, key) result(pos)
    class(imap), intent(in) :: this
    type(str), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_imap

  function GetLength_cmap(this) result(n)
    class(cmap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_cmap

  subroutine append_cmap(this, key, val)
    class(cmap), intent(inout) :: this
    type(str), intent(in) :: key
    type(str), intent(in) :: val
    integer :: n, pos
    type(str), allocatable :: keys(:)
    type(str), allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key(1) = key
      this%val(1) = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_cmap

  subroutine print_cmap(this, iunit)
    class(cmap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(4a)',advance='no') trim(this%key(i)%val), ':', trim(this%val(i)%val), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_cmap

  function Search_cmap(this, key) result(pos)
    class(cmap), intent(in) :: this
    type(str), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_cmap

  function Get_cmap(this, key) result(r)
    class(cmap), intent(in) :: this
    type(str), intent(in) :: key
    type(str) :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key%val)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_cmap

  function GetLength_dmap(this) result(n)
    class(dmap), intent(in) :: this
    integer :: n
    if(allocated(this%key)) then
      n = size(this%key)
    else
      n = 0
    end if
  end function GetLength_dmap

  subroutine append_dmap(this, key, val)
    class(dmap), intent(inout) :: this
    type(str), intent(in) :: key
    real(8), intent(in) :: val
    integer :: n, pos
    type(str), allocatable :: keys(:)
    real(8), allocatable :: vals(:)
    n = this%GetLen()
    if(n == 0) then
      allocate(this%key(1))
      allocate(this%val(1))
      this%key(1) = key
      this%val(1) = val
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        allocate(keys(n+1))
        allocate(vals(n+1))
        keys(1:n) = this%key(:)
        vals(1:n) = this%val(:)
        keys(n+1) = key
        vals(n+1) = val
        deallocate(this%key, this%val)
        allocate(this%key(n+1))
        allocate(this%val(n+1))
        this%key = keys
        this%val = vals
        deallocate(keys)
        deallocate(vals)
      else
        this%val(pos) = val
      end if
    end if
  end subroutine append_dmap

  subroutine print_dmap(this, iunit)
    class(dmap), intent(inout) :: this
    integer, optional, intent(in) :: iunit
    integer :: n, unt = 6, i
    if(present(iunit)) unt = iunit
    n = this%GetLen()
    if(n == 0) then
      write(unt,*) ''
    else
      write(unt,'(a)',advance='no') '{'
      do i = 1, n
        write(unt,'(2a,es12.4,a)',advance='no') trim(this%key(i)%val), ':', this%val(i), ', '
      end do
      write(unt,'(a)') '}'
    end if
  end subroutine print_dmap

  function Get_dmap(this, key) result(r)
    class(dmap), intent(in) :: this
    type(str), intent(in) :: key
    real(8) :: r
    integer :: n, pos
    n = this%GetLen()
    if(n == 0) then
      write(*,*) 'imap is empty'
      stop
    else
      pos = search_key(this%key, key)
      if(pos == 0) then
        write(*,'(2a)') 'key is not found, key:', trim(key%val)
        stop
      else
        r = this%val(pos)
      end if
    end if
  end function Get_dmap

  function Search_dmap(this, key) result(pos)
    class(dmap), intent(in) :: this
    type(str), intent(in) :: key
    logical :: pos
    if(search_key(this%key, key) == 0) pos = .false.
    if(search_key(this%key, key) /= 0) pos = .true.
  end function Search_dmap

  function search_key(keys, key) result(n)
    type(str), intent(in) :: keys(:)
    type(str), intent(in) :: key
    integer :: i, n
    n = 0
    do i = 1, size(keys)
      if(keys(i)%val == key%val) then
        n = i
        return
      end if
    end do
  end function search_key

  function sort(this, a, direction) result(r)
    class(sys), intent(in) :: this
    real(8), intent(in) :: a(:)
    real(8), allocatable :: r(:)
    character(*), optional, intent(in) :: direction
    character(:), allocatable :: direc
    logical, allocatable :: is(:)
    integer :: i, n
    direc = "ascending"
    if(present(direction)) direc = direction
    n = size(a)
    allocate(r(n), is(n))
    is(:) = .true.
    if(direc=="ascending") then
      do i = 1, n
        r(i) = minval(a, dim=1, mask=is)
        is(minloc(a, dim=1, mask=is)) = .false.
      end do
    else if(direc=="descending") then
      do i = 1, n
        r(i) = maxval(a, dim=1, mask=is)
        is(maxloc(a, dim=1, mask=is)) = .false.
      end do
    end if
    deallocate(is)
  end function sort

  function sorted_index(this, a, direction) result(r)
    class(sys), intent(in) :: this
    real(8), intent(in) :: a(:)
    character(*), optional, intent(in) :: direction
    integer, allocatable :: r(:)
    character(:), allocatable :: direc
    logical, allocatable :: is(:)
    integer :: i, n
    direc = "ascending"
    if(present(direction)) direc = direction
    n = size(a)
    allocate(r(n), is(n))
    is(:) = .true.
    if(direc=="ascending") then
      do i = 1, n
        r(i) = minloc(a, dim=1, mask=is)
        is(minloc(a, dim=1, mask=is)) = .false.
      end do
    else if(direc=="decending") then
      do i = 1, n
        r(i) = maxloc(a, dim=1, mask=is)
        is(maxloc(a, dim=1, mask=is)) = .false.
      end do
    end if
    deallocate(is)
  end function sorted_index

end module ClassSys
!
!program test
!  use ClassSys
!  type(sys) :: s
!  type(str) :: a, b, c
!  type(cList) :: cl
!  a = "aaaa"
!  b = "b"
!  c = a + b + s%str(20) + s%str(2.5d0)
!  write(*,*) a%val
!  write(*,*) b%val
!  write(*,*) c%val
!  call cl%append(a%val)
!  call cl%append(b%val)
!  call cl%append(c%val)
!  call cl%prt()
!end program test
