module input
  use item_mod
  use netcdf

  implicit none
  private
  public type_input

  type type_input
    private
    class(item), pointer:: first_item => null()
    class(item), pointer:: current_item => null()
    class(item), pointer:: last_item => null()
  contains
    private
    procedure:: initialize
    procedure:: add_item
    procedure, public:: get
  end type

  interface type_input
    module procedure type_input_constructor
  end interface

contains

  function type_input_constructor(infile)
    character(len=*), intent(in):: infile
    type(type_input):: type_input_constructor

    call type_input_constructor%initialize(infile)

  end function

  subroutine initialize(self, infile)
    class(type_input):: self
    character(len=*), intent(in):: infile

    !class(*):: value
    character(len=64)::vname!, xname, yname
    integer:: i
    integer:: ncid, xtype, ndims, varid
    integer:: dimids(2)
    !integer:: nx, ny

    !call check(nf90_open(infile, nf90_nowrite, ncid))
    !!call check(nf90_inquire_dimension(ncid, 1, xname, nx)
    !!call check(nf90_inquire_dimension(ncid, 2, xname, ny)
    !
    !i=1
    !do while (check_while(nf90_inquire_variable(&
    !          ncid,i,vname,xtype,ndims,dimids)))
    !  call check(nf90_inq_varid(ncid, vname, varid))
    !  call check(nf90_get_var(ncid, varid, value))
    !  
    !  call self%add_item(vname, value)
    !  i=i+1
    !end do
    !
    !call check(nf90_close(ncid))
  end subroutine

  subroutine add_item(self, name, value)
    class(type_input):: self
    character(len=*):: name
    class(*):: value
    class(item), pointer:: new_item

    if (.not. associated(self%first_item)) then
      self%first_item => item(name, value, null())
      self%current_item => self%first_item
      self%last_item => self%first_item
    else
      new_item => item(name, value, null())
      call self%last_item%set_next_item(new_item)
      self%last_item => new_item
    end if
  end subroutine

  function get(self, infile)
    use types
    
    class(type_input):: self
    character(len=*), intent(in):: infile
    class(variable),pointer:: get  
    
  end function

  subroutine check(status)
    integer, intent(in):: status

    if (status .ne. NF90_NOERR) then
      print *, trim(nf90_strerror(status))
      stop
    end if
  end subroutine

  function check_while(status)
    logical check_while
    integer, intent(in):: status

    check_while = .true.
    if (status .ne. NF90_NOERR) then
      check_while = .false.
    end if
  end function

end module
