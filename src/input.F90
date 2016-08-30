module input
  use types
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
    !procedure, public:: get
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

  subroutine initialize(self,infile)
    class(type_input):: self
    character(len=*),intent(in):: infile

    type(alone_variable):: alone_var
    type(variable_1d):: var_1d
    type(variable_2d):: var_2d
    real(rk):: value
    real(rk),allocatable,dimension(:):: value_1d
    real(rk),allocatable,dimension(:,:):: value_2d

    character(len=64):: vname,xname,yname
    integer:: i
    integer:: ncid,xtype,ndims,varid
    integer:: dimids(2)
    integer:: nx,ny

    call check(nf90_open(infile,nf90_nowrite,ncid))
    call check(nf90_inquire_dimension(ncid,1,xname,nx))
    call check(nf90_inquire_dimension(ncid,2,xname,ny))

    allocate(value_1d(nx))
    allocate(value_2d(nx,ny))

    i=1
    do while (check_while(nf90_inquire_variable(&
              ncid,i,vname,xtype,ndims,dimids)))
    !  call check(nf90_inq_varid(ncid, vname, varid))
    !  call check(nf90_get_var(ncid, varid, value))
    !
    !  alone_var = variable()
    !  1d_var = variable()
    !
    !  call self%add_item(vname, var)
      i=i+1
    end do

    deallocate(value_1d)
    deallocate(value_2d)

    call check(nf90_close(ncid))
  end subroutine

  subroutine add_item(self, var)
    class(type_input):: self
    class(variable):: var
    class(item), pointer:: new_item

    if (.not. associated(self%first_item)) then
      self%first_item => item(var, null())
      self%current_item => self%first_item
      self%last_item => self%first_item
    else
      new_item => item(var, null())
      call self%last_item%set_next_item(new_item)
      self%last_item => new_item
    end if
  end subroutine

  !function get(self, inname)
  !  class(type_input):: self
  !  character(len=*), intent(in):: inname
  ! next_item
  ! get_item

  !end function

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
