module input
  use types
  use list_mod
  use netcdf

  implicit none
  private
  public type_input

  type,extends(list):: type_input
  contains
    private
    procedure:: initialize
    procedure:: add_input
    !procedure, public:: get
  end type

  type,extends(list):: type_netcdf_dimension
  contains
    private
    procedure:: add_netcdf_dimension
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

    type(netcdf_dimension):: alone_dim
    type(type_netcdf_dimension):: list_dim
    
    type(alone_variable):: alone_var
    type(variable_1d):: var_1d
    type(variable_2d):: var_2d

    integer:: ncid        !NetCDF ID
    integer:: ndims       !number of dimensions
    integer:: nvars       !number of variables
    integer:: nglobalatts !number of global attributes
    integer:: unlimdimid  !ID of the unlimited dimension
    integer:: i           !Dimension ID, Variable ID
    integer:: len_dim     !Returned length of dimension
    integer:: xtype       !Returned variable type
    !dimids - Returned vector of *ndimsp dimension 
    !IDs corresponding to the variable dimensions
    integer,dimension(NF90_MAX_VAR_DIMS):: dimids
    character(len=NF90_MAX_NAME):: name  !Returned dimension name
    character(len=NF90_MAX_NAME):: vname !Returned variable name

    call check(nf90_open(infile,nf90_nowrite,ncid))
    call check(nf90_inquire(ncid,ndims,nvars,nglobalatts,unlimdimid))
    
    !Dimension ID
    do i=1,ndims
      call check(nf90_inquire_dimension(ncid,i,name,len_dim))
      alone_dim = netcdf_dimension(name,ncid,len_dim)
      call list_dim%add_netcdf_dimension(alone_dim)
    end do

    !Variable ID
    do i=1,nvars 
      call check(nf90_inquire_variable(ncid,i,vname,xtype,ndims,dimids))
    !  call check(nf90_inq_varid(ncid, vname, varid))
    !  call check(nf90_get_var(ncid, varid, value))
    !
    !  alone_var = variable()
    !  1d_var = variable()
    !
    !  call self%add_item(vname, var)
    end do

    call check(nf90_close(ncid))
  end subroutine

  subroutine add_input(self, var)
    class(type_input):: self
    class(variable):: var
    class(*),allocatable:: temp

    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  !function get(self, inname)
  !  class(type_input):: self
  !  character(len=*), intent(in):: inname
  ! next_item
  ! get_item

  !end function
  
  subroutine add_netcdf_dimension(self, var)
    class(type_netcdf_dimension):: self
    class(netcdf_dimension):: var
    class(*),allocatable:: temp

    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  subroutine check(status)
    integer, intent(in):: status

    if (status .ne. NF90_NOERR) then
      print *, trim(nf90_strerror(status))
      stop
    end if
  end subroutine

  !function check_while(status)
  !  logical check_while
  !  integer, intent(in):: status

  !  check_while = .true.
  !  if (status .ne. NF90_NOERR) then
  !    check_while = .false.
  !  end if
  !end function
end module
