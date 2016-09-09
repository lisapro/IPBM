module input
  use types
  use list_mod
  use netcdf
  use fabm_driver

  implicit none
  private
  public type_input

  type,extends(list):: type_input
  contains
    private
    procedure:: initialize
    procedure:: add_input
    procedure:: get_one_item
    procedure:: get_input_alone_variable
    procedure:: get_input_variable_1d
    procedure:: get_input_variable_2d
    generic,public:: get_input => get_input_alone_variable,&
      get_input_variable_1d, get_input_variable_2d
  end type

  type,extends(list):: type_netcdf_dimension
  contains
    private
    procedure:: add_netcdf_dimension
    procedure:: get_netcdf_dimension
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
    real(rk):: value
    real(rk),dimension(:),allocatable:: value_1d
    real(rk),dimension(:,:),allocatable:: value_2d
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
    integer:: varid
    integer,dimension(2):: len_dim_2d  !Returned length of dimensions
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
      alone_dim = netcdf_dimension(trim(name),i,len_dim)
      call list_dim%add_netcdf_dimension(alone_dim)
    end do

    !Variable ID
    do i=1,nvars
      call check(nf90_inquire_variable(ncid,i,vname,xtype,ndims,dimids))
      if (ndims==0) then
        call check(nf90_get_var(ncid,i,value))
        alone_var = alone_variable(vname,'',value)
        call self%add_input(alone_var)
      else if (ndims==1) then
        len_dim_2d = list_dim%get_netcdf_dimension(dimids(1))
        allocate(value_1d(len_dim_2d(1)))
        call check(nf90_get_var(ncid,i,value_1d))
        var_1d = variable_1d(vname,'',value_1d)
        call self%add_input(var_1d)
        deallocate(value_1d)
      else if (ndims==2) then
        len_dim_2d = list_dim%get_netcdf_dimension(dimids(1),dimids(2))
        allocate(value_2d(len_dim_2d(1),len_dim_2d(2)))
        call check(nf90_get_var(ncid,i,value_2d))
        var_2d = variable_2d(vname,'',value_2d)
        call self%add_input(var_2d)
        deallocate(value_2d)
      end if
    end do

    call check(nf90_close(ncid))
  end subroutine

  subroutine add_input(self, var)
    class(type_input)   :: self
    class(variable)     :: var
    class(*),allocatable:: temp

    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  function get_one_item(self,inname)
    class(type_input):: self
    character(len=*),intent(in):: inname
    class(*),pointer:: get_one_item

    call self%reset()
    do while(self%moreitems())
      get_one_item => self%get_item()
      select type(get_one_item)
      class is(variable)
        if (trim(get_one_item%name)==trim(inname)) return
      end select
      call self%next()
    end do
    call fatal_error("Getting variables from NetCDF file",&
                     "can't find '"//inname//&
                     "' variable in the NetCDF file")
  end function

  subroutine get_input_alone_variable(self,inname,get_alone_variable)
    class(type_input):: self
    character(len=*),intent(in):: inname
    type(alone_variable),allocatable,intent(out):: get_alone_variable
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(alone_variable)
      allocate(get_alone_variable,source=curr)
      return
    end select
  end subroutine

  subroutine get_input_variable_1d(self,inname,get_variable_1d)
    class(type_input):: self
    character(len=*),intent(in):: inname
    type(variable_1d),allocatable,intent(out):: get_variable_1d
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(variable_1d)
      allocate(get_variable_1d,source=curr)
      return
    end select
  end subroutine

  subroutine get_input_variable_2d(self,inname,get_variable_2d)
    class(type_input):: self
    character(len=*),intent(in):: inname
    type(variable_2d),allocatable,intent(out):: get_variable_2d
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(variable_2d)
      allocate(get_variable_2d,source=curr)
      return
    end select
  end subroutine

  subroutine add_netcdf_dimension(self, var)
    class(type_netcdf_dimension):: self
    class(netcdf_dimension):: var
    class(*),allocatable:: temp

    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  function get_netcdf_dimension(self, indim_id_1, indim_id_2)
    class(type_netcdf_dimension):: self
    integer,intent(in)          :: indim_id_1
    integer,optional,intent(in) :: indim_id_2
    integer,dimension(2)        :: get_netcdf_dimension
    class(*),pointer            :: curr

    call self%reset()
    do while(self%moreitems())
      curr => self%get_item()
      select type(curr)
      type is(netcdf_dimension)
        if (curr%dim_id==indim_id_1) then
          get_netcdf_dimension(1) = curr%dim_len
          if (.not.present(indim_id_2)) return
        end if
        if (present(indim_id_2) .and. curr%dim_id==indim_id_2) then
          get_netcdf_dimension(2) = curr%dim_len
        end if
      end select
      call self%next()
    end do
  end function

  subroutine check(status)
    integer, intent(in):: status

    if (status .ne. NF90_NOERR) then
      call fatal_error("Netcdf internal",&
                        nf90_strerror(status))
    end if
  end subroutine
end module
