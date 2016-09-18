module ice_host
  use types
  use list_mod
  use input
  use fabm_driver

  implicit none
  private
  public type_host_model

  type,extends(list):: type_host_model
  contains
    private
    procedure:: initialize
    procedure:: add_var
    procedure:: get_one_var
    procedure:: get_variable
  end type

  interface type_host_model
    module procedure type_host_constructor
  end interface
contains
  function type_host_constructor()
    type(type_host_model):: type_host_constructor

    call type_host_constructor%initialize()
  end function

  subroutine initialize(self)
    class(type_host_model):: self
    type(type_input):: kara_input

    kara_input = type_input('KaraSea.nc')
    !vertical variables
    call self%add_var('depth',kara_input)
    call self%add_var('depth_w',kara_input)
    !horizontal variables
    call self%add_var('ocean_time',kara_input)
    call self%add_var('Pair',kara_input)
    call self%add_var('shflux',kara_input)
    !2d variables
    call self%add_var('temp',kara_input)
    call self%add_var('salt',kara_input)
    call self%add_var('rho',kara_input)
    call self%add_var('AKv',kara_input)
  end subroutine

  subroutine add_var(self,inname,name_input)
    class(type_host_model)   :: self
    character(len=*),intent(in):: inname
    type(type_input),intent(in):: name_input
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_input(inname,var)
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  function get_one_var(self,inname)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    class(*),pointer:: get_one_var

    call self%reset()
    do while(self%moreitems())
      get_one_var => self%get_item()
      select type(get_one_var)
      class is(variable)
        if (trim(get_one_var%name)==trim(inname)) return
      end select
      call self%next()
    end do
    call fatal_error("Getting variables from NetCDF file",&
                     "can't find '"//inname//&
                     "' variable in the NetCDF file")
  end function

  subroutine get_variable(self,inname,get_var)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    class(variable),allocatable,intent(out):: get_var
    class(*),pointer:: curr

    curr => self%get_one_var(inname)
    select type(curr)
    type is(alone_variable)
      allocate(get_var,source=curr)
    type is(variable_1d)
      allocate(get_var,source=curr)
    type is(variable_2d)
      allocate(get_var,source=curr)
    end select
  end subroutine
end module
