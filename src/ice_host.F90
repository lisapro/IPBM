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
    procedure:: get_alone_variable
    procedure:: get_variable_1d
    procedure:: get_variable_2d
    generic:: get_var => get_alone_variable,&
      get_variable_1d, get_variable_2d
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

    kara_input = type_input('KaraSea.nc')

    !vertical variables
    call self%add_var('depth')
    call self%add_var('depth_w')
    !horizontal variables
    call self%add_var('ocean_time')
    call self%add_var('Pair')
    call self%add_var('shflux')
    !2d variables
    call self%add_var('temp')
    call self%add_var('salt')
    call self%add_var('rho')
    call self%add_var('AKv')
  end subroutine

  subroutine add_var(self,inname)
    class(type_host_model)   :: self
    character(len=*),intent(in):: inname
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call kara_input%get_input(inname,var)
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  function get_one_var(self,inname)
    class(type_host_model):: self
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

  subroutine get_alone_variable(self,inname,get_alone_var)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    type(alone_variable),allocatable,intent(out):: get_alone_var
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(alone_variable)
      allocate(get_alone_var,source=curr)
      return
    end select
  end subroutine

  subroutine get_variable_1d(self,inname,get_var_1d)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    type(variable_1d),allocatable,intent(out):: get_var_1d
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(variable_1d)
      allocate(get_var_1d,source=curr)
      return
    end select
  end subroutine

  subroutine get_variable_2d(self,inname,get_var_2d)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    type(variable_2d),allocatable,intent(out):: get_var_2d
    class(*),pointer:: curr

    curr => self%get_one_item(inname)
    select type(curr)
    type is(variable_2d)
      allocate(get_var_2d,source=curr)
      return
    end select
  end subroutine
end module
