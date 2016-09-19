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
    procedure:: get_var
    procedure:: set_var
    procedure:: inv_var
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
    call self%inv_var('depth')
    call self%add_var('depth_w',kara_input)
    call self%inv_var('depth_w')
    !horizontal variables
    call self%add_var('ocean_time',kara_input)
    call self%add_var('Pair',kara_input)
    call self%add_var('shflux',kara_input)
    !2d variables
    call self%add_var('temp',kara_input)
    call self%inv_var('temp')
    call self%add_var('salt',kara_input)
    call self%inv_var('salt')
    call self%add_var('rho',kara_input)
    call self%inv_var('rho')
    call self%add_var('AKv',kara_input)
    call self%inv_var('AKv')
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

  subroutine get_var(self,inname,get_variable)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    class(variable),allocatable,intent(out):: get_variable
    class(*),pointer:: curr

    call self%reset()
    do
      curr=>self%get_item()
      select type(curr)
      class is(variable)
        if (trim(curr%name)==trim(inname)) then
          allocate(get_variable,source=curr)
          return
        end if
      end select
      call self%next()
      if (.not.self%moreitems()) then
        call fatal_error("Getting variables",&
                         "can't find '"//inname//&
                         " variable")
      end if
    end do
  end subroutine

  subroutine set_var(self,inname,new_var)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    class(variable),allocatable:: new_var
    class(*),pointer:: curr

    call self%reset()
    do
      curr=>self%get_item()
      select type(curr)
      class is(variable)
        if (trim(curr%name)==trim(inname)) then
          call self%set_item(new_var)
          return
        end if
      end select
      call self%next()
      if (.not.self%moreitems()) then
        call fatal_error("Setting variables",&
                         "can't find '"//inname//&
                         "' variable")
      end if
    end do
  end subroutine

  subroutine inv_var(self,inname)
    class(type_host_model):: self
    character(len=*),intent(in):: inname
    class(variable),allocatable:: var

    call self%get_var(inname,var)
    call var%print_value()
    call var%inverse()
    call self%set_var(inname,var)
    call self%get_var(inname,var)
    call var%print_value()
  end subroutine
end module
