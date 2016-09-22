module variables_mod
  use types_mod
  use input_mod
  use fabm_driver

  implicit none
  private
  public brom_standard_variables

  type,extends(list_variables):: brom_standard_variables
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
    procedure:: add_day_number
  end type

  interface brom_standard_variables
    module procedure brom_standard_variables_constructor
  end interface
contains
  function brom_standard_variables_constructor()
    type(brom_standard_variables):: brom_standard_variables_constructor

    call brom_standard_variables_constructor%initialize()
  end function

  subroutine initialize_standard_variables(self)
    class(brom_standard_variables):: self
    type(type_input):: kara_input

    kara_input = type_input('KaraSea.nc')
    !vertical variables
    call self%add_var('depth',kara_input)
    call self%inv_var('depth')
    call self%print_var('depth')
    call self%add_var('depth_w',kara_input)
    call self%inv_var('depth_w')
    !horizontal variables
    call self%add_var('ocean_time',kara_input)
    call self%print_var('ocean_time')
    call self%add_day_number('day_number')
    call self%print_var('day_number')
    !2d variables
    call self%add_var('temp',kara_input)
    call self%inv_var('temp')
    call self%add_var('salt',kara_input)
    call self%inv_var('salt')
    call self%add_var('AKv',kara_input)
    call self%inv_var('AKv')
    call self%print_list()
  end subroutine

  subroutine add_standard_var(self,inname,name_input)
    class(brom_standard_variables):: self
    character(len=*),intent(in):: inname
    type(type_input),intent(in):: name_input
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_input(inname,var)
    !memory allocation problems occur without it
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  subroutine add_day_number(self,inname)
    class(brom_standard_variables),intent(inout):: self
    character(len=*),intent(in):: inname
    class(variable),allocatable:: var

    call self%get_var('ocean_time',var)
    select type(var)
    type is(variable_1d)
      var%value = var%value/86400._rk
      var%name = inname
      call self%add_item(var)
    class default
      call fatal_error("Add day number",&
        "Wrong type")
    end select
  end subroutine
end module
