#include "../include/brom.h"
#include "../include/parameters.h"

module variables_mod
  use types_mod
  use input_mod
  use fabm_driver

  implicit none

  type,extends(list_variables):: brom_standard_variables
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
    procedure:: add_day_number
    procedure,public:: first_day
  end type

  type,extends(variable_1d):: brom_state_variable
    real(rk) bound_up
    real(rk) bound_low
    logical:: use_bound_up = .false.
    logical:: use_bound_low = .false.
    real(rk):: sinking_velocity = 0._rk
  contains
    procedure:: set_brom_state_variable
    procedure:: print_state_variable
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
    class(brom_standard_variables),intent(inout):: self
    type(type_input):: kara_input
    logical first

    kara_input = type_input(_FILE_NAME_)
    !vertical variables
    call self%add_var(_MIDDLE_LAYER_DEPTH_,kara_input)
    call self%inv_var(_MIDDLE_LAYER_DEPTH_)
    call self%add_var(_DEPTH_ON_BOUNDARY_,kara_input)
    call self%inv_var(_DEPTH_ON_BOUNDARY_)
    !horizontal variables
    call self%add_var(_OCEAN_TIME_,kara_input)
    call self%add_day_number('day_number')
    !2d variables
    call self%add_var(_TEMPERATURE_,kara_input)
    !call self%print_var(_TEMPERATURE_)
    call self%inv_var(_TEMPERATURE_)
    !call self%print_var(_TEMPERATURE_)
    call self%add_var(_SALINITY_,kara_input)
    call self%inv_var(_SALINITY_)
    call self%add_var(_TURBULENCE_,kara_input)
    call self%inv_var(_TURBULENCE_)

    !call self%print_list('Allocated brom_standard_variables:')
  end subroutine

  subroutine add_standard_var(self,inname,name_input)
    class(brom_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    type(type_input)              ,intent(in)   :: name_input
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_var(inname,var)
    !memory allocation problems occur without it
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

  subroutine add_day_number(self,inname)
    class(brom_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var

    call self%get_var(_OCEAN_TIME_,var)
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

  integer function first_day(self)
    class(brom_standard_variables),intent(in):: self
    class(variable),allocatable:: var

    call self%get_var("day_number",var)
    select type(var)
    type is(variable_1d)
      first_day = int(var%value(1))
    end select
  end function

  subroutine set_brom_state_variable(self,use_bound_up,&
      use_bound_low,bound_up,bound_low,sinking_velocity)
    class(brom_state_variable),intent(inout):: self
    logical,optional          ,intent(in)   :: use_bound_up
    logical,optional          ,intent(in)   :: use_bound_low
    real(rk),optional         ,intent(in)   :: bound_up
    real(rk),optional         ,intent(in)   :: bound_low
    real(rk),optional         ,intent(in)   :: sinking_velocity

    if(present(use_bound_up)) self%use_bound_up = use_bound_up
    if(present(use_bound_low)) self%use_bound_low = use_bound_low
    if(present(bound_up)) self%bound_up = bound_up
    if(present(bound_low)) self%bound_low = bound_low
    if(present(sinking_velocity)) self%sinking_velocity = sinking_velocity
  end subroutine

  subroutine print_state_variable(self)
    class(brom_state_variable),intent(in):: self

    _LINE_
    write(*,*) self%name
    write(*,*) self%use_bound_up,self%bound_up
    write(*,*) self%use_bound_low,self%bound_up
    write(*,*) self%sinking_velocity
    _LINE_
  end subroutine
end module
