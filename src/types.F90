module types
  use fabm_types, only: attribute_length, rk

  implicit none

  type:: variable
    character(len=attribute_length):: name  = ''
    character(len=attribute_length):: units = ''
  end type

  type,extends(variable):: variable_1d
    real(rk),allocatable,dimension(:):: value
  end type

  type,extends(variable):: alone_variable
    real(rk):: value = 0._rk
  end type

  type,extends(variable_1d):: horizontal_variable
  end type

  type,extends(variable_1d):: vertical_variable
  end type

  type,extends(variable_1d):: standard_variable
  end type

  type,extends(variable_1d):: state_variable
    logical:: use_bound_up = .false.
    logical:: use_bound_low = .false.
    real(rk):: sinking_velocity = 0._rk
  end type
end module
