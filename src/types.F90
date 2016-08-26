module types

  use fabm_types, only: attribute_length, rk

  implicit none
  private
  public model

  type,abstract:: variable
    character(len=attribute_length):: name  = ''
    character(len=attribute_length):: units = ''
    real(rk),dimension(:):: value = 0._rk
  end type

  type,extends(variable):: horizontal_variable
  end type

  type,extends(variable):: vertical_variable
  end type

  type,extends(variable):: standart_variable
  end type

  type,extends(variable):: state_variable
    logical:: use_bound_up = .false.
    logical:: use_bound_low = .false.
    real(rk):: sinking_velocity = 0._rk
  end type

end module
