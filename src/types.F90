module types
  use fabm_types, only: rk

  type:: variable
    character(len=64):: name  = ''
    character(len=64):: units = ''
  end type

  type,extends(variable):: alone_variable
    real(rk):: value = 0._rk
  end type

  type,extends(variable):: variable_1d
    real(rk),allocatable,dimension(:):: value
  end type

  type,extends(variable_1d):: state_variable
    logical:: use_bound_up = .false.
    logical:: use_bound_low = .false.
    real(rk):: sinking_velocity = 0._rk
  end type

  type,extends(variable):: variable_2d
    real(rk),allocatable,dimension(:,:):: value
  end type
  
  type:: netcdf_dimension
    character(len=64):: name = ''
    integer:: dim_id
    integer:: dim_len
  end type
  
end module
