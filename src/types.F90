module types
  use fabm_types, only: rk

  type,abstract:: variable
    character(len=64):: name  = ''
    character(len=64):: units = ''
  contains
    procedure,non_overridable:: inverse
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
contains
  subroutine inverse(self)
    class(variable),intent(inout):: self
    integer:: temp,temp2,i

    select type(self)
    type is(alone_variable)
      return
    type is(variable_1d)
      self%value = self%value(size(self%value):1:-1)
    type is(variable_2d)
      temp = size(self%value,2)
      temp2 = size(self%value,1)
      do i=1,temp
        self%value(:,i) = &
        self%value(temp2:1:-1,i)
      end do
    end select
  end subroutine
end module
