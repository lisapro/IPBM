module ice_mod
  use fabm_types, only: rk

  implicit none
  private
  public :: ice

  type ice
    private
    !depth of layers(z), upper faces
    real(rk),allocatable,dimension(:):: z
  !contains
  !  private
  end type ice

  interface ice
    procedure constructor_ice
  end interface
contains
  function constructor_ice(number_of_layers)
    type(ice):: constructor_ice
    integer,intent(in):: number_of_layers

    allocate(constructor_ice%z(number_of_layers))
  end function constructor_ice
end module ice_mod
