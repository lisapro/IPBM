module ice_host
  use types
  use input
  use fabm_driver

  implicit none
  private
  public standard_variables

  type,extends(list_variables):: standard_variables
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
  end type

  interface standard_variables
    module procedure standard_variables_constructor
  end interface
contains
  function standard_variables_constructor()
    type(standard_variables):: standard_variables_constructor

    call standard_variables_constructor%initialize()
  end function

  subroutine initialize_standard_variables(self)
    class(standard_variables):: self
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

  subroutine add_standard_var(self,inname,name_input)
    class(standard_variables)  :: self
    character(len=*),intent(in):: inname
    type(type_input),intent(in):: name_input
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_input(inname,var)
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine

end module
