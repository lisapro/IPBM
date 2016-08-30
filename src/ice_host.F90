module ice_host
  use types

  implicit none
  private
  public model

  type model
    type(alone_variable):: dt
    type(alone_variable):: wind_speed
    type(alone_variable):: mole_fraction_of_carbon_dioxide_in_air

    type(horizontal_variable):: ice_thickness
    type(horizontal_variable):: snow_thickness
    type(horizontal_variable):: ice_temperature

    type(vertical_variable):: z
    type(vertical_variable):: dz

    type(standard_variable):: temperature
    type(standard_variable):: practical_salinity
    type(standard_variable):: density
    type(standard_variable):: pressure
    type(standard_variable):: downwelling_photosynthetic_radiative_flux
    type(standard_variable):: kz
    type(standard_variable):: kz_bio

    type(state_variable),allocatable,dimension(:):: state_variables

    integer:: year
    integer:: day
    integer,pointer:: number_of_layers
    integer:: number_of_ice_layers
    integer:: number_of_state_variables
  contains
    private
    procedure:: initialize
  end type

  interface model
    module procedure model_constructor
  end interface

contains

  function model_constructor()
    type(model):: model_constructor

    call model_constructor%initialize()
  end function

  subroutine initialize(self)
    use input

    class(model):: self
    type(type_input):: kara_input

    kara_input = type_input('KaraSea.nc')
    self%number_of_layers => kara_input%get('')
    self%number_of_ice_layers = kara_input%get('')
    self%number_of_state_variables = kara_input%get('')
    self%ice_thickness = kara_input%get('')
    
    allocate(self%temperature%value, source=kara_imput%get(''))
    
  end subroutine

end module
