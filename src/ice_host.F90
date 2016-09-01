module ice_host
  use types

  implicit none
  private
  public type_host_model

  type type_host_model
    type(alone_variable):: dt
    type(alone_variable):: wind_speed
    type(alone_variable):: mole_fraction_of_carbon_dioxide_in_air
    type(alone_variable):: year
    type(alone_variable):: day
    type(alone_variable):: number_of_layers
    type(alone_variable):: number_of_ice_layers
    type(alone_variable):: number_of_state_variables

    type(variable_1d):: ice_thickness
    type(variable_1d):: snow_thickness
    type(variable_1d):: ice_temperature

    type(variable_1d):: z
    type(variable_1d):: dz
    type(variable_1d):: pressure
    type(variable_1d):: downwelling_photosynthetic_radiative_flux
    type(variable_1d):: kz_bio

    type(variable_2d):: temperature
    type(variable_2d):: practical_salinity
    type(variable_2d):: density
    type(variable_2d):: kz

    type(state_variable),allocatable,dimension(:):: state_variables

  contains
    private
    procedure:: initialize
  end type

  interface type_host_model
    module procedure model_constructor
  end interface

contains

  function model_constructor()
    type(type_host_model):: model_constructor

    call model_constructor%initialize()
  end function

  subroutine initialize(self)
    use input

    class(type_host_model):: self
    type(type_input):: kara_input

    kara_input = type_input('KaraSea.nc')
    !self%number_of_layers = kara_input%get_input('h')
    !allocate(self%number_of_layers,source=kara_input%get_input('h'))
    !self%number_of_ice_layers => kara_input%get('')
    !self%number_of_state_variables => kara_input%get('')
    !self%ice_thickness = kara_input%get('')
    !
    !allocate(self%temperature%value, source=kara_imput%get(''))
  end subroutine
end module
