module ice_host
  use types

  implicit none
  private
  public type_host_model

  type type_host_model
    type(alone_variable),allocatable:: dt
    type(alone_variable),allocatable:: wind_speed
    type(alone_variable),allocatable:: mole_fraction_of_carbon_dioxide_in_air
    type(alone_variable),allocatable:: year
    type(alone_variable),allocatable:: day
    type(alone_variable),allocatable:: number_of_layers
    type(alone_variable),allocatable:: number_of_ice_layers
    type(alone_variable),allocatable:: number_of_state_variables

    type(variable_1d),allocatable:: ice_thickness
    type(variable_1d),allocatable:: snow_thickness
    type(variable_1d),allocatable:: ice_temperature

    type(variable_1d),allocatable:: z
    type(variable_1d),allocatable:: z_boundary
    type(variable_1d),allocatable:: time
    type(variable_1d),allocatable:: dz
    type(variable_1d),allocatable:: air_pressure
    type(variable_1d),allocatable:: downwelling_photosynthetic_radiative_flux
    type(variable_1d),allocatable:: kz_bio

    type(variable_2d),allocatable:: temperature
    type(variable_2d),allocatable:: practical_salinity
    type(variable_2d),allocatable:: density_anomaly
    type(variable_2d),allocatable:: kz

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
    !vertical variables
    call kara_input%get_input('depth',self%z)
    call kara_input%get_input('depth_w',self%z_boundary)
    !horizontal variables
    call kara_input%get_input('ocean_time',self%time)
    call kara_input%get_input('Pair',self%air_pressure)
    call kara_input%get_input(&
    'shflux',self%downwelling_photosynthetic_radiative_flux)
    !2d variables
    call kara_input%get_input('temp',self%temperature)
    call kara_input%get_input('salt',self%practical_salinity)
    call kara_input%get_input('rho',self%density_anomaly)
    call kara_input%get_input('AKv',self%kz)
  end subroutine
end module
