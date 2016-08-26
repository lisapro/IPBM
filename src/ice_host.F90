module ice_host

  use types
  use input
  use fabm_types, only: attribute_length, rk

  implicit none
  private
  public model

  type model
    type(horizontal_variable):: ice_thickness
    type(horizontal_variable):: snow_thickness
    type(horizontal_variable):: ice_temperature

    type(vertical_variable):: z
    type(vertical_variable):: dz
    type(vertical_variable):: kz_bio
    type(vertical_variable):: density
    type(vertical_variable):: pressure
    type(vertical_variable):: irradiance

    type(standart_variable):: temperature
    type(standart_variable):: salinity
    type(standart_variable):: kz

    type(state_variable),dimension(:):: state_variables

    integer:: year
    integer:: day
    integer:: number_of_layers
    integer:: number_of_ice_layers
    integer:: number_of_state_variables

    real(rk):: dt
    real(rk):: surface_irradiance
    real(rk):: wind_speed
    real(rk):: pco2_atm
  contains
    private
    procedure:: initialize
  end type

contains

  subroutine initialize(self)
    class(model):: self
    type(type_input):: kara_input

    kara_input = type_input('KaraSea.nc')
    self%ice_thickness = kara_input%get('')
  end subroutine

end module
