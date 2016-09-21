program main
  use types_mod
  use variables_mod
  use fabm
  use fabm_config

  implicit none
  integer i
  integer number_of_parameters
  integer number_of_layers
  real(rk),allocatable,dimension(:):: temp
  real(rk),allocatable,dimension(:):: salt
  real(rk),allocatable,dimension(:):: radiative_flux
  real(rk),allocatable,dimension(:):: pressure
  !fabm model
  type(type_model) fabm_model
  !standatd variables for model
  type(brom_standard_variables) standard_vars
  type(brom_state_variable),allocatable,dimension(:):: state_vars

  !initializing fabm
  call fabm_create_model_from_yaml_file(fabm_model)
  !initializing standard_variables
  standard_vars = brom_standard_variables()

  call standard_vars%get_z_length('depth',number_of_layers)

  call fabm_set_domain(fabm_model,number_of_layers)
  call fabm_model%set_surface_index(1)
  call fabm_model%set_bottom_index(number_of_layers)

  number_of_parameters = size(fabm_model%state_variables)
  allocate(state_vars(number_of_parameters))

  do i = 1,number_of_parameters
    allocate(state_vars(i)%value(number_of_layers))
    call fabm_link_bulk_state_data(&
      fabm_model,i,state_vars(i)%value)
    state_vars(i)%name = fabm_model%state_variables(i)%name
  end do

  !linking bulk variables
  allocate(temp(number_of_layers))
  !call standard_vars%get_column('temp',1,temp)
  call fabm_link_bulk_data(&
    fabm_model,standard_variables%temperature,temp)
  allocate(salt(number_of_layers))
  !call standard_vars%get_column('salt',1,salt)
  call fabm_link_bulk_data(&
    fabm_model,standard_variables%practical_salinity,salt)
  allocate(radiative_flux(number_of_layers))
  call fabm_link_bulk_data(&
    fabm_model,&
    standard_variables%downwelling_photosynthetic_radiative_flux,&
    radiative_flux)
  allocate(pressure(number_of_layers))
  !call standard_vars%get_column('depth',1,pressure)
  !pressure = pressure + 10._rk
  !call fabm_link_bulk_data(&
  !  fabm_model,standard_variables%density,host_model%density)
  call fabm_link_bulk_data(&
    fabm_model,standard_variables%pressure,pressure)
  !
  !linking horizontal variables
  call fabm_link_horizontal_data(&
    fabm_model,standard_variables%wind_speed,5._rk)
  call fabm_link_horizontal_data(&
    fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
    380._rk)

  call fabm_check_ready(fabm_model)
end program
