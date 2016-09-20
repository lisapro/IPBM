program main
  use ice_host
  use fabm
  use fabm_config

  implicit none
  integer:: i
  !fabm model
  type(type_model):: fabm_model
  !ice_host model
  type(standard_variables):: standard_vars

  !initializing fabm
  call fabm_create_model_from_yaml_file(fabm_model)
  !initializing host model
  standard_vars = standard_variables()
  !call fabm_set_domain(fabm_model,host_model%number_of_layers)
  !call fabm_model%set_surface_index(1)
  !call fabm_model%set_bottom_index(host_model%number_of_layers)

  !forall(i = 1:host_model%number_of_layers)
  !  call fabm_link_bulk_state_data(&
  !    fabm_model,i,host_model%state_variables(i))
  !end forall

  !!linking bulk variables
  !call fabm_link_bulk_data(&
  !  fabm_model,standard_variables%temperature,host_model%temperature)
  !call fabm_link_bulk_data(&
  !  fabm_model,standard_variables%practical_salinity,&
  !  host_model%practical_salinity)
  !call fabm_link_bulk_data(&
  !  fabm_model,&
  !  standard_variables%downwelling_photosynthetic_radiative_flux,&
  !  host_model%downwelling_photosynthetic_radiative_flux)
  !call fabm_link_bulk_data(&
  !  fabm_model,standard_variables%density,host_model%density)
  !call fabm_link_bulk_data(&
  !  fabm_model,standard_variables%pressure,host_model%pressure)
  !
  !!linking horizontal variables
  !call fabm_link_horizontal_data(&
  !  fabm_model,standard_variables%wind_speed,host_model%wind_speed)
  !call fabm_link_horizontal_data(&
  !  fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
  !  host_model%mole_fraction_of_carbon_dioxide_in_air)

  call fabm_check_ready(fabm_model)
end program
