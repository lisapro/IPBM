#include "../include/brom.h"

module transport
  use types_mod
  use variables_mod
  use fabm
  use fabm_config

  implicit none
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
contains
  subroutine initialize_brom()
    integer i

    !initializing fabm
    call fabm_create_model_from_yaml_file(fabm_model)
    !initializing standard_variables
    standard_vars = brom_standard_variables()
    number_of_layers = standard_vars%&
      get_1st_dim_length(_MIDDLE_LAYER_DEPTH_)
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
      call fabm_initialize_state(fabm_model,1,i)
    end do
    !linking bulk variables
    allocate(temp(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%temperature,temp)
    allocate(salt(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%practical_salinity,salt)
    allocate(radiative_flux(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,&
      standard_variables%downwelling_photosynthetic_radiative_flux,&
      radiative_flux)
    allocate(pressure(number_of_layers))
    call standard_vars%get_column(_MIDDLE_LAYER_DEPTH_,1,pressure)
    !convert depth to pressure
    pressure = pressure + 10._rk
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%pressure,pressure)
    !linking horizontal variables
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%wind_speed,5._rk)
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
      380._rk)

    call fabm_check_ready(fabm_model)
  end subroutine

  subroutine sarafan()
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer number_of_days
    integer day
    integer i

    number_of_days = standard_vars%get_1st_dim_length('day_number')
    day = standard_vars%first_day()
    call initial_date(day,year)

    do i = 1,number_of_days
      !call standard_vars%get_column('temp',i,temp)
      !call standard_vars%get_column('salt',i,salt)
      call date(day,year)
      write(*,*) i,day,year
      day = day+1
    end do
  end subroutine

  subroutine initial_date(day,first_year)
    integer,intent(inout):: day
    integer,intent(inout):: first_year
    integer nydays

    nydays = 365+merge(1,0,(mod(first_year,4).eq.0))
    do
      if (day<=nydays) exit
      nydays = 365+merge(1,0,(mod(first_year,4).eq.0))
      day = day-nydays
      first_year = first_year+1
    end do
  end subroutine

  subroutine date(day,year)
    integer,intent(inout):: day
    integer,intent(inout):: year
    integer days

    days = 365+merge(1,0,(mod(year,4).eq.0))
    if (day==days+1) then
      year = year+1
      day = day-days
    end if
  end subroutine
end module

program main
  use transport

  call initialize_brom()
  call sarafan()
end program
