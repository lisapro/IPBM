#include "../include/brom.h"
#include "../include/parameters.h"

module transport
  use fabm
  use fabm_config
  use fabm_driver
  use fabm_types,only: rk
  use variables_mod

  implicit none
  integer number_of_parameters
  integer number_of_layers
  real(rk),allocatable,dimension(:):: temp
  real(rk),allocatable,dimension(:):: salt
  real(rk),allocatable,dimension(:):: radiative_flux
  real(rk),allocatable,dimension(:):: depth
  !fabm model
  type(type_model) fabm_model
  !standard variables for model
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
      call state_vars(i)%print_name()
    end do
    call fabm_initialize_state(fabm_model,1,number_of_layers)
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
    allocate(depth(number_of_layers))
    call standard_vars%get_column(_MIDDLE_LAYER_DEPTH_,1,depth)
    !convert depth to pressure
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%pressure,(depth+10._rk))
    !linking horizontal variables
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%wind_speed,5._rk)
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
      380._rk)
    call fabm_check_ready(fabm_model)
    call set_state_variable(state_vars,"niva_brom_redox_SO4",&
      .true.,.true.,25000._rk,25000._rk)
    call set_state_variable(state_vars,"niva_brom_redox_Mn4",&
      .true.,.false.,0.5e-4_rk)
    call set_state_variable(state_vars,"niva_brom_redox_Fe3",&
      .true.,.false.,0.4e-4_rk)
    call set_state_variable(state_vars,"niva_brom_carb_Alk",&
      use_bound_up = .true.,bound_up = 2250._rk)
  end subroutine

  subroutine sarafan()
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer number_of_days
    integer day
    integer i

    number_of_days = standard_vars%get_1st_dim_length("day_number")
    day = standard_vars%first_day()
    call initial_date(day,year)

    do i = 1,1!number_of_days
      call date(day,year)
      write(*,*) i,day,year
      radiative_flux = calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),depth)
      call standard_vars%get_column(_TEMPERATURE_,i,temp)
      call standard_vars%get_column(_SALINITY_,i,salt)
      call set_state_variable(state_vars,"niva_brom_bio_PO4",&
        use_bound_up = .true.,bound_up = sinusoidal(day,0.45_rk))
      call set_state_variable(state_vars,"niva_brom_bio_NO3",&
        use_bound_up = .true.,bound_up = sinusoidal(day,3.8_rk))
      call set_state_variable(state_vars,"niva_brom_redox_Si",&
        use_bound_up = .true.,bound_up = sinusoidal(day,2._rk))
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%temperature,temp)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%practical_salinity,salt)
      call fabm_link_bulk_data(&
        fabm_model,&
        standard_variables%downwelling_photosynthetic_radiative_flux,&
        radiative_flux)
      call day_circle()
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

  pure function surface_radiative_flux(latitude,julian_day)
    real(rk),intent(in):: latitude
    integer, intent(in):: julian_day
    real(rk) surface_radiative_flux

    surface_radiative_flux = 80._rk*cos((latitude-(23.5_rk*&
      sin(2._rk*_PI_*(julian_day-81._rk)/365._rk)))*3.14_rk/180._rk)
    if (surface_radiative_flux<0._rk) surface_radiative_flux = 0._rk
  end function

  pure function calculate_radiative_flux(surface_flux,depth)
    real(rk),intent(in)             :: surface_flux
    real(rk),dimension(:),intent(in):: depth
    real(rk),dimension(size(depth)):: calculate_radiative_flux

    calculate_radiative_flux = surface_flux*exp(-_ERLOV_*depth)
  end function

  pure function sinusoidal(julian_day,multiplier)
    integer, intent(in):: julian_day
    real(rk),intent(in):: multiplier
    real(rk) sinusoidal

    sinusoidal = (1._rk+sin(2._rk*_PI_*(&
                  julian_day-40._rk)/365._rk))*multiplier
  end function

  subroutine day_circle()
    integer i,j
    integer number_of_circles
    type(brom_state_variable):: temporary_variable
    real(rk),dimension(number_of_parameters):: temporary_vector
    real(rk),dimension(number_of_layers,number_of_parameters):: temporary_matrix

    if (mod(60*60*24,_SECONDS_PER_CIRCLE_)/=0) then
      call fatal_error("Check _SECONDS_PER_CIRCLE_",&
                       "It should be multiple of 86400")
    else
      number_of_circles = int(60*60*24/_SECONDS_PER_CIRCLE_)
    end if

    do i = 1,number_of_circles
      !biogeochemistry
      temporary_matrix = 0._rk
      call fabm_do(fabm_model,1,number_of_layers,temporary_matrix)
      temporary_vector = 0._rk
      call fabm_do_surface(fabm_model,temporary_vector)
      temporary_matrix(:,1) = temporary_matrix(:,1)+temporary_vector
      temporary_matrix = _SECONDS_PER_CIRCLE_*temporary_matrix
      forall(j = 1:number_of_parameters)&
        state_vars(j)%value = state_vars(j)%value+temporary_matrix(:,j)
    end do
    temporary_variable = find_state_variable(state_vars,&
                        "niva_brom_bio_O2")
    call temporary_variable%print_state_variable()
  end subroutine

  subroutine set_state_variable(state_vars,inname,use_bound_up,&
      use_bound_low,bound_up,bound_low,sinking_velocity)
    type(brom_state_variable),dimension(:),intent(inout):: state_vars
    character(len=*),                      intent(in):: inname
    logical,optional,                      intent(in):: use_bound_up
    logical,optional,                      intent(in):: use_bound_low
    real(rk),optional,                     intent(in):: bound_up
    real(rk),optional,                     intent(in):: bound_low
    real(rk),optional,                     intent(in):: sinking_velocity
    integer number_of_vars
    integer i

    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq.inname) then
        call state_vars(i)%set_brom_state_variable(use_bound_up,&
          use_bound_low,bound_up,bound_low,sinking_velocity)
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end subroutine

  function find_state_variable(state_vars,inname)
    type(brom_state_variable),dimension(:),intent(in):: state_vars
    character(len=*),                      intent(in):: inname
    type(brom_state_variable):: find_state_variable
    integer number_of_vars
    integer i

    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq.inname) then
        find_state_variable = state_vars(i)
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end function
end module

program main
  use transport

  call initialize_brom()
  call sarafan()
end program
