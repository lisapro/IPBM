#include "../include/brom.h"
#include "../include/parameters.h"

module transport
  use fabm
  use fabm_config
  use fabm_driver
  use fabm_types,only: rk
  use variables_mod

  implicit none
  private
  public initialize_brom,sarafan

  integer number_of_parameters
  integer number_of_layers
  real(rk),allocatable,dimension(:):: temp
  real(rk),allocatable,dimension(:):: salt
  real(rk),allocatable,dimension(:):: turb
  real(rk),allocatable,dimension(:):: radiative_flux
  real(rk),allocatable,dimension(:):: depth
  real(rk),allocatable,dimension(:):: pressure
  real(rk),allocatable,dimension(:):: layer_thicknesses
  !fabm model
  type(type_model) fabm_model
  !standard variables for model
  type(brom_standard_variables) standard_vars
  type(brom_state_variable),allocatable,dimension(:):: state_vars
contains
  subroutine initialize_brom()
    integer i

    !initializing fabm
    _LINE_
    call fabm_create_model_from_yaml_file(fabm_model)
    _LINE_
    !initializing standard_variables
    standard_vars = brom_standard_variables()
    number_of_layers = standard_vars%&
      get_1st_dim_length(_MIDDLE_LAYER_DEPTH_)
    call fabm_set_domain(fabm_model,number_of_layers)
    call fabm_model%set_surface_index(number_of_layers)
    call fabm_model%set_bottom_index(1)
    number_of_parameters = size(fabm_model%state_variables)
    allocate(state_vars(number_of_parameters))
    do i = 1,number_of_parameters
      allocate(state_vars(i)%value(number_of_layers))
      call fabm_link_bulk_state_data(&
        fabm_model,i,state_vars(i)%value)
      state_vars(i)%name = fabm_model%state_variables(i)%name
      call state_vars(i)%set_brom_state_variable(_NEUMANN_,&
        _NEUMANN_,0._rk,0._rk,0._rk)
      call state_vars(i)%print_name()
    end do
    _LINE_
    call fabm_initialize_state(fabm_model,1,number_of_layers)
    !linking bulk variables
    allocate(temp(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%temperature,temp)
    allocate(salt(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%practical_salinity,salt)
    allocate(turb(number_of_layers+1))
    allocate(radiative_flux(number_of_layers))
    call fabm_link_bulk_data(&
      fabm_model,&
      standard_variables%downwelling_photosynthetic_radiative_flux,&
      radiative_flux)
    allocate(depth(number_of_layers))
    depth = standard_vars%get_column(_MIDDLE_LAYER_DEPTH_)
    allocate(layer_thicknesses(number_of_layers))
    layer_thicknesses = standard_vars%get_column(&
                        "layer_thicknesses")
    !convert depth to pressure
    !total=water+atmosphere [dbar]
    pressure = depth + 10._rk
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%pressure,pressure)
    !linking horizontal variables
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%wind_speed,5._rk)
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
      380._rk)
    call fabm_check_ready(fabm_model)
    call find_set_state_variable(state_vars,"niva_brom_redox_SO4",&
      use_bound_up = _DIRICHLET_,use_bound_low = _DIRICHLET_,&
      bound_up = 25000._rk,bound_low = 25000._rk)
    call find_set_state_variable(state_vars,inname = "niva_brom_redox_Mn4",&
      use_bound_up = _DIRICHLET_,bound_up = 0.5e-4_rk)
    call find_set_state_variable(state_vars,"niva_brom_redox_Fe3",&
      use_bound_up = _DIRICHLET_,bound_up = 0.4e-4_rk)
    call find_set_state_variable(state_vars,"niva_brom_carb_Alk",&
      use_bound_up = _DIRICHLET_,bound_up = 2250._rk)
  end subroutine

  subroutine sarafan()
    use output_mod

    type(brom_state_variable):: temporary_variable
    type(type_output):: netcdf
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer number_of_days
    integer day
    integer i
    !cpu time
    real(rk) t1,t2

    netcdf = type_output(fabm_model,_FILE_NAME_WATER_,&
                         1,number_of_layers,&
                         number_of_layers)
    number_of_days = standard_vars%get_1st_dim_length("day_number")
    day = standard_vars%first_day()
    call initial_date(day,year)
    call stabilize(day,year)

    do i = 1,number_of_days
      call date(day,year)
      radiative_flux = calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),depth)
      temp = standard_vars%get_column(_TEMPERATURE_,i)
      salt = standard_vars%get_column(_SALINITY_,i)
      turb = standard_vars%get_column(_TURBULENCE_,i)
      call find_set_state_variable(state_vars,"niva_brom_bio_PO4",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      call find_set_state_variable(state_vars,"niva_brom_bio_NO3",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      call find_set_state_variable(state_vars,"niva_brom_redox_Si",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,2._rk))
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%temperature,temp)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%practical_salinity,salt)
      call fabm_link_bulk_data(&
        fabm_model,&
        standard_variables%downwelling_photosynthetic_radiative_flux,&
        radiative_flux)
      call cpu_time(t1)
      call day_circle()
      call netcdf%save(fabm_model,state_vars,i,&
                  temp,salt,turb,radiative_flux,depth)
      call cpu_time(t2)
      write(*,*) "number / ","julianday / ","year",i,day,year
      write(*,*) "Time taken by day circle:",t2-t1," seconds"
      day = day+1
      temporary_variable = find_state_variable(state_vars,&
                          "niva_brom_bio_O2")
      call temporary_variable%print_state_variable()
    end do
    write(*,*) "Finish"
    _LINE_
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
    real(rk),dimension(number_of_layers,number_of_parameters):: increment

    if (mod(60*60*24,_SECONDS_PER_CIRCLE_)/=0) then
      call fatal_error("Check _SECONDS_PER_CIRCLE_",&
                       "It should be multiple of 86400")
    else
      number_of_circles = int(60*60*24/_SECONDS_PER_CIRCLE_)
    end if

    do i = 1,number_of_circles
      !biogeochemistry
      increment = 0._rk
      call fabm_do(fabm_model,1,number_of_layers,increment)
      increment = _SECONDS_PER_CIRCLE_*increment
      forall(j = 1:number_of_parameters)&
        state_vars(j)%value = state_vars(j)%value+increment(:,j)

      !diffusion
      call brom_do_diffusion()
    end do
  end subroutine

  subroutine brom_do_diffusion()
    use diff_mod

    real(rk),dimension(number_of_layers+1):: ones
    real(rk),dimension(number_of_layers+1):: zeros
    real(rk),dimension(number_of_layers+1):: taur_r
    real(rk),dimension(number_of_parameters):: surface_flux
    real(rk),dimension(0:number_of_layers,&
                         number_of_parameters):: temporary
    integer i

    ones=1._rk
    zeros=0._rk
    taur_r=1.d20
    surface_flux = 0._rk
    call fabm_do_surface(fabm_model,surface_flux)

    do i = 1,number_of_parameters
      if (surface_flux(i)/=0._rk) then
        call state_vars(i)%set_brom_state_variable(&
          use_bound_up = _NEUMANN_,bound_up = surface_flux(i))
      end if
    end do

    forall (i = 1:number_of_parameters)
      temporary(:,i) = do_diffusive(&
          N       = number_of_layers,&
          dt      = _SECONDS_PER_CIRCLE_,&
          cnpar   = 0.6_rk,&
          posconc = 1,&
          h  = (/ 0._rk,layer_thicknesses /),&
          Bcup = state_vars(i)%use_bound_up,&
          Bcdw = state_vars(i)%use_bound_low,&
          Yup  = state_vars(i)%bound_up,&
          Ydw  = state_vars(i)%bound_low,&
          nuY  = turb,&
          Lsour = zeros,&
          Qsour = zeros,&
          Taur  = taur_r,&
          Yobs  = zeros,&
          Y     = (/ 0._rk,state_vars(i)%value /))
      state_vars(i)%value = temporary(1:number_of_layers,i)
    end forall
  end subroutine

  subroutine stabilize(inday,inyear)
    integer,intent(in):: inday,inyear

    type(brom_state_variable):: temporary_variable
    integer day,year
    integer pseudo_day,days_in_year,counter
    integer i

    day = inday; year = inyear;
    days_in_year = 365+merge(1,0,(mod(year,4).eq.0))
    counter = days_in_year*10
    do i = 1,counter
      call date(day,year)
      radiative_flux = calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),depth)
      !day from 1 to 365 or 366
      pseudo_day = i-int(i/days_in_year)*&
                   days_in_year+1
      temp = standard_vars%get_column(_TEMPERATURE_,pseudo_day)
      salt = standard_vars%get_column(_SALINITY_,pseudo_day)
      turb = standard_vars%get_column(_TURBULENCE_,pseudo_day)

      call find_set_state_variable(state_vars,"niva_brom_bio_PO4",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      call find_set_state_variable(state_vars,"niva_brom_bio_NO3",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      call find_set_state_variable(state_vars,"niva_brom_redox_Si",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,2._rk))
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%temperature,temp)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%practical_salinity,salt)
      call fabm_link_bulk_data(&
        fabm_model,&
        standard_variables%downwelling_photosynthetic_radiative_flux,&
        radiative_flux)
      call day_circle()
      write(*,*) "Stabilizing initial array of values, in progress ..."
      write(*,*) "number / ","julianday / ","pseudo day",&
                 i,day,pseudo_day
      day = day+1
      temporary_variable = find_state_variable(state_vars,&
                          "niva_brom_bio_O2")
      call temporary_variable%print_state_variable()
    end do
  end subroutine

  subroutine find_set_state_variable(state_vars,inname,use_bound_up,&
      use_bound_low,bound_up,bound_low,sinking_velocity)
    type(brom_state_variable),dimension(:),intent(inout):: state_vars
    character(len=*),                      intent(in):: inname
    integer,optional,                      intent(in):: use_bound_up
    integer,optional,                      intent(in):: use_bound_low
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
