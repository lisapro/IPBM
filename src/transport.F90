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
  real(rk),allocatable,dimension(:),target:: temp
  real(rk),allocatable,dimension(:),target:: salt
  real(rk),allocatable,dimension(:),target:: radiative_flux
  real(rk),allocatable,dimension(:),target:: pressure
  real(rk),allocatable,dimension(:):: turb
  real(rk),allocatable,dimension(:):: depth
  !fabm model
  type(type_model) fabm_model
  !standard variables for model
  type(brom_standard_variables) standard_vars
  type(brom_state_variable),allocatable,&
                       dimension(:),target:: state_vars
contains
  subroutine initialize_brom()
    integer i

    !initializing fabm
    _LINE_
    call fabm_create_model_from_yaml_file(fabm_model)
    _LINE_
    !initializing standard_variables
    standard_vars = brom_standard_variables()
    number_of_layers = standard_vars%get_value(&
                           "number_of_layers")
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
      call state_vars(i)%set_brom_state_variable(.false.,_NEUMANN_,&
        _NEUMANN_,0._rk,0._rk,0._rk,0._rk)
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
    depth = standard_vars%get_column("middle_layer_depths",1)
    !convert depth to pressure
    !total=water+atmosphere [dbar]
    allocate(pressure(number_of_layers))
    pressure = depth + 10._rk
    pressure(standard_vars%get_value("ice_water_index"):) = 10._rk
    call fabm_link_bulk_data(&
      fabm_model,standard_variables%pressure,pressure)
    !linking horizontal variables
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%wind_speed,5._rk)
    call fabm_link_horizontal_data(&
      fabm_model,standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
      380._rk)
    call fabm_check_ready(fabm_model)
    call configurate_state_variables()
  end subroutine

  subroutine sarafan()
    use output_mod

    !type(brom_state_variable):: temporary_variable
    type(type_output):: netcdf_water
    type(type_output):: netcdf_sediments
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer ice_water_index,water_bbl_index,number_of_days
    integer day
    integer i
    !cpu time
    real(rk) t1,t2

    !temp
    !type(brom_state_variable):: oxygen
    !real(rk),dimension(number_of_layers+1):: kz_mol
    !real(rk),dimension(number_of_layers+1):: kz_bio
    !real(rk),dimension(number_of_layers+1):: kz_tot
    !integer bbl_sed_index
    !real(rk) O2stat
    !temp

    ice_water_index = standard_vars%get_value("ice_water_index")
    water_bbl_index = standard_vars%get_value("water_bbl_index")
    netcdf_water = type_output(fabm_model,_FILE_NAME_WATER_,&
                         water_bbl_index,ice_water_index-1,&
                         number_of_layers)
    netcdf_sediments = type_output(fabm_model,_FILE_NAME_SEDIMENTS_,&
                         1,water_bbl_index-1,&
                         number_of_layers)
    number_of_days = standard_vars%get_1st_dim_length("day_number")
    day = standard_vars%first_day()
    call initial_date(day,year)
    !call stabilize(day,year)

    do i = 1,number_of_days
      call date(day,year)
      call calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),&
        standard_vars%get_value(_SNOW_THICKNESS_,i),&
        standard_vars%get_value(_ICE_THICKNESS_ ,i))
      temp  = standard_vars%get_column(_TEMPERATURE_,i)
      salt  = standard_vars%get_column(_SALINITY_,i)
      turb  = standard_vars%get_column(_TURBULENCE_,i)
      call find_set_state_variable("niva_brom_bio_PO4",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      call find_set_state_variable("niva_brom_bio_NO3",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      call find_set_state_variable("niva_brom_redox_Si",&
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
      depth = standard_vars%get_column("middle_layer_depths",i)
      call day_circle(i)
      call netcdf_water%save(fabm_model,state_vars,i,&
                  temp,salt,turb,radiative_flux,depth)
      call netcdf_sediments%save(fabm_model,state_vars,i,&
                  temp,salt,turb,radiative_flux,depth)
      call cpu_time(t2)
      write(*,*) "number / ","julianday / ","year",i,day,year
      write(*,*) "Time taken by day circle:",t2-t1," seconds"
      day = day+1
      !temporary_variable = find_state_variable("niva_brom_bio_O2")
      !call temporary_variable%print_state_variable()

      !temp
      !kz_mol = standard_vars%get_column("molecular_diffusivity")
      !kz_bio = standard_vars%get_column("bioturbation_diffusivity")
      !oxygen = find_state_variable("niva_brom_bio_O2")
      !bbl_sed_index = standard_vars%get_value("bbl_sediments_index")
      !!oxygen status of sediments
      !O2stat = oxygen%value(bbl_sed_index)/&
      !(oxygen%value(bbl_sed_index)+_KO2_)
      !kz_tot = turb+kz_mol+kz_bio*O2stat
      !write(*,'(f20.15)') (/ kz_tot(size(turb,1):1:-1) /)
      !temp
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

  subroutine calculate_radiative_flux(surface_flux,snow_thick,ice_thick)
    real(rk),intent(in):: surface_flux
    real(rk),intent(in):: snow_thick
    real(rk),intent(in):: ice_thick
    integer  ice_water_index
    real(rk) par_alb,par_scat
    real(rk),allocatable:: ice_depths(:)
    
    ice_water_index = standard_vars%get_value("ice_water_index")
    !ice_column
    !surface_flux in Watts, to calculate it in micromoles photons per m2*s =>
    !=> [w] = 4.6*[micromole photons]
    !Grenfell and Maykutt 1977 indicate that the magnitude and shape
    !of the albedo curves depend strongly on the amount of liquid
    !water present in the upper part of the ice, so it fluctuates
    !throught year (true also for extinction coefficient)
    par_alb = surface_flux*(1._rk-_ICE_ALBEDO_)
    !after scattered surface of ice
    par_scat = par_alb*_ICE_SCATTERED_
    allocate(ice_depths,source=ice_thick-depth(ice_water_index:))
    radiative_flux(ice_water_index:) = &
      par_scat*exp(-_ICE_EXTINCTION_*ice_depths)
    !water column calculations
    radiative_flux(:ice_water_index-1) = &
      surface_flux*exp(-_ERLOV_*depth(:ice_water_index-1))
  end subroutine

  pure function sinusoidal(julian_day,multiplier)
    integer, intent(in):: julian_day
    real(rk),intent(in):: multiplier
    real(rk) sinusoidal

    sinusoidal = (1._rk+sin(2._rk*_PI_*(&
                  julian_day-40._rk)/365._rk))*multiplier
  end function

  subroutine day_circle(id)
    integer,intent(in):: id !number of the count
    
    real(rk),dimension(number_of_layers+1):: pF1_solutes
    real(rk),dimension(number_of_layers+1):: pF2_solutes
    real(rk),dimension(number_of_layers+1):: pF1_solids
    real(rk),dimension(number_of_layers+1):: pF2_solids
    real(rk),dimension(number_of_layers+1):: kz_mol
    real(rk),dimension(number_of_layers+1):: kz_bio
    real(rk),dimension(number_of_layers+1):: layer_thicknesses
    integer i,j,bbl_sed_index
    integer number_of_circles
    real(rk),dimension(number_of_layers,number_of_parameters):: increment

    bbl_sed_index = standard_vars%get_value("bbl_sediments_index")
    pF1_solutes = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solutes_1",id) /)
    pF2_solutes = standard_vars%get_column("porosity_factor_solutes_2",id)
    pF1_solids = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solids_1",id) /)
    pF2_solids = standard_vars%get_column("porosity_factor_solids_2",id)
    kz_mol = standard_vars%get_column("molecular_diffusivity",id)
    kz_bio = standard_vars%get_column("bioturbation_diffusivity",id)
    layer_thicknesses = &
    (/ 0._rk,standard_vars%get_column("layer_thicknesses",id) /)

    if (mod(60*60*24,_SECONDS_PER_CIRCLE_)/=0) then
      call fatal_error("Check _SECONDS_PER_CIRCLE_",&
                       "Wrong value")
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
      call brom_do_diffusion(bbl_sed_index,&
                             pF1_solutes,pF2_solutes,pF1_solids,&
                             pF2_solids,kz_mol,kz_bio,&
                             layer_thicknesses)
    end do
  end subroutine

  subroutine brom_do_diffusion(bbl_sed_index,&
                               pF1_solutes,pF2_solutes,pF1_solids,&
                               pF2_solids,kz_mol,kz_bio,&
                               layer_thicknesses)
    use diff_mod
    integer,intent(in):: bbl_sed_index
    real(rk),dimension(number_of_layers+1),intent(in):: pF1_solutes
    real(rk),dimension(number_of_layers+1),intent(in):: pF2_solutes
    real(rk),dimension(number_of_layers+1),intent(in):: pF1_solids
    real(rk),dimension(number_of_layers+1),intent(in):: pF2_solids
    real(rk),dimension(number_of_layers+1),intent(in):: kz_mol
    real(rk),dimension(number_of_layers+1),intent(in):: kz_bio
    real(rk),dimension(number_of_layers+1),intent(in):: layer_thicknesses
    
    type(brom_state_variable):: oxygen
    real(rk),dimension(number_of_layers+1):: ones
    real(rk),dimension(number_of_layers+1):: zeros
    real(rk),dimension(number_of_layers+1):: taur_r
    real(rk),dimension(number_of_layers+1):: kz_tot
    real(rk),dimension(number_of_parameters):: surface_flux
    real(rk),dimension(0:number_of_layers,&
                         number_of_parameters):: temporary
    integer i
    real(rk) O2stat
    real(rk) pFSWIup_solutes, pFSWIdw_solutes
    real(rk) pFSWIup_solids , pFSWIdw_solids

    ones  = 1._rk
    zeros = 0._rk
    taur_r= 1.e20_rk

    oxygen = find_state_variable("niva_brom_bio_O2")
    !oxygen status of sediments
    O2stat = oxygen%value(bbl_sed_index)/&
      (oxygen%value(bbl_sed_index)+_KO2_)
    kz_tot = turb+kz_mol+kz_bio*O2stat

    surface_flux = 0._rk
    call fabm_do_surface(fabm_model,surface_flux)
    do i = 1,number_of_parameters
      if (surface_flux(i)/=0._rk) then
        call state_vars(i)%set_brom_state_variable(&
          use_bound_up = _NEUMANN_,bound_up = surface_flux(i))
      end if
    end do
    !
    !adopted from Phil Wallhead (PW):
    !solutes:
    !fick_SWI = -phi_0*Km/dz*(C_1/phi_1 - C_-1) - Kb/dz*(C_1 - C_-1)
    !         = -phi_0*(Km+Kb)/dz * (p_1*C_1 - p_-1*C_-1)
    !         [intraphase molecular + interphase bioturb]
    !where:
    !       p_-1 = (phi_0*Km + Kb) / (phi_0*(Km+Kb))
    !       p_1  = (phi_0/phi_1*Km + Kb) / (phi_0*(Km+Kb))
    !and subscripts refer to SWI (0), below (1), and above (-1)
    !
    pFSWIup_solutes = (pF2_solutes(bbl_sed_index)*&
      kz_mol(bbl_sed_index)+kz_bio(bbl_sed_index)*O2stat)/&
      (pF2_solutes(bbl_sed_index)*(kz_mol(bbl_sed_index)+&
      kz_bio(bbl_sed_index)*O2stat))
    pFSWIdw_solutes = (pF1_solutes(bbl_sed_index)*&
      pF2_solutes(bbl_sed_index)*kz_mol(bbl_sed_index)+&
      kz_bio(bbl_sed_index)*O2stat)/(pF2_solutes(bbl_sed_index)*&
      (kz_mol(bbl_sed_index)+kz_bio(bbl_sed_index)*O2stat))
    !
    !(PW) solids:
    !fick_SWI = -Kb/dz*(C_1 - C_-1)   [interphase bioturb]
    !         = -(1-phi_0)*Kb/dz * (p_1*C_1 - p_-1*C_-1)
    !where: p_-1 = 1 / (1 - phi_0)
    !       p_1  = 1 / (1 - phi_0)
    !and subscripts refer to SWI (0), below (1), and above (-1)
    !
    pFSWIup_solids = 1.0_rk/(1.0_rk-pF2_solutes(bbl_sed_index))
    pFSWIdw_solids = pFSWIup_solids

    forall (i = 1:number_of_parameters)
      temporary(:,i) = do_diffusive(&
          N       = number_of_layers,&
          dt      = _SECONDS_PER_CIRCLE_,&
          cnpar   = 0.6_rk,&
          posconc = 1,&
          h    = layer_thicknesses,&
          Bcup = state_vars(i)%use_bound_up,&
          Bcdw = state_vars(i)%use_bound_low,&
          Yup  = state_vars(i)%bound_up,&
          Ydw  = state_vars(i)%bound_low,&
          nuY_in  = kz_tot,&
          Lsour = zeros,&
          Qsour = zeros,&
          Taur  = taur_r,&
          Yobs  = zeros,&
          Y     = (/ 0._rk,state_vars(i)%value /),&
          i_sed_top = bbl_sed_index-1,&
          is_solid = state_vars(i)%is_solid,&
          pF1_solutes = pF1_solutes,&
          pF2_solutes = pF2_solutes,&
          pF1_solids = pF1_solids,&
          pF2_solids = pF2_solids,&
          pFSWIup_solutes = pFSWIup_solutes,&
          pFSWIdw_solutes = pFSWIdw_solutes,&
          pFSWIup_solids = pFSWIup_solids,&
          pFSWIdw_solids = pFSWIdw_solids)
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
      call calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),&
        standard_vars%get_value(_SNOW_THICKNESS_,i),&
        standard_vars%get_value(_ICE_THICKNESS_ ,i))
      !day from 1 to 365 or 366
      pseudo_day = i-int(i/days_in_year)*&
                   days_in_year+1
      temp = standard_vars%get_column(_TEMPERATURE_,pseudo_day)
      salt = standard_vars%get_column(_SALINITY_,pseudo_day)
      turb = standard_vars%get_column(_TURBULENCE_,pseudo_day)

      call find_set_state_variable("niva_brom_bio_PO4",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      call find_set_state_variable("niva_brom_bio_NO3",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      call find_set_state_variable("niva_brom_redox_Si",&
        use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,2._rk))
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%temperature,temp)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%practical_salinity,salt)
      call fabm_link_bulk_data(&
        fabm_model,&
        standard_variables%downwelling_photosynthetic_radiative_flux,&
        radiative_flux)
      call day_circle(i)
      write(*,*) "Stabilizing initial array of values, in progress ..."
      write(*,*) "number / ","julianday / ","pseudo day",&
                 i,day,pseudo_day
      day = day+1
      !temporary_variable = find_state_variable("niva_brom_bio_O2")
      !call temporary_variable%print_state_variable()
    end do
  end subroutine

  subroutine find_set_state_variable(inname,is_solid,&
      use_bound_up,use_bound_low,bound_up,bound_low,&
      density,sinking_velocity)
    character(len=*),                      intent(in):: inname
    logical,optional,                      intent(in):: is_solid
    integer,optional,                      intent(in):: use_bound_up
    integer,optional,                      intent(in):: use_bound_low
    real(rk),optional,                     intent(in):: bound_up
    real(rk),optional,                     intent(in):: bound_low
    real(rk),optional,                     intent(in):: density
    real(rk),optional,                     intent(in):: sinking_velocity
    integer number_of_vars
    integer i

    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq.inname) then
        call state_vars(i)%set_brom_state_variable(is_solid,&
          use_bound_up,use_bound_low,bound_up,bound_low,&
          density,sinking_velocity)
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end subroutine

  subroutine configurate_state_variables()
    !call find_set_state_variable("niva_brom_redox_SO4",&
    !  use_bound_up = _DIRICHLET_,use_bound_low = _DIRICHLET_,&
    !  bound_up = 25000._rk,bound_low = 25000._rk)
    call find_set_state_variable(inname = "niva_brom_redox_Mn4",&
      use_bound_up = _DIRICHLET_,bound_up = 0.5e-4_rk)
    call find_set_state_variable("niva_brom_redox_Fe3",&
      use_bound_up = _DIRICHLET_,bound_up = 0.4e-4_rk)
    !call find_set_state_variable("niva_brom_carb_Alk",&
    !  use_bound_up = _DIRICHLET_,bound_up = 2250._rk)

    call find_set_state_variable("niva_brom_bio_Phy",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_bio_PON",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_bio_Het",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_redox_Baae",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_redox_Bhae",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_redox_Baan",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_redox_Bhan",&
      is_solid = .true.,density = 1.5E7_rk)
    call find_set_state_variable("niva_brom_redox_CaCO3",&
      is_solid = .true.,density = 2.80E7_rk)
    call find_set_state_variable("niva_brom_redox_Fe3",&
      is_solid = .true.,density = 3.27E7_rk)
    call find_set_state_variable("niva_brom_redox_FeCO3",&
      is_solid = .true.,density = 2.93E7_rk)
    call find_set_state_variable("niva_brom_redox_FeS",&
      is_solid = .true.,density = 5.90E7_rk)
    call find_set_state_variable("niva_brom_redox_FeS2",&
      is_solid = .true.,density = 4.17E7_rk)
    call find_set_state_variable("niva_brom_redox_Mn4",&
      is_solid = .true.,density = 5.78E7_rk)
    call find_set_state_variable("niva_brom_redox_MnCO3",&
      is_solid = .true.,density = 3.20E7_rk)
    call find_set_state_variable("niva_brom_redox_MnS",&
      is_solid = .true.,density = 4.60E7_rk)
    call find_set_state_variable("niva_brom_redox_S0",&
      is_solid = .true.,density = 6.56E7_rk)
    call find_set_state_variable("niva_brom_redox_Sipart",&
      is_solid = .true.,density = 4.40E7_rk)
  end subroutine

  function find_state_variable(inname)
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
