!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------
! Original author(s): Shamil Yakubov
!-----------------------------------------------------------------------

#include "../include/brom.h"
#include "../include/parameters.h"

module transport
  use fabm
  use fabm_config
  use fabm_driver
  use fabm_types
  use variables_mod,only: brom_standard_variables,&
                          brom_state_variable

  implicit none
  private
  public initialize_brom,sarafan

  integer previous_ice_index
  integer number_of_parameters
  integer number_of_layers
  real(rk),allocatable,dimension(:):: depth
  real(rk),allocatable,dimension(:):: porosity
  !absorption_of_silt value - ersem light
  !real(rk),allocatable,dimension(:):: aos_value
  !for coupling with fabm
  !realday for link scalar subroutine, bdepth - bottom depth
  real(rk),                         target:: bdepth!,realday
  real(rk),allocatable,dimension(:),target:: cell! - dz
  real(rk),allocatable,dimension(:),target:: temp
  real(rk),allocatable,dimension(:),target:: salt
  real(rk),allocatable,dimension(:),target:: density
  real(rk),allocatable,dimension(:),target:: pressure
  real(rk),allocatable,dimension(:),target:: radiative_flux
  !taub - bottom stress
  real(rk),allocatable,             target:: taub
  !bcc and scc - arrays for bottom and surface fabm variables
  real(rk),allocatable,dimension(:),target:: bcc,scc
  !variables for model
  type(brom_standard_variables) standard_vars
  type(brom_state_variable),allocatable,&
                       dimension(:),target:: state_vars
  !fabm model
  type(type_model) fabm_model
  !absorption_of_silt - ersem light
  !type(type_bulk_standard_variable) aos
  !ids for fabm
  !type(type_scalar_variable_id),save:: id_yearday !- ersem zenith
  type(type_bulk_variable_id),  save:: temp_id,salt_id,h_id
  type(type_bulk_variable_id),  save:: pres_id,rho_id,par_id
  type(type_horizontal_variable_id),save:: lon_id,lat_id,ws_id
  type(type_horizontal_variable_id),save:: taub_id,bdepth_id
  !type(type_horizontal_variable_id),save:: ssf_id !- ersem light
contains
  !
  !initialize brom
  !
  subroutine initialize_brom()
    real(rk),allocatable,dimension(:):: air_ice_indexes
    !NaN value
    REAL(rk), PARAMETER :: D_QNAN = &
              TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
    integer water_sediments_index
    integer i

    !initializing fabm from fabm.yaml file
    _LINE_
    call fabm_create_model_from_yaml_file(fabm_model)
    _LINE_
    !initializing brom standard_variables
    !makes grid, it starts from bottom (1) to surface (end point)
    standard_vars = brom_standard_variables()
    number_of_layers = standard_vars%get_value(&
                           "number_of_layers")
    !fabm grid
    call fabm_set_domain(fabm_model,number_of_layers)
    call fabm_model%set_surface_index(number_of_layers)
    call fabm_model%set_bottom_index(1)
    !linking state variables to fabm
    number_of_parameters = size(fabm_model%state_variables)
    allocate(state_vars(number_of_parameters))
    do i = 1,number_of_parameters
      call state_vars(i)%set_brom_state_variable(.false.,.false.,&
        _NEUMANN_,_NEUMANN_,0._rk,0._rk,0._rk,0._rk)
      state_vars(i)%name = fabm_model%state_variables(i)%name
      allocate(state_vars(i)%value(number_of_layers))
      state_vars(i)%value = fabm_model%state_variables(i)%initial_value
      !vertical movement rates (m/s, positive for upwards),
      !and set these to the values provided by the model.
      state_vars(i)%sinking_velocity = &
                    fabm_model%state_variables(i)%vertical_movement
      call fabm_link_bulk_state_data(fabm_model,i,state_vars(i)%value)
      call state_vars(i)%print_name()
    end do
    !in case of existing bottom and surface variables
    allocate(bcc(size(fabm_model%bottom_state_variables)))
    do i = 1,size(fabm_model%bottom_state_variables)
      bcc(i) = fabm_model%bottom_state_variables(i)%initial_value
      call fabm_link_bottom_state_data(fabm_model,i,bcc(i))
    end do
    allocate(scc(size(fabm_model%surface_state_variables)))
    do i = 1,size(fabm_model%surface_state_variables)
      scc(i) = fabm_model%surface_state_variables(i)%initial_value
      call fabm_link_surface_state_data(fabm_model,i,scc(i))
    end do
    _LINE_
    !initializing values
    call fabm_initialize_state(fabm_model,1,number_of_layers)
    allocate(air_ice_indexes,source=&
             standard_vars%get_column("air_ice_indexes"))
    allocate(porosity(number_of_layers))
    porosity=standard_vars%get_column("porosity",1)
    !recalculating concentrations in the sediments
    !concentrations are per total volume of sediment
    water_sediments_index = standard_vars%get_value("bbl_sediments_index")
    forall (i = 1:number_of_parameters)
      state_vars(i)%value(:water_sediments_index)&
        = state_vars(i)%value(:water_sediments_index)*&
          porosity(:water_sediments_index)
      state_vars(i)%value(air_ice_indexes(1):) = D_QNAN
    end forall
    !linking variables
    !yearday - ersem zenith_angle
    !id_yearday = fabm_model%get_scalar_variable_id(&
    !  standard_variables%number_of_days_since_start_of_the_year)
    !call fabm_model%link_scalar(id_yearday,realday)
    !cell thickness - ersem
    h_id    = fabm_model%get_bulk_variable_id(&
              standard_variables%cell_thickness)
    allocate(cell(number_of_layers))
    !cell = standard_vars%get_column("layer_thicknesses",1)
    call fabm_link_bulk_data(fabm_model,h_id,cell)
    !temperature
    temp_id = fabm_model%get_bulk_variable_id(&
              standard_variables%temperature)
    allocate(temp(number_of_layers))
    call fabm_link_bulk_data(fabm_model,temp_id,temp)
    !salinity
    salt_id = fabm_model%get_bulk_variable_id(&
              standard_variables%practical_salinity)
    allocate(salt(number_of_layers))
    call fabm_link_bulk_data(fabm_model,salt_id,salt)
    !density anomaly
    rho_id  = fabm_model%get_bulk_variable_id(standard_variables%density)
    allocate(density(number_of_layers))
    !density = standard_vars%get_column(_RHO_,1)
    call fabm_link_bulk_data(fabm_model,rho_id,density)
    !pressure - does fabm able to calculate it from density?
    pres_id = fabm_model%get_bulk_variable_id(&
              standard_variables%pressure)
    allocate(depth(number_of_layers))
    allocate(pressure(number_of_layers))
    depth = standard_vars%get_column("middle_layer_depths",1)
    !convert depth to pressure
    !total=water+atmosphere [dbar]
    pressure = depth+10._rk
    pressure(int(standard_vars%get_value("ice_water_index")):) = 10._rk
    call fabm_link_bulk_data(fabm_model,pres_id,pressure)
    !surface shortwave flux - ersem light
    !ssf_id  = fabm_model%get_horizontal_variable_id(&
    !          standard_variables%surface_downwelling_shortwave_flux)
    !allocate(radiative_flux)
    !call fabm_link_horizontal_data(fabm_model,ssf_id,radiative_flux)
    !photosynthetic_radiative_flux
    par_id  = fabm_model%get_bulk_variable_id(&
              standard_variables%downwelling_photosynthetic_radiative_flux)
    allocate(radiative_flux(number_of_layers))
    call fabm_link_bulk_data(fabm_model,par_id,radiative_flux)
    !longtitude
    lon_id  = fabm_model%get_horizontal_variable_id(&
              standard_variables%longitude)
    call fabm_link_horizontal_data(fabm_model,lon_id,_LONGITUDE_)
    !latitude
    lat_id  = fabm_model%get_horizontal_variable_id(&
              standard_variables%latitude)
    call fabm_link_horizontal_data(fabm_model,lat_id,_LATITUDE_)
    !wind speed
    ws_id   = fabm_model%get_horizontal_variable_id(&
              standard_variables%wind_speed)
    call fabm_link_horizontal_data(fabm_model,ws_id,5._rk)
    !bottom stress - ersem
    taub_id = fabm_model%get_horizontal_variable_id(&
              standard_variables%bottom_stress)
    allocate(taub)
    taub = 0._rk
    call fabm_link_horizontal_data(fabm_model,taub_id,taub)
    !bottom depth - fix to depth on boundary - ersem
    bdepth_id = fabm_model%get_horizontal_variable_id(&
                standard_variables%bottom_depth_below_geoid)
    bdepth = depth(1)
    call fabm_link_horizontal_data(fabm_model,bdepth_id,bdepth)
    !carbon dioxide in air
    call fabm_link_horizontal_data(fabm_model,&
         standard_variables%mole_fraction_of_carbon_dioxide_in_air,&
         380._rk)
    !absorption of silt - ersem light
    !aos = type_bulk_standard_variable(name="absorption_of_silt",units="1")
    !allocate(aos_value(number_of_layers))
    !aos_value = 4.e-5_rk
    !call fabm_link_bulk_data(fabm_model,aos,aos_value)

    !check all needed by fabm model variables
    call fabm_check_ready(fabm_model)
    !brom needs to know is variable a solid or gas
    !call configurate_state_variables()

    previous_ice_index=0
  end subroutine

  subroutine sarafan()
    use output_mod

    type(type_output):: netcdf_ice
    type(type_output):: netcdf_water
    type(type_output):: netcdf_sediments
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer ice_water_index,water_bbl_index,number_of_days
    integer surface_index
    integer day,i
    !ice thickness
    !real(rk) ice
    !cpu time
    real(rk) t1,t2
    real(rk),allocatable,dimension(:):: indices
    real(rk),allocatable,dimension(:):: air_ice_indexes

    allocate(indices(number_of_layers))
    indices = (/(i,i=number_of_layers,1,-1)/)

    ice_water_index = standard_vars%get_value ("ice_water_index")
    water_bbl_index = standard_vars%get_value ("water_bbl_index")
    number_of_days = standard_vars%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(number_of_days))
    air_ice_indexes = standard_vars%get_column("air_ice_indexes")

    netcdf_ice = type_output(fabm_model,_FILE_NAME_ICE_,&
                         ice_water_index,number_of_layers,&
                         number_of_layers)
    netcdf_water = type_output(fabm_model,_FILE_NAME_WATER_,&
                         water_bbl_index,ice_water_index-1,&
                         number_of_layers)
    netcdf_sediments = type_output(fabm_model,_FILE_NAME_SEDIMENTS_,&
                         1,water_bbl_index-1,&
                         number_of_layers)

    day = standard_vars%first_day()
    call initial_date(day,year)
    !call stabilize(day,year)

    do i = 1,number_of_days
      call date(day,year)

      !ice   = standard_vars%get_value(_ICE_THICKNESS_,i)
      !porosity = standard_vars%get_column("porosity",i)
      depth = standard_vars%get_column("middle_layer_depths",i)
      temp  = standard_vars%get_column(_TEMPERATURE_,i)
      salt  = standard_vars%get_column(_SALINITY_,i)

      !change surface index due to ice depth
      !index for boundaries so for layers it should be -1
      surface_index = air_ice_indexes(i)
      call fabm_model%set_surface_index(surface_index-1)

      !update links
      !realday = day !to convert integer to real - ersem zenith_angle
      !call fabm_model%link_scalar(id_yearday,realday)
      !cell thickness - ersem
      cell = standard_vars%get_column("layer_thicknesses",i)
      call fabm_link_bulk_data(fabm_model,h_id,cell)
      !temperature
      call fabm_link_bulk_data(&
           fabm_model,standard_variables%temperature,temp)
      !salinity
      call fabm_link_bulk_data(&
           fabm_model,standard_variables%practical_salinity,salt)
      !density
      density = standard_vars%get_column(_RHO_,i)+1000._rk
      call fabm_link_bulk_data(fabm_model,rho_id,density)
      !par
      call calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),&
        standard_vars%get_value(_SNOW_THICKNESS_,i),&
        standard_vars%get_value(_ICE_THICKNESS_ ,i))
      call fabm_link_bulk_data(fabm_model,par_id,radiative_flux)
      !bottom stress - ersem
      !call fabm_link_horizontal_data(fabm_model,taub_id,taub)

      !
      !have to change to flux on ice_water_index layer
      !
      !call find_set_state_variable("niva_brom_bio_PO4",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      !call find_set_state_variable("niva_brom_bio_NO3",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      !call find_set_state_variable("niva_brom_redox_Si",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,2._rk))

      call cpu_time(t1)
      call day_circle(i,surface_index)
      call netcdf_ice%save(fabm_model,state_vars,indices,i,&
                           temp,salt,depth,radiative_flux,&
                           int(air_ice_indexes(i)))
      call netcdf_water%save(fabm_model,state_vars,depth,i,&
                             temp,salt,depth,radiative_flux,&
                             int(air_ice_indexes(i)))
      call netcdf_sediments%save(fabm_model,state_vars,depth,i,&
                                 temp,salt,depth,radiative_flux,&
                                 int(air_ice_indexes(i)))
      call cpu_time(t2)

      write(*,*) "number / ","julianday / ","year",i,day,year
      write(*,*) "Time taken by day circle:",t2-t1," seconds"
      day = day+1
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
    !surface_flux in Watts,
    !to calculate it in micromoles photons per m2*s =>
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
    if (ice_thick==0._rk) then
      radiative_flux(:ice_water_index-1) = &
        surface_flux*exp(-_ERLOV_*depth(:ice_water_index-1))
    else
      radiative_flux(:ice_water_index-1) = &
        radiative_flux(ice_water_index)*&
        exp(-_ERLOV_*depth(:ice_water_index-1))
    end if
  end subroutine

  pure function sinusoidal(julian_day,multiplier)
    integer, intent(in):: julian_day
    real(rk),intent(in):: multiplier
    real(rk) sinusoidal

    sinusoidal = (1._rk+sin(2._rk*_PI_*(&
                  julian_day-40._rk)/365._rk))*multiplier
  end function

  subroutine day_circle(id,surface_index)
    integer,intent(in):: id !number of the count
    integer,intent(in):: surface_index

    real(rk),dimension(number_of_layers+1):: pF1_solutes
    real(rk),dimension(number_of_layers+1):: pF2_solutes
    real(rk),dimension(number_of_layers+1):: pF1_solids
    real(rk),dimension(number_of_layers+1):: pF2_solids
    real(rk),dimension(number_of_layers+1):: kz_mol
    real(rk),dimension(number_of_layers+1):: kz_bio
    real(rk),dimension(number_of_layers+1):: kz_turb
    real(rk),dimension(number_of_layers+1):: layer_thicknesses
    integer i,j,bbl_sed_index,ice_water_index
    integer number_of_circles
    !fabm logical parameters
    logical:: repair=.true.,valid
    !index for boundaries so for layers it should be -1
    real(rk),dimension(surface_index-1,number_of_parameters):: increment

    bbl_sed_index = standard_vars%get_value("bbl_sediments_index")
    ice_water_index = standard_vars%get_value("ice_water_index")
    pF1_solutes = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solutes_1",id) /)
    pF2_solutes = standard_vars%get_column("porosity_factor_solutes_2",id)
    pF1_solids = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solids_1",id) /)
    pF2_solids = standard_vars%get_column("porosity_factor_solids_2",id)
    kz_mol = standard_vars%get_column("molecular_diffusivity",id)
    kz_bio = standard_vars%get_column("bioturbation_diffusivity",id)
    kz_turb = standard_vars%get_column(_TURBULENCE_,id)
    layer_thicknesses = &
    (/ 0._rk,standard_vars%get_column("layer_thicknesses",id) /)

    if (mod(60*60*24,_SECONDS_PER_CIRCLE_)/=0) then
      call fatal_error("Check _SECONDS_PER_CIRCLE_",&
                       "should fit 86400/_SECONDS_PER_CIRCLE_=integer")
    else
      number_of_circles = int(60*60*24/_SECONDS_PER_CIRCLE_)
    end if
    call recalculate_ice(id)
    do i = 1,number_of_circles
      !diffusion
      call brom_do_diffusion(surface_index,bbl_sed_index,ice_water_index,&
                             pF1_solutes,pF2_solutes,pF1_solids,&
                             pF2_solids,kz_mol,kz_bio,kz_turb,&
                             layer_thicknesses)

      !biogeochemistry
      call fabm_check_state(fabm_model,1,surface_index-1,repair,valid)
      increment = 0._rk
      call fabm_do(fabm_model,1,surface_index-1,increment)
      increment = _SECONDS_PER_CIRCLE_*increment
      forall(j = 1:number_of_parameters)&
        state_vars(j)%value(:surface_index-1) = &
          state_vars(j)%value(:surface_index-1)+increment(:,j)
    end do
  end subroutine

  subroutine brom_do_diffusion(surface_index,bbl_sed_index,ice_water_index,&
                               pF1_solutes,pF2_solutes,pF1_solids,&
                               pF2_solids,kz_mol,kz_bio,kz_turb,&
                               layer_thicknesses)
    use diff_mod
    integer,intent(in):: surface_index
    integer,intent(in):: bbl_sed_index
    integer,intent(in):: ice_water_index
    real(rk),dimension(number_of_layers+1),intent(in):: pF1_solutes
    real(rk),dimension(number_of_layers+1),intent(in):: pF2_solutes
    real(rk),dimension(number_of_layers+1),intent(in):: pF1_solids
    real(rk),dimension(number_of_layers+1),intent(in):: pF2_solids
    real(rk),dimension(number_of_layers+1),intent(in):: kz_mol
    real(rk),dimension(number_of_layers+1),intent(in):: kz_bio
    real(rk),dimension(number_of_layers+1),intent(in):: kz_turb
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

    oxygen = find_state_variable("O2_o")
    !oxygen status of sediments
    O2stat = oxygen%value(bbl_sed_index)/&
      (oxygen%value(bbl_sed_index)+_KO2_)
    kz_tot = kz_turb+kz_mol+kz_bio*O2stat

    !calculate surface fluxes only for ice free periods
    surface_flux = 0._rk
    if (surface_index == ice_water_index) then
      call fabm_do_surface(fabm_model,surface_flux)
    end if
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
      temporary(:surface_index-1,i) = do_diffusive(&
          N       = surface_index-1,&
          dt      = _SECONDS_PER_CIRCLE_,&
          cnpar   = 0.6_rk,&
          posconc = 1,&
          h    = layer_thicknesses(:surface_index),&
          Bcup = state_vars(i)%use_bound_up,&
          Bcdw = state_vars(i)%use_bound_low,&
          Yup  = state_vars(i)%bound_up,&
          Ydw  = state_vars(i)%bound_low,&
          nuY_in  = kz_tot(:surface_index),&
          Lsour = zeros(:surface_index),&
          Qsour = zeros(:surface_index),&
          Taur  = taur_r(:surface_index),&
          Yobs  = zeros(:surface_index),&
          Y     = (/ 0._rk,state_vars(i)%value(:surface_index-1) /),&
          i_sed_top = bbl_sed_index-1,&
          is_solid = state_vars(i)%is_solid,&
          i_ice_water = ice_water_index,&
          is_gas = state_vars(i)%is_gas,&
          pF1_solutes = pF1_solutes(:surface_index),&
          pF2_solutes = pF2_solutes(:surface_index),&
          pF1_solids = pF1_solids(:surface_index),&
          pF2_solids = pF2_solids(:surface_index),&
          pFSWIup_solutes = pFSWIup_solutes,&
          pFSWIdw_solutes = pFSWIdw_solutes,&
          pFSWIup_solids = pFSWIup_solids,&
          pFSWIdw_solids = pFSWIdw_solids)
      state_vars(i)%value(1:surface_index-1) = temporary(1:surface_index-1,i)
    end forall
  end subroutine

  subroutine stabilize(inday,inyear)
    integer,intent(in):: inday,inyear

    type(brom_state_variable):: temporary_variable
    integer day,year
    integer pseudo_day,days_in_year,counter
    integer i
    integer surface_index
    real(rk),allocatable,dimension(:):: air_ice_indexes

    allocate(air_ice_indexes,source = &
             standard_vars%get_column("air_ice_indexes"))
    day = inday; year = inyear;
    days_in_year = 365+merge(1,0,(mod(year,4).eq.0))
    counter = days_in_year*5
    do i = 1,counter
      !day from 1 to 365 or 366
      if (mod(i,days_in_year).eq.0) then
        pseudo_day = days_in_year
      else
        pseudo_day = i-int(i/days_in_year)*&
                     days_in_year
      end if
      !ice   = standard_vars%get_value(_ICE_THICKNESS_,pseudo_day)
      !porosity = standard_vars%get_column("porosity",pseudo_day)
      depth = standard_vars%get_column("middle_layer_depths",pseudo_day)
      temp  = standard_vars%get_column(_TEMPERATURE_,pseudo_day)
      salt  = standard_vars%get_column(_SALINITY_,pseudo_day)

      surface_index = air_ice_indexes(pseudo_day)
      call date(day,year)
      !call calculate_radiative_flux(&
      !  surface_radiative_flux(_LATITUDE_,pseudo_day),&
      !  standard_vars%get_value(_SNOW_THICKNESS_,pseudo_day),&
      !  standard_vars%get_value(_ICE_THICKNESS_ ,pseudo_day))

      !call find_set_state_variable("niva_brom_bio_PO4",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,0.45_rk))
      !call find_set_state_variable("niva_brom_bio_NO3",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,3.8_rk))
      !call find_set_state_variable("niva_brom_redox_Si",&
      !  use_bound_up = _DIRICHLET_,bound_up = sinusoidal(day,2._rk))

      !change surface index due to ice depth
      !index for boundaries so for layers it should be -1
      call fabm_model%set_surface_index(surface_index-1)
      !call fabm_link_bulk_data(fabm_model,fabm_porosity,porosity)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%temperature,temp)
      call fabm_link_bulk_data(&
        fabm_model,standard_variables%practical_salinity,salt)
      !call fabm_link_bulk_data(&
      !  fabm_model,&
      !  standard_variables%downwelling_photosynthetic_radiative_flux,&
      !  radiative_flux)
      call day_circle(pseudo_day,surface_index)
      write(*,*) "Stabilizing initial array of values, in progress ..."
      write(*,*) "number / ","julianday / ","pseudo day",&
                 i,day,pseudo_day
      day = day+1
      !temporary_variable = find_state_variable("niva_brom_bio_O2")
      !call temporary_variable%print_state_variable()
    end do
  end subroutine

  subroutine recalculate_ice(id)
    integer,intent(in)   :: id
    real(rk),dimension(:),allocatable:: air_ice_indexes
    real(rk),dimension(:),allocatable:: layer_thicknesses
    !real(rk) ice_water_layer_thickness
    integer i,j
    integer ice_growth_save,ice_growth,ice_water_index

    allocate(air_ice_indexes,&
             source=standard_vars%get_column("air_ice_indexes"))
    allocate(layer_thicknesses,&
             source=standard_vars%get_column("layer_thicknesses",id))
    ice_water_index = standard_vars%get_value("ice_water_index")
    !ice_water_layer_thickness = layer_thicknesses(ice_water_index-1)
    !calculates number of layers freezed or melted
    if (previous_ice_index==0) then
      ice_growth_save = 0
    else
      ice_growth_save = int(air_ice_indexes(id))-previous_ice_index
    end if
    previous_ice_index = int(air_ice_indexes(id))

    do j = 1,number_of_parameters
      !recalculating
      ice_growth = ice_growth_save
      !melting
      if (ice_growth<0) then
        ice_growth = abs(ice_growth)
        do i=1,ice_growth
          state_vars(j)%value(ice_water_index-1) =&
            state_vars(j)%value(ice_water_index-1)+&
            state_vars(j)%value(ice_water_index-1+i)*&
            _ICE_LAYERS_RESOLUTION_/&
            layer_thicknesses(ice_water_index-1)
            !ice_water_layer_thickness
          state_vars(j)%value(ice_water_index-1+i) = 0._rk
        end do
        do i=ice_water_index,int(air_ice_indexes(id))+ice_growth-1
          state_vars(j)%value(i) = state_vars(j)%value(i+ice_growth)
        end do
      !freezing
      else if (ice_growth>0) then
        do i=int(air_ice_indexes(id)-1),ice_water_index+ice_growth,-1
          state_vars(j)%value(i)=state_vars(j)%value(i-ice_growth)
        end do
        do i=ice_water_index,ice_water_index+ice_growth-1
          state_vars(j)%value(i) = state_vars(j)%value(ice_water_index-1)
          state_vars(j)%value(ice_water_index-1) = &
            state_vars(j)%value(ice_water_index-1)-&
            state_vars(j)%value(ice_water_index-1)*&
            _ICE_LAYERS_RESOLUTION_/&
            layer_thicknesses(ice_water_index-1)
        end do
      end if
    end do
  end subroutine recalculate_ice

  subroutine find_set_state_variable(inname,is_solid,&
      is_gas,use_bound_up,use_bound_low,bound_up,&
      bound_low,density,sinking_velocity)
    character(len=*),                      intent(in):: inname
    logical,optional,                      intent(in):: is_solid
    logical,optional,                      intent(in):: is_gas
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
          is_gas,use_bound_up,use_bound_low,bound_up,&
          bound_low,density,sinking_velocity)
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end subroutine

  subroutine configurate_state_variables()
    !
    !have to change to flux on ice_water_index layer
    !
    !call find_set_state_variable("niva_brom_redox_SO4",&
    !  use_bound_up = _DIRICHLET_,use_bound_low = _DIRICHLET_,&
    !  bound_up = 25000._rk,bound_low = 25000._rk)
    !call find_set_state_variable(inname = "niva_brom_redox_Mn4",&
    !  use_bound_up = _DIRICHLET_,bound_up = 0.5e-4_rk)
    !call find_set_state_variable("niva_brom_redox_Fe3",&
    !  use_bound_up = _DIRICHLET_,bound_up = 0.4e-4_rk)
    !call find_set_state_variable("niva_brom_carb_Alk",&
    !  use_bound_up = _DIRICHLET_,bound_up = 2250._rk)

    call find_set_state_variable("niva_brom_bio_NH4"  ,is_gas = .true.)
    call find_set_state_variable("niva_brom_bio_O2"   ,is_gas = .true.)
    call find_set_state_variable("niva_brom_redox_H2S",is_gas = .true.)
    call find_set_state_variable("niva_brom_redox_CH4",is_gas = .true.)

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
    character(len=*),intent(in):: inname
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
