!-----------------------------------------------------------------------
! IPBM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the IPBM distribution.
!-----------------------------------------------------------------------
! Original author(s): Shamil Yakubov
!-----------------------------------------------------------------------

#include "../include/ipbm.h"
#include "../include/parameters.h"
!
!main module for transport calculations
!
module transport
  use fabm
  use fabm_config
  use fabm_driver
  use fabm_types
  use variables_mod,only: ipbm_standard_variables,&
                          ipbm_state_variable

  implicit none
  private
  public initialize_ipbm,sarafan

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
  type(ipbm_standard_variables) standard_vars
  type(ipbm_state_variable),allocatable,&
                       dimension(:),target:: state_vars
  !fabm model
  type(type_model) fabm_model
  !absorption_of_silt - ersem light
  !type(type_bulk_standard_variable) aos
  !ids for fabm
  !type(type_scalar_variable_id),save:: id_yearday !- ersem zenith
  type(type_bulk_variable_id),      save:: temp_id,salt_id,h_id
  type(type_bulk_variable_id),      save:: pres_id,rho_id,par_id
  type(type_horizontal_variable_id),save:: lon_id,lat_id,ws_id
  type(type_horizontal_variable_id),save:: taub_id,bdepth_id
  !type(type_horizontal_variable_id),save:: ssf_id !- ersem light
  
  interface do_relaxation
    module procedure do_relaxation_single
    module procedure do_relaxation_array
  end interface
  
contains
  !
  !initialize ipbm
  !
  subroutine initialize_ipbm()
    real(rk),allocatable,dimension(:):: air_ice_indexes
    real(rk),allocatable,dimension(:):: zeros
    !NaN value
    !REAL(rk), PARAMETER :: D_QNAN = &
    !          TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
    real(rk) D_QNAN

    integer water_sediments_index
    integer i

    !NaN
    D_QNAN = 0._rk
    D_QNAN = D_QNAN / D_QNAN

    !initializing fabm from fabm.yaml file
    _LINE_
    call fabm_create_model_from_yaml_file(fabm_model)
    _LINE_
    !
    !initializing ipbm standard_variables
    !makes grid, it starts from bottom (1) to surface (end point)
    standard_vars = ipbm_standard_variables()
    number_of_layers = standard_vars%get_value(&
                           "number_of_layers")
    !
    !fabm grid
    call fabm_set_domain(fabm_model,number_of_layers)
    call fabm_model%set_surface_index(number_of_layers)
    call fabm_model%set_bottom_index(1)
    !
    !linking state variables to fabm
    number_of_parameters = size(fabm_model%state_variables)
    allocate(state_vars(number_of_parameters))
    allocate(zeros(number_of_layers))
    zeros = 0._rk
    do i = 1,number_of_parameters
      state_vars(i)%name = fabm_model%state_variables(i)%name
      allocate(state_vars(i)%value(number_of_layers))
      !allocate(state_vars(i)%fabm_value(number_of_layers))
      allocate(state_vars(i)%sinking_velocity(number_of_layers))
      state_vars(i)%value = fabm_model%state_variables(i)%initial_value
      call state_vars(i)%set_ipbm_state_variable(.false.,.false.,&
        _NEUMANN_,_NEUMANN_,0._rk,0._rk,0._rk,zeros)
      !vertical movement rates (m/s, positive for upwards),
      !and set these to the values provided by the model.
      !state_vars(i)%sinking_velocity = &
      !              fabm_model%state_variables(i)%vertical_movement
      !if (fabm_model%state_variables(i)%vertical_movement/=0._rk) then
      !  state_vars(i)%is_solid=.true.
      !  state_vars(i)%density=1.5E7_rk
      !end if
      call fabm_link_bulk_state_data(fabm_model,i,state_vars(i)%value)
      call state_vars(i)%print_name()
    end do
    !
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
    !
    !initializing values
    !ipbm needs to know is variable a solid or gas
    call configurate_state_variables()
    call fabm_initialize_state(fabm_model,1,number_of_layers)
    allocate(air_ice_indexes,source=&
             standard_vars%get_column("air_ice_indexes"))
    allocate(porosity(number_of_layers))
    porosity=standard_vars%get_column("porosity",1)
    !recalculating concentrations in the sediments
    !concentrations are per total volume of sediment
    water_sediments_index = standard_vars%get_value("bbl_sediments_index")
    do i = 1,number_of_parameters
      if (state_vars(i)%is_solid.neqv..true.) then
        state_vars(i)%value(:air_ice_indexes(1)-1)&
          = state_vars(i)%value(:air_ice_indexes(1)-1)*&
            porosity(:air_ice_indexes(1)-1)
      else
        state_vars(i)%value(:water_sediments_index-1)&
          = state_vars(i)%value(:water_sediments_index-1)*&
            (1-porosity(:water_sediments_index-1))
      end if
      state_vars(i)%value(air_ice_indexes(1):) = D_QNAN
    end do
    !
    !linking variables
    !yearday - ersem zenith_angle
    !id_yearday = fabm_model%get_scalar_variable_id(&
    !  standard_variables%number_of_days_since_start_of_the_year)
    !call fabm_model%link_scalar(id_yearday,realday)
    !cell thickness - ersem
    h_id    = fabm_model%get_bulk_variable_id(&
              standard_variables%cell_thickness)
    allocate(cell(number_of_layers))
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
    !bottom stress - ersem - is needed by do_bottom
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

    previous_ice_index=0
  contains
    subroutine configurate_state_variables()
      !ammonium NH4+
      call find_set_state_variable(_NH4_,is_gas = .true.)
      !oxygen O2
      call find_set_state_variable(_O2_,is_gas = .true.)
      call find_set_state_variable(_H2S_,is_gas = .true.)
      call find_set_state_variable(_CH4_,is_gas = .true.)

      !calcite - CaCO3
      call find_set_state_variable(_CaCO3_,&
        is_solid = .true.,density = 2.80E7_rk)
      !S0
      call find_set_state_variable(_S0_,&
        is_solid = .true.,density = 6.56E7_rk)
      !Fe
      call find_set_state_variable(_Fe3_,&
        is_solid = .true.,density = 3.27E7_rk)
      call find_set_state_variable(_FeCO3_,&
        is_solid = .true.,density = 2.93E7_rk)
      call find_set_state_variable(_FeS_,&
        is_solid = .true.,density = 5.90E7_rk)
      call find_set_state_variable(_FeS2_,&
        is_solid = .true.,density = 4.17E7_rk)
      !Mn
      call find_set_state_variable(_Mn4_,&
        is_solid = .true.,density = 5.78E7_rk)
      call find_set_state_variable(_MnCO3_,&
        is_solid = .true.,density = 3.20E7_rk)
      call find_set_state_variable(_MnS_,&
        is_solid = .true.,density = 4.60E7_rk)
      !Silicon particulate
      call find_set_state_variable(_Sipart_,&
        is_solid = .true.,density = 4.40E7_rk)
      !organic compounds
      call find_set_state_variable(_Phy_,&
        is_solid = .true.,density = 1.5E7_rk)
      call find_set_state_variable(_PON_,&
        is_solid = .true.,density = 1.5E7_rk)
      !!small-size POM
      !call find_set_state_variable("R4_c",&
      !  is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      !call find_set_state_variable("R4_n",&
      !  is_solid = .true.,density = 1.5E7_rk)
      !call find_set_state_variable("R4_p",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("R4_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      !!medium-size POM
      !call find_set_state_variable("R6_c",&
      !  is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      !call find_set_state_variable("R6_n",&
      !  is_solid = .true.,density = 1.5E7_rk)
      !call find_set_state_variable("R6_p",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("R6_s",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/1000._rk)
      !call find_set_state_variable("R6_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      !!large-size POM
      !call find_set_state_variable("R8_c",&
      !  is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      !call find_set_state_variable("R8_n",&
      !  is_solid = .true.,density = 1.5E7_rk)
      !call find_set_state_variable("R8_p",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("R8_s",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/1000._rk)
      !diatoms
      call find_set_state_variable("P1_c",&
        is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      call find_set_state_variable("P1_n",&
        is_solid = .true.,density = 1.5E7_rk)
      call find_set_state_variable("P1_p",&
        is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("P1_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      call find_set_state_variable("P1_s",& !rios 1998
        is_solid = .true.,density = 1.5E7_rk*6.2_rk/15.7_rk)
      call find_set_state_variable("P1_Chl",&
        is_solid = .true.,density = 1.2E7_rk)!1200e6 from wiki
      !microzooplankton
      call find_set_state_variable("Z5_c",&
        is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      call find_set_state_variable("Z5_n",&
        is_solid = .true.,density = 1.5E7_rk)
      call find_set_state_variable("Z5_p",&
        is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !!nanophytoplankton
      !call find_set_state_variable("P2_c",&
      !  is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      !call find_set_state_variable("P2_n",&
      !  is_solid = .true.,density = 1.5E7_rk)
      !call find_set_state_variable("P2_p",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("P2_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      !call find_set_state_variable("P2_Chl",&
      !  is_solid = .true.,density = 1.2E7_rk)!1200e6 from wiki
      !!picophytoplankton
      !call find_set_state_variable("P3_c",&
      !  is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      !call find_set_state_variable("P3_n",&
      !  is_solid = .true.,density = 1.5E7_rk)
      !call find_set_state_variable("P3_p",&
      !  is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("P3_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      !call find_set_state_variable("P3_Chl",&
      !  is_solid = .true.,density = 1.2E7_rk)!1200e6 from wiki
      !microphytoplankton
      call find_set_state_variable("P4_c",&
        is_solid = .true.,density = 1.5E7_rk*106._rk/16._rk)
      call find_set_state_variable("P4_n",&
        is_solid = .true.,density = 1.5E7_rk)
      call find_set_state_variable("P4_p",&
        is_solid = .true.,density = 1.5E7_rk*1._rk/16._rk)
      !call find_set_state_variable("P4_f",&
      !  is_solid = .true.,density = 1.5E7_rk*5._rk/1260._rk)
      call find_set_state_variable("P4_Chl",&
        is_solid = .true.,density = 1.2E7_rk)!1200e6 from wiki
    end subroutine
  end subroutine
  !
  !procedure with the main cycle
  !
  subroutine sarafan()
    use output_mod

    type(type_output):: netcdf_ice
    type(type_output):: netcdf_water
    type(type_output):: netcdf_sediments
    integer:: year = _INITIALIZATION_SINCE_YEAR_
    integer ice_water_index,water_bbl_index,bbl_sediments_index
    integer number_of_days
    integer surface_index
    integer day,i
    !ice thickness
    !real(rk) ice
    !cpu time
    real(rk) t1,t2
    real(rk),allocatable,dimension(:):: depth_faces
    real(rk),allocatable,dimension(:):: indices,indices_faces
    real(rk),allocatable,dimension(:):: air_ice_indexes

    allocate(indices(number_of_layers))
    indices = (/(i,i=number_of_layers,1,-1)/)
    allocate(indices_faces(number_of_layers+1))
    indices_faces = (/(i,i=number_of_layers+1,1,-1)/)

    ice_water_index = standard_vars%get_value ("ice_water_index")
    water_bbl_index = standard_vars%get_value ("water_bbl_index")
    bbl_sediments_index = standard_vars%get_value ("bbl_sediments_index")
    number_of_days = standard_vars%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(number_of_days))
    air_ice_indexes = standard_vars%get_column("air_ice_indexes")

    day = standard_vars%first_day()
    call initial_date(day,year)
    !first day cycle
    call first_day_circle(100,ice_water_index,&
                          water_bbl_index,bbl_sediments_index,&
                          indices,indices_faces)
    !cycle first year 10 times
    call first_year_circle(day,year,ice_water_index,&
                           water_bbl_index,bbl_sediments_index,&
                           indices,indices_faces)

    netcdf_ice = type_output(fabm_model,standard_vars,_FILE_NAME_ICE_,&
                             ice_water_index,ice_water_index+20,&
                             number_of_layers)
    netcdf_water = type_output(fabm_model,standard_vars,_FILE_NAME_WATER_,&
                               bbl_sediments_index,ice_water_index-1,&
                               number_of_layers)
    netcdf_sediments = type_output(fabm_model,standard_vars,_FILE_NAME_SEDIMENTS_,&
                                   1,water_bbl_index-1,&
                                   number_of_layers)
    do i = 1,number_of_days
      call date(day,year)
      !ice   = standard_vars%get_value(_ICE_THICKNESS_,i)
      !porosity = standard_vars%get_column("porosity",i)

      !for netcdf output
      depth = standard_vars%get_column("middle_layer_depths",i)
      allocate(depth_faces,source=&
               standard_vars%get_column(_DEPTH_ON_BOUNDARY_,i))

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
      temp  = standard_vars%get_column(_TEMPERATURE_,i)
      call fabm_link_bulk_data(fabm_model,temp_id,temp)
      !salinity
      salt  = standard_vars%get_column(_SALINITY_,i)
      call fabm_link_bulk_data(fabm_model,salt_id,salt)
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

      call cpu_time(t1)
      call day_circle(i,surface_index,day)
      call netcdf_ice%save(fabm_model,standard_vars,state_vars,&
                           indices,indices_faces,i,&
                           int(air_ice_indexes(i)))
      call netcdf_water%save(fabm_model,standard_vars,state_vars,&
                             depth,depth_faces,i,&
                             int(air_ice_indexes(i)))
      call netcdf_sediments%save(fabm_model,standard_vars,state_vars,&
                                 depth,depth_faces,i,&
                                 int(air_ice_indexes(i)))
      call cpu_time(t2)

      write(*,*) "number / ","day / ","year",i,day,year
      write(*,*) "Time taken by day circle:",t2-t1," seconds"
      day = day+1
      deallocate(depth_faces)
    end do
    call netcdf_ice%close()
    call netcdf_water%close()
    call netcdf_sediments%close()
    write(*,*) "Finish"
    _LINE_
  contains
    subroutine initial_date(day,first_year)
      integer,intent(inout):: day
      integer,intent(inout):: first_year
      integer nydays

      nydays = 365+is_leap(first_year)
      do
        if (day<=nydays) exit
        nydays = 365+is_leap(first_year)
        day = day-nydays
        first_year = first_year+1
      end do
    end subroutine
  end subroutine
  !
  !first day iteration
  !
  subroutine first_day_circle(counter,ice_water_index,&
                              water_bbl_index,bbl_sediments_index,&
                              indices,indices_faces)
    use output_mod

    integer,intent(in):: counter
    integer,intent(in):: ice_water_index,water_bbl_index
    integer,intent(in):: bbl_sediments_index
    real(rk),allocatable,intent(in):: indices(:)
    real(rk),allocatable,intent(in):: indices_faces(:)

    type(type_output):: netcdf_ice
    type(type_output):: netcdf_water
    type(type_output):: netcdf_sediments
    integer i
    integer surface_index
    real(rk),allocatable,dimension(:):: air_ice_indexes
    real(rk),allocatable,dimension(:):: depth_faces

    allocate(air_ice_indexes,source = &
             standard_vars%get_column("air_ice_indexes"))
    !
    netcdf_ice = type_output(fabm_model,standard_vars,'ice_day.nc',&
                         ice_water_index,number_of_layers,&
                         number_of_layers)
    netcdf_water = type_output(fabm_model,standard_vars,'water_day.nc',&
                         water_bbl_index,ice_water_index-1,&
                         number_of_layers)
    netcdf_sediments = type_output(fabm_model,standard_vars,'sediments_day.nc',&
                         1,water_bbl_index-1,&
                         number_of_layers)
    !
    !change surface index due to ice depth
    !index for boundaries so for layers it should be -1
    surface_index = air_ice_indexes(1)
    call fabm_model%set_surface_index(surface_index-1)
    !
    !for netcdf output
    depth = standard_vars%get_column("middle_layer_depths",1)
    !
    !update links
    !realday = day !to convert integer to real - ersem zenith_angle
    !call fabm_model%link_scalar(id_yearday,realday)
    !cell thickness - ersem
    cell = standard_vars%get_column("layer_thicknesses",1)
    call fabm_link_bulk_data(fabm_model,h_id,cell)
    !temperature
    temp  = standard_vars%get_column(_TEMPERATURE_,1)
    call fabm_link_bulk_data(fabm_model,temp_id,temp)
    !salinity
    salt  = standard_vars%get_column(_SALINITY_,1)
    call fabm_link_bulk_data(fabm_model,salt_id,salt)
    !density
    density = standard_vars%get_column(_RHO_,1)+1000._rk
    call fabm_link_bulk_data(fabm_model,rho_id,density)
    !par
    call calculate_radiative_flux(&
      surface_radiative_flux(_LATITUDE_,14),& !kara case starts 14jan
      standard_vars%get_value(_SNOW_THICKNESS_,1),&
      standard_vars%get_value(_ICE_THICKNESS_ ,1))
    call fabm_link_bulk_data(fabm_model,par_id,radiative_flux)
    !
    do i = 1,counter
      call day_circle(1,surface_index,14)
      allocate(depth_faces,source=&
               standard_vars%get_column(_DEPTH_ON_BOUNDARY_,1))

      call netcdf_ice%save(fabm_model,standard_vars,state_vars,&
                           indices,indices_faces,i,&
                           int(air_ice_indexes(1)))
      call netcdf_water%save(fabm_model,standard_vars,state_vars,&
                             depth,depth_faces,i,&
                             int(air_ice_indexes(1)))
      call netcdf_sediments%save(fabm_model,standard_vars,state_vars,&
                                 depth,depth_faces,i,&
                                 int(air_ice_indexes(1)))

      write(*,*) "Stabilizing initial array of values, in progress ..."
      write(*,*) "number / ",i
      deallocate(depth_faces)
    end do
    call netcdf_ice%close()
    call netcdf_water%close()
    call netcdf_sediments%close()
  end subroutine
  !
  !first year iteration
  !
  subroutine first_year_circle(inday,inyear,ice_water_index,&
                               water_bbl_index,bbl_sediments_index,&
                               indices,indices_faces)
    use output_mod

    integer,intent(in):: inday,inyear
    integer,intent(in):: ice_water_index,water_bbl_index
    integer,intent(in):: bbl_sediments_index
    real(rk),allocatable,intent(in):: indices(:)
    real(rk),allocatable,intent(in):: indices_faces(:)

    type(type_output):: netcdf_ice
    type(type_output):: netcdf_water
    type(type_output):: netcdf_sediments
    integer day,year
    integer pseudo_day,days_in_year,counter
    integer i
    integer surface_index
    real(rk),allocatable,dimension(:):: air_ice_indexes
    real(rk),allocatable,dimension(:):: depth_faces

    allocate(air_ice_indexes,source = &
             standard_vars%get_column("air_ice_indexes"))
    day = inday; year = inyear;
    days_in_year = 365+merge(1,0,(mod(year,4).eq.0))
    counter = days_in_year*10

    netcdf_ice = type_output(fabm_model,standard_vars,'ice_year.nc',&
                         ice_water_index,ice_water_index+20,&
                         number_of_layers)
    netcdf_water = type_output(fabm_model,standard_vars,'water_year.nc',&
                         bbl_sediments_index,ice_water_index-1,&
                         number_of_layers)
    netcdf_sediments = type_output(fabm_model,standard_vars,'sediments_year.nc',&
                         1,water_bbl_index-1,&
                         number_of_layers)
    do i = 1,counter
      !day from 1 to 365 or 366
      if (mod(i,days_in_year).eq.0) then
        pseudo_day = days_in_year
      else
        pseudo_day = i-int(i/days_in_year)*&
                     days_in_year
      end if

      !for netcdf output
      depth = standard_vars%get_column("middle_layer_depths",pseudo_day)
      allocate(depth_faces,source=&
               standard_vars%get_column(_DEPTH_ON_BOUNDARY_,pseudo_day))

      call date(day,year)
      year = inyear
      !change surface index due to ice depth
      !index for boundaries so for layers it should be -1
      surface_index = air_ice_indexes(pseudo_day)
      call fabm_model%set_surface_index(surface_index-1)

      !update links
      !realday = day !to convert integer to real - ersem zenith_angle
      !call fabm_model%link_scalar(id_yearday,realday)
      !cell thickness - ersem
      cell = standard_vars%get_column("layer_thicknesses",pseudo_day)
      call fabm_link_bulk_data(fabm_model,h_id,cell)
      !temperature
      temp  = standard_vars%get_column(_TEMPERATURE_,pseudo_day)
      call fabm_link_bulk_data(fabm_model,temp_id,temp)
      !salinity
      salt  = standard_vars%get_column(_SALINITY_,pseudo_day)
      call fabm_link_bulk_data(fabm_model,salt_id,salt)
      !density
      density = standard_vars%get_column(_RHO_,pseudo_day)+1000._rk
      call fabm_link_bulk_data(fabm_model,rho_id,density)
      !par
      call calculate_radiative_flux(&
        surface_radiative_flux(_LATITUDE_,day),&
        standard_vars%get_value(_SNOW_THICKNESS_,pseudo_day),&
        standard_vars%get_value(_ICE_THICKNESS_ ,pseudo_day))
      call fabm_link_bulk_data(fabm_model,par_id,radiative_flux)
      !bottom stress - ersem
      !call fabm_link_horizontal_data(fabm_model,taub_id,taub)

      call day_circle(pseudo_day,surface_index,day)

      call netcdf_ice%save(fabm_model,standard_vars,state_vars,&
                           indices,indices_faces,pseudo_day,&
                           int(air_ice_indexes(pseudo_day)))
      call netcdf_water%save(fabm_model,standard_vars,state_vars,&
                             depth,depth_faces,pseudo_day,&
                             int(air_ice_indexes(pseudo_day)))
      call netcdf_sediments%save(fabm_model,standard_vars,state_vars,&
                                 depth,depth_faces,pseudo_day,&
                                 int(air_ice_indexes(pseudo_day)))

      write(*,*) "Stabilizing initial array of values, in progress ..."
      write(*,*) "number / ","day / ","pseudo day",&
                 i,day,pseudo_day
      day = day+1
      deallocate(depth_faces)
    end do
    call netcdf_ice%close()
    call netcdf_water%close()
    call netcdf_sediments%close()
  end subroutine
  !
  !keeps actual day and year
  !
  subroutine date(day,year)
    integer,intent(inout):: day
    integer,intent(inout):: year
    integer days

    days = 365+is_leap(year)
    if (day==days+1) then
      year = year+1
      day = day-days
    end if
  end subroutine
  !
  !return 1 (integer) in case of leap year
  !
  pure function is_leap(year)
    integer,intent(in):: year
    integer is_leap
    
    if (mod(year,400).eq.0) then
      is_leap = 1
      return
    else if (mod(year,100).eq.0) then
      is_leap = 0
      return
    else if (mod(year,4).eq.0) then
      is_leap = 1
      return
    else
      is_leap = 0
    end if
  end function is_leap
  !
  !returns surface PAR
  !
  pure function surface_radiative_flux(latitude,day)
    real(rk),intent(in):: latitude
    integer, intent(in):: day
    real(rk) surface_radiative_flux

    real(rk) Io,decl
    !Theoretical maximum 24-hr average surface downwelling
    !shortwave irradiance in air [W/m2] (default = 180 W/m2)
    !http://www.soda-pro.com
    !This should include that effect of average cloud cover (local)
    Io = 180._rk
    !Compute surface shortwave downwelling irradiance [W m^-2, 24-hr average]
    !Solar declination in degrees
    decl = 23.5_rk*sin(2.0_rk*_PI_*(real(day,rk)-81.0_rk)/365.0_rk)
    !This is the approximation used in Yakushev and Sorensen (2013) for OXYDEP
    surface_radiative_flux = max(0.0_rk, Io*cos((latitude-decl)*_PI_/180.0_rk))
    surface_radiative_flux = _PAR_PART_*surface_radiative_flux
  end function
  !
  !calculates PAR in the ice and water columns
  !
  subroutine calculate_radiative_flux(surface_flux,snow_thick,ice_thick)
    real(rk),intent(in):: surface_flux
    real(rk),intent(in):: snow_thick
    real(rk),intent(in):: ice_thick
    integer  ice_water_index
    real(rk) par_alb,par_scat
    real(rk),allocatable:: ice_depths(:)

    ice_water_index = standard_vars%get_value("ice_water_index")
    !ice_column
    !Grenfell and Maykutt 1977 indicate that the magnitude and shape
    !of the albedo curves depend strongly on the amount of liquid
    !water present in the upper part of the ice, so it fluctuates
    !throught year (true also for extinction coefficient)
    !if (snow_thick <= 0.005) then
      par_alb = surface_flux*(1._rk-_ICE_ALBEDO_)
    !else
    !  par_alb = surface_flux*(1._rk-_SNOW_ALBEDO_)*&
    !            exp(-_SNOW_EXTINCTION_*snow_thick)
    !end if
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
  !
  !calculate iterations within a day
  !
  subroutine day_circle(id,surface_index,day)
    integer,intent(in):: id !number of the count
    integer,intent(in):: surface_index
    integer,intent(in):: day !day

    real(rk),dimension(number_of_layers+1):: face_porosity
    real(rk),dimension(number_of_layers+1):: pF1_solutes
    real(rk),dimension(number_of_layers+1):: pF2_solutes
    real(rk),dimension(number_of_layers+1):: pF1_solids
    real(rk),dimension(number_of_layers+1):: pF2_solids
    real(rk),dimension(number_of_layers+1):: kz_mol
    real(rk),dimension(number_of_layers+1):: kz_bio
    real(rk),dimension(number_of_layers+1):: kz_turb
    real(rk),dimension(number_of_layers+1):: kz_ice_gravity
    real(rk),dimension(number_of_layers+1):: layer_thicknesses
    real(rk),dimension(number_of_layers)  :: porosity
    real(rk),dimension(number_of_layers-1):: dz
    !indices of layer interfaces in the sediments including the SWI
    real(rk),dimension(:),allocatable:: k_sed1
    real(rk),dimension(:),allocatable:: w_b
    real(rk),dimension(:),allocatable:: u_b
    real(rk) brine_release
    real(rk) dphidz_SWI
    integer i,j,bbl_sed_index,ice_water_index
    integer number_of_circles
    !fabm logical parameters
    logical:: repair=.true.,valid
    !index for boundaries so for layers it should be -1
    !increment for fabm do
    real(rk),dimension(surface_index-1,number_of_parameters):: increment
    !increment for diffusion
    real(rk),dimension(surface_index-1,number_of_parameters):: dcc

    bbl_sed_index   = standard_vars%get_value("bbl_sediments_index")
    ice_water_index = standard_vars%get_value("ice_water_index")
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume pore water] for solutes in sediments
    !first zero for gotm diff solver
    pF1_solutes = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solutes_1",id) /)
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    pF2_solutes = standard_vars%get_column("porosity_factor_solutes_2",id)
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume solids] for solids in sediments
    !first zero for gotm diff solver
    pF1_solids = &
    (/ 0._rk,standard_vars%get_column("porosity_factor_solids_1",id) /)
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    pF2_solids = standard_vars%get_column("porosity_factor_solids_2",id)
    kz_mol = standard_vars%get_column("molecular_diffusivity",id)
    kz_bio = standard_vars%get_column("bioturbation_diffusivity",id)
    kz_turb = standard_vars%get_column(_TURBULENCE_,id)
    kz_ice_gravity = standard_vars%get_column("ice_gravity_drainage",id)
    !first zero for gotm diff solver
    layer_thicknesses = &
    (/ 0._rk,standard_vars%get_column("layer_thicknesses",id) /)
    !
    !is needed by ipbm_do_sedimentation
    dphidz_SWI    = standard_vars%get_value("dphidz_SWI")
    !is needed by ipbm_do_sedimentation
    face_porosity = &
    standard_vars%get_column("porosity_on_interfaces",id)
    !is needed to constrain production in ice and seds
    porosity = &
    standard_vars%get_column("porosity",id)
    !indices of sediments
    allocate(k_sed1(bbl_sed_index))
    k_sed1 = standard_vars%get_column("k_sed1")
    !vertical velocities in the sediments, depends on porosity, constant
    allocate(w_b(bbl_sed_index))
    w_b = standard_vars%get_column("w_b")
    allocate(u_b(bbl_sed_index))
    u_b = standard_vars%get_column("u_b")
    !distance between centers of the layers
    dz = standard_vars%get_column("dz",id)

    if (mod(60*60*24,_SECONDS_PER_CIRCLE_)/=0) then
      call fatal_error("Check _SECONDS_PER_CIRCLE_",&
                       "should fit 86400/_SECONDS_PER_CIRCLE_=integer")
    else
      number_of_circles = int(60*60*24/_SECONDS_PER_CIRCLE_)
    end if
    call recalculate_ice(id,brine_release)
    do i = 1,number_of_circles
      !
      call relaxation(ice_water_index,bbl_sed_index,day)
      
      !diffusion
      !dcc = 0._rk
      call ipbm_do_diffusion(surface_index,bbl_sed_index,ice_water_index,&
                             pF1_solutes,pF2_solutes,pF1_solids,&
                             pF2_solids,kz_mol,kz_bio,kz_turb,kz_ice_gravity,&
                             layer_thicknesses,brine_release,dcc)
      !call check_array("after_diffusion",surface_index,id,i)
      call fabm_check_state(fabm_model,1,surface_index-1,repair,valid)
      
      !biogeochemistry
      increment = 0._rk
      do j = 1,number_of_parameters
        !if (state_vars(j)%is_solid.neqv..true. .and. &
        !    state_vars(j)%is_gas  .neqv..true.) then
        !for solutes:
        if (state_vars(j)%is_solid.neqv..true.) then
          !sedimants domain
          state_vars(j)%value(:bbl_sed_index-1)&
            = state_vars(j)%value(:bbl_sed_index-1)&
            / porosity(:bbl_sed_index-1)
          !ice domain
          state_vars(j)%value(ice_water_index:surface_index-1)&
            = state_vars(j)%value(ice_water_index:surface_index-1)&
            / porosity(ice_water_index:surface_index-1)
        !for solids:
        else
          !sediments domain
          state_vars(j)%value(:bbl_sed_index-1)&
            = state_vars(j)%value(:bbl_sed_index-1)&
            * pF1_solids(2:bbl_sed_index)
          !
          !ice domain          
          !all solids in the ice domain have similar
          !both brine and total volume concentrations
        end if
      end do
      !call check_array("after_concentrations_recalculation",surface_index,id,i)
      call fabm_do(fabm_model,1,surface_index-1,increment)
      do j = 1,number_of_parameters
        !if (state_vars(j)%is_solid.neqv..true. .and. &
        !    state_vars(j)%is_gas  .neqv..true.) then
        !for solutes:
        if (state_vars(j)%is_solid.neqv..true.) then
          !sediments_domain
          increment(:bbl_sed_index-1,j) &
            = _SECONDS_PER_CIRCLE_*increment(:bbl_sed_index-1,j)
          state_vars(j)%value(:bbl_sed_index-1)&
            = (state_vars(j)%value(:bbl_sed_index-1)&
            + increment(:bbl_sed_index-1,j))&
            * porosity(:bbl_sed_index-1)
          !water domain
          state_vars(j)%value(bbl_sed_index:ice_water_index-1)&
            = state_vars(j)%value(bbl_sed_index:ice_water_index-1)&
            + _SECONDS_PER_CIRCLE_&
            * increment(bbl_sed_index:ice_water_index-1,j)
          !ice_domain
          increment(ice_water_index:surface_index-1,j) &
            = _SECONDS_PER_CIRCLE_*increment(ice_water_index:surface_index-1,j)
          state_vars(j)%value(ice_water_index:surface_index-1)&
            = (state_vars(j)%value(ice_water_index:surface_index-1)&
            + increment(ice_water_index:surface_index-1,j))&
            * porosity(ice_water_index:surface_index-1)
        !for solids:
        else
          !sediments domain
          increment(:bbl_sed_index-1,j) &
            = _SECONDS_PER_CIRCLE_*increment(:bbl_sed_index-1,j)
          state_vars(j)%value(:bbl_sed_index-1)&
            = (state_vars(j)%value(:bbl_sed_index-1)&
            + increment(:bbl_sed_index-1,j))&
            / pF1_solids(2:bbl_sed_index)
          !water domain
          state_vars(j)%value(bbl_sed_index:ice_water_index-1)&
            = state_vars(j)%value(bbl_sed_index:ice_water_index-1)&
            + _SECONDS_PER_CIRCLE_&
            * increment(bbl_sed_index:ice_water_index-1,j)
          !ice_domain
          !porosity is needed to constrain production
          increment(ice_water_index:surface_index-1,j)&
            = _SECONDS_PER_CIRCLE_&
            * increment(ice_water_index:surface_index-1,j)&
            * porosity(ice_water_index:surface_index-1)
          state_vars(j)%value(ice_water_index:surface_index-1)&
            = state_vars(j)%value(ice_water_index:surface_index-1)&
            + increment(ice_water_index:surface_index-1,j)
          
        end if
      end do
      !do j = 1,number_of_parameters
      !    state_vars(j)%value(:surface_index-1) = &
      !      state_vars(j)%value(:surface_index-1)+&
      !      _SECONDS_PER_CIRCLE_*increment(:,j)
      !end do
      !call check_array("after_fabm_do",surface_index,id,i)
      call fabm_check_state(fabm_model,1,surface_index-1,repair,valid)
      
      !sedimentation
      call ipbm_do_sedimentation(surface_index,bbl_sed_index,&
                                 ice_water_index,k_sed1,w_b,u_b,&
                                 dphidz_SWI,&
                                 increment,&
                                 face_porosity(:surface_index),&
                                 kz_bio(:surface_index),&
                                 layer_thicknesses(2:surface_index),&
                                 dz(:surface_index-2))
      !call check_array("after_sedimentation",surface_index,id,i)
      call fabm_check_state(fabm_model,1,surface_index-1,repair,valid)
    end do
  end subroutine
  !
  !diffusion part
  !
  subroutine ipbm_do_diffusion(&
             surface_index,bbl_sed_index,ice_water_index,&
             pF1_solutes,pF2_solutes,pF1_solids,&
             pF2_solids,kz_mol,kz_bio,kz_turb,kz_ice_gravity,&
             layer_thicknesses,brine_release,increment)
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
    real(rk),dimension(number_of_layers+1),intent(in):: kz_ice_gravity
    real(rk),dimension(number_of_layers+1),intent(in):: layer_thicknesses
    real(rk)                              ,intent(in):: brine_release

    real(rk),dimension(surface_index-1,number_of_parameters),&
                                           intent(out):: increment

    type(ipbm_state_variable):: oxygen
    real(rk),dimension(number_of_layers+1):: ones
    real(rk),dimension(number_of_layers+1):: zeros
    real(rk),dimension(number_of_layers+1):: taur_r
    real(rk),dimension(number_of_layers+1):: brine_flux
    real(rk),dimension(number_of_parameters):: surface_flux
    real(rk),dimension(number_of_layers+1,&
                         number_of_parameters):: kz_tot
    real(rk),dimension(0:number_of_layers,&
                         number_of_parameters):: temporary
    integer i
    real(rk) O2stat
    real(rk) pFSWIup_solutes, pFSWIdw_solutes
    real(rk) pFSWIup_solids , pFSWIdw_solids

    ones  = 1._rk
    zeros = 0._rk
    taur_r= 1.e20_rk

    !calculating kz in the whole column
    oxygen = find_state_variable(_O2_)
    !oxygen status of sediments
    O2stat = oxygen%value(bbl_sed_index)/&
      (oxygen%value(bbl_sed_index)+_KO2_)
    brine_flux = 0._rk
    brine_flux(ice_water_index) = brine_release
    !add turb to ice water boundary
    if (surface_index/=ice_water_index) then
      brine_flux(ice_water_index) = &
        brine_flux(ice_water_index)+1.e-5_rk
    end if
    do i = 1,number_of_parameters
      kz_tot(:,i) = brine_flux+kz_ice_gravity+&
                    kz_turb+kz_mol+kz_bio*O2stat
    end do

    !calculate surface fluxes only for ice free periods
    surface_flux = 0._rk
    if (surface_index == ice_water_index) then
      call fabm_do_surface(fabm_model,surface_flux)
    end if
    !
    do i = 1,number_of_parameters
      call state_vars(i)%set_ipbm_state_variable(&
        use_bound_up = _NEUMANN_,bound_up = surface_flux(i))
    end do

    !
    !adapted from Phil Wallhead (PW):
    !solutes:
    !fick_SWI = -phi_0*Km/dz*(C_1/phi_1 - C_-1) - Kb/dz*(C_1 - C_-1)
    !         = -phi_0*(Km+Kb)/dz * (p_1*C_1 - p_-1*C_-1)
    !         [intraphase molecular + interphase bioturb]
    !where:
    !       p_-1 = (phi_0*Km + Kb) / (phi_0*(Km+Kb))
    !       p_1  = (phi_0/phi_1*Km + Kb) / (phi_0*(Km+Kb))
    !and subscripts refer to SWI (0), below (1), and above (-1)
    !
    !top cell of sediments
    pFSWIup_solutes = (pF1_solutes(bbl_sed_index+1)*&
      pF2_solutes(bbl_sed_index)*kz_mol(bbl_sed_index)+&
      kz_bio(bbl_sed_index)*O2stat)/(pF2_solutes(bbl_sed_index)*&
      (kz_mol(bbl_sed_index)+kz_bio(bbl_sed_index)*O2stat))
    !bottom cell of water column
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
    pFSWIup_solids = 1.0_rk/(1.0_rk-pF2_solids(bbl_sed_index))
    pFSWIdw_solids = pFSWIup_solids

    !forall (i = 1:number_of_parameters)
    do i = 1,number_of_parameters
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
          nuY_in  = kz_tot(:surface_index,i),&
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
      increment(:surface_index-1,i) = temporary(1:surface_index-1,i)-&
                                      state_vars(i)%value(:surface_index-1)/&
                                      _SECONDS_PER_CIRCLE_
      state_vars(i)%value(1:surface_index-1) = &
                          temporary(1:surface_index-1,i)
    !end forall
    end do
  end subroutine
  !
  !adapted from Phil Wallhead code
  !calculates vertical advection (sedimentation)
  !in the water column and sediments
  !
  subroutine ipbm_do_sedimentation(surface_index,bbl_sed_index,&
                                   ice_water_index,k_sed1,w_b,u_b,&
                                   dphidz_SWI,increment,&
                                   face_porosity,kz_bio,&
                                   hz,dz)
    integer ,intent(in):: surface_index
    integer ,intent(in):: bbl_sed_index
    integer ,intent(in):: ice_water_index

    real(rk),dimension(bbl_sed_index),intent(in):: k_sed1
    real(rk),dimension(bbl_sed_index),intent(in):: w_b,u_b

    real(rk)                         ,intent(in):: dphidz_SWI
    real(rk),dimension(surface_index-1,number_of_parameters),&
                                      intent(in):: increment
    real(rk),dimension(surface_index),intent(in):: face_porosity
    real(rk),dimension(surface_index),intent(in):: kz_bio
    ! layer thickness (m)
    real(rk),intent(in):: hz(surface_index-1)
    ! distance between centers of the layers
    real(rk),intent(in):: dz(surface_index-2)

    !Local variables
    type(ipbm_state_variable):: oxygen

    !NaN value
    REAL(rk), PARAMETER :: D_QNAN = &
              TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)

    real(rk):: dcc (surface_index-1,number_of_parameters)
    real(rk):: wbi (surface_index-1,number_of_parameters)
    real(rk):: wti (surface_index  ,number_of_parameters)
    real(rk):: sink(surface_index  ,number_of_parameters)
    real(rk):: w_1m(surface_index  ,number_of_parameters)
    real(rk):: w_1 (surface_index)
    real(rk):: u_1 (surface_index)
    real(rk):: w_1c(surface_index)
    real(rk):: u_1c(surface_index)

    real(rk) O2stat
    integer i,k,ip

    dcc  = 0.0_rk
    wbi  = 0.0_rk
    wti  = 0.0_rk
    w_1m = 0.0_rk
    w_1  = 0.0_rk
    u_1  = 0.0_rk
    w_1c = 0.0_rk
    u_1c = 0.0_rk

    !Compute vertical velocity in water column (sinking/floating)
    !using the FABM.
    call fabm_get_vertical_movement(&
         fabm_model,1,surface_index-1,wbi)
    do i=1,number_of_parameters
      !Note: sinking velocity is on layer midpoints
      !call fabm_get_vertical_movement(&
      !     fabm_model,1,number_of_layers,&
      !     state_vars(i)%sinking_velocity(:))
      state_vars(i)%sinking_velocity(:) = D_QNAN
      state_vars(i)%sinking_velocity(:surface_index-1) = wbi(:,i)
    end do

    !Phil Wallhead comments:
    !Compute vertical velocity components due to modelled (reactive)
    !particles, if required
    !
    !Assuming no externally impressed porewater flow, the equations
    !for liquid and solid volume fractions are:
    !
    !dphi/dt = -d/dz(phi*u - Dbip*dphi/dz) - sum_i Rp_i/rhop_i (1)
    !d(1-phi)/dt =
    ! -d/dz((1-phi)*w - Dbip*d/dz(1-phi)) + sum_i Rp_i/rhop_i  (2)
    !
    !where u is the solute advective velocity,
    !      w is the particulate advective velocity,
    !Dbip is the interphase component of bioturbation diffusivity
    !      ( = Db at SWI, 0 elsewhere)
    !Rp_i is the net reaction term for the i^th
    !      particulate substance [mmol/m3/s]
    !rhop_i is the molar density of the i^th
    !      particulate substance [mmol/m3]
    !
    !(1)+(2) => phi*u + (1-phi)*w = const      (3)
    !
    !Given u = w at some (possibly infinite) depth where compaction
    !      ceases and phi = phi_inf, w = w_inf:
    !phi*u + (1-phi)*w = phi_inf*w_inf + (1-phi_inf)*w_inf
    !        => phi*u = w_inf - (1-phi)*w      (4)
    !
    !We calculate w(z) by assuming steady state compaction (dphi/dt = 0)
    !and matching the solid volume flux across the SWI with the
    !(approximated) sinking flux of suspended particles in the fluff layer:
    !
    !(2) -> (1-phi)*w + Dbip*dphi/dz =
    !       (1-phi_0)*w_0 + Dbip*dphi/dz_0 +
    !       sum_i (1/rhop_i)*int_z0^z Rp_i(z') dz' (5)
    !
    !Matching the volume flux across the SWI gives:
    !       (1-phi_0)*w_0 + Dbip*dphi/dz_0 =
    !       F_b0 + sum_i (1/rhop_i)*wbio_f(i)*Cp_sf(i)
    !
    !where F_b0 is the background volume flux [m/s]
    !      wbio_f(i) is the sinking speed in the fluff layer [m/s]
    !      Cp_sf(i) is the suspended particle concentration
    !      in the fluff layer [mmol/m3]
    !      (approximated by the minimum of the particle concentration in
    !       the fluff layer and layer above)
    !
    !For dynamic_w_sed = 0, we set the reactions and modelled volume
    !      fluxes in (5) to zero:
    !
    !(1-phi)*w + Dbip*dphi/dz = F_b0 = (1-phi_inf)*w_binf
    !
    !Then using (4):
    !
    !phi*u = w_inf - (1-phi)*w = phi_inf*w_binf + Dbip*dphi/dz
    !
    !Decomposing the velocities as (w,u) = (w,u)_b + (w,u)_1, we get:
    !
    !(1-phi)*w_1 = -Dbip*dphi/dz
    !phi*u_1     = Dbip*dphi/dz
    !
    !where (1-phi)*w_b = (1-phi_inf)*w_binf and phi*u_b = phi_inf*w_binf
    !
    !For dynamic_w_sed = 1, we have further corrections (w,u)_1c, where:
    !
    !(1-phi)*w_1c =
    !   sum_i [ (1/rhop_i)*(wbio_f(i)*Cp_sf(i) + int_zTF^z Rp_i(z') dz') ]
    !
    !and using (4):
    !
    !phi*u_1c     = w_1cinf - (1-phi)*w_1c
    !
    !where w_1cinf can be approximated by the deepest value of w_1c
    !

    oxygen = find_state_variable(_O2_)
    !oxygen status of sediments set by O2
    !level just above sediment surface
    O2stat = oxygen%value(bbl_sed_index)/&
      (oxygen%value(bbl_sed_index)+_KO2_)

    w_1(bbl_sed_index) = -1.0_rk*O2stat*kz_bio(bbl_sed_index)*&
                       dphidz_SWI/(1.0_rk-face_porosity(bbl_sed_index))
    u_1(bbl_sed_index) = O2stat*kz_bio(bbl_sed_index)*dphidz_SWI/&
                       face_porosity(bbl_sed_index)

    !!Sum over contributions from each particulate variable
    !do i=1,number_of_parameters
    !  if (state_vars(i)%is_solid) then
    !    !First set rhop_i*(1-phi)*w_1i at the SWI
    !    w_1m(bbl_sed_index,i) = &
    !      state_vars(i)%sinking_velocity(bbl_sed_index)*&
    !      min(state_vars(i)%value(bbl_sed_index),&
    !      state_vars(i)%value(bbl_sed_index+1))
    !      !Now set rhop_i*(1-phi)*w_1i in the sediments
    !      !  by integrating the reaction terms
    !    do k=bbl_sed_index-1,1,-1
    !      w_1m(k,i) = w_1m(k+1,i)+increment(k,i)*hz(k)
    !    end do
    !    !Divide by rhop_i*(1-phi) to get w_1c(i)
    !    w_1m(k_sed1,i) = w_1m(k_sed1,i)/&
    !      (state_vars(i)%density*(1.0_rk-face_porosity(k_sed1)))
    !    !Add to total w_1c
    !    w_1c(k_sed1) = w_1c(k_sed1) + w_1m(k_sed1,i)
    !  end if
    !end do
    !!Now calculate u from w using (4) above
    !u_1c(k_sed1) = (w_1c(1)-&
    !               (1.0_rk-face_porosity(k_sed1))*&
    !               w_1(k_sed1))/face_porosity(k_sed1)

    !find phy to specify its sinking behavior
    !i = find_index_of_state_variable(_Phy_)
    !Interpolate velocities from FABM (defined on layer midpoints, as for
    !  concentrations) to wti on the layer interfaces
    do ip=1,number_of_parameters
      !Air-sea interface (unused)
      wti(surface_index,ip) = &
        state_vars(ip)%sinking_velocity(surface_index-1)
      !Water column layer interfaces, not including SWI
      wti(bbl_sed_index+1:surface_index-1,ip) = &
        state_vars(ip)%sinking_velocity(bbl_sed_index:surface_index-2)+&
        0.5_rk*hz(bbl_sed_index:surface_index-2)*&
        (state_vars(ip)%sinking_velocity(bbl_sed_index+1:surface_index-1)-&
        state_vars(ip)%sinking_velocity(bbl_sed_index:surface_index-2))/&
        dz(bbl_sed_index:surface_index-2)
      !Sediment layer interfaces, including SWI
      if (state_vars(ip)%is_solid) then
        wti(k_sed1,ip) = w_b(k_sed1)+w_1(k_sed1)+w_1c(k_sed1)
      else
        wti(k_sed1,ip) = u_b(k_sed1)+u_1(k_sed1)+u_1c(k_sed1)
      end if
      wti(1,ip) = wti(2,ip)
      !special vertical sedimentation for diatoms in the ice
      if (ip==find_index_of_state_variable('P1_c').or.&
          ip==find_index_of_state_variable('P1_n').or.&
          ip==find_index_of_state_variable('P1_p').or.&
          ip==find_index_of_state_variable('P1_s').or.&
          ip==find_index_of_state_variable('P1_Chl')&
          ) then
        !set velocity 3 cm/day
        wti(ice_water_index+1:surface_index-1,ip) = -0.03_rk/86400._rk
        wti(ice_water_index,ip) = 0._rk
      end if
    end do

    !Perform advective flux calculation and concentrations update
    !This uses a simple first order upwind differencing scheme (FUDM)
    !It uses the fluxes sink in a consistent manner and therefore
    !conserves mass
    !Calculate sinking fluxes at layer interfaces
    !(sink, units strictly [mass/unit total area/second])
    !Air-sea interface
    sink = 0.0_rk
    !Water column and sediment layer interfaces
    do i = 1,number_of_parameters
      do k = 1,surface_index-1
        sink(k,i) = wti(k,i)*state_vars(i)%value(k)
        !This is an upwind differencing approx., hence the use of cc(k-1)
        !Note: factors phi, (1-phi) are not needed in the sediments
        !because cc is in units [mass per unit total volume]
      end do
    end do
    !Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints
    !(top/bottom not used where Dirichlet bc imposed)
    do k=1,surface_index-1
      dcc(k,:) = -1.0_rk * (sink(k+1,:)-sink(k,:)) / hz(k)
    end do

    !Time integration
    do i = 1,number_of_parameters
      do k = 1,surface_index-1
          state_vars(i)%value(k) = state_vars(i)%value(k)+&
                                   _SECONDS_PER_CIRCLE_*dcc(k,i)
      end do
    end do
    !call state_vars(1)%print_state_variable()
  end subroutine ipbm_do_sedimentation
  !
  !recalculates state variables concentrations due to freezing or melting
  !
  subroutine recalculate_ice(id,brine_release)
    integer ,intent(in) :: id
    real(rk),intent(out):: brine_release
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
        if ((j==find_index_of_state_variable('P1_c').or.&
            j==find_index_of_state_variable('P1_n').or.&
            j==find_index_of_state_variable('P1_p').or.&
            j==find_index_of_state_variable('P1_s').or.&
            j==find_index_of_state_variable('P1_Chl'))&
            .and.air_ice_indexes(id)/=ice_water_index) then
          do i=1,ice_growth
            state_vars(j)%value(ice_water_index) =&
              state_vars(j)%value(ice_water_index)+&
              state_vars(j)%value(ice_water_index+i)
              !ice_water_layer_thickness
            state_vars(j)%value(ice_water_index+i) = 0._rk
          end do
          do i=ice_water_index+1,int(air_ice_indexes(id))+ice_growth-1
            state_vars(j)%value(i) = state_vars(j)%value(i+ice_growth)
          end do
          cycle
        end if
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
    if (ice_growth_save>0) then
      brine_release = calc_brine_release(ice_growth)
    else
      brine_release = 0._rk
    end if
  contains
    !
    !adapted from Arrigo 1993, Duarte 2007
    !
    real(rk) function calc_brine_release(ice_growth)
      integer, intent(in):: ice_growth

      real(rk) growth_rate

      growth_rate = ice_growth*_ICE_LAYERS_RESOLUTION_
      growth_rate = growth_rate*100._rk !m to cm
      growth_rate = growth_rate/86400._rk !per day to per second
      calc_brine_release = (9.667e-9_rk+4.49e-6_rk*growth_rate &
                                       -1.39e-7_rk*growth_rate**2)&
                                       *_ICE_LAYERS_RESOLUTION_
      calc_brine_release = calc_brine_release/100._rk
    end function
  end subroutine recalculate_ice
  !
  !sets state variable value and other parameters by name
  !
  subroutine find_set_state_variable(inname,is_solid,&
      is_gas,use_bound_up,use_bound_low,bound_up,&
      bound_low,density,sinking_velocity,value,layer)
    character(len=*),                      intent(in):: inname
    logical,optional,                      intent(in):: is_solid
    logical,optional,                      intent(in):: is_gas
    integer,optional,                      intent(in):: use_bound_up
    integer,optional,                      intent(in):: use_bound_low
    real(rk),optional,                     intent(in):: bound_up
    real(rk),optional,                     intent(in):: bound_low
    real(rk),optional,                     intent(in):: density
    real(rk),allocatable,dimension(:),optional,intent(in)&
                                                     :: sinking_velocity
    real(rk),optional,                         intent(in):: value
    integer ,optional,                         intent(in):: layer

    integer number_of_vars
    integer i

    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq.inname) then
        call state_vars(i)%set_ipbm_state_variable(is_solid,&
          is_gas,use_bound_up,use_bound_low,bound_up,&
          bound_low,density,sinking_velocity,value,layer)
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end subroutine
  !
  !
  !
  subroutine relaxation(ice_water_index,bbl_sed_index,day)
    integer,intent(in):: ice_water_index,bbl_sed_index
    integer,intent(in):: day

    integer number_of_vars
    integer i
    
    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq._DIC_) then
        call do_relaxation(1930._rk,ice_water_index-1,i)
        call do_relaxation(2280._rk,bbl_sed_index,i)
      else if (state_vars(i)%name.eq._Alk_) then
        call do_relaxation(2000._rk,ice_water_index-1,i)
        call do_relaxation(2350._rk,bbl_sed_index,i)
      !else if (state_vars(i)%name.eq._NH4_) then
      !  call do_relaxation(1._rk,bbl_sed_index-1,i)
      else if (state_vars(i)%name.eq._PO4_) then
        call do_relaxation(sinusoidal(day,0.5_rk),ice_water_index-1,i)
      else if (state_vars(i)%name.eq._NO3_) then
        call do_relaxation(sinusoidal(day,0.5_rk),ice_water_index-1,i)
      else if (state_vars(i)%name.eq._Si_) then
        call do_relaxation(sinusoidal(day,1.0_rk),ice_water_index-1,i)
      end if
    end do
  contains
    pure function sinusoidal(day,multiplier)
      integer, intent(in):: day
      real(rk),intent(in):: multiplier
      real(rk) sinusoidal

      sinusoidal = (1._rk+sin(2._rk*_PI_*(&
                    day-40._rk)/365._rk))*multiplier
    end function sinusoidal
  end subroutine
  !
  !
  !
  subroutine do_relaxation_single(value,index,i)
    real(rk),intent(in):: value
    integer ,intent(in):: index
    integer, intent(in):: i
    
    real(rk) dcc
      
    dcc = _HMIX_RATE_ &
      * (value-state_vars(i)%value(index))
    state_vars(i)%value(index) = state_vars(i)%value(index) &
      + _SECONDS_PER_CIRCLE_*dcc

  end subroutine do_relaxation_single
  !
  !indexes from bottom upwards
  !
  subroutine do_relaxation_array(from,till,i,value)
    integer ,intent(in):: from
    integer ,intent(in):: till
    integer, intent(in):: i
    real(rk),dimension(from:till),intent(in):: value
    
    integer j
    real(rk) dcc
      
    do j = from,till
      dcc = _HMIX_RATE_ &
        * (value(j)-state_vars(i)%value(j))
      state_vars(i)%value(j) = state_vars(i)%value(j) &
        + _SECONDS_PER_CIRCLE_*dcc
    end do
  end subroutine do_relaxation_array
  !
  !returns type ipbm_state_variable by name
  !
  function find_state_variable(inname)
    character(len=*),intent(in):: inname
    type(ipbm_state_variable),target:: find_state_variable
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
  !
  !returns index of variable
  !
  function find_index_of_state_variable(inname)
    character(len=*),intent(in):: inname
    integer:: find_index_of_state_variable
    integer number_of_vars
    integer i

    number_of_vars = size(state_vars)
    do i = 1,number_of_vars
      if (state_vars(i)%name.eq.inname) then
        find_index_of_state_variable = i
        return
      end if
    end do
    call fatal_error("Search state variable",&
                     "No such variable")
  end function
  !
  !checking for negative and NaNs
  !
  subroutine check_array(location,surface_index,day,i)
    character(len=*),intent(in):: location
    integer         ,intent(in):: surface_index
    integer         ,intent(in):: day
    integer         ,intent(in):: i

    character(10)::sday,si
    integer j

    write(si,'(i10)') i
    write(sday,'(i10)') day

    do j = 1,number_of_parameters
      if (any(state_vars(j)%value(1:surface_index-1)<0._rk)) then
        _LINE_
        write(*,*) "value=", state_vars(j)%value(1:surface_index-1)
        _LINE_
        call fatal_error("Negative value: ",&
                         location//": day="//trim(sday)//&
                         " iteration="//trim(si)//&
                         " name="//state_vars(j)%name)
      end if
      if (any(isnan(state_vars(j)%value(1:surface_index-1)))) then
        call fatal_error("NaN: ",&
                         location//": "//state_vars(j)%name)
      end if
    end do
  end subroutine check_array
end module

program main
  use transport

  call initialize_ipbm()
  call sarafan()
end program
