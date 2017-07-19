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

module variables_mod
  use types_mod
  use input_mod
  use ice_mod
  use fabm_driver

  implicit none
  !NaN value
  !REAL(rk), PARAMETER :: D_QNAN = &
  !          TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
  real(rk) D_QNAN

  type,extends(list_variables):: ipbm_standard_variables
    type(ice) type_ice
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
    procedure:: add_ice_thickness
    procedure:: add_grid_on_faces
    procedure:: add_grid_on_centers
    procedure:: add_day_number
    procedure:: add_layer_thicknesses
    procedure:: add_constant_in_sed
    procedure:: add_porosity
    procedure:: add_diffusivity
    procedure,public:: first_day
  end type

  type,extends(variable_1d):: ipbm_state_variable
    logical  is_solid
    logical  is_gas
    integer  use_bound_up
    integer  use_bound_low
    real(rk) bound_up
    real(rk) bound_low
    real(rk) density
    real(rk),allocatable,dimension(:):: sinking_velocity
    !real(rk),allocatable,dimension(:):: fabm_value
  contains
    procedure:: set_ipbm_state_variable
    procedure:: print_state_variable
  end type

  interface ipbm_standard_variables
    module procedure ipbm_standard_variables_constructor
  end interface
contains
  function ipbm_standard_variables_constructor()
    type(ipbm_standard_variables):: ipbm_standard_variables_constructor

    !NaN
    D_QNAN = 0._rk
    D_QNAN = D_QNAN / D_QNAN

    call ipbm_standard_variables_constructor%initialize()
  end function
  !
  !Initialize standard variables list
  !
  subroutine initialize_standard_variables(self)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input):: kara_input

    !open input netcdf and make list with all variables
    kara_input = type_input(_FILE_NAME_)
    !horizontal variables
    call self%add_var(kara_input,_OCEAN_TIME_)
    call self%add_day_number("day_number")
    !ice variables
    call self%add_ice_thickness(kara_input)
    call self%add_var(kara_input,_SNOW_THICKNESS_)
    call self%add_var(kara_input,_ICE_SURFACE_TEMPERATURE_)
    self%type_ice = ice(self%get_1st_dim_length("day_number"),&
      self%get_column(_ICE_THICKNESS_))
    !vertical variables
    call self%add_grid_on_faces(kara_input,&
      _DEPTH_ON_BOUNDARY_,self%type_ice%get_number_of_layers(),&
      "ice_water_index","water_bbl_index","bbl_sediments_index",&
      "number_of_boundaries")
    call self%add_grid_on_centers("middle_layer_depths","dz",&
                                  "number_of_layers","air_ice_indexes")
    call self%add_layer_thicknesses("layer_thicknesses")
    !2d variables
    !Add variables which are constants in sediments
    call self%add_constant_in_sed(kara_input,_TEMPERATURE_)
    call self%add_constant_in_sed(kara_input,_SALINITY_)
    call self%add_constant_in_sed(kara_input,_RHO_)
    !other variables
    call self%add_porosity("porosity","porosity_on_interfaces",&
      "porosity_factor_solutes_1","porosity_factor_solutes_2",&
      "porosity_factor_solids_1","porosity_factor_solids_2",&
      "tortuosity_on_interfaces")
    call self%add_diffusivity(kara_input,_TURBULENCE_,&
                              "molecular_diffusivity",&
                              "bioturbation_diffusivity")
    call self%print_list_variables('Allocated ipbm_standard_variables:')
    !call self%print_var("porosity")
    !call self%print_var("porosity_on_interfaces")
    !delete unneeded list
    call kara_input%delete_list()
  end subroutine

  subroutine add_standard_var(self,name_input,inname)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var
    class(*),allocatable:: temp

    call name_input%get_var(inname,var)
    !memory allocation problems occur without it
    allocate(temp,source=var)
    call self%add_item(temp)
  end subroutine
  !
  !Add discretized ice_thicknesses
  !
  subroutine add_ice_thickness(self,name_input)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    class(variable),allocatable:: var
    type(variable_1d) new_var_1d
    real(rk),dimension(:),allocatable:: value_1d

    call name_input%get_var(_ICE_THICKNESS_,var)
    select type(var)
    class is(variable_1d)
      allocate(value_1d(size(var%value,1)))
      value_1d = var%value-mod(var%value,_ICE_LAYERS_RESOLUTION_)
      new_var_1d = variable_1d(var%name,var%long_name,var%units,value_1d)
      call self%add_item(new_var_1d)
    end select
  end subroutine
  !
  !Adds bbl and sediments to depths of layers faces
  !
  subroutine add_grid_on_faces(self,name_input,&
             inname,ice_layers,ice_water_index,&
             water_bbl_index,bbl_sediments_index,&
             number_of_boundaries)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input),intent(in):: name_input
    character(len=*),intent(in):: inname
    integer         ,intent(in):: ice_layers
    character(len=*),intent(in):: ice_water_index
    character(len=*),intent(in):: water_bbl_index
    character(len=*),intent(in):: bbl_sediments_index
    character(len=*),intent(in):: number_of_boundaries

    class(variable)        ,allocatable:: var
    real(rk),dimension(:)  ,allocatable:: ice_thickness
    real(rk),dimension(:,:),allocatable:: value_2d
    type(alone_variable) new_var
    type(variable_2d) new_var_2d
    integer length,bbl_count,sediments_count
    integer ice_water,water_bbl,bbl_sediments,total_boundaries
    integer time
    integer i

    real(rk):: width_bbl = _WIDTH_BBL_
    real(rk):: resolution_bbl = _RESOLUTION_BBL_
    real(rk):: width_sediments = _WIDTH_SEDIMENTS_
    real(rk):: resolution_sediments = _RESOLUTION_SEDIMENTS_

    call name_input%get_var(inname,var)
    select type(var)
    class is(variable_1d)
      length=size(var%value,1)
    end select
    time = self%get_1st_dim_length("day_number")
    allocate(ice_thickness(time))
    ice_thickness = self%get_column(_ICE_THICKNESS_)

    bbl_count = width_bbl/resolution_bbl
    sediments_count = width_sediments/resolution_sediments
    allocate(value_2d(ice_layers+length+bbl_count+sediments_count,time))
    water_bbl = 1+bbl_count+sediments_count
    bbl_sediments = 1+sediments_count
    ice_water = length+bbl_count+sediments_count
    total_boundaries = ice_layers+length+bbl_count+sediments_count

    !adding indexes of inner boundaries
    new_var = alone_variable(ice_water_index,'','',ice_water)
    call self%add_item(new_var)
    new_var = alone_variable(water_bbl_index,'','',water_bbl)
    call self%add_item(new_var)
    new_var = alone_variable(bbl_sediments_index,'','',bbl_sediments)
    call self%add_item(new_var)
    new_var = alone_variable(number_of_boundaries,'','',total_boundaries)
    call self%add_item(new_var)

    select type(var)
    class is(variable_1d)
      value_2d = 0._rk
      forall (i=1:time) value_2d(water_bbl:ice_water,i) = var%value
      value_2d(ice_water+1:,:) = self%type_ice%do_grid(ice_thickness)
      value_2d(bbl_sediments,:) = var%value(1)
      do i = bbl_sediments+1,water_bbl
        value_2d(i,:) = value_2d(i-1,:)-resolution_bbl
      end do
      if (value_2d(water_bbl,1)<=value_2d(water_bbl+1,1)) then
        call fatal_error("BBL configurating","Wrong _BBL_WIDTH_")
      end if
      value_2d(1,:) = value_2d(bbl_sediments,:)+width_sediments
      do i = 2,(bbl_sediments-1)
        value_2d(i,:) = value_2d(i-1,:)-resolution_sediments
      end do
      new_var_2d = variable_2d(var%name,var%long_name,var%units,value_2d)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Adding layers","Wrong type")
    end select
  end subroutine

  subroutine add_grid_on_centers(self,inname,inname_increments,&
                                 number_of_layers,air_ice_indexes)
    class(ipbm_standard_variables),intent(inout):: self

    character(len=*),intent(in):: inname
    character(len=*),intent(in):: inname_increments
    character(len=*),intent(in):: number_of_layers
    character(len=*),intent(in):: air_ice_indexes

    class(variable)        ,allocatable:: var
    real(rk),dimension(:,:),allocatable:: mid_depths
    !increments between layer midpoints
    real(rk),dimension(:,:),allocatable:: dz

    type(alone_variable):: new_var
    type(variable_1d):: new_var_1d
    type(variable_2d):: new_var_2d

    integer i,j,length,time
    integer ice_water_i
    integer,dimension(:),allocatable:: active_layers

    time = self%get_1st_dim_length("day_number")
    allocate(active_layers(time))
    ice_water_i = self%get_value("ice_water_index")
    active_layers = ice_water_i+self%type_ice%get_active_layers()
    new_var_1d = variable_1d(air_ice_indexes,'','',active_layers)
    call self%add_item(new_var_1d)

    length = self%get_value("number_of_boundaries")-1._rk
    new_var = alone_variable(number_of_layers,'','',length)
    call self%add_item(new_var)

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    select type(var)
    class is(variable_2d)
      allocate(mid_depths(length,time))
      mid_depths = D_QNAN
      do j = 1,time
        do i = 1,(active_layers(j)-1)
          mid_depths(i,j) = abs((var%value(i+1,j)+&
                          var%value(i,j))/2._rk)
        end do
      end do
      new_var_2d = variable_2d(inname,'',var%units,mid_depths)
      call self%add_item(new_var_2d)
      !for example: for 4 boundaries exists 3 layers
      !and 3 midpoints of these layers and 2 dz
      allocate(dz(length-1,time))
      dz = D_QNAN
      do j = 1,time
        do i = 2,(active_layers(j)-1)
          dz(i-1,j) = abs(mid_depths(i,j)-mid_depths(i-1,j))
        end do
      end do
      new_var_2d = variable_2d(inname_increments,'',var%units,dz)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Add grid on centers",&
        "Wrong type")
    end select
  end subroutine

  subroutine add_day_number(self,inname)
    class(ipbm_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var

    call self%get_var(_OCEAN_TIME_,var)
    select type(var)
    class is(variable_1d)
      var%value = var%value/86400._rk
      var%name = inname
      call self%add_item(var)
    class default
      call fatal_error("Add day number",&
        "Wrong type")
    end select
  end subroutine

  integer function first_day(self)
    class(ipbm_standard_variables),intent(in):: self
    class(variable),allocatable:: var

    call self%get_var("day_number",var)
    select type(var)
    class is(variable_1d)
      first_day = int(var%value(1))
    end select
  end function

  subroutine add_layer_thicknesses(self,inname)
    class(ipbm_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable)        ,allocatable:: var
    real(rk),dimension(:)  ,allocatable:: air_ice_indexes
    real(rk),dimension(:,:),allocatable:: value_2d
    type(variable_2d):: new_var
    integer i,j,length,time

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    length = self%get_value("number_of_layers")
    time = self%get_1st_dim_length("day_number")
    allocate(value_2d(length,time))
    value_2d = D_QNAN
    allocate(air_ice_indexes(time))
    air_ice_indexes = self%get_column("air_ice_indexes")

    select type(var)
    class is(variable_2d)
      !forall(i = 1:length)&
      !  value_2d(i,:) = abs(var%value(i+1,:)-var%value(i,:))
      do j = 1,time
        do i = 1,(air_ice_indexes(j)-1)
          value_2d(i,j) = abs(var%value(i+1,j)-&
                          var%value(i,j))
        end do
      end do
      new_var = variable_2d(inname,'',var%units,value_2d)
      call self%add_item(new_var)
    class default
      call fatal_error("Add layer_thicknesses",&
        "Wrong type")
    end select
  end subroutine

  subroutine add_constant_in_sed(self,name_input,inname)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable        :: var
    real(rk),dimension(:,:),allocatable:: value_2d
    real(rk),dimension(:,:),allocatable:: salinity
    real(rk),dimension(:)  ,allocatable:: air_temp,ice_thickness
    type(variable_2d):: new_var
    integer i,length,time
    integer water_bbl_index,ice_water_index

    length          = self%get_value("number_of_layers")
    water_bbl_index = self%get_value("water_bbl_index")
    ice_water_index = self%get_value("ice_water_index")
    time            = self%get_1st_dim_length("day_number")
    allocate(value_2d(length,time))
    value_2d = 0._rk

    call name_input%get_var(inname,var)
    select type(var)
    class is(variable_2d)
      value_2d(water_bbl_index:ice_water_index-1,:time) = var%value
      select case(inname)
      case(_TEMPERATURE_)
        allocate(air_temp(time))
        allocate(ice_thickness(time))
        air_temp      = self%get_column(_ICE_SURFACE_TEMPERATURE_)
        ice_thickness = self%get_column(_ICE_THICKNESS_)
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_temperature(&
          air_temp,value_2d(ice_water_index-1,:time),ice_thickness)
      case(_SALINITY_)
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_brine_salinity(&
          value_2d(ice_water_index-1,:time))
      case(_RHO_)
        allocate(salinity(length,time))
        salinity = self%get_array(_SALINITY_)
        value_2d(ice_water_index:,:time) = self%type_ice%do_ice_brine_density(&
          salinity(ice_water_index:,:time))/1000._rk
      case default
        forall (i = 1:time)&
          value_2d(ice_water_index:,i) = &
          value_2d(ice_water_index-1,i)
      end select
      forall (i = 1:time)&
        value_2d(:water_bbl_index-1,i) =&
        value_2d(water_bbl_index,i)
      new_var = variable_2d(inname,var%long_name,var%units,value_2d)
      call self%add_item(new_var)
    class default
      call fatal_error("Adding constant in sediments variable",&
                       "Wrong type")
    end select
  end subroutine
  !
  !adopted from Phil Wallhead (PW)
  !Adds porosity,tortuosity and porosity factors, PW:
  !"These allow us to use a single equation to model diffusivity
  !updates in the water column and sediments, for both solutes and solids:
  !dC/dt = d/dz(pF2*kzti*d/dz(pF1*C))
  !where C has units [mass per unit total volume (water+sediments)]"
  !
  subroutine add_porosity(self,name_porosity,name_porosity_faces,&
      name_pf_solutes_1,name_pf_solutes_2,&
      name_pf_solids_1,name_pf_solids_2,name_tortuosity)
    !
    ! equations for solutes, PW:
    ! dC/dt = d/dz(kzti*dC/dz)            in the water column
    ! dC/dt = d/dz(phi*kzti*d/dz(C/phi))  in the sediments
    ! equations for solids, PW:
    ! dC/dt = d/dz(kzti*dC/dz)                   in the water column
    ! dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi))) in the sediments
    !
    class(ipbm_standard_variables),intent(inout):: self
    character(len=*),intent(in):: name_porosity
    character(len=*),intent(in):: name_porosity_faces
    character(len=*),intent(in):: name_pf_solutes_1
    character(len=*),intent(in):: name_pf_solutes_2
    character(len=*),intent(in):: name_pf_solids_1
    character(len=*),intent(in):: name_pf_solids_2
    character(len=*),intent(in):: name_tortuosity

    real(rk),dimension(:)  ,allocatable:: w_b,u_b
    real(rk),dimension(:,:),allocatable:: porosity
    real(rk),dimension(:,:),allocatable:: tortuosity
    real(rk),dimension(:,:),allocatable:: porosity_factor
    type(alone_variable):: dphidz_SWI
    type(variable_1d):: var_1d
    type(variable_2d):: new_var

    integer ice_water_index,swi_index,length,time,i
    real(rk) max_porosity,min_porosity,porosity_decay

    real(rk),dimension(:)  ,allocatable:: air_ice_indexes
    real(rk),dimension(:)  ,allocatable:: swi_depth
    real(rk),dimension(:,:),allocatable:: depth_center
    real(rk),dimension(:,:),allocatable:: depth_boundary

    max_porosity   = _MAX_POROSITY_
    min_porosity   = _MIN_POROSITY_
    porosity_decay = _POROSITY_DECAY_

    ice_water_index = self%get_value("ice_water_index")
    swi_index       = self%get_value("bbl_sediments_index")
    length          = self%get_value("number_of_layers")
    time            = self%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(time))
    air_ice_indexes = self%get_column("air_ice_indexes")
    allocate(depth_center,source=self%get_array("middle_layer_depths"))
    allocate(depth_boundary,source=self%get_array(_DEPTH_ON_BOUNDARY_))
    allocate(swi_depth,source = depth_boundary(swi_index,:))

    !Indices of layer interfaces in the sediments
    var_1d = variable_1d("k_sed","Indices of layer interfaces",&
                         "",(/(i,i=1,swi_index-1)/))
    call self%add_item(var_1d)
    !Indices of layer interfaces in the sediments including the SWI
    var_1d = variable_1d("k_sed1","Indices of layer interfaces with SWI",&
                         "",(/(i,i=1,swi_index)/))
    call self%add_item(var_1d)

    !alone variable
    !dphi/dz on SWI for sedimentation velocity caclulations
    dphidz_SWI = alone_variable("dphidz_SWI",&
                 "dphi/dz on SWI for sedimentation velocity caclulations",&
                 "m^-1",1.0_rk*(max_porosity-min_porosity)/porosity_decay)
    call self%add_item(dphidz_SWI)

    !for layers
    allocate(porosity(length,time))
    porosity = 1._rk
    forall (i = 1:time)
      porosity(1:swi_index-1,i) = min_porosity+(&
        max_porosity-min_porosity)*exp(-1._rk*(&
        depth_center(1:swi_index-1,i)-swi_depth)/&
        porosity_decay)
    end forall
    porosity(ice_water_index:,:) = self%type_ice%do_brine_relative_volume(&
          .true.,self%get_column(_ICE_THICKNESS_))
    new_var = variable_2d(name_porosity,'','',porosity)
    call self%add_item(new_var)

    !porosity factor 1 for solutes
    allocate(porosity_factor(length,time))
    porosity_factor = 1._rk/porosity
    !PW:
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume pore water] for solutes in sediments
    new_var = variable_2d(name_pf_solutes_1,&
                          "Total vol. to pore water vol.",&
                          '',porosity_factor)
    call self%add_item(new_var)

    !porosity factor 1 for solids
    porosity_factor = 1._rk
    do i = 1,time
      porosity_factor(air_ice_indexes(i):,i) = D_QNAN
    end do
    porosity_factor(1:swi_index-1,:) = 1._rk/&
      (1._rk-porosity(1:swi_index-1,:))
    !PW:
    !Factor to convert [mass per unit total volume]
    !to [mass per unit volume solids] for solids in sediments
    new_var = variable_2d(name_pf_solids_1,&
                          "Total vol. to vol. solids",&
                          '',porosity_factor)
    call self%add_item(new_var)

    deallocate(porosity)
    deallocate(porosity_factor)

    !for boundaries
    allocate(porosity(length+1,time))

    !background vertical advective velocities of particulates and
    !solutes on layer interfaces in the sediments (w_b, u_b)
    !these assume steady state compaction and neglect reaction terms
    allocate(w_b(swi_index))
    allocate(u_b(swi_index))

    porosity = 1._rk
    forall (i = 1:self%get_1st_dim_length("day_number"))
      porosity(1:swi_index,i) = min_porosity+(&
        max_porosity-min_porosity)*exp(-1._rk*(&
        depth_boundary(1:swi_index,i)-swi_depth)/&
        porosity_decay)
    end forall
    porosity(ice_water_index,:) = 0.5_rk
    porosity(ice_water_index+1:,:) = &
          self%type_ice%do_brine_relative_volume(&
          .false.,self%get_column(_ICE_THICKNESS_))
    new_var = variable_2d(name_porosity_faces,'','',porosity)
    call self%add_item(new_var)
    w_b = 0._rk
    u_b = 0._rk
    w_b(1:swi_index) = -((1._rk-min_porosity)/&
                       (1._rk-porosity(1:swi_index,1)))*_BURIAL_VELOCITY_
    u_b(1:swi_index) = -(min_porosity)/&
                       (porosity(1:swi_index,1))*_BURIAL_VELOCITY_
    var_1d = variable_1d("w_b","Particulates velocities","m s^-1",w_b)
    call self%add_item(var_1d)
    var_1d = variable_1d("u_b","Solutes velocities","m s^-1",u_b)
    call self%add_item(var_1d)

    !porosity factor 2 for solutes
    allocate(porosity_factor(length+1,time))
    porosity_factor = porosity
    !PW:
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    new_var = variable_2d(name_pf_solutes_2,&
                          "Area restriction factor - solutes",&
                          '',porosity_factor)
    call self%add_item(new_var)

    !porosity factor 2 for solids
    porosity_factor = 1._rk
    do i = 1,time
      porosity_factor(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    porosity_factor(1:swi_index,:) = &
      1._rk-porosity(1:swi_index,:)
    !PW:
    !Porosity-related area restriction factor for fluxes
    !across layer interfaces
    new_var = variable_2d(name_pf_solids_2,&
                          "Area restriction factor - solids",&
                          '',porosity_factor)
    call self%add_item(new_var)

    !tortuosity on layer interfaces
    !Boudreau 1996, eq. 4.120
    allocate(tortuosity(length+1,time))
    tortuosity = sqrt(1._rk-2._rk*log(porosity))
    new_var = variable_2d(name_tortuosity,&
                          "Tortousity on layer interfaces",&
                          '',tortuosity)
    call self%add_item(new_var)
  end subroutine

  subroutine add_diffusivity(self,name_input,name_eddy_diffusivity,&
                            name_molecular_diffusivity,&
                            name_bioturbation_diffusivity)
    class(ipbm_standard_variables),intent(inout):: self
    type(type_input),intent(in):: name_input
    character(len=*),intent(in):: name_eddy_diffusivity
    character(len=*),intent(in):: name_molecular_diffusivity
    character(len=*),intent(in):: name_bioturbation_diffusivity

    class(variable)        ,allocatable:: var
    !real(rk),dimension(:)  ,allocatable:: value_1d
    real(rk),dimension(:,:),allocatable:: eddy_kz
    real(rk),dimension(:,:),allocatable:: value_2d
    real(rk),dimension(:,:),allocatable:: tortuosity
    real(rk),dimension(:,:),allocatable:: porosity
    real(rk),dimension(:,:),allocatable:: depth_boundary
    real(rk),dimension(:)  ,allocatable:: air_ice_indexes
    type(variable_2d):: new_var_2d
    !type(variable_1d):: new_var_1d
    real    z_conv
    integer i,j,time
    integer number_of_boundaries
    integer water_bbl_index,ice_water_index
    integer bbl_sediments_index

    number_of_boundaries &
      = self%get_value("number_of_boundaries")
    water_bbl_index = self%get_value("water_bbl_index")
    ice_water_index = self%get_value("ice_water_index")
    bbl_sediments_index &
      = self%get_value("bbl_sediments_index")
    time = self%get_1st_dim_length("day_number")
    allocate(air_ice_indexes(time))
    air_ice_indexes = self%get_column("air_ice_indexes")
    allocate(eddy_kz(number_of_boundaries,time))
    eddy_kz = 0._rk

    !add eddy diffusivity
    allocate(depth_boundary(number_of_boundaries,time))
    depth_boundary = self%get_array(_DEPTH_ON_BOUNDARY_)
    call name_input%get_var(name_eddy_diffusivity,var)
    select type(var)
    class is(variable_2d)
      eddy_kz(water_bbl_index:ice_water_index-1,:) = var%value
      !linear interpolation in bbl
      forall (i = bbl_sediments_index+1:water_bbl_index)
        eddy_kz(i,:) = &
        !(i-bbl_sediments_index)*eddy_kz(water_bbl_index+1,:)/(&
        !water_bbl_index+1-bbl_sediments_index)
        eddy_kz(bbl_sediments_index,:)+&
        ((eddy_kz(water_bbl_index+1,:)-eddy_kz(bbl_sediments_index,:))/&
        (depth_boundary(water_bbl_index+1,:)-&
         depth_boundary(bbl_sediments_index,:)))*&
        (depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:))
      end forall
      do i = 1,time
        eddy_kz(air_ice_indexes(i)+1:,i) = D_QNAN
      end do
      new_var_2d = variable_2d(name_eddy_diffusivity,&
                               'Eddy diffusivity',var%units,eddy_kz)
      call self%add_item(new_var_2d)
    class default
      call fatal_error("Adding turbulence",&
        "Wrong type")
    end select

    !add molecular diffusivity
    allocate(tortuosity(number_of_boundaries,time))
    tortuosity = self%get_array("tortuosity_on_interfaces")
    allocate(value_2d(number_of_boundaries,time))
    value_2d = _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_
    value_2d(:bbl_sediments_index,:) = &
      _RELATIVE_DYNAMIC_VISCOSITY_*&
      _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_/&
      tortuosity(:bbl_sediments_index,:)**2
    do i = 1,time
      value_2d(ice_water_index:air_ice_indexes(i),i) = 0._rk
      value_2d(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    new_var_2d = variable_2d(name_molecular_diffusivity,&
                             "Molecular diffusivity",var%units,value_2d)
    call self%add_item(new_var_2d)

    !add bioturbation diffusivity
    value_2d = 0._rk
    forall (i = 1:bbl_sediments_index)
    !accuracy problem
      where (depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:)&
          <_MIXED_LAYER_DEPTH_)
        value_2d(i,:) = _MAX_BIOTURBATION_DIFFUSIVITY_
      else where
        value_2d(i,:) = _MAX_BIOTURBATION_DIFFUSIVITY_*exp(-1._rk*(&
          depth_boundary(i,:)-depth_boundary(bbl_sediments_index,:)-&
          _MIXED_LAYER_DEPTH_)/_DECAY_BIOTURBATION_SCALE_)
      end where
    end forall
    do i = 1,time
      value_2d(air_ice_indexes(i)+1:,i) = D_QNAN
    end do
    new_var_2d = variable_2d(name_bioturbation_diffusivity,&
                             "Bioturbation diffusivity",var%units,value_2d)
    call self%add_item(new_var_2d)

    !add brine channel gravity drainage diffusivity
    allocate(porosity(number_of_boundaries,time))
    value_2d = 0._rk
    porosity = self%get_array("porosity_on_interfaces")
    do i = 1,time
      z_conv = 0._rk
      !if no ice, go to the next day
      if (air_ice_indexes(i)/=ice_water_index) then
        do j = ice_water_index,air_ice_indexes(i)
          if (porosity(j,i)>0.072_rk) then
            z_conv = z_conv+_ICE_LAYERS_RESOLUTION_
          else
            exit
          end if
        end do
        value_2d(ice_water_index:j-1,i) = &
          z_conv*_GRAVITY_DRAINAGE_
      end if
    end do
    new_var_2d = variable_2d("ice_gravity_drainage",&
                             "Ice gravity drainage diffusivity",&
                              var%units,value_2d)
    call self%add_item(new_var_2d)
  end subroutine

  pure subroutine set_ipbm_state_variable(self,is_solid,&
      is_gas,use_bound_up,use_bound_low,bound_up,&
      bound_low,density,sinking_velocity,value,layer)
    class(ipbm_state_variable),intent(inout):: self
    logical,optional          ,intent(in)   :: is_solid
    logical,optional          ,intent(in)   :: is_gas
    integer,optional          ,intent(in)   :: use_bound_up
    integer,optional          ,intent(in)   :: use_bound_low
    real(rk),optional         ,intent(in)   :: bound_up
    real(rk),optional         ,intent(in)   :: bound_low
    real(rk),optional         ,intent(in)   :: density
    real(rk),allocatable,dimension(:),optional,intent(in)&
                                            :: sinking_velocity
    real(rk),optional         ,intent(in)   :: value
    integer ,optional         ,intent(in)   :: layer

    if(present(is_solid)) self%is_solid = is_solid
    if(present(is_gas)) self%is_gas = is_gas
    if(present(use_bound_up)) self%use_bound_up = use_bound_up
    if(present(use_bound_low)) self%use_bound_low = use_bound_low
    if(present(bound_up)) self%bound_up = bound_up
    if(present(bound_low)) self%bound_low = bound_low
    if(present(density)) self%density = density
    if(present(sinking_velocity)) self%sinking_velocity = sinking_velocity
    if(present(value).and.present(layer)) self%value(layer) = value
  end subroutine

  subroutine print_state_variable(self)
    class(ipbm_state_variable),intent(in):: self

    integer i

    do i = 1,72!size(self%value,1)
      write(*,'(f9.7,2x,i3)') self%value(i),i
    end do
    write(*,*) self%name
    !write(*,*) "is_solid:",self%is_solid
    !write(*,*) "density:",self%density
    !write(*,'(f9.3)') (/ self%value(size(self%value,1):1:-1) /)
    !write(*,'(x,a,2x,i2,f9.6)') 'up',self%use_bound_up,self%bound_up
    !write(*,'(x,a,i2,f9.6)') 'down',self%use_bound_low,self%bound_low
    !write(*,*) 'sinking',self%sinking_velocity
    _LINE_
  end subroutine
end module
