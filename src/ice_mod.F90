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

module ice_mod
  use fabm_types, only: rk

  implicit none
  private
  !NaN value
  !REAL(rk), PARAMETER :: D_QNAN = &
  !          TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
  real(rk) D_QNAN

  public:: ice

  type ice
    private
    integer number_of_layers
    integer number_of_days
    integer,allocatable,dimension(:):: air_ice_index
    !depth, upper faces, centers
    real(rk),allocatable,dimension(:,:):: depth_face
    real(rk),allocatable,dimension(:,:):: depth_center
    !brine ice temperature
    real(rk),allocatable,dimension(:,:):: t_face
    real(rk),allocatable,dimension(:,:):: t_center
    !brine ice salinity
    real(rk),allocatable,dimension(:,:):: s_brine_face
    real(rk),allocatable,dimension(:,:):: s_brine_center
  contains
    private
    procedure,public:: get_active_layers
    procedure,public:: get_number_of_layers
    procedure,public:: do_grid
    procedure,public:: do_ice_temperature
    procedure,public:: do_ice_brine_salinity
    procedure,public:: do_brine_relative_volume
    procedure:: do_depths
    procedure:: do_ice_bulk_salinity
    procedure,public:: do_ice_brine_density
    procedure:: do_ice_bulk_density
  end type ice

  interface ice
    procedure constructor_ice
  end interface
contains
  function constructor_ice(number_of_days,ice_thickness)
    type(ice):: constructor_ice
    integer,intent(in):: number_of_days
    real(rk),dimension(:),intent(in):: ice_thickness

    !NaN
    D_QNAN = 0._rk
    D_QNAN = D_QNAN / D_QNAN

    constructor_ice%number_of_layers = &
      int(maxval(ice_thickness)/_ICE_LAYERS_RESOLUTION_)+2
    constructor_ice%number_of_days   = number_of_days
    allocate(constructor_ice%&
    air_ice_index (constructor_ice%number_of_days))
    allocate(constructor_ice%&
    depth_face    (constructor_ice%number_of_layers,number_of_days))
    allocate(constructor_ice%&
    depth_center  (constructor_ice%number_of_layers,number_of_days))
    allocate(constructor_ice%&
    t_face        (constructor_ice%number_of_layers,number_of_days))
    allocate(constructor_ice%&
    t_center      (constructor_ice%number_of_layers,number_of_days))
    allocate(constructor_ice%&
    s_brine_face  (constructor_ice%number_of_layers,number_of_days))
    allocate(constructor_ice%&
    s_brine_center(constructor_ice%number_of_layers,number_of_days))
    call constructor_ice%do_depths(ice_thickness)
  end function constructor_ice
  !
  !Private, saves depths for upper faces and layer centers
  !saves air_ice_index
  !
  subroutine do_depths(self,ice_thickness)
    class(ice)           ,intent(inout):: self
    real(rk),dimension(:),intent(in):: ice_thickness
    integer i

    forall (i = 1:self%number_of_layers)
      self%depth_face(i,:) = ice_thickness-i*_ICE_LAYERS_RESOLUTION_
    end forall
    where (self%depth_face<0._rk)
      self%depth_face = D_QNAN
    end where
    self%air_ice_index = minloc(self%depth_face,1)
    where (ice_thickness < 0.03_rk)
      self%air_ice_index = 0
    end where
    forall (i = 2:self%number_of_layers)
      self%depth_center(i,:) = &
        (self%depth_face(i,:)+self%depth_face(i-1,:))/2._rk
    end forall
    self%depth_center(1,:) = &
      (self%depth_face(1,:)+ice_thickness)/2._rk
  end subroutine do_depths
  !
  !Returns number of ice layers
  !
  integer function get_number_of_layers(self)
    class(ice)           ,intent(inout):: self

    get_number_of_layers = self%number_of_layers
  end function get_number_of_layers
  !
  !Returns number of ice layers
  !
  function get_active_layers(self)
    integer,allocatable  ,dimension(:) :: get_active_layers
    class(ice)           ,intent(inout):: self

    allocate(get_active_layers,source=self%air_ice_index)
  end function get_active_layers
  !
  !Returns range from ice-water, upper faces - z [m]
  !
  function do_grid(self,ice_thickness)
    real(rk),allocatable,dimension(:,:):: do_grid
    class(ice)           ,intent(inout):: self
    real(rk),dimension(:),intent(in)   :: ice_thickness

    integer i,j

    allocate(do_grid(self%number_of_layers,self%number_of_days))
    do_grid = D_QNAN
    !do i = 1,self%number_of_layers
    !  do_grid(i,:) = i*_ICE_LAYERS_RESOLUTION_
    !  where (do_grid(i,:)>ice_thickness)
    !    do_grid(i,:) = 0._rk
    !  end where
    !end do
    do j = 1,self%number_of_days
      do i = 1,self%air_ice_index(j)
        do_grid(i,j) = i*_ICE_LAYERS_RESOLUTION_
      end do
    end do
  end function do_grid
  !
  !Saves brine temperature for both faces and centers of layers [C]
  !brine temperature equals bulk temperature
  !Returns ice temperature column for centers of layers
  !
  function do_ice_temperature(self,air_temp,water_temp,ice_thickness)
    real(rk),allocatable,dimension(:,:):: do_ice_temperature
    class(ice),intent(inout)           :: self
    real(rk),dimension(:),intent(in)   :: air_temp,water_temp,ice_thickness
    integer i

    forall (i = 1:self%number_of_layers)
      self%t_face(i,:) = air_temp+((water_temp-air_temp)*&
        self%depth_face(i,:))/ice_thickness
      self%t_center(i,:) = air_temp+((water_temp-air_temp)*&
        self%depth_center(i,:))/ice_thickness
    end forall
    !self%t_face(1,:) = water_temp
    !self%t_center(1,:) = water_temp
    where (self%t_face>-2._rk) self%t_face = -2._rk
    where (self%t_center>-2._rk) self%t_center = -2._rk
    !where (self%t_face/=self%t_face) self%t_face = 0._rk
    !where (self%t_center/=self%t_center) self%t_center = 0._rk
    allocate(do_ice_temperature,source=self%t_center)
  end function do_ice_temperature
  !
  !Saves brine salinity for both faces and centers of layers [ppt]
  !Returns ice salinity column for centers of layers
  !
  function do_ice_brine_salinity(self,water_salt)
    real(rk),allocatable,dimension(:,:):: do_ice_brine_salinity
    class(ice),intent(inout)           :: self
    real(rk),dimension(:),intent(in)   :: water_salt

    integer i,j

    self%s_brine_face = D_QNAN
    do j = 1,self%number_of_days
      do i = 1,self%air_ice_index(j)
        if (self%t_face(i,j)<0._rk.and.self%t_face(i,j)>=-22.9_rk) then
          self%s_brine_face(i,j) = -3.9921_rk+(-22.700_rk*self%t_face(i,j))+(&
                   -1.0015_rk*self%t_face(i,j)**2)+(&
                   -0.019956_rk*self%t_face(i,j)**3)
        else if (self%t_face(i,j)<-22.9_rk.and.self%t_face(i,j)>=-44._rk) then
          self%s_brine_face(i,j) = 206.24_rk+(-1.8907_rk*self%t_face(i,j))+(&
                   -0.060868_rk*self%t_face(i,j)**2)+(&
                   -0.0010247_rk*self%t_face(i,j)**3)
        else if (self%t_face(i,j)<-44._rk) then
          self%s_brine_face(i,j) = -4442.1_rk+(-277.86_rk*self%t_face(i,j))+(&
                   -5.501_rk*self%t_face(i,j)**2)+(&
                   -0.03669_rk*self%t_face(i,j)**3)
        !else
        !  self%s_brine_face(i,j) = 0._rk
        end if
      end do
    end do
    !where (self%t_face<0._rk.and.self%t_face>=-22.9_rk)
    !  self%s_brine_face = -3.9921_rk+(-22.700_rk*self%t_face)+(&
    !           -1.0015_rk*self%t_face**2)+(&
    !           -0.019956_rk*self%t_face**3)
    !else where (self%t_face<-22.9_rk.and.self%t_face>=-44._rk)
    !  self%s_brine_face = 206.24_rk+(-1.8907_rk*self%t_face)+(&
    !           -0.060868_rk*self%t_face**2)+(&
    !           -0.0010247_rk*self%t_face**3)
    !else where (self%t_face<-44._rk)
    !  self%s_brine_face = -4442.1_rk+(-277.86_rk*self%t_face)+(&
    !           -5.501_rk*self%t_face**2)+(&
    !           -0.03669_rk*self%t_face**3)
    !else where
    !  self%s_brine_face = 0._rk
    !end where

    self%s_brine_center = D_QNAN
    do j = 1,self%number_of_days
      do i = 1,self%air_ice_index(j)
        if (self%t_center(i,j)<0._rk.and.self%t_center(i,j)>=-22.9_rk) then
          self%s_brine_center(i,j) = -3.9921_rk+(-22.700_rk*self%t_center(i,j))+(&
                   -1.0015_rk*self%t_center(i,j)**2)+(&
                   -0.019956_rk*self%t_center(i,j)**3)
        else if (self%t_center(i,j)<-22.9_rk.and.self%t_center(i,j)>=-44._rk) then
          self%s_brine_center(i,j) = 206.24_rk+(-1.8907_rk*self%t_center(i,j))+(&
                   -0.060868_rk*self%t_center(i,j)**2)+(&
                   -0.0010247_rk*self%t_center(i,j)**3)
        else if (self%t_center(i,j)<-44._rk) then
          self%s_brine_center(i,j) = -4442.1_rk+(-277.86_rk*self%t_center(i,j))+(&
                   -5.501_rk*self%t_center(i,j)**2)+(&
                   -0.03669_rk*self%t_center(i,j)**3)
        !else
        !  self%s_brine_center(i,j) = 0._rk
        end if
      end do
    end do
    !where (self%t_center<0._rk.and.self%t_center>=-22.9_rk)
    !  self%s_brine_center = -3.9921_rk+(-22.700_rk*self%t_center)+(&
    !           -1.0015_rk*self%t_center**2)+(&
    !           -0.019956_rk*self%t_center**3)
    !else where (self%t_center<-22.9_rk.and.self%t_center>=-44._rk)
    !  self%s_brine_center = 206.24_rk+(-1.8907_rk*self%t_center)+(&
    !           -0.060868_rk*self%t_center**2)+(&
    !           -0.0010247_rk*self%t_center**3)
    !else where (self%t_center<-44._rk)
    !  self%s_brine_center = -4442.1_rk+(-277.86_rk*self%t_center)+(&
    !           -5.501_rk*self%t_center**2)+(&
    !           -0.03669_rk*self%t_center**3)
    !else where
    !  self%s_brine_center = 0._rk
    !end where
    !self%s_brine_face(1,:) = water_salt
    !self%s_brine_center(1,:) = water_salt
    allocate(do_ice_brine_salinity,source=self%s_brine_center)
  end function do_ice_brine_salinity
  !
  !Private procedure for calculating ice_bulk_salinity [ppt]
  !
  function do_ice_bulk_salinity(self,depth,ice_thickness)
    real(rk),allocatable,dimension(:,:):: do_ice_bulk_salinity
    class(ice)           ,intent(inout):: self
    real(rk),dimension(:,:) ,intent(in):: depth
    real(rk),dimension(:)   ,intent(in):: ice_thickness
    !ratio btwn the distance from the ice surface and ice thickness
    real(rk),allocatable,dimension(:,:):: z_p
    integer i

    allocate(do_ice_bulk_salinity(self%number_of_layers,self%number_of_days))
    allocate(z_p(self%number_of_layers,self%number_of_days))
    forall (i = 1:self%number_of_layers)
      z_p(i,:) = depth(i,:)/ice_thickness
      do_ice_bulk_salinity(i,:) = 19.539_rk*(z_p(i,:)**2)&
                                 -19.93_rk*z_p(i,:)+8.913_rk
    end forall
    !do_ice_bulk_salinity(1,:) =  self%s_brine_center(1,:)
  end function do_ice_bulk_salinity
  !
  !Private procedure, brine density through ice [g*m-3]
  !
  function do_ice_brine_density(self,brine_salinity)
    real(rk),allocatable,dimension(:,:):: do_ice_brine_density
    class(ice)           ,intent(inout):: self
    real(rk),dimension(:,:) ,intent(in):: brine_salinity
    real(rk):: c ![g*m-3*ppt-1]

    c = 8.e-4_rk
    allocate(do_ice_brine_density(self%number_of_layers,self%number_of_days))
    do_ice_brine_density = (1._rk+c*brine_salinity)*1.e6_rk
  end function do_ice_brine_density
  !
  !Private procedure, bulk density through ice [g*m-3]
  !
  function do_ice_bulk_density(self,brine_salinity,bulk_salinity)
    real(rk),allocatable,dimension(:,:):: do_ice_bulk_density
    class(ice)           ,intent(inout):: self
    real(rk),dimension(:,:) ,intent(in):: brine_salinity
    real(rk),dimension(:,:) ,intent(in):: bulk_salinity
    real(rk),allocatable,dimension(:,:):: bs
    real(rk),allocatable,dimension(:,:):: brine_density
    real(rk):: dens_pure !density of pure ice

    dens_pure = 912000._rk ![g*m-3]
    allocate(bs,source=brine_salinity)
    allocate(brine_density,source=&
      self%do_ice_brine_density(brine_salinity))
    allocate(do_ice_bulk_density(self%number_of_layers,self%number_of_days))
    where (bs < 2._rk) bs = 2._rk! function has spikes below x=1
    do_ice_bulk_density = dens_pure*brine_density*bs/&
    (brine_density*bs-bulk_salinity*(brine_density-dens_pure))
  end function do_ice_bulk_density
  !
  !Returns ice porosity
  !
  function do_brine_relative_volume(self,is_center,ice_thickness)
    real(rk),allocatable,dimension(:,:):: do_brine_relative_volume
    !real(rk),allocatable,dimension(:,:):: temporary
    class(ice)           ,intent(inout):: self
    logical              ,intent(in)   :: is_center
    real(rk),dimension(:),intent(in)   :: ice_thickness

    real(rk),allocatable,dimension(:,:):: brine_salinity
    real(rk),allocatable,dimension(:,:):: bulk_salinity
    real(rk),allocatable,dimension(:,:):: depth

    allocate(do_brine_relative_volume(&
      self%number_of_layers,self%number_of_days))
    !allocate(temporary(&
    !  self%number_of_layers,self%number_of_days))
    select case(is_center)
    case(.true.)
      allocate(brine_salinity,source=self%s_brine_center)
      allocate(depth,source=self%depth_center)
    case(.false.)
      allocate(brine_salinity,source=self%s_brine_face)
      allocate(depth,source=self%depth_face)
    end select
    where (brine_salinity<10._rk) brine_salinity = 10._rk!to fix above 1 values
    allocate(bulk_salinity,&
      source=self%do_ice_bulk_salinity(depth,ice_thickness))
    do_brine_relative_volume = &
      (self%do_ice_bulk_density(brine_salinity,bulk_salinity)*&
       bulk_salinity)/(self%do_ice_brine_density(brine_salinity)*&
       brine_salinity)
    !do_brine_relative_volume(1,:) = 0.5_rk
    where (do_brine_relative_volume > 0.9_rk) &
      do_brine_relative_volume = 0.9_rk
    !temporary = do_brine_relative_volume
  end function do_brine_relative_volume
end module ice_mod
