module ice_mod
  use fabm_types, only: rk

  implicit none
  private
  public:: ice

  type ice
    private
    integer number_of_layers
    integer number_of_days
    !range from ice-water, upper faces
    real(rk),allocatable,dimension(:):: z
    !depth, upper faces
    real(rk),allocatable,dimension(:):: depth
    !brine ice temperature
    real(rk),allocatable,dimension(:,:):: t
    !brine ice salinity
    real(rk),allocatable,dimension(:,:):: s
  contains
    private
    procedure,public:: do_grid_faces
    procedure,public:: do_ice_temperature
    procedure,public:: do_ice_salinity
  end type ice

  interface ice
    procedure constructor_ice
  end interface
contains
  function constructor_ice(number_of_layers,number_of_days)
    type(ice):: constructor_ice
    integer,intent(in):: number_of_layers
    integer,intent(in):: number_of_days

    constructor_ice%number_of_layers = number_of_layers
    constructor_ice%number_of_days = number_of_days
    allocate(constructor_ice%z(number_of_layers))
    allocate(constructor_ice%depth(number_of_layers))
    allocate(constructor_ice%t(number_of_layers,number_of_days))
    allocate(constructor_ice%s(number_of_layers,number_of_days))
  end function constructor_ice
  
  function do_grid_faces(self,ice_thickness)
    real(rk),allocatable,dimension(:):: do_grid_faces
    class(ice),intent(inout):: self
    real(rk)  ,intent(in)   :: ice_thickness
    integer i
    real(rk) delta
    
    allocate(do_grid_faces(self%number_of_layers))
    self%z(1) = 0.03_rk
    self%depth(1) = ice_thickness-0.03_rk
    delta = (ice_thickness-0.03)&
        /(self%number_of_layers-1)
    do i = 2,self%number_of_layers
      self%z(i) = self%z(i-1)+delta
      self%depth(i) = self%depth(i-1)-delta
    end do
    self%depth(self%number_of_layers) = 0._rk
    do_grid_faces = self%z
  end function do_grid_faces
  
  function do_ice_temperature(self,air_temp,water_temp,ice_thickness)
    real(rk),allocatable,dimension(:,:):: do_ice_temperature
    class(ice),intent(inout)           :: self
    real(rk),dimension(:),intent(in)   :: air_temp,water_temp,ice_thickness
    integer i
    
    do i = 1,self%number_of_layers
      self%t(i,:) = air_temp+((water_temp-air_temp)*self%depth(i))/ice_thickness
    end do
    self%t(1,:) = water_temp
    where (self%t>-0.2_rk) self%t = -0.2_rk
    allocate(do_ice_temperature,source=self%t)
  end function do_ice_temperature
  
  function do_ice_salinity(self,water_salt)
    real(rk),allocatable,dimension(:,:):: do_ice_salinity
    class(ice),intent(inout)           :: self
    real(rk),dimension(:),intent(in)   :: water_salt
    
    where (self%t<0._rk.and.self%t>=-22.9_rk)
      self%s = -3.9921_rk+(-22.700_rk*self%t)+(&
               -1.0015_rk*self%t**2)+(-0.019956_rk*self%t**3)
    else where (self%t<-22.9_rk.and.self%t>=-44._rk)
      self%s = 206.24_rk+(-1.8907_rk*self%t)+(&
               -0.060868_rk*self%t**2)+(&
               -0.0010247_rk*self%t**3)
    else where (self%t<-44._rk)
      self%s = -4442.1_rk+(-277.86_rk*self%t)+(&
               -5.501_rk*self%t**2)+(&
               -0.03669_rk*self%t**3)
    end where
    self%s(1,:) = water_salt
    allocate(do_ice_salinity,source=self%s)
  end function do_ice_salinity
end module ice_mod
