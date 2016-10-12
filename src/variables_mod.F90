#include "../include/brom.h"
#include "../include/parameters.h"

module variables_mod
  use types_mod
  use input_mod
  use fabm_driver

  implicit none

  type,extends(list_variables):: brom_standard_variables
  contains
    private
    procedure:: initialize=>initialize_standard_variables
    procedure:: add_var=>add_standard_var
    procedure:: add_grid_on_faces
    procedure:: add_grid_on_centers
    procedure:: add_day_number
    procedure:: add_layer_thicknesses
    procedure:: add_constant_in_sed
    procedure,public:: first_day
  end type

  type,extends(variable_1d):: brom_state_variable
    integer  use_bound_up
    integer  use_bound_low
    real(rk) bound_up
    real(rk) bound_low
    real(rk) sinking_velocity
  contains
    procedure:: set_brom_state_variable
    procedure:: print_state_variable
  end type

  interface brom_standard_variables
    module procedure brom_standard_variables_constructor
  end interface
contains
  function brom_standard_variables_constructor()
    type(brom_standard_variables):: brom_standard_variables_constructor

    call brom_standard_variables_constructor%initialize()
  end function

  !
  !Initialize standard variables list
  !
  subroutine initialize_standard_variables(self)
    class(brom_standard_variables),intent(inout):: self
    type(type_input):: kara_input

    !open input netcdf and make list with all variables
    kara_input = type_input(_FILE_NAME_)
    !vertical variables
    call self%add_grid_on_faces(kara_input,&
      _DEPTH_ON_BOUNDARY_,"water_bbl_index",&
      "bbl_sediments_index","number_of_boundaries")
    call self%add_grid_on_centers("middle_layer_depths",&
                                  "number_of_layers")
    call self%add_layer_thicknesses("layer_thicknesses")
    !horizontal variables
    call self%add_var(kara_input,_OCEAN_TIME_)
    call self%add_day_number("day_number")
    !2d variables
    call self%add_constant_in_sed(kara_input,_TEMPERATURE_)
    call self%add_constant_in_sed(kara_input,_SALINITY_)
    !call self%add_var(kara_input,_TURBULENCE_)

    call self%print_var(_TEMPERATURE_)
    call self%print_var(_SALINITY_)
    !call self%print_var(_DEPTH_ON_BOUNDARY_)
    !call self%print_var("middle_layer_depths")
    !call self%print_var("layer_thicknesses")
    !call self%print_var("water_bbl_index")
    !call self%print_var("bbl_sediments_index")
    !call self%print_var("number_of_boundaries")
    !call self%print_var("number_of_layers")
    !call self%print_list_variables('Allocated brom_standard_variables:')
    !delete unneeded list
    call kara_input%delete_list()
    stop
  end subroutine

  subroutine add_standard_var(self,name_input,inname)
    class(brom_standard_variables),intent(inout):: self
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
  !Adds bbl and sediments to depths of layers faces
  !
  subroutine add_grid_on_faces(self,name_input,&
             inname,water_bbl_index,bbl_sediments_index,&
             number_of_boundaries)
    class(brom_standard_variables),intent(inout):: self
    type(type_input),intent(in):: name_input
    character(len=*),intent(in):: inname
    character(len=*),intent(in):: water_bbl_index
    character(len=*),intent(in):: bbl_sediments_index
    character(len=*),intent(in):: number_of_boundaries

    class(variable)      ,allocatable:: var
    real(rk),dimension(:),allocatable:: value_1d
    type(alone_variable) new_var
    type(variable_1d) new_var_1d
    integer length,bbl_count,sediments_count
    integer water_bbl,bbl_sediments,total_boundaries
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

    bbl_count = int(width_bbl/resolution_bbl)
    sediments_count = int(width_sediments/resolution_sediments)
    allocate(value_1d(length+bbl_count+sediments_count))
    water_bbl = 1+bbl_count+sediments_count
    bbl_sediments = 1+sediments_count
    total_boundaries = length+bbl_count+sediments_count

    !adding indexes of inner boundaries
    new_var = alone_variable(water_bbl_index,'',water_bbl)
    call self%add_item(new_var)
    new_var = alone_variable(bbl_sediments_index,'',bbl_sediments)
    call self%add_item(new_var)
    new_var = alone_variable(number_of_boundaries,'',total_boundaries)
    call self%add_item(new_var)

    select type(var)
    class is(variable_1d)
      value_1d = 0._rk
      value_1d(water_bbl:) = var%value
      value_1d(bbl_sediments) = var%value(1)
      do i = bbl_sediments+1,water_bbl
        value_1d(i) = value_1d(i-1)-resolution_bbl
      end do
      if (value_1d(water_bbl)<=value_1d(water_bbl+1)) then
        call fatal_error("BBL configurating","Wrong _BBL_WIDTH_")
      end if
      value_1d(1) = value_1d(bbl_sediments)+width_sediments
      do i = 2,(bbl_sediments-1)
        value_1d(i) = value_1d(i-1)-resolution_sediments
      end do
      new_var_1d = variable_1d(var%name,'',value_1d)
      call self%add_item(new_var_1d)
    class default
      call fatal_error("Adding layers","Wrong type")
    end select
  end subroutine

  subroutine add_grid_on_centers(self,inname,number_of_layers)
    class(brom_standard_variables),intent(inout):: self
    character(len=*),intent(in):: inname
    character(len=*),intent(in):: number_of_layers
    class(variable)      ,allocatable:: var
    real(rk),dimension(:),allocatable:: value_1d
    type(alone_variable):: new_var
    type(variable_1d):: new_var_1d
    integer i,length

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    length = int(self%get_value("number_of_boundaries")-1._rk)
    allocate(value_1d(length))
    new_var = alone_variable(number_of_layers,'',length)
    call self%add_item(new_var)

    select type(var)
    type is(variable_1d)
      forall(i = 1:length)&
        value_1d(i) = abs((var%value(i+1)+var%value(i))/2._rk)
      new_var_1d = variable_1d(inname,'',value_1d)
      call self%add_item(new_var_1d)
    class default
      call fatal_error("Add grid on centers",&
        "Wrong type")
    end select
  end subroutine

  subroutine add_day_number(self,inname)
    class(brom_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var

    call self%get_var(_OCEAN_TIME_,var)
    select type(var)
    type is(variable_1d)
      var%value = var%value/86400._rk
      var%name = inname
      call self%add_item(var)
    class default
      call fatal_error("Add day number",&
        "Wrong type")
    end select
  end subroutine

  integer function first_day(self)
    class(brom_standard_variables),intent(in):: self
    class(variable),allocatable:: var

    call self%get_var("day_number",var)
    select type(var)
    type is(variable_1d)
      first_day = int(var%value(1))
    end select
  end function

  subroutine add_layer_thicknesses(self,inname)
    class(brom_standard_variables),intent(inout):: self
    character(len=*)              ,intent(in)   :: inname
    class(variable)      ,allocatable:: var
    real(rk),dimension(:),allocatable:: value_1d
    type(variable_1d):: new_var
    integer i,length

    call self%get_var(_DEPTH_ON_BOUNDARY_,var)
    length = int(self%get_value("number_of_layers"))
    allocate(value_1d(length))

    select type(var)
    type is(variable_1d)
      forall(i = 1:length)&
        value_1d(i) = abs(var%value(i+1)-var%value(i))
      new_var = variable_1d(inname,'',value_1d)
      call self%add_item(new_var)
    class default
      call fatal_error("Add layer_thicknesses",&
        "Wrong type")
    end select
  end subroutine

  subroutine add_constant_in_sed(self,name_input,inname)
    class(brom_standard_variables),intent(inout):: self
    type(type_input)              ,intent(in)   :: name_input
    character(len=*)              ,intent(in)   :: inname
    class(variable),allocatable:: var
    real(rk),dimension(:,:),allocatable:: value_2d
    type(variable_2d):: new_var
    integer i,length,time
    integer water_bbl_index

    length = int(self%get_value("number_of_layers"))
    water_bbl_index = int(self%get_value("water_bbl_index"))
    time = int(self%get_1st_dim_length("day_number"))
    allocate(value_2d(length,time))
    value_2d = 0._rk

    call name_input%get_var(inname,var)
    select type(var)
    type is(variable_2d)
      value_2d(water_bbl_index:,:time) = var%value
      forall (i = 1:time)&
        value_2d(:water_bbl_index-1,i) =&
        value_2d(water_bbl_index,i)
      new_var = variable_2d(inname,'',value_2d)
      call self%add_item(new_var)
    class default
      call fatal_error("Add temperature",&
        "Wrong type")
    end select
  end subroutine

  pure subroutine set_brom_state_variable(self,use_bound_up,&
      use_bound_low,bound_up,bound_low,sinking_velocity)
    class(brom_state_variable),intent(inout):: self
    integer,optional          ,intent(in)   :: use_bound_up
    integer,optional          ,intent(in)   :: use_bound_low
    real(rk),optional         ,intent(in)   :: bound_up
    real(rk),optional         ,intent(in)   :: bound_low
    real(rk),optional         ,intent(in)   :: sinking_velocity

    if(present(use_bound_up)) self%use_bound_up = use_bound_up
    if(present(use_bound_low)) self%use_bound_low = use_bound_low
    if(present(bound_up)) self%bound_up = bound_up
    if(present(bound_low)) self%bound_low = bound_low
    if(present(sinking_velocity)) self%sinking_velocity = sinking_velocity
  end subroutine

  subroutine print_state_variable(self)
    class(brom_state_variable),intent(in):: self

    _LINE_
    write(*,*) self%name
    write(*,'(f9.3)') (/ self%value(size(self%value,1):1:-1) /)
    write(*,'(x,a,2x,i2,f9.6)') 'up',self%use_bound_up,self%bound_up
    write(*,'(x,a,i2,f9.6)') 'down',self%use_bound_low,self%bound_low
    write(*,*) 'sinking',self%sinking_velocity
    _LINE_
  end subroutine
end module
