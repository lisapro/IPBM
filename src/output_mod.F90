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

module output_mod
  use types_mod
  use variables_mod,only: ipbm_standard_variables,&
                          ipbm_state_variable
  use fabm_driver
  use fabm      , only: type_model,fabm_get_bulk_diagnostic_data,&
                        type_bulk_variable_id,standard_variables
  use fabm_types, only: rk
  use netcdf

  implicit none
  private
  public type_output

  type:: type_output
    private
    logical first
    integer first_layer
    integer last_layer
    integer number_of_layers
    !netCDF file id
    integer            :: nc_id
    !parameter_ids
    integer            :: z_id,z_id_faces
    integer            :: par_id
    integer,allocatable:: standard_id(:)
    integer,allocatable:: parameter_id(:)
    integer,allocatable:: parameter_id_diag(:)
  contains
    private
    procedure:: initialize
    procedure,public:: save
    procedure,public:: close
  end type

  interface type_output
    module procedure type_output_constructor
  end interface
contains
  function type_output_constructor(model,standard_vars,infile,&
                                   first_layer,last_layer,&
                                   number_of_layers)
    type(type_output):: type_output_constructor
    type(type_model)              ,intent(in):: model
    type(ipbm_standard_variables) ,intent(in):: standard_vars
    character(len=*)              ,intent(in):: infile
    integer                       ,intent(in):: first_layer
    integer                       ,intent(in):: last_layer
    integer                       ,intent(in):: number_of_layers

    call type_output_constructor%initialize(model,standard_vars,infile,&
                                            first_layer,last_layer,&
                                            number_of_layers)
  end function

  subroutine initialize(self,model,standard_vars,infile,&
                        first_layer,last_layer,number_of_layers)
    class(type_output)                  ,intent(inout):: self
    type(type_model)                     ,intent(in):: model
    class(ipbm_standard_variables),target,intent(in):: standard_vars
    character(len=*)                     ,intent(in):: infile
    integer                              ,intent(in):: first_layer
    integer                              ,intent(in):: last_layer
    integer                              ,intent(in):: number_of_layers

    type(type_bulk_variable_id):: fabm_standard_id
  
    class(variable)               ,allocatable:: curr
    class(ipbm_standard_variables),pointer    :: temporary

    !dimension lengths
    integer,parameter:: time_len = NF90_UNLIMITED
    integer nlev
    integer dim_ids(2), dim_ids_faces(2)

    integer z_dim_id, z_dim_id_faces, time_dim_id
    integer ip, number_of_variables

    self%first = .true.
    self%first_layer = first_layer
    self%last_layer = last_layer
    self%number_of_layers = number_of_layers

    nlev = last_layer-first_layer+1
    self%nc_id = -1
    call check(nf90_create(infile,NF90_CLOBBER,self%nc_id))
    !define the dimensions
    call check(nf90_def_dim(self%nc_id,_OCEAN_TIME_,time_len,time_dim_id))
    call check(nf90_def_dim(self%nc_id,"z",nlev,z_dim_id))
    call check(nf90_def_dim(self%nc_id,"z_faces",&
                            nlev+1,z_dim_id_faces))
    dim_ids(1) = z_dim_id
    dim_ids(2) = time_dim_id
    dim_ids_faces(1) = z_dim_id_faces
    dim_ids_faces(2) = time_dim_id
    
    !define variables
    call check(nf90_def_var(self%nc_id,"z",&
                NF90_REAL,z_dim_id,self%z_id))
    call check(nf90_def_var(self%nc_id,"z_faces",&
                NF90_REAL,z_dim_id_faces,self%z_id_faces))
    
    !to make standard_vars intent(in)
    temporary => standard_vars
    number_of_variables = temporary%count_list()
    allocate(self%standard_id(number_of_variables))
    do ip = 1,number_of_variables
      call temporary%get_var_by_number(ip,curr)
      select type(curr)
      class is(variable_1d)
        if(curr%name == _OCEAN_TIME_) then
          call check(nf90_def_var(self%nc_id,curr%name,&
                      NF90_REAL,time_dim_id,self%standard_id(ip)))
          call check(set_attributes(ncid=self%nc_id,id=self%standard_id(ip),&
                      units=curr%units,&
                      long_name=curr%long_name))
        end if
      class is(variable_2d)
        !for variables on layers
        if (curr%variable_1d_size() == number_of_layers) then
          call check(nf90_def_var(self%nc_id,curr%name,&
                      NF90_REAL,dim_ids,self%standard_id(ip)))
          call check(set_attributes(ncid=self%nc_id,id=self%standard_id(ip),&
                      units=curr%units,&
                      long_name=curr%long_name))
        !for variables on interfaces
        else if (curr%variable_1d_size() == number_of_layers+1) then
          call check(nf90_def_var(self%nc_id,curr%name,&
                      NF90_REAL,dim_ids_faces,self%standard_id(ip)))
          call check(set_attributes(ncid=self%nc_id,id=self%standard_id(ip),&
                      units=curr%units,&
                      long_name=curr%long_name))
        end if
      end select
      deallocate(curr)
    end do
    
    allocate(self%parameter_id(size(model%state_variables)))
    do ip = 1,size(model%state_variables)
      call check(nf90_def_var(self%nc_id,&
                 model%state_variables(ip)%name,&
                 NF90_REAL,dim_ids,self%parameter_id(ip)))
      call check(set_attributes(ncid=self%nc_id,id=self%parameter_id(ip),&
                 units=model%state_variables(ip)%units,&
                 long_name=model%state_variables(ip)%long_name,&
                 missing_value=model%state_variables(ip)%missing_value))
    end do
    
    allocate(self%parameter_id_diag(size(model%diagnostic_variables)))
    do ip = 1,size(model%diagnostic_variables)
      if (model%diagnostic_variables(ip)%save) then
        call check(nf90_def_var(&
                   self%nc_id,model%diagnostic_variables(ip)%name,&
                   NF90_REAL,dim_ids,self%parameter_id_diag(ip)))
        call check(set_attributes(ncid=self%nc_id,&
                   id=self%parameter_id_diag(ip),&
                   units=model%diagnostic_variables(ip)%units,&
                   long_name=model%diagnostic_variables(ip)%long_name,&
                   missing_value = &
                   model%diagnostic_variables(ip)%missing_value))
      end if
    end do

    !fabm standard variables
    fabm_standard_id = model%get_bulk_variable_id(&
               standard_variables%downwelling_photosynthetic_radiative_flux)
    call check(nf90_def_var(&
               self%nc_id,fabm_standard_id%variable%name,&
               NF90_REAL,dim_ids,self%par_id))
    call check(set_attributes(ncid=self%nc_id,&
               id=self%par_id,&
               units=fabm_standard_id%variable%units,&
               long_name=fabm_standard_id%variable%long_name,&
               missing_value = &
               fabm_standard_id%variable%missing_value))
    !end define
    call check(nf90_enddef(self%nc_id))
  end subroutine initialize

  subroutine save(self,model,standard_vars,state_vars,&
                  z,z_faces,day,air_ice_index)

    class(type_output),intent(inout):: self
    type (type_model)                    ,intent(in):: model
    type(ipbm_standard_variables),target ,intent(in):: standard_vars
    type(ipbm_state_variable),allocatable,intent(in):: state_vars(:)
    real(rk),allocatable,dimension(:)    ,intent(in):: z
    real(rk),allocatable,dimension(:)    ,intent(in):: z_faces
    integer                              ,intent(in):: day
    integer                              ,intent(in):: air_ice_index

    class(variable)               ,allocatable:: curr
    class(ipbm_standard_variables),pointer    :: temporary
    
    !NaN value
    !REAL(rk), PARAMETER :: D_QNAN = &
    !          TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_rk)
    real(rk) D_QNAN

    integer ip,number_of_variables
    integer edges(2),edges_faces(2),start(2),start_time(1),edges_time(1)
    real(rk) temp_matrix(self%number_of_layers)
    real(rk) temp_matrix_faces(self%number_of_layers+1)
    real(rk) foo(1)

    type(type_bulk_variable_id)  :: fabm_standard_id
    real(rk),dimension(:),pointer:: dat
    
    !NaN
    D_QNAN = 0._rk
    D_QNAN = D_QNAN / D_QNAN

    !write data
    edges(1) = self%last_layer-self%first_layer+1
    edges(2) = 1
    edges_faces(1) = edges(1)+1
    edges_faces(2) = 1
    start(1) = 1
    start(2) = day
    start_time(1) = day
    edges_time(1) = 1
    
    if (self%first) then
      call check(nf90_put_var(self%nc_id,self%z_id,&
                 z(self%first_layer:self%last_layer),start,edges))
      call check(nf90_put_var(self%nc_id,self%z_id_faces,&
                 z_faces(self%first_layer:self%last_layer+1),start,edges_faces))
      self%first = .false.
    end if

    if (self%nc_id.ne.-1) then
      !to make standard_vars intent(in)
      temporary => standard_vars
      number_of_variables = temporary%count_list()
      do ip = 1,number_of_variables
        call temporary%get_var_by_number(ip,curr)
        select type(curr)
        class is(variable_1d)
          if(curr%name == _OCEAN_TIME_) then
            foo(1) = curr%value(day)
            call check(nf90_put_var(self%nc_id,self%standard_id(ip),foo,&
                                    start_time,edges_time))
          end if
        class is(variable_2d)
        !for variables on layers
          if (curr%variable_1d_size() == self%number_of_layers) then
            temp_matrix = curr%value(:,day)
            temp_matrix(air_ice_index:) = D_QNAN
            call check(nf90_put_var(self%nc_id,self%standard_id(ip),&
                        temp_matrix(self%first_layer:self%last_layer),&
                        start,edges))
          else if (curr%variable_1d_size() == self%number_of_layers+1) then
            temp_matrix_faces = curr%value(:,day)
            temp_matrix_faces(air_ice_index+1:) = D_QNAN
            call check(nf90_put_var(self%nc_id,self%standard_id(ip),&
                        temp_matrix_faces(self%first_layer:self%last_layer+1),&
                        start,edges_faces))
          end if
        end select
      deallocate(curr)
      end do
      
      do ip = 1,size(model%state_variables)
        call check(nf90_put_var(self%nc_id,self%parameter_id(ip),&
                   state_vars(ip)%value(self%first_layer:self%last_layer),&
                   start,edges))
      end do
      do ip = 1,size(model%diagnostic_variables)
        if (model%diagnostic_variables(ip)%save) then
          temp_matrix = fabm_get_bulk_diagnostic_data(model,ip)
          temp_matrix(air_ice_index:) = D_QNAN
          if (maxval(abs(temp_matrix)).lt.1.0E37) then
            call check(nf90_put_var(self%nc_id,self%parameter_id_diag(ip),&
                       temp_matrix(self%first_layer:self%last_layer),&
                       start,edges))
          end if
        end if
      end do
      
      !fabm standard variables
      fabm_standard_id = model%get_bulk_variable_id(&
                standard_variables%downwelling_photosynthetic_radiative_flux)
      dat => model%get_data(fabm_standard_id)
      call check(nf90_put_var(self%nc_id,self%par_id,&
                  dat(self%first_layer:self%last_layer),&
                  start,edges))
    
      call check(nf90_sync(self%nc_id))
    end if
  end subroutine save

  subroutine close(self)
    class(type_output),intent(inout):: self

    if (self%nc_id.ne.-1) then
      call check(nf90_close(self%nc_id))
      deallocate(self%parameter_id)
      deallocate(self%parameter_id_diag)
    end if
    self%nc_id = -1
  end subroutine close

  integer function set_attributes(ncid,id,&
                                  units,long_name,                 &
                                  valid_min,valid_max,valid_range, &
                                  scale_factor,add_offset,         &
                                  FillValue,missing_value,         &
                                  C_format,FORTRAN_format)
    !
    ! !DESCRIPTION:
    !  This routine is used to set a number of attributes for
    !  variables. The routine makes heavy use of the {\tt optional} keyword.
    !  The list of recognized keywords is very easy to extend. We have
    !  included a sub-set of the COARDS conventions.
    !
    ! !USES:
    !  IMPLICIT NONE
    !
    ! !INPUT PARAMETERS:
    integer, intent(in)                     :: ncid,id
    character(len=*), optional              :: units,long_name
    real, optional                          :: valid_min,valid_max
    real, optional                          :: valid_range(2)
    real, optional                          :: scale_factor,add_offset
    double precision, optional              :: FillValue,missing_value
    character(len=*), optional              :: C_format,FORTRAN_format
    !
    ! !REVISION HISTORY:
    !  Original author(s): Karsten Bolding & Hans Burchard
    !
    ! !LOCAL VARIABLES:
    integer                                 :: iret
    real                                    :: vals(2)
    !
    !
    !-----------------------------------------------------------------------
    !
    if(present(units)) then
      iret = nf90_put_att(ncid,id,'units',trim(units))
    end if

    if(present(long_name)) then
      iret = nf90_put_att(ncid,id,'long_name',trim(long_name))
    end if

    if(present(C_format)) then
      iret = nf90_put_att(ncid,id,'C_format',trim(C_format))
    end if

    if(present(FORTRAN_format)) then
      iret = nf90_put_att(ncid,id,'FORTRAN_format',trim(FORTRAN_format))
    end if

    if(present(valid_min)) then
      vals(1) = valid_min
      iret = nf90_put_att(ncid,id,'valid_min',vals(1:1))
    end if

    if(present(valid_max)) then
      vals(1) = valid_max
      iret = nf90_put_att(ncid,id,'valid_max',vals(1:1))
    end if

    if(present(valid_range)) then
      vals(1) = valid_range(1)
      vals(2) = valid_range(2)
      iret = nf90_put_att(ncid,id,'valid_range',vals(1:2))
    end if

    if(present(scale_factor)) then
      vals(1) = scale_factor
      iret = nf90_put_att(ncid,id,'scale_factor',vals(1:1))
    end if

    if(present(add_offset)) then
      vals(1) = add_offset
      iret = nf90_put_att(ncid,id,'add_offset',vals(1:1))
    end if

    if(present(FillValue)) then
      vals(1) = FillValue
      iret = nf90_put_att(ncid,id,'_FillValue',vals(1:1))
    end if

    if(present(missing_value)) then
      vals(1) = missing_value
      iret = nf90_put_att(ncid,id,'missing_value',vals(1:1))
    end if

    set_attributes = 0
  end function set_attributes

  subroutine check(status)
    integer, intent(in):: status

    if (status .ne. NF90_NOERR) then
      call fatal_error("Netcdf internal",&
                        nf90_strerror(status))
    end if
  end subroutine check
end module
