module output_mod
  use variables_mod
  use fabm_driver
  use fabm, only: type_model,fabm_get_bulk_diagnostic_data
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
    integer            :: z_id,time_id,iz_id
    integer            :: t_id,s_id,kz2_id
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
  function type_output_constructor(model,infile,&
                                   first_layer,last_layer,&
                                   number_of_layers)
    type(type_output):: type_output_constructor
    type(type_model),intent(in):: model
    character(len=*),intent(in):: infile
    integer         ,intent(in):: first_layer
    integer         ,intent(in):: last_layer
    integer         ,intent(in):: number_of_layers

    call type_output_constructor%initialize(model,infile,&
                                            first_layer,last_layer,&
                                            number_of_layers)
  end function

  subroutine initialize(self,model,infile,&
                        first_layer,last_layer,number_of_layers)
    class(type_output),intent(inout):: self
    type(type_model)  ,intent(in):: model
    character(len=*)  ,intent(in):: infile
    integer           ,intent(in):: first_layer
    integer           ,intent(in):: last_layer
    integer           ,intent(in):: number_of_layers

    !dimension lengths
    integer,parameter:: time_len = NF90_UNLIMITED
    integer nlev
    integer dim1d
    integer dim_ids(2)

    integer z_dim_id, time_dim_id
    integer ip,ilast

    self%first = .true.
    self%first_layer = first_layer
    self%last_layer = last_layer
    self%number_of_layers = number_of_layers

    nlev = last_layer-first_layer+1
    self%nc_id = -1
    call check(nf90_create(infile,NF90_CLOBBER,self%nc_id))
    !define the dimensions
    call check(nf90_def_dim(self%nc_id,"depth",nlev,z_dim_id))
    call check(nf90_def_dim(self%nc_id,"time",time_len,time_dim_id))
    !define coordinates
    dim1d = z_dim_id
    call check(nf90_def_var(self%nc_id,"depth",NF90_REAL,dim1d,self%z_id))
    dim1d = time_dim_id
    call check(nf90_def_var(self%nc_id,"time",NF90_REAL,dim1d,self%time_id))
    !define variables
    dim_ids(1) = z_dim_id
    dim_ids(2) = time_dim_id
    allocate(self%parameter_id(size(model%state_variables)))
    do ip = 1,size(model%state_variables)
      ilast = index(model%state_variables(ip)%path,'/',.true.)
      call check(nf90_def_var(self%nc_id,&
                 model%state_variables(ip)%path(ilast+1:),&
                 NF90_REAL,dim_ids,self%parameter_id(ip)))
      call check(set_attributes(ncid=self%nc_id,id=self%parameter_id(ip),&
                 units=model%state_variables(ip)%units,&
                 long_name=model%state_variables(ip)%long_name,&
                 missing_value=model%state_variables(ip)%missing_value))
    end do
    allocate(self%parameter_id_diag(size(model%diagnostic_variables)))
    do ip = 1,size(model%diagnostic_variables)
      if (model%diagnostic_variables(ip)%save) then
        ilast = index(model%diagnostic_variables(ip)%path,'/',.true.)
        call check(nf90_def_var(&
                   self%nc_id,model%diagnostic_variables(ip)%&
                   path(ilast+1:),&
                   NF90_REAL,dim_ids,self%parameter_id_diag(ip)))
        call check(set_attributes(ncid=self%nc_id,&
                   id=self%parameter_id_diag(ip),&
                   units=model%diagnostic_variables(ip)%units,&
                   long_name=model%diagnostic_variables(ip)%long_name,&
                   missing_value = &
                   model%diagnostic_variables(ip)%missing_value))
      end if
    end do
    call check(nf90_def_var(self%nc_id,"temp",NF90_REAL,dim_ids,self%t_id))
    call check(nf90_def_var(self%nc_id,"salt",NF90_REAL,dim_ids,self%s_id))
    call check(nf90_def_var(self%nc_id,"turb",NF90_REAL,dim_ids,&
                            self%kz2_id))
    call check(nf90_def_var(self%nc_id,"radiative_flux",&
                            NF90_REAL,dim_ids,self%iz_id))
    !end define
    call check(nf90_enddef(self%nc_id))
  end subroutine initialize

  subroutine save(self,model,state_vars,day,&
                  temp,salt,turb,radiative_flux,depth)
    class(type_output),intent(inout):: self
    type (type_model)                    ,intent(in):: model
    type(brom_state_variable),allocatable,intent(in):: state_vars(:)
    integer                              ,intent(in):: day
    real(rk),allocatable,dimension(:)    ,intent(in):: temp
    real(rk),allocatable,dimension(:)    ,intent(in):: salt
    real(rk),allocatable,dimension(:)    ,intent(in):: turb
    real(rk),allocatable,dimension(:)    ,intent(in):: radiative_flux
    real(rk),allocatable,dimension(:)    ,intent(in):: depth

    integer ip,i
    integer edges(2),start(2),start_time(1),edges_time(1)
    real(rk) temp_matrix(self%number_of_layers)
    real(rk) foo(1)

    !write data
    edges(1) = self%last_layer-self%first_layer+1
    edges(2) = 1
    start(1) = 1
    start(2) = day
    start_time(1) = day
    edges_time(1) = 1

    if (self%first) then
      call check(nf90_put_var(self%nc_id,self%z_id,&
                 depth(self%first_layer:self%last_layer),start,edges))
      self%first = .false.
    end if

    foo(1) = real(day)
    if (self%nc_id.ne.-1) then
      call check(nf90_put_var(self%nc_id,self%time_id,foo,&
                              start_time,edges_time))
      do ip = 1,size(model%state_variables)
        call check(nf90_put_var(self%nc_id,self%parameter_id(ip),&
                   state_vars(ip)%value(self%first_layer:self%last_layer),&
                   start,edges))
      end do
      do ip = 1,size(model%diagnostic_variables)
        if (model%diagnostic_variables(ip)%save) then
          temp_matrix = fabm_get_bulk_diagnostic_data(model,ip)
          if (maxval(abs(temp_matrix)).lt.1.0E37) then
            call check(nf90_put_var(self%nc_id,self%parameter_id_diag(ip),&
                       temp_matrix(self%first_layer:self%last_layer),&
                       start,edges))
          end if
        end if
      end do
      call check(nf90_put_var(self%nc_id,self%t_id,&
                              temp(self%first_layer:self%last_layer),&
                              start,edges))
      call check(nf90_put_var(self%nc_id,self%s_id,&
                              salt(self%first_layer:self%last_layer),&
                              start,edges))
      call check(nf90_put_var(self%nc_id,self%kz2_id,&
                              turb(self%first_layer:self%last_layer),&
                              start,edges))
      call check(nf90_put_var(self%nc_id,self%iz_id,&
                              radiative_flux(&
                              self%first_layer:self%last_layer),&
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