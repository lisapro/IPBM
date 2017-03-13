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

module sed_mod
  use fabm, only: type_model,fabm_get_vertical_movement,fabm_do_bottom
  use fabm_types, only: rk
  use io_ascii, only: get_brom_par

  implicit none
  private
  public calculate_sed
contains
  !
  !adapted from Phil Wallhead code
  !Calculates vertical advection (sedimentation)
  !in the water column and sediments
  !
  subroutine calculate_sed(model,cc,&
                           i,kmax,par_max,id_O2,julianday,&
                           k_bbl_sed,dynamic_w_sed,&
                           bctype_top,bctype_bottom,&
                           is_solid,k_sed1,&
                           dcc_R,fick,bc_top,bc_bottom,phi1,&
                           w_b,u_b,kz_bio,hz,dz,rho,dt,K_O2s,&
                           dphidz_SWI,cc0,&
                           wbio,wti,sink,dcc,bott_flux,bott_source)
    !REVISION HISTORY:
    !Original author(s): Phil Wallhead

    !Input/output variables
    type (type_model), intent(inout)         :: model
    real(rk), dimension(:,:,:), intent(inout):: cc

    !Input variables
    integer,               intent(in):: i
    integer,               intent(in):: k_max
    integer,               intent(in):: par_max
    integer,               intent(in):: id_O2
    integer,               intent(in):: julianday
    integer,               intent(in):: k_bbl_sed
    integer,               intent(in):: dynamic_w_sed
    integer,dimension(:,:),intent(in):: bctype_top
    integer,dimension(:,:),intent(in):: bctype_bottom
    integer,dimension(:),  intent(in):: is_solid
    integer,dimension(:),  intent(in):: k_sed1

    real(rk),dimension(:,:,:),intent(in):: dcc_R
    real(rk),dimension(:,:,:),intent(in):: fick
    real(rk),dimension(:,:),  intent(in):: bc_top
    real(rk),dimension(:,:),  intent(in):: bc_bottom
    real(rk),dimension(:,:),  intent(in):: phi1
    real(rk),dimension(:,:),  intent(in):: w_b
    real(rk),dimension(:,:),  intent(in):: u_b
    real(rk),dimension(:,:),  intent(in):: kz_bio
    real(rk),dimension(:),    intent(in):: hz
    real(rk),dimension(:),    intent(in):: dz
    real(rk),dimension(:),    intent(in):: rho
    !Sedimentation model time step [seconds]
    real(rk),                 intent(in):: dt
    real(rk),                 intent(in):: K_O2s
    real(rk),                 intent(in):: dphidz_SWI
    real(rk),                 intent(in):: cc0

    !Output variables
    real(rk), dimension(:,:,:),intent(out):: wbio
    real(rk), dimension(:,:,:),intent(out):: wti
    real(rk), dimension(:,:,:),intent(out):: sink
    real(rk), dimension(:,:,:),intent(out):: dcc
    real(rk), dimension(:,:),  intent(out):: bott_flux
    real(rk), dimension(:,:),  intent(out):: bott_source

    !Local variables
    real(rk):: O2stat
    real(rk):: w_1m(k_max+1,par_max)
    real(rk):: w_1(k_max+1)
    real(rk):: u_1(k_max+1)
    real(rk):: w_1c(k_max+1)
    real(rk):: u_1c(k_max+1)

    integer k, ip, idf

    dcc = 0.0_rk
    w_1 = 0.0_rk
    u_1 = 0.0_rk
    w_1c = 0.0_rk
    u_1c = 0.0_rk
    sink(i,:,:) = 0.0_rk

    !Compute vertical velocity in water column (sinking/floating)
    !using the FABM.
    wbio = 0.0_rk
    do k=1,k_max
      !Note: wbio is on layer midpoints
      call fabm_get_vertical_movement(model,i,i,k,wbio(i:i,k,:))
    end do

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

    if(k_bbl_sed.ne.k_max) then
        O2stat = cc(i,k_bbl_sed,id_O2) / (cc(i,k_bbl_sed,id_O2) + K_O2s)   !Oxygen status of sediments set by O2 level just above sediment surface
    else
        O2stat = 0.0_rk
    endif
    w_1(k_bbl_sed+1) = -1.0_rk*O2stat*kz_bio(i,k_bbl_sed+1)*dphidz_SWI / (1.0_rk-phi1(i,k_bbl_sed+1))
    u_1(k_bbl_sed+1) = O2stat*kz_bio(i,k_bbl_sed+1)*dphidz_SWI / phi1(i,k_bbl_sed+1)
    if (dynamic_w_sed.eq.1) then
        w_1m = 0.0_rk
        do ip=1,par_max !Sum over contributions from each particulate variable
            if (is_solid(ip).eq.1) then
                !First set rhop_i*(1-phi)*w_1i at the SWI
                w_1m(k_bbl_sed+1,ip) = wbio(i,k_bbl_sed,ip)*min(cc(i,k_bbl_sed,ip),cc(i,k_bbl_sed-1,ip))
                !Now set rhop_i*(1-phi)*w_1i in the sediments by integrating the reaction terms
                do k=k_bbl_sed+2,k_max+1
                    w_1m(k,ip) = w_1m(k-1,ip) + dcc_R(i,k-1,ip)*hz(k-1)
                end do
                w_1m(k_sed1,ip) = w_1m(k_sed1,ip) / (rho(ip)*(1.0_rk-phi1(i,k_sed1))) !Divide by rhop_i*(1-phi) to get w_1c(i)
                w_1c(k_sed1) = w_1c(k_sed1) + w_1m(k_sed1,ip)                         !Add to total w_1c
            end if
        end do
        !Now calculate u from w using (4) above
        u_1c(k_sed1) = (w_1c(k_max+1) - (1.0_rk-phi1(i,k_sed1))*w_1(k_sed1)) / phi1(i,k_sed1)
    end if


    !Interpolate wbio from FABM (defined on layer midpoints, as for concentrations) to wti on the layer interfaces
    wti(i,1,:) = 0.0_rk
    do ip=1,par_max
        !Air-sea interface (unused)
        wti(i,1,ip) = wbio(i,1,ip)
        !Water column layer interfaces, not including SWI
        wti(i,2:k_bbl_sed,ip) = wbio(i,1:k_bbl_sed-1,ip) + 0.5_rk*hz(1:k_bbl_sed-1)*(wbio(i,2:k_bbl_sed,ip) - wbio(i,1:k_bbl_sed-1,ip)) / dz(1:k_bbl_sed-1)
        !Sediment layer interfaces, including SWI
        if (is_solid(ip).eq.0) then
            wti(i,k_sed1,ip) = u_b(i,k_sed1) + u_1(k_sed1) + u_1c(k_sed1)
        else
            wti(i,k_sed1,ip) = w_b(i,k_sed1) + w_1(k_sed1) + w_1c(k_sed1)
        end if
        wti(i,k_max+1,ip) = wti(i,k_max,ip)
    end do

    !! Perform advective flux calculation and cc update
    !This uses a simple first order upwind differencing scheme (FUDM)
    !It uses the fluxes sink in a consistent manner and therefore conserves mass
    do idf=1,freq_sed
        !!Calculate sinking fluxes at layer interfaces (sink, units strictly [mass/unit total area/second])
        !Air-sea interface
        sink(i,1,:) = 0.0_rk
        !Water column and sediment layer interfaces
        do k=2,k_max+1
            sink(i,k,:) = wti(i,k,:)*cc(i,k-1,:)
            !This is an upwind differencing approx., hence the use of cc(k-1)
            !Note: factors phi, (1-phi) are not needed in the sediments because cc is in units [mass per unit total volume]
        end do


        !!Calculate tendencies dcc = dcc/dt = -dF/dz on layer midpoints (top/bottom not used where Dirichlet bc imposed)
        do k=1,k_max
            dcc(i,k,:) = -1.0_rk * (sink(i,k+1,:)-sink(i,k,:)) / hz(k)
        end do


        !!Time integration
        !cc surface layer
        do ip=1,par_max
            if (bctype_top(i,ip).gt.0) then
                cc(i,1,ip) = bc_top(i,ip)
            else
                cc(i,1,ip) = cc(i,1,ip) + dt*dcc(i,1,ip) !Simple Euler time step of size dt = 86400.*dt/freq_sed [s]
            end if
        end do
        !cc intermediate layers
        do k=2,(k_max-1)
            cc(i,k,:) = cc(i,k,:) + dt*dcc(i,k,:)
        end do
        !cc bottom layer
        bott_flux = 0.0_rk
        bott_source = 0.0_rk
        call fabm_do_bottom(model, i, i, bott_flux(i:i,:),bott_source(i:i,:))
        do ip=1,par_max
            if (is_solid(ip).eq.1) then
                dcc(i,k_max,ip) = dcc(i,k_max,ip) + bott_flux(i,ip) / hz(k_max)
                sink(i,k_max+1,:) = bott_flux(i,:)
            endif
        end do
!            cc(i,k_max,ip) = cc(i,k_max,ip) + dt*dcc(i,k_max,ip)
        do ip=1,par_max
            if (bctype_bottom(i,ip).gt.0) then
                cc(i,k_max,ip) = bc_bottom(i,ip)
            else
                cc(i,k_max,ip) = cc(i,k_max,ip) + dt*dcc(i,k_max,ip)
            end if
        end do
        cc(i,:,:) = max(cc0, cc(i,:,:)) !Impose resilient concentration
    end do

    end subroutine calculate_sed
!=======================================================================================================================


end module calculate
