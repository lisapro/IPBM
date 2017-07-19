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

#include "../include/parameters.h"

module diff_mod
  use fabm_types,only: rk

  implicit none
  private
  public do_diffusive
contains
  pure function do_diffusive(N,dt,cnpar,posconc,h,Bcup,Bcdw,&
                Yup,Ydw,nuY_in,Lsour,Qsour,Taur,Yobs,Y,&
                i_sed_top,is_solid,i_ice_water,is_gas,&
                pF1_solutes,pF2_solutes,&
                pF1_solids,pF2_solids,pFSWIup_solutes,&
                pFSWIdw_solutes,pFSWIup_solids,pFSWIdw_solids)

    !REVISION HISTORY:
    !Original author(s): Lars Umlauf
    !-----------------------------------------------------------------------
    ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
    !-----------------------------------------------------------------------
    !
    ! !DESCRIPTION:
    ! This subroutine solves the one-dimensional diffusion equation
    ! including source terms,
    !  \begin{equation}
    !   \label{YdiffCenter}
    !    \partder{Y}{t}
    !    = \partder{}{z} \left( \nu_Y \partder{pY}{z} \right)
    !    - \frac{1}{\tau_R}(Y-Y_{obs})
    !    + Y L_{\text{sour}} + Q_{\text{sour}}
    !    \comma
    !  \end{equation}
    ! for al variables defined at the centers of the grid cells, and
    ! a diffusion coefficient $\nu_Y$ defined at the faces.
    ! Relaxation with time scale $\tau_R$ towards observed values
    ! $Y_{\text{obs}}$ is possible. $L_{\text{sour}}$ specifies a
    ! linear source term, and $Q_{\text{sour}}$ a constant source term.
    ! Central differences are used to discretize the problem
    ! as discussed in \sect{SectionNumericsMean}. The diffusion term,
    ! the linear source term, and the linear part arising from the
    ! relaxation term are treated
    ! with an implicit method, whereas the constant source term is treated
    ! fully explicit.
    ! p is a porosity factor to account any necessary changes of units from
    ! those of Y (assumed to be [mass per unit total volume]) to either
    ! [mass per unit volume pore water] or [mass per unit volume solids]
    ! when calculating Fickian gradients in the sediments
    ! (see Berner, 1980; Boudreau, 1997).
    ! p = 1 in the water column
    ! p = 1/phi for solutes in the sediments
    ! p = 1/(1-phi) for solids in the sediments
    ! A further volume factor is not needed to account for porosity
    ! because this subroutine assumes that [mass per unit total volume]
    ! is the modelling currency for all concentrations.
    ! This routine also makes special cases of the cells neighbouring the
    ! sediment-water interface, allowing for interphase mixing.
    !
    ! The input parameters {\tt Bcup} and {\tt Bcdw} specify the type
    ! of the upper and lower boundary conditions, which can be either
    ! Dirichlet or Neumann-type. {\tt Bcup} and {\tt Bcdw} must have integer
    ! values corresponding to the parameters {\tt Dirichlet}
    ! and {\tt Neumann}
    ! defined in the module {\tt util}, see \sect{sec:utils}.
    ! {\tt Yup} and {\tt Ydw} are the values of the boundary conditions at
    ! the surface and the bottom. Depending on the values of {\tt Bcup} and
    ! {\tt Bcdw}, they represent either fluxes or prescribed values.
    ! The integer {\tt posconc} indicates if a quantity is
    ! non-negative by definition ({\tt posconc}=1, such as for
    ! concentrations)
    ! or not ({\tt posconc}=0). For {\tt posconc}=1 and negative
    ! boundary fluxes, the source term linearisation according to
    ! \cite{Patankar80} is applied.
    !
    ! Note that fluxes \emph{entering} a boundary cell are counted positive
    ! by convention. The lower and upper position for prescribing these
    ! fluxes
    ! are located at the lowest und uppermost grid faces with index "0" and
    ! index "N", respectively. If values are prescribed, they are located at
    ! the centers with index "1" and index "N", respectively.
    !
    !
    !Phil Wallhead -
    ! 18/05/2016: Introduced porosity factors pF to convert
    !units for calculating Fickian gradients in the sediments
    ! 13/06/2016: Introduced special case porosity factors
    !pFSWIup, pFSWIdw and index i_sed_top to allow for interphase mixing
    !across SWI
    !
    !SY - minor changes to implement pure function
    !     choose btw solute, solid and gas
    !
    real(rk),dimension(0:N):: do_diffusive

    ! number of vertical layers
    integer,  intent(in)                :: N
    ! time step (s)
    integer,  intent(in)                :: dt
    ! "implicitness" parameter
    real(rk), intent(in)                :: cnpar
    ! 1: non-negative concentration, 0: else
    integer, intent(in)                 :: posconc
    ! layer thickness (m)
    real(rk), intent(in)                :: h(0:N)
    ! type of upper BC
    integer,  intent(in)                :: Bcup
    ! type of lower BC
    integer,  intent(in)                :: Bcdw
    ! value of upper BC
    real(rk), intent(in)                :: Yup
    ! value of lower BC
    real(rk), intent(in)                :: Ydw
    ! diffusivity of Y
    real(rk), intent(in)                :: nuY_in(0:N)
    ! linear source term
    ! (treated implicitly)
    real(rk), intent(in)                :: Lsour(0:N)
    ! constant source term
    ! (treated explicitly)
    real(rk), intent(in)                :: Qsour(0:N)
    ! relaxation time (s)
    real(rk), intent(in)                :: Taur(0:N)
    ! observed value of Y
    real(rk), intent(in)                :: Yobs(0:N)
    ! concentranions
    real(rk), intent(in)                :: Y(0:N)
    ! index to indentify the cell just below the SWI
    integer,  intent(in)                :: i_sed_top
    logical,  intent(in)                :: is_solid
    ! for gases handling in the ice core
    integer,  intent(in)                :: i_ice_water
    logical,  intent(in)                :: is_gas
    ! porosity factors
    real(rk), intent(in)                :: pF1_solutes(0:N)
    real(rk), intent(in)                :: pF2_solutes(0:N)
    real(rk), intent(in)                :: pF1_solids(0:N)
    real(rk), intent(in)                :: pF2_solids(0:N)
    ! porosity factors for sediment water interface (SWI)
    real(rk), intent(in)                :: pFSWIup_solutes
    real(rk), intent(in)                :: pFSWIdw_solutes
    real(rk), intent(in)                :: pFSWIup_solids
    real(rk), intent(in)                :: pFSWIdw_solids

    ! diffusivity * porosity
    real(rk):: nuY(0:N)
    ! porosity factor for Fickian gradient
    real(rk):: pF(0:N)
    ! porosity factor for cell just above the sediment-water interface
    real(rk):: pFSWIup
    ! porosity factor for cell just below the sediment-water interface
    real(rk):: pFSWIdw

    integer                   :: i
    real(rk)                  :: a,c,l
    real(rk), dimension(0:N)  :: au,bu,cu,du

    !choose btw solute, solid and gas
    if (.not.is_solid) then
      if (.not.is_gas) then
        !liquid case
        !dC/dt = d/dz(kzti*dC/dz)            in the water column
        !dC/dt = d/dz(phi*kzti*d/dz(C/phi))  in the ice and the sediments
        pF = pF1_solutes
        nuY = nuY_in*pF2_solutes
      else
        !gas case - currently like liquid
        !dC/dt = d/dz(phi*kzti*dC/dz)        in the ice
        !dC/dt = d/dz(kzti*dC/dz)            in the water column
        !dC/dt = d/dz(phi*kzti*d/dz(C/phi))  in the sediments
        !pF(:i_ice_water-1) = pF1_solutes(:i_ice_water-1)
        !pF(i_ice_water:) = 1._rk
        pF = pF1_solutes
        nuY = nuY_in*pF2_solutes
        !nuY(:i_ice_water-1) = nuY_in(:i_ice_water-1)*&
        !                      pF2_solutes(:i_ice_water-1)
        !nuY(i_ice_water:) = nuY_in(i_ice_water:)
      end if
      pFSWIup = pFSWIup_solutes
      pFSWIdw = pFSWIdw_solutes
    else
      !solid case
      !dC/dt = d/dz(phi*kzti*dC/dz)        in the ice
      !dC/dt = d/dz(kzti*dC/dz)            in the water column
      !dC/dt = d/dz((1-phi)*kzti*d/dz(C/(1-phi)))  in the sediments
      pF = pF1_solids
      !nuY = nuY_in*pF2_solids
      nuY(:i_ice_water-1) = nuY_in(:i_ice_water-1)*&
                            pF2_solids(:i_ice_water-1)
      nuY(i_ice_water:)   = nuY_in(i_ice_water:)*&
                            pF2_solutes(i_ice_water:)
      pFSWIup = pFSWIup_solids
      pFSWIdw = pFSWIdw_solids
    end if

    ! set up matrix
    do i=2,N-1
      c     = 2.0_rk*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0_rk*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     = dt*Lsour(i)
      !Special treatment for top cell of sediments (PWA)
      if (i.eq.i_sed_top) then
        !pF(i+1) -> pFSWIup, pF(i)*c -> pFSWIdw*c
        cu(i) =-cnpar*pFSWIup*c
        au(i) =-cnpar*pF(i-1)*a
        bu(i) = 1.0_rk + cnpar*(pF(i)*a + pFSWIdw*c) - l
        du(i) = (1.0_rk-(1.0_rk-cnpar)*(pF(i)*a+pFSWIdw*c))*Y(i)+&
                (1.0_rk-cnpar)*(pF(i-1)*a*Y(i-1)+pFSWIup*c*Y(i+1) )+dt*Qsour(i)
      !Special treatment for bottom cell of water column (PWA)
      else if (i.eq.i_sed_top+1) then
        !pF(i-1) -> pFSWIdw, pF(i)*a -> pFSWIup*a
        cu(i) =-cnpar*pF(i+1)*c
        au(i) =-cnpar*pFSWIdw*a
        bu(i) = 1.0_rk+cnpar*(pFSWIup*a + pF(i)*c)-l
        du(i) = (1.0_rk-(1.0_rk-cnpar)*(pFSWIup*a+pF(i)*c))*Y(i)+&
                (1.0_rk-cnpar)*(pFSWIdw*a*Y(i-1)+pF(i+1)*c*Y(i+1) )+dt*Qsour(i)
      else
        cu(i) =-cnpar*pF(i+1)*c
        au(i) =-cnpar*pF(i-1)*a
        bu(i) = 1.0_rk+cnpar*pF(i)*(a + c)-l
        du(i) = (1.0_rk-(1.0_rk-cnpar)*pF(i)*(a+c))*Y(i)+&
                (1.0_rk-cnpar)*(pF(i-1)*a*Y(i-1)+pF(i+1)*c*Y(i+1))+dt*Qsour(i)
      end if
      !Note: all terms of coefficients involving cnpar or (1-cnpar)
      ! derive from the Fickian gradient and therefore need to be
      ! multiplied by the appropriate porosity factors (pF)
    end do

    ! set up upper boundary condition
    select case(Bcup)
    case(_NEUMANN_)
      a     = 2.0d0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*pF(N-1)*a
      if (posconc.eq.1.and.Yup.lt.0.0_rk) then ! Patankar (1980) trick
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l  - dt*Yup/Y(N)/h(N)
         du(N) = Y(N) + dt*Qsour(N)   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      else
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l
         du(N) = Y(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      end if
    case(_DIRICHLET_)
      au(N) = 0.0_rk
      bu(N) = 1.0_rk
      du(N) = Yup
    end select

    ! set up lower boundary condition
    select case(Bcdw)
    case(_NEUMANN_)
      c     = 2.0d0*dt*nuY(1)/(h(1)+h(2))/h(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*pF(2)*c
      if (posconc.eq.1.and.Ydw.lt.0.0_rk) then !Patankar (1980) trick
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l - dt*Ydw/Y(1)/h(1)
         du(1) = Y(1) + dt*(Qsour(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      else
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l
         du(1) = Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      end if
    case(_DIRICHLET_)
      cu(1) = 0.0_rk
      bu(1) = 1.0_rk
      du(1) = Ydw
    end select

    ! relaxation to observed value
    if (minval(Taur).lt.1.d10) then
      do i=1,N
        bu(i)=bu(i)+dt/Taur(i)
        du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
    end if

    ! solve linear system
    do_diffusive = tridiagonal(bu,au,cu,du,&
                               N,1,N)
  end function

  pure function tridiagonal(bu,au,cu,du,&
                            N,fi,lt)
    ! adopted from
    ! !DESCRIPTION:
    ! A linear equation with tridiagonal matrix structure is solved here. The main
    ! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
    ! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
    ! The method used here is the simplified Gauss elimination, also called
    ! \emph{Thomas algorithm}.
    !
    ! !REVISION HISTORY:
    !  Original author(s): Hans Burchard & Karsten Bolding
    !  SY - minor changes to implement pure function
    !-----------------------------------------------------------------------
    ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
    !-----------------------------------------------------------------------

    real(rk)                          :: tridiagonal(0:N)
    ! !INPUT PARAMETERS:
    real(rk), dimension(0:N),intent(in) :: au,bu,cu,du
    integer, intent(in)                 :: N,fi,lt
    !
    ! !LOCAL VARIABLES:
    real(rk), dimension(0:N):: ru,qu
    integer                 :: i

    ru(lt)=au(lt)/bu(lt)
    qu(lt)=du(lt)/bu(lt)

    do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
    end do

    qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

    tridiagonal(fi)=qu(fi)
    do i=fi+1,lt
      tridiagonal(i)=qu(i)-ru(i)*tridiagonal(i-1)
    end do
  end function
end module
