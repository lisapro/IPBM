#include "../include/parameters.h"

module diff_mod
  use fabm_types,only: rk
  implicit none

  private
  public do_diffusive
contains
  pure function do_diffusive(N,dt,cnpar,posconc,h,Bcup,Bcdw,&
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Y)
    !adopted from
    !Original author(s): Lars Umlauf
    !SY - minor changes to implement pure function
    !-----------------------------------------------------------------------
    ! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
    !-----------------------------------------------------------------------
    real(rk),dimension(1:N):: do_diffusive

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
    real(rk), intent(in)                :: nuY(0:N)
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

    integer                   :: i
    real(rk)                  :: a,c,l
    real(rk), dimension(0:N)  :: au,bu,cu,du
    real(rk), dimension(0:N)  :: result

    ! set up matrix
    do i=2,N-1
      c     = 2.0_rk*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0_rk*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     = dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
      bu(i) = 1.0_rk + cnpar*(a + c) - l
      du(i) = (1.0_rk - (1.0_rk-cnpar)*(a + c))*Y(i)                  &
            + (1.0_rk - cnpar)*( a*Y(i-1) + c*Y(i+1) ) + dt*Qsour(i)
    end do

    ! set up upper boundary condition
    select case(Bcup)
    case(_NEUMANN_)
      a     = 2.0d0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*a
      if (posconc .eq. 1 .and. Yup.lt.0.0_rk) then ! Patankar (1980) trick
         bu(N) =  1.0_rk - au(N) - l  - dt*Yup/Y(N)/h(N)
         du(N) = Y(N) + dt*Qsour(N)   &
               + (1.0_rk - cnpar)*a*(Y(N-1)-Y(N))
      else
         bu(N) =  1.0_rk - au(N) - l
         du(N) = Y(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1.0_rk - cnpar)*a*(Y(N-1)-Y(N))
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

      cu(1) =-cnpar*c
      if (posconc.eq.1 .and. Ydw.lt.0.0_rk) then ! Patankar (1980) trick
        bu(1) = 1.0_rk - cu(1) - l - dt*Ydw/Y(1)/h(1)
        du(1) = Y(1) + dt*(Qsour(1))   &
              + (1.0_rk - cnpar)*c*(Y(2)-Y(1))
      else
        bu(1) = 1.0_rk - cu(1) - l
        du(1) = Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
              + (1.0_rk - cnpar)*c*(Y(2)-Y(1))
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
    result = tridiagonal(bu,au,cu,du,&
                               N,1,N)
    do_diffusive = result(1:N)
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
