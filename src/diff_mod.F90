module diff_mod
  use fabm_types,only: rk
  implicit none

  private
  public do_diffusive
contains
  !
  !adopted from
  !Original author(s): Lars Umlauf
  !
  pure function do_diffusive(N,dt,cnpar,posconc,h,Vc,Af,Bcup,Bcdw, &
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Y)
!  number of vertical layers
   integer,  intent(in)                :: N
!  time step (s)
   real(rk), intent(in)                :: dt
!  "implicitness" parameter
   real(rk), intent(in)                :: cnpar
!  1: non-negative concentration, 0: else
   integer, intent(in)                 :: posconc
!  layer thickness (m)
   real(rk), intent(in)                :: h(0:N)
!  hypsograph at grid centre
   real(rk), intent(in)                :: Vc(0:N)
!  hypsograph at grid face
   real(rk), intent(in)                :: Af(0:N)
!  type of upper BC
   integer,  intent(in)                :: Bcup
!  type of lower BC
   integer,  intent(in)                :: Bcdw
!  value of upper BC
   real(rk), intent(in)                :: Yup
!  value of lower BC
   real(rk), intent(in)                :: Ydw
!  diffusivity of Y
   real(rk), intent(in)                :: nuY(0:N)
!  linear source term
!  (treated implicitly)
   real(rk), intent(in)                :: Lsour(0:N)
!  constant source term
!  (treated explicitly)
   real(rk), intent(in)                :: Qsour(0:N)
!  relaxation time (s)
   real(rk), intent(in)                :: Taur(0:N)
!  observed value of Y
   real(rk), intent(in)                :: Yobs(0:N)
!  concentranions
   real(rk), intent(in)                :: Y(0:N)

   real(rk),dimension(N):: do_diffusive
  end function
end module
