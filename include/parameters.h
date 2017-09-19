#define _LINE_      write(*,*) repeat('*',79)

#define _PI_        3.141592653589793_rk
#define _DIRICHLET_ 0
#define _NEUMANN_   1

!light parameters:
#define _ICE_ALBEDO_       0.744_rk
#define _SNOW_ALBEDO_      0.9_rk
#define _ICE_SCATTERED_    0.97_rk
#define _ICE_EXTINCTION_   0.93_rk
#define _SNOW_EXTINCTION_  4.3_rk
#define _ERLOV_            0.05_rk
!-Factor to convert downwelling shortwave in air to scalar PAR in water(default = 0.5)
!-Radiative transfer models suggest an average value ~0.5 but with ~10 % variability
!-at mid / high latitudes depending on season, latitude, and wind speed
!-see Mobley and Boss(2012), Figs. 5b, 8b.
#define _PAR_PART_        0.5_rk

!porosity parameters:
!(Soetaert etc., 1996)
#define _MAX_POROSITY_    0.95_rk
#define _MIN_POROSITY_    0.80_rk
#define _POROSITY_DECAY_  0.04_rk
#define _BURIAL_VELOCITY_ 1.e-10_rk

!molecular diffusivity parameters:
#define _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_ 1.e-9_rk 
!(Boudreau, 1997)
#define _RELATIVE_DYNAMIC_VISCOSITY_               0.94_rk
!bioturbation diffusivity parameters:
#define _MIXED_LAYER_DEPTH_                        0.02_rk
#define _MAX_BIOTURBATION_DIFFUSIVITY_             1.e-11_rk
#define _DECAY_BIOTURBATION_SCALE_                 0.01_rk
!ice gravity drainage
#define _GRAVITY_DRAINAGE_                         1.e-8_rk
!Half-saturation constant for the effect of oxygen on bioturbation
!and bioirrigation [uM] (default = 5.0 uM)
#define _KO2_                                      5._rk

!horizontal flux rate
#define _HMIX_RATE_								   0.00333_rk
