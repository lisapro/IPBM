#define _LINE_      write(*,*) repeat('*',79)

#define _PI_        3.141592653589793_rk
#define _DIRICHLET_ 0
#define _NEUMANN_   1

!light parameters:
#define _ERLOV_     0.05_rk

!porosity parameters:
!(Soetaert etc., 1996)
#define _MAX_POROSITY_    0.95_rk
#define _MIN_POROSITY_    0.80_rk
#define _POROSITY_DECAY_  0.04_rk

!molecular diffusivity parameters:
#define _INFINITE_DILLUTION_MOLECULAR_DIFFUSIVITY_ 1.e-9_rk 
!(Boudreau, 1997)
#define _RELATIVE_DYNAMIC_VISCOSITY_               0.94_rk
!bioturbation diffusivity parameters:
#define _MIXED_LAYER_DEPTH_                        0.02_rk
#define _MAX_BIOTURBATION_DIFFUSIVITY_             1.e-11_rk
#define _DECAY_BIOTURBATION_SCALE_                 0.01_rk
