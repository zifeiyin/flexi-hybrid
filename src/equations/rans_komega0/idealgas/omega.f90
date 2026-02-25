#include "flexi.h"
#include "eos.h"

MODULE MOD_Omega
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ComputeOmega
  MODULE PROCEDURE ComputeOmega
END INTERFACE

INTERFACE ComputeOmega2
  MODULE PROCEDURE ComputeOmega2
END INTERFACE

INTERFACE SolveOmega0
  MODULE PROCEDURE SolveOmega0
END INTERFACE

PUBLIC::ComputeOmega
PUBLIC::ComputeOmega2
PUBLIC::SolveOmega0

CONTAINS

SUBROUTINE ComputeOmega()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY: nElems
USE MOD_DG_Vars, ONLY: UPrim
USE MOD_Equation_Vars, ONLY: Comega2, Cexp
USE MOD_EddyVisc_Vars, ONLY: ywall,omega
USE MOD_Lifting_Vars     ,ONLY: gradUx,gradUy,gradUz
USE MOD_Viscosity

INTEGER :: iElem,i,j,k
REAL    :: nuS, xi, omega0Pos, magS

DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    nuS = VISCOSITY_TEMPERATURE(UPrim(TEMP,i,j,k,iElem)) / UPrim(DENS,i,j,k,iElem)
    IF (UPrim(OMG,i,j,k,iElem) .GE. 0.0) THEN
      omega0Pos = MAX(UPrim(OMG,i,j,k,iElem),0.0)
      xi = omega0Pos * ywall(i,j,k,0,iElem)**2 / nuS
      omega(i,j,k,iElem) = omega0Pos + (6.0 / Comega2) * nuS / ywall(i,j,k,0,iElem)**2 * EXP(-Cexp * xi)
    ELSE 
      omega0Pos = SQRT( &
        2. * (gradUx(LIFT_VEL1,i,j,k,iElem)**2 + gradUy(LIFT_VEL2,i,j,k,iElem)**2 + gradUz(LIFT_VEL3,i,j,k,iElem)**2) + &
             (gradUy(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL2,i,j,k,iElem))**2 + &
             (gradUz(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL3,i,j,k,iElem))**2 + &
             (gradUz(LIFT_VEL2,i,j,k,iElem) + gradUy(LIFT_VEL3,i,j,k,iElem))**2) * SQRT(6.0)
      xi = omega0Pos * ywall(i,j,k,0,iElem)**2 / nuS
      omega(i,j,k,iElem) = omega0Pos + (6.0 / Comega2) * nuS / ywall(i,j,k,0,iElem)**2 * EXP(-Cexp * xi)
    END IF
    IF (ywall(i,j,k,0,iElem).LE.1.0E-16) THEN
      ! Here we assume y1 >= 1.0E-6
      omega(i,j,k,iElem) = (60.0 / Comega2) * nuS / 1.0E-6**2
    END IF
  END DO; END DO; END DO
END DO

END SUBROUTINE ComputeOmega

SUBROUTINE ComputeOmega2(UPrim, ywall, omega)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: Comega2, Cexp
USE MOD_Viscosity

REAL,INTENT(IN)  :: UPrim(PP_nVarPrim)
REAL,INTENT(IN)  :: ywall
REAL,INTENT(OUT) :: omega

REAL    :: nuS, xi, omega0Pos

nuS = VISCOSITY_TEMPERATURE(UPrim(TEMP)) / UPrim(DENS)
omega0Pos = MAX(UPrim(OMG), 0.0)
xi = omega0Pos * ywall**2 / nuS
omega = omega0Pos + (6.0 / Comega2) * nuS / ywall**2 * EXP(-Cexp * xi)

END SUBROUTINE ComputeOmega2

SUBROUTINE SolveOmega0(UPrim, ywall, omega0)
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars, ONLY: Comega2, Cexp
USE MOD_Viscosity

REAL,INTENT(IN)  :: UPrim(PP_nVarPrim)
REAL,INTENT(IN)  :: ywall
REAL,INTENT(OUT) :: omega0

REAL    :: nuS, xi, omega, omega_asymptotic
REAL    :: upper, lower, mid

! UPrim = (rho, u, v, w, p, T, k, omega) instead of omega0

nuS = VISCOSITY_TEMPERATURE(UPrim(TEMP)) / UPrim(DENS)

omega_asymptotic = (6.0 / Comega2) * nuS / ywall**2

IF (UPrim(OMG) .LE. omega_asymptotic) THEN
  omega0 = 0.0
  RETURN
END IF

! lower = 0.0
! upper = UPrim(OMG)
! DO
!   mid = 0.5 * (lower + upper)
!   xi = mid * ywall**2 / nuS
!   omega = mid + omega_asymptotic * EXP(-Cexp * xi)
!   IF ((ABS(omega - UPrim(OMG)) / UPrim(OMG) .LE. 1.0E-2).OR.&
!       ((upper - lower) / mid .LE. 1.0E-2)) THEN
!     EXIT
!   ELSEIF (omega .GT. UPrim(OMG)) THEN
!     upper = mid
!   ELSE
!     lower = mid
!   END IF
! END DO

! omega0 = mid

! omega0 = MAX(0.0, UPrim(OMG) - omega_asymptotic)
omega0 = UPrim(OMG)
DO
  IF ((ABS(omega - UPrim(OMG)) / UPrim(OMG) .LE. 1.0E-2)) THEN
    EXIT
  END IF
  xi = omega0 * ywall**2 / nuS
  omega = omega0 + omega_asymptotic * EXP(-Cexp * xi)
  omega0 = omega0 - &
      (omega - UPrim(OMG)) / &
      (1.0 - (6.0 / Comega2 * Cexp) * EXP(-Cexp * xi))
END DO

END SUBROUTINE SolveOmega0

END MODULE MOD_Omega