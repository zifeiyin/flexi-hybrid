#include "flexi.h"
#include "eos.h"

MODULE MOD_TKE_Reconstruction
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::tke_factor,ComputeTKER

CONTAINS

PPURE FUNCTION tke_factor(nutplus)
USE MOD_Equation_Vars, ONLY: rhokContribution, kReconstruction

REAL,INTENT(IN) :: nutplus

REAL            :: tke_factor

REAL :: f, phi

IF (kReconstruction) THEN
  f   = 1.0 / MAX(1.0E-4, (1.0 - EXP(-MAX(nutplus, 0.0) / 56.0))**0.25)
  phi = TANH(nutplus / 70.0)
  tke_factor = (1.0 - phi) * f + phi * 1.0
ELSE
  tke_factor = 1.0
END IF

END FUNCTION tke_factor

PPURE SUBROUTINE ComputeTKER(UPrim, omega, y)
USE MOD_Equation_Vars, ONLY: kReconstruction
USE MOD_Omega,         ONLY: ComputeOmega2
USE MOD_Viscosity
! USE MOD_EOS_Vars,      ONLY: kappaM1

REAL,INTENT(INOUT)       :: UPrim(PP_nVarPrim)
REAL,INTENT(IN),OPTIONAL :: omega
REAL,INTENT(IN),OPTIONAL :: y

REAL :: omega2,muS,muT

IF (.NOT. kReconstruction) THEN
  RETURN
END IF

IF (PRESENT(omega)) THEN
  omega2 = omega
ELSE
  CALL ComputeOmega2(UPrim, y, omega2)
END IF

muS = VISCOSITY_TEMPERATURE(UPrim(TEMP))
muT = MAX(UPrim(DENS) * UPrim(TKE) / omega2, 0.0)
UPrim(TKER) = MAX(tke_factor(muT / muS) * UPrim(TKE), 0.0)
! UPrim(TKER) = MIN(UPrim(TKER), )

IF (PRESENT(y)) THEN
  IF (y .LE. 1.0E-6) THEN
    UPrim(TKER) = UPrim(TKE)
  END IF
END IF

END SUBROUTINE ComputeTKER

END MODULE MOD_TKE_Reconstruction