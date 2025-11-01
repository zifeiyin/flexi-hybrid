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

PUBLIC::ComputeOmega

CONTAINS

SUBROUTINE ComputeOmega()
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY: nElems
USE MOD_DG_Vars, ONLY: UPrim
USE MOD_Equation_Vars, ONLY: Comega2
USE MOD_EddyVisc_Vars, ONLY: ywall,omega
USE MOD_Viscosity

INTEGER :: iElem,i,j,k
REAL    :: nuS, xi

DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    nuS = VISCOSITY_TEMPERATURE(UPrim(TEMP,i,j,k,iElem)) / UPrim(DENS,i,j,k,iElem)
    xi = MAX(UPrim(OMG,i,j,k,iElem), 0.0) * ywall(i,j,k,0,iElem)**2 / nuS
    ! omega(i,j,k,iElem) = (6.0 / Comega2) * nuS / ywall(i,j,k,0,iElem)**2 * (xi / (6.0 / Comega2) + EXP(-xi / 95.0))
    omega(i,j,k,iElem) = UPrim(OMG,i,j,k,iElem) + (6.0 / Comega2) * nuS / ywall(i,j,k,0,iElem)**2 * EXP(-xi / 95.0)
  END DO; END DO; END DO
END DO

END SUBROUTINE ComputeOmega

END MODULE MOD_Omega