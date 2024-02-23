!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"
#include "eos.h"

!==================================================================================================================================
!> Subroutines needed by the eddy viscosity model. This is the default  model, which means no eddy viscosity.
!> We need to initialize some arrays that are used in the interfaces to the rounties. The eddy viscosity model itself will
!> return 0 as turbulent viscosity.
!==================================================================================================================================
MODULE MOD_DefaultEddyVisc
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitDefaultEddyVisc
   MODULE PROCEDURE InitDefaultEddyVisc
END INTERFACE

INTERFACE DefaultEddyVisc
   MODULE PROCEDURE DefaultEddyVisc_Point
   MODULE PROCEDURE DefaultEddyVisc_Volume
END INTERFACE

PUBLIC::InitDefaultEddyVisc, DefaultEddyVisc_Volume, FinalizeDefaultEddyViscosity
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Initialize kg model
!===================================================================================================================================
SUBROUTINE InitDefaultEddyVisc()
! MODULES
USE MOD_PreProc,           ONLY: PP_N
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars
IMPLICIT NONE

ALLOCATE(ProdK(0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(DissK(0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(ProdG(0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(DissG(0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(CrossG(0:PP_N,0:PP_N,0:PP_NZ,nElems))
END SUBROUTINE InitDefaultEddyVisc

!===================================================================================================================================
!> Compute the default turbulent eddy viscosity based on k-g RANS model
!===================================================================================================================================
PPURE SUBROUTINE DefaultEddyVisc_Point(gradUx, gradUy, gradUz, U, muSGS)
! MODULES
USE MOD_Equation_Vars,  ONLY: Cmu, epsOMG, epsTKE
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVar)       ,INTENT(IN)  :: U          !> conserved variables
REAL,DIMENSION(2)             ,INTENT(OUT) :: muSGS      !> pointwise eddyviscosity
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL       :: sRho
REAL       :: magS, mutLim, muOrig
!===================================================================================================================================
sRho  = 1.0 / U(DENS)

#if PP_dim==3
magS  = SQRT ( 2.*(gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2) &
             + ( gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2) )**2                    &
             + ( gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3) )**2                    &
             + ( gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3) )**2 )
#else
magS  = SQRT ( 2.*(gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2) &
             + ( gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2) )**2    )
#endif

mutLim = ABS(U(RHOK))/ MAX( SQRT(6.0) * magS, 1.0e-16)

muOrig = Cmu * sRho * MAX(U(RHOK)*sRho, epsTKE) * U(RHOG)**2

muSGS(1) = MIN( mutLim , muOrig ) ! muSGS(1) = rho Cmu k g^2 with limiter
muSGS(2) = muSGS(1) ! muSGS(2) = muSGS(1)

END SUBROUTINE DefaultEddyVisc_Point

!===================================================================================================================================
!> Compute the default turbulent eddy viscosity based on k-g RANS model
!===================================================================================================================================
SUBROUTINE DefaultEddyVisc_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_EddyVisc_Vars,      ONLY: muSGS
USE MOD_Lifting_Vars,       ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,            ONLY: U
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!===================================================================================================================================
DO iElem = 1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      CALL DefaultEddyVisc_Point(gradUx(:,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                                 U(:,i,j,k,iElem),      muSGS(:,i,j,k,iElem) )
    END DO; END DO; END DO ! i,j,k
END DO

END SUBROUTINE DefaultEddyVisc_Volume

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by the default eddy viscosity
!===============================================================================================================================
SUBROUTINE FinalizeDefaultEddyviscosity()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!===============================================================================================================================
SDEALLOCATE(ProdK)
SDEALLOCATE(DissK)
SDEALLOCATE(ProdG)
SDEALLOCATE(DissG)
SDEALLOCATE(CrossG)
END SUBROUTINE FinalizeDefaultEddyViscosity

END MODULE MOD_DefaultEddyVisc
