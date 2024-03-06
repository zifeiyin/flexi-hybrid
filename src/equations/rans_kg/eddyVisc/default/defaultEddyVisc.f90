!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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

PUBLIC::DefaultEddyVisc, FinalizeDefaultEddyViscosity
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Dummy for default eddy viscosity (meaning no eddy viscosity), do nothing since the muSGS arrays will be passed here and they
!> are zero all the time.
!===================================================================================================================================
SUBROUTINE DefaultEddyVisc()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems
USE MOD_EddyVisc_Vars,     ONLY: muSGS
USE MOD_Lifting_Vars,      ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,           ONLY: U
USE MOD_EOS_Vars,          ONLY: mu0
USE MOD_Equation_Vars,     ONLY: s43,s23,epsTKE,epsOMG,Cmu
IMPLICIT NONE
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: magS
! REAL    :: Sxx, Sxy, Syy
! #if PP_dim == 3
! REAL    :: Sxz, Syz, Szz
! #endif
REAL    :: rho, kPos, gPos, muu
!===================================================================================================================================

DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N

! #if PP_dim == 2
!     Sxx = 0.5 * ( s43 * gradUx(LIFT_VEL1,i,j,k,iElem) - s23 * gradUy(LIFT_VEL2,i,j,k,iElem) )
!     Syy = 0.5 * (-s23 * gradUx(LIFT_VEL1,i,j,k,iElem) + s43 * gradUy(LIFT_VEL2,i,j,k,iElem) )
!     Sxy = 0.5 * (gradUy(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL2,i,j,k,iElem))
! #else
!     Sxx = 0.5 * ( s43 * gradUx(LIFT_VEL1,i,j,k,iElem) - s23 * gradUy(LIFT_VEL2,i,j,k,iElem) - s23 * gradUz(LIFT_VEL3,i,j,k,iElem))
!     Syy = 0.5 * (-s23 * gradUx(LIFT_VEL1,i,j,k,iElem) + s43 * gradUy(LIFT_VEL2,i,j,k,iElem) - s23 * gradUz(LIFT_VEL3,i,j,k,iElem))
!     Szz = 0.5 * (-s23 * gradUx(LIFT_VEL1,i,j,k,iElem) - s23 * gradUy(LIFT_VEL2,i,j,k,iElem) + s43 * gradUz(LIFT_VEL3,i,j,k,iElem))
!     Sxy = 0.5 * (gradUy(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL2,i,j,k,iElem))
!     Sxz = 0.5 * (gradUz(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL3,i,j,k,iElem))
!     Syz = 0.5 * (gradUz(LIFT_VEL2,i,j,k,iElem) + gradUy(LIFT_VEL3,i,j,k,iElem))
! #endif

#if PP_dim==3
    magS = SQRT( &
        2. * (gradUx(LIFT_VEL1,i,j,k,iElem)**2 + gradUy(LIFT_VEL2,i,j,k,iElem)**2 + gradUz(LIFT_VEL3,i,j,k,iElem)**2) + &
        (gradUy(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL2,i,j,k,iElem))**2 + &
        (gradUz(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL3,i,j,k,iElem))**2 + &
        (gradUz(LIFT_VEL2,i,j,k,iElem) + gradUy(LIFT_VEL3,i,j,k,iElem))**2)
#else
    magS = SQRT( &
        2. * (gradUx(LIFT_VEL1,i,j,k,iElem)**2 + gradUy(LIFT_VEL2,i,j,k,iElem)**2) + &
        (gradUy(LIFT_VEL1,i,j,k,iElem) + gradUx(LIFT_VEL2,i,j,k,iElem))**2)
#endif

    rho = U(DENS,i,j,k,iElem)
    kPos = MAX(U(RHOK,i,j,k,iElem) / rho, epsTKE)
    gPos = MAX(U(RHOG,i,j,k,iElem) / rho, epsOMG)
    muu = Cmu * rho * kPos * gPos**2
    muSGS(1,i,j,k,iElem) = MIN(muu, rho * kPos / MAX(SQRT(6.) * magS, epsOMG))
  END DO; END DO; END DO
END DO

END SUBROUTINE DefaultEddyVisc

!===============================================================================================================================
!> Deallocate arrays and finalize variables used by the default eddy viscosity
!===============================================================================================================================
SUBROUTINE FinalizeDefaultEddyviscosity()
! MODULES
IMPLICIT NONE
!===============================================================================================================================
END SUBROUTINE FinalizeDefaultEddyViscosity

END MODULE MOD_DefaultEddyVisc
