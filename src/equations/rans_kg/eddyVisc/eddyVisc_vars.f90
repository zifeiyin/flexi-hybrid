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

!===================================================================================================================================
!> Contains the parameters needed for eddy viscosity models
!===================================================================================================================================
MODULE MOD_EddyVisc_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

ABSTRACT INTERFACE
  SUBROUTINE EddyViscInt()
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE FinalizeEddyViscosityInt()
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                                      :: eddyViscType          !< type of eddy viscosity
PROCEDURE(EddyViscInt),POINTER               :: ComputeEddyViscosity  !< pointer to routine for computing volume eddy viscosity
PROCEDURE(FinalizeEddyViscosityInt),POINTER  :: FinalizeEddyViscosity !< pointer tofinalize routine

INTEGER,ALLOCATABLE  :: averageType(:)                !< type of averaging for dynamic Smagorinsky model

LOGICAL,ALLOCATABLE  :: doFilterDir(:,:)              !< specifies in which directions the test filter should be applied in each element

REAL,ALLOCATABLE  :: yWall(:,:,:,:,:)                 !< wall distance
REAL,ALLOCATABLE  :: gradyWall(:,:,:,:,:,:)           !< gradient of wall distance
REAL,ALLOCATABLE  :: damp(:,:,:,:,:)                  !< damping factor
REAL,ALLOCATABLE  :: Cdes2(:,:,:,:)                   !< precomputed (model constant)**2 => dynamic Smagorinsky DDES
REAL,ALLOCATABLE  :: IntElem(:,:,:,:)                 !< integration weights for dynamic Smagorinsky model
                                                      !< for dynamic Smagorinsky model

REAL,ALLOCATABLE  :: DeltaS(:)                        !< filter width
REAL,ALLOCATABLE  :: DeltaRatio(:)                    !< test-filter width / filter width
REAL,ALLOCATABLE  :: CSdeltaS2(:)                     !< precomputed (model constant*filter width)**2 => Vreman,Sigma model
REAL,ALLOCATABLE  :: muSGS(:,:,:,:,:)                 !< Sub-grid eddy viscosity
REAL,ALLOCATABLE  :: muSGS_master(:,:,:,:)            !< Sub-grid eddy viscosity on master sides
REAL,ALLOCATABLE  :: muSGS_slave (:,:,:,:)            !< Sub-grid eddy viscosity on slave sides
REAL,ALLOCATABLE  :: muSGSmax(:)                      !< maximum eddy viscosity per element
REAL,ALLOCATABLE  :: FilterMat_TestFilter(:,:)        !< filter matrix for test filter for dynamic Smagorinsky model
REAL              :: muSGS_limits(2)                  !< allowed range of eddy viscosity as multiple of physical viscosit
REAL              :: CS                               !< Model coefficient for eddy viscosity models
REAL,PARAMETER    :: CDES0 = 0.12                     !< Model coefficient for DDES
REAL              :: PrSGS                            !< turbulent Prandtl number for the sub-grid scales

LOGICAL           :: VanDriest=.FALSE.                !< Logical indicating if Van Driest damping is activated (only use for channel flow)
LOGICAL           :: SmagorinskyInitIsDone=.FALSE.    !< Logical indicating if Smagorinsky model has been initialized
LOGICAL           :: DynSmagorinskyInitIsDone=.FALSE. !< Logical indicating if Smagorinsky model has been initialized
LOGICAL           :: VremanInitIsDone=.FALSE.         !< Logical indicating if Vreman model has been initialized
LOGICAL           :: SigmaModelInitIsDone=.FALSE.     !< Logical indicating if sigma model has been initialized

REAL,ALLOCATABLE  :: prodK (:,:,:,:,:)                !< production K
REAL,ALLOCATABLE  :: dissK (:,:,:,:,:)                !< dissipation K
REAL,ALLOCATABLE  :: prodG (:,:,:,:,:)                !< production G
REAL,ALLOCATABLE  :: dissG (:,:,:,:,:)                !< dissipation G
REAL,ALLOCATABLE  :: crossG(:,:,:,:,:)                !< cross G
REAL,ALLOCATABLE  :: fd(:,:,:,:)                      !< fd in DDES

REAL,ALLOCATABLE  :: SijUij(:,:,:,:)                  !< magS
REAL,ALLOCATABLE  :: dGidGi(:,:,:,:)                  !< dGdG

END MODULE MOD_EddyVisc_Vars
