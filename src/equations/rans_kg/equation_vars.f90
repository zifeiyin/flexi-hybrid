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

!==================================================================================================================================
!> Contains the (physical) parameters needed for the Navier Stokes calculation
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL           :: doCalcSource      !< automatically set by calcsource itself
INTEGER           :: IniExactFunc      !< number identifying the used exact function
INTEGER           :: IniRefState       !< RefState for initialization (case IniExactFunc=1 only)
INTEGER           :: nRefState         !< number of refstates defined in parameter file
REAL,ALLOCATABLE  :: RefStatePrim(:,:) !< refstates in primitive variables (as read from ini file)
REAL,ALLOCATABLE  :: RefStateCons(:,:) !< refstates in conservative variables
CHARACTER(LEN=255):: BCStateFile       !< file containing the reference solution on the boundary to be used as BC

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:) !< array with precomputed BC values (conservative)
REAL,ALLOCATABLE     :: BCDataPrim(:,:,:,:) !< array with precomputed BC values (primitive)
INTEGER,ALLOCATABLE  :: nBCByType(:)   !< number of sides with specific BC type
INTEGER,ALLOCATABLE  :: BCSideID(:,:)  !< array storing side IDs of sides with different BCs

REAL                 :: s43            !< precomputed 4./3.
REAL                 :: s23            !< precomputed 2./3.

! k-omega specific variables and parameters
REAL              :: PrTurb  = 0.9     !< Turbulent Prandtl number
REAL, PARAMETER   :: sigmaK  = 0.5     !< diffusivity coefficient for k equation
REAL, PARAMETER   :: sigmaG  = 0.5     !< diffusivity coefficient for omega equation
REAL, PARAMETER   :: Comega1 = 5./9.    !< alpha in k omega model 
REAL, PARAMETER   :: Comega2 = 0.075   !< beta in k omega model 
REAL, PARAMETER   :: Cmu     = 0.09    !< Cmu in k omega model
REAL, PARAMETER   :: epsTKE  = 1.e-16
REAL, PARAMETER   :: epsOMG  = 1.e-16

CHARACTER(LEN=255),DIMENSION(7),PARAMETER :: StrVarNames =&
  (/ CHARACTER(LEN=255) :: 'Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity','RhoK','RhoG'/) !< conservative variable names
CHARACTER(LEN=255),DIMENSION(9),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Temperature','TurbK','TurbG','TurbMut'/) !< primitive variable names

LOGICAL           :: EquationInitIsDone=.FALSE.
!==================================================================================================================================

END MODULE MOD_Equation_Vars
