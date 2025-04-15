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

!==================================================================================================================================
!> Contains the (physical) parameters needed for the NS-kg calculation
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

INTEGER           :: IniSourceTerm     !< number identifying the used source term
REAL              :: ConstantBodyForce(3) !< Constant body force to be added, IniSourceTerm==1
REAL              :: ConstantBodyHeat  !< Constant body heat to be added, IniSourceTerm==1
REAL              :: Fluctuation       !< Fluctuation of the constant body force

! Boundary condition arrays
REAL,ALLOCATABLE     :: BCData(:,:,:,:) !< array with precomputed BC values (conservative)
REAL,ALLOCATABLE     :: BCDataPrim(:,:,:,:) !< array with precomputed BC values (primitive)
INTEGER,ALLOCATABLE  :: nBCByType(:)   !< number of sides with specific BC type
INTEGER,ALLOCATABLE  :: BCSideID(:,:)  !< array storing side IDs of sides with different BCs

REAL                 :: s43            !< precomputed 4./3.
REAL                 :: s23            !< precomputed 2./3.

REAL                 :: sqrt6          !< precomputed sqrt(6)

REAL,PARAMETER       :: Cmu = 0.09      !< aka beta*
REAL,PARAMETER       :: COmega1 = 5./9. !< aka alpha
REAL,PARAMETER       :: COmega2 = 0.075 !< aka beta
REAL,PARAMETER       :: invSigmaK = 0.5 !< inverse of sigmaK = 2
REAL,PARAMETER       :: invSigmaG = 0.5 !< inverse of sigmaG = 2

LOGICAL              :: rhokContribution      = .TRUE.!< rhok contribution in Reynolds stress, might be needed in compressible flow
LOGICAL              :: danisDurbinCorrection = .FALSE.!< danis durbin correction
LOGICAL              :: crossDiffusionTerm    = .FALSE.!< cross diffusion term 
LOGICAL              :: RiemannInvariantBC    = .FALSE.!< Riemann invariant BC


CHARACTER(LEN=255),DIMENSION(7),PARAMETER :: StrVarNames =&
  (/ CHARACTER(LEN=255) :: 'Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity','DensityK','DensityG'/) !< conservative variable names
CHARACTER(LEN=255),DIMENSION(8),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Temperature','turbK','turbOmega'/) !< primitive variable names

LOGICAL           :: EquationInitIsDone=.FALSE.
!==================================================================================================================================

END MODULE MOD_Equation_Vars
