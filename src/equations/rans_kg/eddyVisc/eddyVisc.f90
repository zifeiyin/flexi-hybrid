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

!===================================================================================================================================
!> In this module the parameters for all the eddy viscosity models will be set, so we only need to call this single routine
!> from flexi.f90.
!===================================================================================================================================
MODULE MOD_EddyVisc
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER      :: EDDYVISCTYPE_NONE     = 0
INTEGER,PARAMETER      :: EDDYVISCTYPE_SMAGO    = 1
INTEGER,PARAMETER      :: EDDYVISCTYPE_VREMAN   = 2
INTEGER,PARAMETER      :: EDDYVISCTYPE_SIGMA    = 3
INTEGER,PARAMETER      :: EDDYVISCTYPE_DYNSMAGO = 4

INTERFACE DefineParametersEddyVisc
  MODULE PROCEDURE DefineParametersEddyVisc
END INTERFACE

PUBLIC:: DefineParametersEddyVisc
PUBLIC:: InitEddyVisc
PUBLIC:: FinalizeEddyVisc
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Define the parameters for the eddy viscosity models.
!===================================================================================================================================
SUBROUTINE DefineParametersEddyVisc()
! MODULES
USE MOD_ReadInTools,        ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("EddyViscParameters")
CALL prms%CreateIntFromStringOption('eddyViscType',  'Type of eddy viscosity. None, Smagorinsky, dynSmagorinsky, Vreman, Sigma',&
                                                     'none')
CALL addStrListEntry(               'eddyViscType',  'none',          EDDYVISCTYPE_NONE    )
CALL addStrListEntry(               'eddyViscType',  'smagorinsky',   EDDYVISCTYPE_SMAGO   )
CALL addStrListEntry(               'eddyViscType',  'dynsmagorinsky',EDDYVISCTYPE_DYNSMAGO)
CALL addStrListEntry(               'eddyViscType',  'vreman',        EDDYVISCTYPE_VREMAN  )
CALL addStrListEntry(               'eddyViscType',  'sigma',         EDDYVISCTYPE_SIGMA   )
CALL prms%CreateIntOption(          'N_testFilter',  'Polynomial degree of test filter (modal cutoff filter).','-1')
CALL prms%CreateRealOption(         'CS',            'EddyViscParameters constant')
CALL prms%CreateRealOption(         'PrSGS',         'Turbulent Prandtl number','0.9')
CALL prms%CreateRealArrayOption(    'eddyViscLimits','Limits for the computed eddy viscosity as multiples of the physical &
                                                     &viscosity.','(/0.,100./)')
CALL prms%CreateLogicalOption(      'VanDriest',     'Van Driest damping, only for channel flow!', '.FALSE.')
CALL prms%CreateStringOption(       'WallDistFile',  'File containing the distances to the nearest walls in the domain.')
END SUBROUTINE DefineParametersEddyVisc

!===================================================================================================================================
!> Initialize eddy viscosity routines
!===================================================================================================================================
SUBROUTINE InitEddyVisc()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_EddyVisc_Vars
USE MOD_DefaultEddyVisc
USE MOD_Smagorinsky
USE MOD_DynSmagorinsky
USE MOD_Vreman
USE MOD_SigmaModel
USE MOD_Mesh_Vars  ,ONLY: nElems,nSides
USE MOD_ReadInTools,ONLY: GETINTFROMSTR, GETREAL
USE MOD_IO_HDF5    ,ONLY: AddToFieldData,FieldOut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
eddyViscType = GETINTFROMSTR('eddyViscType')

! Allocate arrays needed by all SGS models
ALLOCATE(DeltaS(nElems))
DeltaS=0.
ALLOCATE(muSGS(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
muSGS = 0.
ALLOCATE(muSGS_master(1,0:PP_N,0:PP_NZ,nSides))
ALLOCATE(muSGS_slave (1,0:PP_N,0:PP_NZ,nSides))
muSGS_master=0.
muSGS_slave =0.

ALLOCATE(prodK (1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(dissK (1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(prodG (1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(dissG (1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(crossG(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(SijUij(  0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(dGidGi(  0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Turbulent Prandtl number
PrSGS  = GETREAL('PrSGS')

SELECT CASE(eddyViscType)
  CASE(EDDYVISCTYPE_NONE) ! k-g model with sqrt(6) limiter
    CALL InitDefaultEddyVisc()
    ComputeEddyViscosity  => DefaultEddyVisc_Volume
    FinalizeEddyViscosity => FinalizeDefaultEddyViscosity
  CASE(EDDYVISCTYPE_SMAGO)  ! Smagorinsky Model with optional Van Driest damping for channel flow
    CALL InitSmagorinsky()
    ComputeEddyViscosity  => Smagorinsky_Volume
    FinalizeEddyViscosity => FinalizeSmagorinsky
  CASE(EDDYVISCTYPE_DYNSMAGO)  ! Smagorinsky Model with dynamic procedure of Lilly
    CALL InitDynSmagorinsky()
    ComputeEddyViscosity  => DynSmagorinsky_Volume
    FinalizeEddyViscosity => FinalizeDynSmagorinsky
  CASE(EDDYVISCTYPE_VREMAN) ! Vreman Model (Vreman, 2004)
    CALL InitVreman()
    ComputeEddyViscosity  => Vreman_Volume
    FinalizeEddyViscosity => FinalizeVreman
  CASE(EDDYVISCTYPE_SIGMA)  ! sigma Model (Nicoud et al., 2011)
    CALL InitSigmaModel()
    ComputeEddyViscosity  => SigmaModel_Volume
    FinalizeEddyViscosity => FinalizeSigmaModel
  CASE DEFAULT
    CALL CollectiveStop(__STAMP__,&
      'Eddy Viscosity Type not specified!')
END SELECT
CALL AddToFieldData(FieldOut,(/1,PP_N+1,PP_N+1,PP_NZ+1/),'muSGS',(/'muSGS'/),RealArray=muSGS)
END SUBROUTINE

!===================================================================================================================================
!> Finalize eddy viscosity routines
!===================================================================================================================================
SUBROUTINE FinalizeEddyVisc()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(DeltaS)
SDEALLOCATE(CSdeltaS2)
SDEALLOCATE(muSGS)
SDEALLOCATE(muSGS_master)
SDEALLOCATE(muSGS_slave)
SDEALLOCATE(prodK)
SDEALLOCATE(dissK)
SDEALLOCATE(prodG)
SDEALLOCATE(dissG)
SDEALLOCATE(crossG)
SDEALLOCATE(SijUij)
SDEALLOCATE(dGidGi)
IF (ASSOCIATED(FinalizeEddyViscosity)) CALL FinalizeEddyViscosity()
END SUBROUTINE

END MODULE MOD_EddyVisc
