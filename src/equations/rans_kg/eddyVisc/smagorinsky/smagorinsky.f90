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

!===================================================================================================================================
!> Subroutines necessary for calculating Smagorinsky Eddy-Viscosity, originally derived in
!>   - Smagorinsky, Joseph. "General circulation experiments with the primitive equations: I. The basic experiment."
!>     Monthly weather review 91.3 (1963): 99-164.
!>
!> The Van-Driest damping for the Smagorinsky model for channel flow is originally used in
!>   - Moin, Parviz, and John Kim. "Numerical investigation of turbulent channel flow."
!>     Journal of fluid mechanics 118 (1982): 341-377.
!===================================================================================================================================
MODULE MOD_Smagorinsky
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitSmagorinsky
   MODULE PROCEDURE InitSmagorinsky
END INTERFACE

INTERFACE Smagorinsky
   MODULE PROCEDURE Smagorinsky_Point
   MODULE PROCEDURE Smagorinsky_Volume
END INTERFACE

INTERFACE FinalizeSmagorinsky
   MODULE PROCEDURE FinalizeSmagorinsky
END INTERFACE

PUBLIC::InitSmagorinsky, Smagorinsky_Volume, FinalizeSmagorinsky
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Get model parameters and initialize Smagorinsky model
!===================================================================================================================================
SUBROUTINE InitSmagorinsky()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL,GETLOGICAL
USE MOD_ReadInTools        ,ONLY: CountOption,GETREALARRAY,GETSTR,GETREAL,GETLOGICAL
USE MOD_HDF5_Input         ,ONLY: ReadArray,OpenDataFile,CloseDataFile,GetDataSize,ReadAttribute
USE MOD_IO_HDF5
USE MOD_2D                 ,ONLY: ExpandArrayTo3D
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,wGP
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,nElems,sJ,Elem_xGP,MeshFile,offsetElem
USE MOD_EOS_Vars           ,ONLY: mu0
USE MOD_DG_Vars            ,ONLY: D
USE MOD_Mesh_Vars          ,ONLY: sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Interpolation      ,ONLY: GetNodesAndWeights
USE MOD_Interpolation_Vars ,ONLY: NodeType
USE MOD_Basis              ,ONLY: PolynomialDerivativeMatrix,LagrangeInterpolationPolys,PolynomialMassMatrix
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: file_exists
INTEGER             :: i,j,k,iElem
INTEGER             :: HSize_proc(4)
CHARACTER(LEN=255)  :: FileName
REAL                :: CellVol
REAL,ALLOCATABLE    :: yWall_local(:,:,:,:)
INTEGER             :: l
REAL,ALLOCATABLE    :: xRef(:), D(:,:)
!===================================================================================================================================
IF (PP_dim.EQ.2) THEN
  CALL CollectiveStop(__STAMP__,"Smagorinsky cannot run in 2D.")
END IF

IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.SmagorinskyInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitSmagorinsky not ready to be called or already called.")
END IF
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SMAGORINSKY...'

ALLOCATE(xRef(0:PP_N))
ALLOCATE(D(0:PP_N,0:PP_N))

CALL GetNodesAndWeights(PP_N,NodeType,xRef)
CALL PolynomialDerivativeMatrix(PP_N,xRef,D)

! allocate memory
ALLOCATE(yWall(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
ALLOCATE(gradyWall(3,0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
ALLOCATE(fd(0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Read in a file that contains the distance to the nearest wall for each solution point (generated by POSTI tool).
yWall = SQRT(HUGE(1.))
FileName = MeshFile(1:INDEX(MeshFile,'_mesh.h5')-1)//'_walldistance.h5'
file_exists = FILEEXISTS(FileName)
IF (file_exists) THEN
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  CALL GetDataSize(File_ID,'walldistance',nDims,HSize)
  IF (HSize(1)-1.NE.PP_N) CALL Abort(__STAMP__,"Polynomial degree of walldistance file does not match!")
  ALLOCATE(yWall_local(0:HSize(1)-1,0:HSize(2)-1,0:HSize(3)-1,nElems))
  HSize_proc = INT(HSize)
  HSize_proc(4) = nElems
  CALL ReadArray('walldistance',4,&
                  HSize_proc,&
                  offsetElem,4,RealArray=yWall_local)

  IF (HSize(3).EQ.1) THEN
    ! Walldistance was created by 2D tool, expand in third dimension
    CALL ExpandArrayTo3D(4,(/PP_N+1,PP_N+1,1,nElems/),3,PP_N+1,yWall_local,yWall(:,:,:,0,:))
  ELSE
    ! 3D walldistance tool and 3D Flexi
    yWall(:,:,:,0,:) = yWall_local
  END IF

  DEALLOCATE(HSize)
  DEALLOCATE(yWall_local)
  CALL CloseDataFile()

  gradyWall = 0.0
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      DO l=0,PP_N
        gradyWall(1,i,j,k,0,iElem) = gradyWall(1,i,j,k,0,iElem) + D(i,l) * ywall(l,j,k,0,iElem)
        gradyWall(2,i,j,k,0,iElem) = gradyWall(2,i,j,k,0,iElem) + D(j,l) * ywall(i,l,k,0,iElem)
        gradyWall(3,i,j,k,0,iElem) = gradyWall(3,i,j,k,0,iElem) + D(k,l) * ywall(i,j,l,0,iElem)
      END DO
      gradyWall(:,i,j,k,0,iElem) = ( &
        gradyWall(1,i,j,k,0,iElem) * Metrics_fTilde(:,i,j,k,iElem,0) + &
        gradyWall(2,i,j,k,0,iElem) * Metrics_gTilde(:,i,j,k,iElem,0) + &
        gradyWall(3,i,j,k,0,iElem) * Metrics_hTilde(:,i,j,k,iElem,0) &
      ) * sJ(i,j,k,iElem,0)
    END DO; END DO; END DO
  END DO
ELSE
  SWRITE(UNIT_stdOut, *) "WARNING: No walldistance file found! DDES shielding not working!"
  CALL CollectiveStop(__STAMP__,'Please use POSTI to compute wall distance first!')
ENDIF

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol + wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
  DeltaS(iElem) = CellVol**(1./3.)  / (REAL(PP_N)+1.)
END DO

SmagorinskyInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SMAGORINSKY k-g DDES DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitSmagorinsky


!===================================================================================================================================
!> Compute Smagorinsky Eddy-Visosity
!===================================================================================================================================
PPURE SUBROUTINE Smagorinsky_Point(gradUx,gradUy,gradUz,UPrim,y,Delta,hmax,muSGS,fd)
! MODULES
USE MOD_Equation_Vars,  ONLY: Cmu, sqrt6
USE MOD_EddyVisc_Vars,  ONLY: CDES0
USE MOD_VISCOSITY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim)   ,INTENT(IN)  :: UPrim   !> pointwise conserved variable
REAL                          ,INTENT(IN)  :: y       !> pointwise wall distance
REAL                          ,INTENT(IN)  :: Delta   !> pointwise cell spacing
REAL                          ,INTENT(IN)  :: hmax    !> local grid size
REAL                          ,INTENT(OUT) :: muSGS   !> pointwise eddyviscosity
REAL                          ,INTENT(OUT) :: fd      !> fd
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: sRho
REAL                                    :: magS
REAL                                    :: kPos, gPos, muTOrig, muTLim
REAL                                    :: lLES, lRANS, lDDES, rd
REAL                                    :: dUdU
REAL                                    :: muS
REAL,PARAMETER                          :: kappa=0.41
!===================================================================================================================================
sRho  = 1.0 / UPrim(DENS)
muS   = VISCOSITY_PRIM(UPrim)

#if PP_dim==2
magS = SQRT( &
    2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2) + &
    (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))**2)
#else
magS = SQRT( &
    2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2) + &
    (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))**2 + &
    (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))**2 + &
    (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))**2)
#endif

kPos    = MAX( UPrim(TKE), 1.e-16 )
gPos    = MAX( UPrim(OMG), 1.e-16 )
muTOrig = Cmu * UPrim(DENS) * kPos * gPos**2

muTLim = MIN(muTOrig,  UPrim(DENS) * kPos / MAX(sqrt6 * magS, 1.e-16))

lRANS = SQRT(mutLim * sRho * Cmu * gPos**2 )

dUdU  = gradUx(LIFT_VEL1)**2 + gradUx(LIFT_VEL2)**2 + gradUx(LIFT_VEL3)**2 &
      + gradUy(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUy(LIFT_VEL3)**2 &
      + gradUz(LIFT_VEL1)**2 + gradUz(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2

rd    = ( Cmu * kPos * gPos**2 + muS * sRho) / ((kappa * y)**2 * SQRT(MAX(1.e-16, dUdU)))
fd    = 1.0 - TANH((8.0 * rd)**3.0)

lLES = CDES0 * (fd * Delta + (1. - fd) * hmax)

lDDES = lRANS - fd * MAX(0., lRANS-lLES)

muSGS = UPrim(DENS) * lDDES**2 / MAX( Cmu * gPos**2, 1.e-16 )

END SUBROUTINE Smagorinsky_Point


!===================================================================================================================================
!> Compute Smagorinsky Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE Smagorinsky_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems, Elem_hmx
USE MOD_EddyVisc_Vars,     ONLY: muSGS, yWall, DeltaS, fd
USE MOD_Lifting_Vars,      ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,           ONLY: UPrim
!USE MOD_FV_Vars,           ONLY: FV_Elems
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!===================================================================================================================================
DO iElem = 1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CALL Smagorinsky_Point( gradUx(:,i,j,k,iElem),    gradUy(:,i,j,k,iElem),              gradUz(:,i,j,k,iElem), &
                            UPrim(:,i,j,k,iElem),     yWall(i,j,k,0,iElem),               DeltaS(iElem),         &
                            Elem_hmx(iElem),          muSGS(1,i,j,k,iElem),               fd(i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE Smagorinsky_Volume


!===============================================================================================================================
!> Deallocate arrays and finalize variables used by Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE FinalizeSmagorinsky()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===============================================================================================================================
SDEALLOCATE(yWall)
SDEALLOCATE(fd)
SmagorinskyInitIsDone = .FALSE.

END SUBROUTINE FinalizeSmagorinsky

END MODULE MOD_Smagorinsky
