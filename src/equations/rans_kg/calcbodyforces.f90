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
!> Compute and integrate force from fluid onto wall surfaces (e.g. adiabatic, isothermal, Euler walls)
!==================================================================================================================================
MODULE MOD_CalcBodyForces
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CalcBodyForces
  MODULE PROCEDURE CalcBodyForces
END INTERFACE

INTERFACE CalcWallFluxes
  MODULE PROCEDURE CalcWallFluxes
END INTERFACE

PUBLIC :: CalcBodyForces, CalcWallFluxes
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Control routine for CalcBodyforces
!==================================================================================================================================
SUBROUTINE CalcBodyForces(BodyForce,Fp,Fv)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,         ONLY:UPrim_master
#if PARABOLIC
USE MOD_Lifting_Vars,    ONLY:gradUx_master,gradUy_master,gradUz_master
#endif
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem,nBCSides,BC,nBCs,Face_xGP
USE MOD_AnalyzeEquation_Vars,ONLY:isWall,doCalcBodyForcesDOF
USE MOD_TimeDisc_Vars,   ONLY:t
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)               :: Fp(3,nBCs)              !< integrated pressure force per wall BC
REAL,INTENT(OUT)               :: Fv(3,nBCs)              !< integrated friction force per wall BC
REAL,INTENT(OUT)               :: BodyForce(3,nBCs)       !< Sum of pressure/friction force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: Fp_loc(3)
REAL,ALLOCATABLE               :: Fpa(:,:,:)
#if PARABOLIC
REAL                           :: Fv_loc(3)
REAL,ALLOCATABLE               :: Fva(:,:,:)
#endif
INTEGER                        :: SideID,iBC
INTEGER                        :: i,j
INTEGER                        :: OutputFileID
CHARACTER(LEN=64)              :: OutputFileName
CHARACTER(LEN=64)              :: FormatStr
#if USE_MPI
REAL                           :: Box(6,nBCs)
#endif /*USE_MPI*/
!==================================================================================================================================
! Calculate body forces  ! Attention: during the initialization phase no face data / gradients available!

Fp=0.
Fv=0.
BodyForce=0.
DO SideID=1,nBCSides
  iBC=BC(SideID)
  IF(.NOT.isWall(iBC)) CYCLE
  ! Calculate pressure force (Euler wall / Navier-Stokes wall)
  CALL CalcPressureForce(Fp_loc,UPrim_master(PRES,:,:,SideID),SurfElem(:,:,0,SideID),NormVec(:,:,:,0,SideID))
  Fp(:,iBC)=Fp(:,iBC)+Fp_loc
#if PARABOLIC
  ! Calculate viscous force (Navier-Stokes wall)
  CALL CalcViscousForce(Fv_loc,                      &
                        UPrim_master(:,:,:,SideID),  &
                        gradUx_master(:,:,:,SideID), &
                        gradUy_master(:,:,:,SideID), &
                        gradUz_master(:,:,:,SideID), &
                        SurfElem(:,:,0,SideID),      &
                        NormVec(:,:,:,0,SideID))
  Fv(:,iBC)=Fv(:,iBC)+Fv_loc
#endif
END DO

#if USE_MPI
Box(1:3,1:nBCs)=Fv; Box(4:6,1:nBCs)=Fp
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,Box,6*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  Fv=Box(1:3,1:nBCs); Fp=Box(4:6,1:nBCs)
  BodyForce=Fv+Fp
ELSE
  CALL MPI_REDUCE(Box         ,0  ,6*nBCs,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#else
BodyForce=Fv+Fp
#endif

IF (doCalcBodyForcesDOF) THEN
  ALLOCATE(Fpa(3,0:PP_N,0:PP_NZ))
  ALLOCATE(Fva(3,0:PP_N,0:PP_NZ))
  Fpa = 0.
  Fva = 0.

  WRITE(OutputFileName, '(A, "_", I4.4, ".csv")') trim(TIMESTAMP("bodyforce", t)), myRank
  OPEN(UNIT=OutputFileID, FILE=trim(OutputFileName), STATUS='UNKNOWN', ACTION='WRITE')

  ! #if !PARABOLIC
  ! WRITE(OutputFileID, *) "x,y,z,px,py,pz"
  ! #else
  ! WRITE(OutputFileID, *) "x,y,z,px,py,pz,tx,ty,tz"
  ! #endif

  DO SideID=1,nBCSides
    iBC=BC(SideID)
    IF(.NOT.isWall(iBC)) CYCLE
    DO j=0,PP_NZ; DO i=0,PP_N
      Fpa(:,i,j)=UPrim_master(PRES,i,j,SideID)*NormVec(:,i,j,0,SideID)
    END DO; END DO
#if PARABOLIC
      CALL CalcViscousForceDOF(Fva(:,:,:),                  &
                               UPrim_master(:,:,:,SideID),  &
                               gradUx_master(:,:,:,SideID), &
                               gradUy_master(:,:,:,SideID), &
                               gradUz_master(:,:,:,SideID), &
                               SurfElem(:,:,0,SideID),      &
                               NormVec(:,:,:,0,SideID))
#endif

    DO j=0,PP_NZ; DO i=0,PP_N
#if !PARABOLIC
    WRITE(formatStr,'(A10,I2,A18)')'(E23.14E5,',5,'(",",1X,E23.14E5))'
    ! WRITE(FormatStr, '(A,5A,A)') "(", 'F17.9, ",",', 'F17.9)'
    WRITE(OutputFileID, *) &
        Face_xGP(:,i,j,0,SideID), Fpa(:,i,j)
        ! Face_xGP(1,i,j,0,SideID), Face_xGP(2,i,j,0,SideID), Face_xGP(3,i,j,0,SideID), &
        ! Fpa(1,i,j), Fpa(2,i,j), Fpa(3,i,j)
#else
    WRITE(formatStr,'(A10,I2,A18)')'(E23.14E5,',8,'(",",1X,E23.14E5))'
    ! WRITE(FormatStr, '(A,8A,A)') "(", 'F17.9, ",",', 'F17.9)'
    WRITE(OutputFileID, *) &
        Face_xGP(:,i,j,0,SideID), Fpa(:,i,j), Fva(:,i,j)
        ! Face_xGP(1,i,j,0,SideID), Face_xGP(2,i,j,0,SideID), Face_xGP(3,i,j,0,SideID), &
        ! Fpa(1,i,j), Fpa(2,i,j), Fpa(3,i,j), Fva(1,i,j), Fva(2,i,j), Fva(3,i,j)
#endif
    END DO; END DO
  END DO

  CLOSE(OutputFileID)

  SDEALLOCATE(Fpa)
  SDEALLOCATE(Fva)
END IF

END SUBROUTINE CalcBodyForces



!==================================================================================================================================
!> Compute integral pressure force per face
!==================================================================================================================================
SUBROUTINE CalcPressureForce(Fp,p_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Analyze_Vars,      ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: p_Face(0:PP_N,0:PP_NZ)        !< (IN) pressure on face
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_NZ)      !< (IN) face surface
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_NZ)     !< (IN) face normal vectors
REAL, INTENT(OUT)              :: Fp(3)                        !< (OUT) integrated pressure force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: dA
INTEGER                        :: i, j
!==================================================================================================================================
Fp=0.
DO j=0,PP_NZ; DO i=0,PP_N
  dA=wGPSurf(i,j)*SurfElem(i,j)
  Fp=Fp+p_Face(i,j)*NormVec(:,i,j)*dA
END DO; END DO
END SUBROUTINE CalcPressureForce


#if PARABOLIC
!==================================================================================================================================
!> Compute integral viscous force per face (only if compiled with parabolic terms)
!==================================================================================================================================
SUBROUTINE CalcViscousForce(Fv,UPrim_Face,gradUx_Face,gradUy_Face,gradUz_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Analyze_Vars, ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: UPrim_Face( PP_nVarPrim,0:PP_N,0:PP_NZ)    !< (IN) primitive solution on face
REAL, INTENT(IN)               :: gradUx_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients x-dir on face
REAL, INTENT(IN)               :: gradUy_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients y-dir on face
REAL, INTENT(IN)               :: gradUz_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients z-dir on face
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_NZ)                   !< (IN) face surface
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_NZ)                  !< (IN) face normal vectors
REAL, INTENT(OUT)              :: Fv(3)                                      !< (OUT) integrated pressure force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: tau(3,3)                  ! Viscous stress tensor
REAL                           :: muS
REAL                           :: GradV(3,3),DivV,prim(PP_nVarPrim)
INTEGER                        :: i, j
!==================================================================================================================================
Fv       =0.

DO j=0,PP_NZ; DO i=0,PP_N
  ! calculate viscosity
  prim = UPrim_Face(:,i,j)
  muS=VISCOSITY_PRIM(prim)

  ! velocity gradients
  GradV(:,1)=gradUx_Face(LIFT_VELV,i,j)
  GradV(:,2)=gradUy_Face(LIFT_VELV,i,j)
#if PP_dim==3
  GradV(:,3)=gradUz_Face(LIFT_VELV,i,j)
#else
  GradV(:,3)=0.
#endif

  ! Velocity divergence
  DivV=GradV(1,1)+GradV(2,2)+GradV(3,3)
  ! Calculate shear stress tensor
  tau=muS*(GradV + TRANSPOSE(GradV))
  tau(1,1)=tau(1,1)-2./3.*muS*DivV
  tau(2,2)=tau(2,2)-2./3.*muS*DivV
#if PP_dim==3
  tau(3,3)=tau(3,3)-2./3.*muS*DivV
#endif
  ! Calculate viscous force vector
  Fv=Fv+MATMUL(tau,NormVec(:,i,j))*wGPSurf(i,j)*SurfElem(i,j)
END DO; END DO

Fv=-Fv  ! Change direction to get the force acting on the wall

END SUBROUTINE CalcViscousForce

!==================================================================================================================================
!> Compute integral viscous force per DOF
!==================================================================================================================================
SUBROUTINE CalcViscousForceDOF(Fv,UPrim_Face,gradUx_Face,gradUy_Face,gradUz_Face,SurfElem,NormVec)
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Analyze_Vars, ONLY:wGPSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN)               :: UPrim_Face( PP_nVarPrim,0:PP_N,0:PP_NZ)    !< (IN) primitive solution on face
REAL, INTENT(IN)               :: gradUx_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients x-dir on face
REAL, INTENT(IN)               :: gradUy_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients y-dir on face
REAL, INTENT(IN)               :: gradUz_Face(PP_nVarLifting,0:PP_N,0:PP_NZ) !< (IN) sln. gradients z-dir on face
REAL, INTENT(IN)               :: SurfElem(0:PP_N,0:PP_NZ)                   !< (IN) face surface
REAL, INTENT(IN)               :: NormVec(3,0:PP_N,0:PP_NZ)                  !< (IN) face normal vectors
REAL, INTENT(OUT)              :: Fv(3,0:PP_N,0:PP_NZ)                       !< (OUT) viscous force
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: tau(3,3)                  ! Viscous stress tensor
REAL                           :: muS
REAL                           :: GradV(3,3),DivV,prim(PP_nVarPrim)
INTEGER                        :: i, j
!==================================================================================================================================
Fv       =0.

DO j=0,PP_NZ; DO i=0,PP_N
  ! calculate viscosity
  prim = UPrim_Face(:,i,j)
  muS=VISCOSITY_PRIM(prim)

  ! velocity gradients
  GradV(:,1)=gradUx_Face(LIFT_VELV,i,j)
  GradV(:,2)=gradUy_Face(LIFT_VELV,i,j)
#if PP_dim==3
  GradV(:,3)=gradUz_Face(LIFT_VELV,i,j)
#else
  GradV(:,3)=0.
#endif

  ! Velocity divergence
  DivV=GradV(1,1)+GradV(2,2)+GradV(3,3)
  ! Calculate shear stress tensor
  tau=muS*(GradV + TRANSPOSE(GradV))
  tau(1,1)=tau(1,1)-2./3.*muS*DivV
  tau(2,2)=tau(2,2)-2./3.*muS*DivV
#if PP_dim==3
  tau(3,3)=tau(3,3)-2./3.*muS*DivV
#endif
  ! Calculate viscous force vector
  Fv(:,i,j) = -MATMUL(tau,NormVec(:,i,j))
END DO; END DO

END SUBROUTINE CalcViscousForceDOF
#endif /*PARABOLIC*/

!==================================================================================================================================
!> Control routine for CalcWallFluxes
!==================================================================================================================================
SUBROUTINE CalcWallFluxes()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars,         ONLY:UPrim_master
USE MOD_Lifting_Vars,    ONLY:gradUx_master,gradUy_master,gradUz_master
USE MOD_Mesh_Vars,       ONLY:NormVec,SurfElem,nBCSides,BC,nBCs,Face_xGP
USE MOD_AnalyzeEquation_Vars,ONLY:isWall
USE MOD_TimeDisc_Vars,   ONLY:t
USE MOD_EOS_Vars,        ONLY: cp,Pr
USE MOD_Viscosity
USE MOD_Output_Vars,     ONLY: ProjectName
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nWallDoFs,iWallDoF
INTEGER                        :: SideID,iBC
REAL,ALLOCATABLE               :: WallFluxes(:,:),AllWallFluxes(:,:)
INTEGER                        :: i,j
REAL                           :: UPrim(PP_nVarPrim),gradU(PP_nVarLifting,3),gradV(3,3),gradT(3),tau(3,3)
REAL                           :: divV,muS,lambda
INTEGER                        :: OutputFileID
CHARACTER(LEN=64)              :: OutputFileName
CHARACTER(LEN=64)              :: FormatStr
INTEGER                        :: nAllWallDoFs
#if USE_MPI
INTEGER,ALLOCATABLE            :: nNodeWallDoFs(:),nOffSetWallDoFs(:)
#endif
!==================================================================================================================================
nWallDoFs = 0
DO SideID=1,nBCSides
  iBC = BC(SideID)
  IF (isWall(iBC)) THEN
    nWallDoFs = nWallDoFs + (1 + PP_N) * (1 + PP_NZ)
  END IF
END DO
! x, y, z, S, nx, ny, nz, p, q, tx, ty, tz
ALLOCATE(WallFluxes(12,nWallDoFs))
iWallDoF = 0
DO SideID=1,nBCSides
  iBC = BC(SideID)
  IF (.NOT.isWall(iBC)) CYCLE
  DO j=0,PP_NZ; DO i=0,PP_N
    iWallDoF = iWallDoF + 1
    UPrim = UPrim_master(:,i,j,SideID)
    gradU(:,1) = gradUx_master(:,i,j,SideID)
    gradU(:,2) = gradUy_master(:,i,j,SideID)
    gradU(:,3) = gradUz_master(:,i,j,SideID)
    gradV = gradU(LIFT_VELV,:)
    gradT = gradU(LIFT_TEMP,:)
    divV = gradV(1,1) + gradV(2,2) + gradV(3,3)
    muS = VISCOSITY_PRIM(UPrim)
    lambda = THERMAL_CONDUCTIVITY_H(muS)
    tau = gradV + TRANSPOSE(gradV)
    tau(1,1) = tau(1,1) - (2.0/3.0) * divV
    tau(2,2) = tau(2,2) - (2.0/3.0) * divV
    tau(3,3) = tau(3,3) - (2.0/3.0) * divV
    WallFluxes(1:3,   iWallDoF) = Face_xGP(:,i,j,0,SideID)
    WallFluxes(4,     iWallDoF) = SurfElem(i,j,0,SideID)
    WallFluxes(5:7,   iWallDoF) = NormVec(:,i,j,0,SideID)
    WallFluxes(8,     iWallDoF) = UPrim(PRES)
    WallFluxes(9,     iWallDoF) = -lambda * DOT_PRODUCT(NormVec(:,i,j,0,SideID), gradT)
    WallFluxes(10:12, iWallDoF) = -muS * MATMUL(NormVec(:,i,j,0,SideID), tau)
  END DO; END DO
END DO
#if USE_MPI
IF (MPIRoot) THEN
  ALLOCATE(nNodeWallDoFs(nProcessors))
  ALLOCATE(nOffsetWallDoFs(nProcessors))
END IF
CALL MPI_GATHER(nWallDoFs,     1, MPI_INTEGER, &
                nNodeWallDoFs, 1, MPI_INTEGER, 0, MPI_COMM_FLEXI, iError)
IF (MPIRoot) THEN
  nAllWallDoFs = SUM(nNodeWallDoFs)
  ALLOCATE(AllWallFluxes(12,nAllWallDoFs))
  nOffSetWallDoFs(1) = 0
  DO i=2,nProcessors
    nOffSetWallDoFs(i) = nOffSetWallDoFs(i-1) + 12 * nNodeWallDoFs(i-1)
  END DO
  ! array is of size (n, 12)
  nNodeWallDoFs = 12 * nNodeWallDoFs
END IF
CALL MPI_GATHERV(WallFluxes,    12 * nWallDoFs,                 MPI_DOUBLE_PRECISION, &
                 AllWallFluxes, nNodeWallDoFs, nOffSetWallDoFs, MPI_DOUBLE_PRECISION, &
                 0, MPI_COMM_FLEXI, iError)
#else
nAllWallDoFs = nWallDoFs
ALLOCATE(AllWallFluxes(12,nAllWallDoFs))
AllWallFluxes = WallFluxes
#endif

IF (MPIRoot) THEN
  OutputFileName = TRIM(TIMESTAMP(TRIM(ProjectName)//'_WallFluxes',t))//'.csv'
  OPEN(UNIT=OutputFileID, FILE=OutputFileName, STATUS='UNKNOWN', ACTION='WRITE')
  WRITE(OutputFileID, '(A)') "x, y, z, S, nx, ny, nz, p, q, tx, ty, tz"
  DO i=1,nAllWallDoFs
    WRITE(OutputFileID, '(E23.14E5, 11(",",1X,E23.14E5))') AllWallFluxes(:,i)
  END DO
  CLOSE(OutputFileID)
END IF

SDEALLOCATE(WallFluxes)
SDEALLOCATE(AllWallFluxes)
#if USE_MPI
SDEALLOCATE(nNodeWallDoFs)
SDEALLOCATE(nOffsetWallDoFs)
#endif
END SUBROUTINE CalcWallFluxes

END MODULE MOD_CalcBodyForces
