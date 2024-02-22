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

!===================================================================================================================================
!> Subroutines used for calculation of dynamic Smagorinksy SGS model.
!> See Germano, Massimo, et al. "A dynamic subgrid‐scale eddy viscosity model." Physics of Fluids A: Fluid Dynamics 3.7 (1991):
!> 1760-1765 for details of model.
!===================================================================================================================================
MODULE MOD_DynSmagorinsky
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE InitDynSmagorinsky
   MODULE PROCEDURE InitDynSmagorinsky
END INTERFACE

INTERFACE DynSmagorinsky
   MODULE PROCEDURE DynSmagorinsky_Point
   MODULE PROCEDURE DynSmagorinsky_Volume
END INTERFACE

INTERFACE FinalizeDynSmagorinsky
   MODULE PROCEDURE FinalizeDynSmagorinsky
END INTERFACE

PUBLIC::InitDynSmagorinsky, DynSmagorinsky_Volume, FinalizeDynSmagorinsky
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Read in user-defined parameters and initialize arrays.
!> We define what directions should be used to average and filter.
!===================================================================================================================================
SUBROUTINE InitDynSmagorinsky()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars
USE MOD_StringTools        ,ONLY: STRICMP
USE MOD_ReadInTools        ,ONLY: GETREAL,GETINT,GETSTR,GETREALARRAY,GETLOGICAL,CountOption
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,Elem_xGP,offsetElem
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType,wGP,NodeTypeCL
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: changeBasis3D
USE MOD_Interpolation_Vars ,ONLY: wGP
USE MOD_Mesh_Vars          ,ONLY: sJ,nElems
USE MOD_Mesh_Vars          ,ONLY: dXCL_N
USE MOD_Mesh_Vars          ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Testcase_Vars      ,ONLY: testcase
USE MOD_HDF5_input         ,ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataSize,ReadAttribute
USE MOD_2D                 ,ONLY: ExpandArrayTo3D
USE MOD_Mesh_Vars          ,ONLY: MeshFile
USE MOD_EOS_Vars           ,ONLY: mu0
USE MOD_IO_HDF5
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL :: isn
INTEGER :: i,j,k,iElem
INTEGER :: N_testFilter
REAL    :: WallDistTrsh, distNorm
REAL    :: vec(3)
REAL    :: WallDist(3,nElems)
REAL    :: Vdm_CLN_N(0:PP_N,0:PP_N)
REAL    :: dX_N(3,0:PP_N,0:PP_N,0:PP_N)
REAL    :: DeltaS_m(3,nElems)
LOGICAL :: doAvgDir(3,nElems)
CHARACTER(LEN=255) :: WallDistFile

LOGICAL             :: file_exists
INTEGER             :: HSize_proc(4)
CHARACTER(LEN=255)  :: FileName
REAL                :: CellVol
REAL,ALLOCATABLE    :: yWall_local(:,:,:,:)   
!===================================================================================================================================
IF(((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone)).OR.DynSmagorinskyInitIsDone)THEN
  CALL CollectiveStop(__STAMP__,&
    "InitDynSmagorinsky not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT Dynamic Smagorinsky...'

#if FV_ENABLED
CALL CollectiveStop(__STAMP__,"The Dynamic Smagorinsky model is not tested for FV yet!.")
#endif

! allocate memory
ALLOCATE(yWall(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
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

! Allocate necessary arrays
ALLOCATE(damp(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(IntElem(0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(doFilterDir(3,nElems))
ALLOCATE(averageType(nElems))

! Build filter matrix for test filter (modal cut-off filter)
N_testFilter = GETINT('N_testFilter')
IF (N_testFilter .LT. 1) N_testFilter = PP_N-2
IF (N_testFilter .LT. 1) CALL CollectiveStop(__STAMP__,&
    "Please use at least N=3 for the dynamic Smagorinsky model.")

ALLOCATE(FilterMat_testFilter(0:PP_N,0:PP_N))
FilterMat_testFilter(:,:) = 0.
DO i=0,N_testFilter
  FilterMat_testFilter(i,i) = 1.
END DO
SWRITE(UNIT_StdOut,'(A)',ADVANCE='NO')'TEST FILTER, FILTER DIAGONAL: '
DO i=0,PP_N
  SWRITE(UNIT_StdOut,'(F7.3)')FilterMat_testFilter(i,i)
END DO
FilterMat_testFilter=MATMUL(MATMUL(Vdm_Leg,FilterMat_testFilter),sVdm_Leg)

! Do filter in all directions
DO iElem=1,nElems
  doFilterDir(:,iElem) = .TRUE.
  doAvgDir(   :,iElem) = .TRUE.
  averageType(iElem)   = 1
END DO

! Build integration weights for each of possibly 7 average types
! sort doAvgDir to averageType for select cases
CALL GetVandermonde(PP_N, NodeTypeCL, PP_N, NodeType, Vdm_CLN_N, modal=.FALSE.)
averageType = 0
DO iElem=1,nElems
IF (doAvgDir(1,iElem).AND.doAvgDir(2,iElem).AND.doAvgDir(3,iElem)) THEN ! Volume
    IntElem(:,:,:,iElem) = 1./sJ(:,:,:,iElem,0)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)*wGP(j)*wGP(k)
    END DO; END DO; END DO
    averageType(iElem) = 1
  ELSEIF ((.NOT.doAvgDir(1,iElem)).AND.doAvgDir(2,iElem).AND.doAvgDir(3,iElem)) THEN ! I-Plane
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(Metrics_ftilde(:,i,j,k,iElem,0))*wGP(j)*wGP(k)
    END DO; END DO; END DO
    averageType(iElem) = 2
  ELSEIF (doAvgDir(1,iElem).AND.(.NOT.doAvgDir(2,iElem)).AND.doAvgDir(3,iElem)) THEN ! J-Plane
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(Metrics_gtilde(:,i,j,k,iElem,0))*wGP(i)*wGP(k)
    END DO; END DO; END DO
    averageType(iElem) = 3
  ELSEIF (doAvgDir(1,iElem).AND.doAvgDir(2,iElem).AND.(.NOT.doAvgDir(3,iElem))) THEN ! K-Plane
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(Metrics_htilde(:,i,j,k,iElem,0))*wGP(i)*wGP(j)
    END DO; END DO; END DO
    averageType(iElem) = 4
  ELSEIF (doAvgDir(1,iElem).AND.(.NOT.doAvgDir(2,iElem)).AND.(.NOT.doAvgDir(3,iElem))) THEN ! I-Line
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(1,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(dX_N(:,i,j,k))*wGP(i)
    END DO; END DO; END DO
    averageType(iElem) = 5
  ELSEIF ((.NOT.doAvgDir(1,iElem)).AND.doAvgDir(2,iElem).AND.(.NOT.doAvgDir(3,iElem))) THEN ! J-Line
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(2,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(dX_N(:,i,j,k))*wGP(j)
    END DO; END DO; END DO
    averageType(iElem) = 6
  ELSEIF ((.NOT.doAvgDir(1,iElem)).AND.(.NOT.doAvgDir(2,iElem)).AND.doAvgDir(3,iElem)) THEN ! K-Line
    CALL ChangeBasis3D(3,PP_N,PP_N,Vdm_CLN_N,dXCL_N(3,1:3,:,:,:,iElem),dX_N(:,:,:,:))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      IntElem(i,j,k,iElem) = NORM2(dX_N(:,i,j,k))*wGP(k)
    END DO; END DO; END DO
    averageType(iElem) = 7
  ELSE
    CALL CollectiveStop(__STAMP__,&
    "average type not defined, some error in wall dist file")
  END IF
END DO

DynSmagorinskyInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Dynamic Smagorinsky DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDynSmagorinsky

!===================================================================================================================================
!> Compute Dynamic Smagorinsky Eddy-Visosity
!===================================================================================================================================
PPURE SUBROUTINE DynSmagorinsky_Point(gradUx,gradUy,gradUz,UPrim,damp,y,delta,hmax,muSGS,fd)
! MODULES
USE MOD_Equation_Vars,  ONLY: Cmu, epsOMG, epsTKE
USE MOD_EddyVisc_Vars,  ONLY: CDES
USE MOD_VISCOSITY 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim)   ,INTENT(IN)  :: UPrim   !> pointwise UPrim
REAL                          ,INTENT(IN)  :: damp    !> constant factor (Cdes * Delta)**2
REAL                          ,INTENT(IN)  :: y       !> wall distance
REAL                          ,INTENT(IN)  :: delta   !> V^(1/3)
REAL                          ,INTENT(IN)  :: hmax    !> local grid size
REAL,DIMENSION(2)             ,INTENT(OUT) :: muSGS   !> pointwise eddyviscosity
REAL                          ,INTENT(OUT) :: fd      !> fd

!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: Cdyn, Clim, nu, Clim_epsilon, Clim_Lk
REAL                                    :: g, omega

REAL                                    :: rho, sRho
REAL                                    :: magS, mutLim, muOrig
REAL                                    :: lLES, lRANS, lDDES, rd
REAL                                    :: nuEff, dUdU
REAL                                    :: muS
REAL,PARAMETER                          :: kappa=0.41

REAL,PARAMETER                          :: C_des0 = 0.12, Clim_alpha = 25., Clim_beta = 0.05
!===================================================================================================================================
#if PP_dim==3
magS = SQRT(2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2) &
               + (gradUy(LIFT_VEL1)    + gradUx(LIFT_VEL2))**2 &
               + (gradUz(LIFT_VEL1)    + gradUx(LIFT_VEL3))**2 &
               + (gradUz(LIFT_VEL2)    + gradUy(LIFT_VEL3))**2 )
#else
magS = SQRT(2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2) &
               + (gradUy(LIFT_VEL1)    + gradUx(LIFT_VEL2))**2)
#endif

omega = 1. / MAX(Cmu * (UPrim(OMG))**2, epsOMG)

Cdyn = SQRT(MAX(damp, 0.)) / delta

nu = VISCOSITY_PRIM(UPrim) / UPrim(DENS)
Clim_epsilon = 2. * (C_des0 * hmax)**2 * omega * magS**2 + Cmu * UPrim(TKE) * omega
Clim_Lk = (nu**3 / Clim_epsilon)**(0.25)
Clim = C_des0 * (1. - TANH(Clim_alpha * EXP(-Clim_beta * hmax / Clim_Lk)))

sRho  = 1.0 / UPrim(DENS)

mutLim = UPrim(DENS) * ABS(UPrim(TKE)) / MAX( SQRT(6.0) * magS, 1.0e-16)

muOrig = Cmu * UPrim(DENS) * MAX(UPrim(TKE), epsTKE) * UPrim(OMG)**2

muSGS(2) = muOrig ! muSGS(2) = rho Cmu k g^2 without limiter

lRANS = SQRT(MIN(mutLim, muOrig) * sRho * Cmu * UPrim(OMG)**2)

muS   = VISCOSITY_PRIM(UPrim)
nuEff = MAX( 0., Cmu * UPrim(TKE) * UPrim(OMG)**2 ) + muS / UPrim(DENS)

dUdU  = gradUx(LIFT_VEL1) * gradUx(LIFT_VEL1)&
      + gradUx(LIFT_VEL2) * gradUx(LIFT_VEL2)&
      + gradUx(LIFT_VEL3) * gradUx(LIFT_VEL3)&
      + gradUy(LIFT_VEL1) * gradUy(LIFT_VEL1)&
      + gradUy(LIFT_VEL2) * gradUy(LIFT_VEL2)&
      + gradUy(LIFT_VEL3) * gradUy(LIFT_VEL3)&
      + gradUz(LIFT_VEL1) * gradUz(LIFT_VEL1)&
      + gradUz(LIFT_VEL2) * gradUz(LIFT_VEL2)&
      + gradUz(LIFT_VEL3) * gradUz(LIFT_VEL3)

rd    = nuEff / ((kappa * y)**2 * SQRT(MAX(1.e-16, dUdU)))
fd    = 1.0 - TANH((8.0 * rd)**3.0)

lLES = MAX(Clim, Cdyn) * (fd * Delta + (1. - fd) * hmax)

lDDES = lRANS - fd * MAX(0., lRANS-lLES)

muSGS(1) = UPrim(DENS) * lDDES**2 / MAX( Cmu * (UPrim(OMG))**2, epsOMG ) ! muSGS(1) = DDES muSGS
END SUBROUTINE DynSmagorinsky_Point

!===================================================================================================================================
!> Compute Dynamic Smagorinsky Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE DynSmagorinsky_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems, Elem_hmx
USE MOD_EddyVisc_Vars,     ONLY: damp, muSGS, muSGS_limits, yWall, DeltaS, fd
USE MOD_Lifting_Vars,      ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,           ONLY: U, UPrim
USE MOD_EOS_Vars,          ONLY: mu0
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!===================================================================================================================================
! Compute Dynamic Smagorinsky Coefficients
CALL Compute_Cd(U)

DO iElem = 1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CALL DynSmagorinsky_Point(gradUx(:   ,i,j,k,iElem), gradUy(:,i,j,k,iElem), gradUz(:,i,j,k,iElem), &
                              UPrim(:,i,j,k,iElem),     damp(1,i,j,k,iElem),   yWall(i,j,k,0,iElem),  &
                              DeltaS(iElem),            Elem_hmx(iElem),       muSGS(:,i,j,k,iElem),  &
                              fd(i,j,k,iElem))
  END DO; END DO; END DO ! i,j,k
END DO
!WRITE(*,*) MAXVAL(damp)
END SUBROUTINE DynSmagorinsky_Volume

!===============================================================================================================================
!> Compute TKE indicator needed for the modification of Smagorinskys viscosity, storing (C_des * delta)^2 in damp
!===============================================================================================================================
SUBROUTINE Compute_Cd(U_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EddyVisc_Vars ,ONLY: damp,DeltaS,FilterMat_testFilter
USE MOD_EddyVisc_Vars ,ONLY: doFilterDir,averageType,IntElem
USE MOD_Lifting_Vars  ,ONLY: gradUx,gradUy,gradUz
USE MOD_Mesh_Vars     ,ONLY: nElems
USE MOD_Equation_Vars ,ONLY: Cmu
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< Conservative volume solution
!-------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                      :: i,j,k,l,m,iElem
REAL,DIMENSION(3,3,0:PP_N,0:PP_N,0:PP_N)     :: S_lm, M_lm, L_lm, gradV
REAL,DIMENSION(  3,0:PP_N,0:PP_N,0:PP_N)     :: V,V_filtered
REAL,DIMENSION(  1,0:PP_N,0:PP_N,0:PP_N)     :: g,g_filtered
REAL,DIMENSION(  1,0:PP_N,0:PP_N,0:PP_N)     :: r,r_filtered ! rho
REAL,DIMENSION(  1,0:PP_N,0:PP_N,0:PP_N)     :: traceL
REAL,DIMENSION(    0:PP_N,0:PP_N,0:PP_N)     :: S_eN,S_eN_filtered
REAL,DIMENSION(    0:PP_N,0:PP_N,0:PP_N)     :: MM ,ML ,divV
REAL                                         :: MMa,MLa
!===============================================================================================================================
DO iElem=1,nElems
  ! TODO: Use UPrim here!!
  V(1,:,:,:) = U_in(MOM1,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  V(2,:,:,:) = U_in(MOM2,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  V(3,:,:,:) = U_in(MOM3,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  g(1,:,:,:) = U_in(RHOG,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  r(1,:,:,:) = U_in(DENS,:,:,:,iElem)

  ! Store gradients in matrix for readability and filtering
  gradV(1:3,1,:,:,:) = gradUx(VELV,:,:,:,iElem)
  gradV(1:3,2,:,:,:) = gradUy(VELV,:,:,:,iElem)
  gradV(1:3,3,:,:,:) = gradUz(VELV,:,:,:,iElem)
  divV(:,:,:) = 1./3.*(gradV(1,1,:,:,:)+gradV(2,2,:,:,:)+gradV(3,3,:,:,:))

  !                             _
  ! Filter velocities to obtain u
  V_Filtered = V
  CALL Filter_Selective(3,FilterMat_testFilter,V_filtered,doFilterDir(:,iElem))

  g_filtered = g
  CALL Filter_Selective(1,FilterMat_testFilter,g_filtered,doFilterDir(:,iElem))

  r_filtered = r
  CALL Filter_Selective(1,FilterMat_testFilter,r_filtered,doFilterDir(:,iElem))

  !              _ __ __ _2  ___________
  ! Compute L = (r ui uj g - r ui uj g^2)
  DO l=1,3
    DO m=1,3
      L_lm(l,m,:,:,:) = -r(1,:,:,:) * V(l,:,:,:) * V(m,:,:,:) * (g(1,:,:,:))**2
    END DO                                                                               !  ___________
    CALL Filter_Selective(3,FilterMat_testFilter,L_lm(l,1:3,:,:,:),doFilterDir(:,iElem)) ! -r ui uj g^2
  END DO
  DO l=1,3
    DO m=1,3
      L_lm(l,m,:,:,:) = L_lm(l,m,:,:,:) + r_filtered(1,:,:,:) * V_filtered(l,:,:,:) * V_filtered(m,:,:,:) * (g_filtered(1,:,:,:))**2
    END DO
  END DO

  ! compute dev(L) = L - 1/3 tr(L) I
  traceL = 0.
  DO l=1,3
    traceL(1,:,:,:) = traceL(1,:,:,:) + L_lm(l,l,:,:,:)
  END DO
  traceL = traceL / 3.0
  DO l=1,3
    L_lm(l,l,:,:,:) = L_lm(l,l,:,:,:) - traceL(1,:,:,:)
  END DO

  !                _____          _ ___   _____
  ! Compute M = { (Delta/Delta)^2 r Sij - r Sij } / Cmu
  !                               _ ___   _____
  !           = { 2^2             r Sij - r Sij } / Cmu
  !                 _ ___   _____
  !           = { 4 r Sij - r Sij } / Cmu

  ! Compute S
  DO m=1,3; DO l=1,3
    S_lm(l,m,:,:,:) = 0.5*(gradV(l,m,:,:,:)+gradV(m,l,:,:,:))
  END DO; END DO

  ! TODO: Do we need to do this?
  ! Correct for compressibility
  S_lm(1,1,:,:,:) = S_lm(1,1,:,:,:) - divV(:,:,:)
  S_lm(2,2,:,:,:) = S_lm(2,2,:,:,:) - divV(:,:,:)
  S_lm(3,3,:,:,:) = S_lm(3,3,:,:,:) - divV(:,:,:)

  !                 _____
  ! Save first term r Sij
  DO m=1,3
    DO l=1,3
      M_lm(l,m,:,:,:) = -r(1,:,:,:) * S_lm(l,m,:,:,:)
    END DO ! l
    CALL Filter_Selective(3,FilterMat_testFilter,M_lm(1:3,m,:,:,:),doFilterDir(:,iElem))
  END DO ! m

  ! Filter gradients
  ! ATTENTION: Overwrite gradients with filtered version
  CALL Filter_Selective(3,FilterMat_testFilter,gradV(:,1,:,:,:),doFilterDir(:,iElem))
  CALL Filter_Selective(3,FilterMat_testFilter,gradV(:,2,:,:,:),doFilterDir(:,iElem))
  CALL Filter_Selective(3,FilterMat_testFilter,gradV(:,3,:,:,:),doFilterDir(:,iElem))

  !         _
  ! Compute S
  DO m=1,3; DO l=1,3
    S_lm(l,m,:,:,:) = 0.5*(gradV(l,m,:,:,:)+gradV(m,l,:,:,:))
  END DO; END DO ! l,m
  divV(:,:,:) = 1./3.*(gradV(1,1,:,:,:)+gradV(2,2,:,:,:)+gradV(3,3,:,:,:))

  ! Correct for compressibility
  S_lm(1,1,:,:,:) = S_lm(1,1,:,:,:) - divV(:,:,:)
  S_lm(2,2,:,:,:) = S_lm(2,2,:,:,:) - divV(:,:,:)
  S_lm(3,3,:,:,:) = S_lm(3,3,:,:,:) - divV(:,:,:)

  !                _____          _ ___   _____
  ! Compute M = { (Delta/Delta)^2 r Sij - r Sij } / Cmu
  !                               _ ___   _____
  !           = { 2^2             r Sij - r Sij } / Cmu
  !                 _ ___   _____
  !           = { 4 r Sij - r Sij } / Cmu
  DO m=1,3; DO l=1,3
    M_lm(l,m,:,:,:) = M_lm(l,m,:,:,:) + 4. * r_filtered(1,:,:,:) * S_lm(l,m,:,:,:)
  END DO; END DO ! l,m
  M_lm = M_lm / Cmu

  ! Compute C_d**2 = 0.5 * (M_lm * L_lm) / (M_lm * M_lm)
  ! ATTENTION: The resulting coeffcient corresponds to C_d**2 = (C_s * Delta)^2, i.e. does already incorporate the filter width!

  ! contract with M_lm according to least square approach of Lilly
  MM=0.
  ML=0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO m=1,3; DO l=1,3
        ML(i,j,k) = ML(i,j,k) + (M_lm(l,m,i,j,k)*L_lm(l,m,i,j,k))
        MM(i,j,k) = MM(i,j,k) + (M_lm(l,m,i,j,k)*M_lm(l,m,i,j,k))
    END DO; END DO ! l,m
  END DO; END DO; END DO ! i,j,k

  ! Selective cell average (on selected directions)
  SELECT CASE (averageType(iElem))
  CASE (1) ! Volume
    MMa = 0.; MLa = 0.
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
      MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
    END DO; END DO; END DO ! i,j,k
    damp(1,:,:,:,iElem) = 0.5*MLa/MMa
  CASE (2) ! I-Plane
    DO i=0,PP_N
      MMa = 0.; MLa = 0.
      DO k=0,PP_N; DO j=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO ! j,k
      damp(1,i,:,:,iElem) = 0.5*MLa/MMa
    END DO !i
  CASE (3) ! J-Plane
    DO j=0,PP_N
      MMa = 0.; MLa = 0.
      DO k=0,PP_N; DO i=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO ! i,k
      damp(1,:,j,:,iElem) = 0.5*MLa/MMa
    END DO !j
  CASE (4) ! K-Plane
    DO k=0,PP_N
      MMa = 0.; MLa = 0.
      DO j=0,PP_N; DO i=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO; END DO !i,j
      damp(1,:,:,k,iElem) = 0.5*MLa/MMa
    END DO !k
  CASE (5) ! I-Line
    DO k=0,PP_N; DO j=0,PP_N
      MMa = 0.; MLa = 0.
      DO i=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! i
      damp(1,:,j,k,iElem) = 0.5*MLa/MMa
    END DO; END DO ! j,k
  CASE (6) ! J-Line
    DO k=0,PP_N; DO i=0,PP_N
      MMa = 0.; MLa = 0.
      DO j=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! j
      damp(1,i,:,k,iElem) = 0.5*MLa/MMa
    END DO; END DO ! i,k
  CASE (7) ! K-Line
    DO j=0,PP_N; DO i=0,PP_N
      MMa = 0.; MLa = 0.
      DO k=0,PP_N
        MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
        MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
      END DO ! k
      damp(1,i,j,:,iElem) = 0.5*MLa/MMa
    END DO; END DO! i,j
  END SELECT
END DO
END SUBROUTINE Compute_Cd

!===============================================================================================================================
!> Is vector vec1 "more normal" to vec2 "than parallel"?
!===============================================================================================================================
FUNCTION ISNORMAL(vec1,vec2) RESULT(norm)
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL, INTENT(IN) :: vec1(:), vec2(:) !< Two vectors that should be compared
LOGICAL          :: norm             !< TRUE meaning that vec1 is normal to vec2
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: vec2Norm(size(vec2,1)), parcomp, normvec(size(vec1,1)), normcomp
!==================================================================================================================================
! norm of vec2
vec2Norm = vec2/NORM2(vec2)
! component of vec1 parallel to vec2
parcomp = DOT_PRODUCT(vec1,vec2Norm)
! part of vec1 normal to vec2
normvec = vec1 - parcomp*vec2Norm
! component of vec1 normal to vec2
normcomp = NORM2(normvec)

! Check if normal component larger than parallel component
norm = (normcomp .GE. ABS(parcomp))

END FUNCTION ISNORMAL

!===============================================================================================================================
!> Filters a volume-sized array with a given filter matrix selectively in xi, eta and zeta direction depending on doFilterDir.
!===============================================================================================================================
SUBROUTINE Filter_Selective(NVar,FilterMat,U_in,doFilterDir)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: NVar                            !< Number of variables in the first dimension
REAL,INTENT(IN)     :: FilterMat(0:PP_N,0:PP_N)        !< filter matrix to be used
REAL,INTENT(INOUT)  :: U_in(NVar,0:PP_N,0:PP_N,0:PP_N) !< solution vector to be filtered
LOGICAL,INTENT(IN)  :: doFilterDir(:)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: i,j,k,l
REAL         :: U_Xi( NVar,0:PP_N,0:PP_N,0:PP_N)
REAL         :: U_Eta(NVar,0:PP_N,0:PP_N,0:PP_N)
!==================================================================================================================================
! Perform filtering

! Xi
IF(doFilterDir(1)) THEN
  U_Xi = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      U_Xi(:,i,j,k) = U_Xi(:,i,j,k) + FilterMat(i,l)*U_in(:,l,j,k)
    END DO !l
  END DO; END DO; END DO ! i,j,k
ELSE
  U_Xi = U_in
END IF

! Eta
IF(doFilterDir(2)) THEN
  U_Eta = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      U_Eta(:,i,j,k) = U_Eta(:,i,j,k) + FilterMat(j,l)*U_Xi(:,i,l,k)
    END DO !l
  END DO; END DO; END DO ! i,j,k
ELSE
  U_Eta = U_Xi
END IF

! Zeta
IF(doFilterDir(3)) THEN
  U_in(:,:,:,:)=0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      U_in(:,i,j,k) = U_in(:,i,j,k) + FilterMat(k,l)*U_Eta(:,i,j,l)
    END DO !l
  END DO; END DO; END DO ! i,j,k
ELSE
  U_in = U_Eta
END IF
END SUBROUTINE Filter_Selective


!===============================================================================================================================
!> Deallocate arrays and finalize variables used by dynamic Smagorinsky SGS model
!===============================================================================================================================
SUBROUTINE FinalizeDynSmagorinsky()
! MODULES
USE MOD_EddyVisc_Vars
IMPLICIT NONE
!===============================================================================================================================
SDEALLOCATE(fd)
SDEALLOCATE(ywall)
SDEALLOCATE(damp)
SDEALLOCATE(IntElem)
SDEALLOCATE(doFilterDir)
SDEALLOCATE(averageType)
DynSmagorinskyInitIsDone = .FALSE.
END SUBROUTINE FinalizeDynSmagorinsky

END MODULE MOD_DynSmagorinsky
