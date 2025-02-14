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
USE MOD_Mesh_Vars          ,ONLY: MeshInitIsDone,nElems,sJ,MeshFile,offsetElem
USE MOD_Interpolation_Vars ,ONLY: InterpolationInitIsDone,Vdm_Leg,sVdm_Leg,NodeType,wGP,NodeTypeCL
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_ChangeBasis        ,ONLY: changeBasis3D
USE MOD_Interpolation_Vars ,ONLY: wGP
USE MOD_Testcase_Vars      ,ONLY: testcase
USE MOD_HDF5_Input        ,ONLY: ReadArray,OpenDataFile,CloseDataFile,GetDataSize,ReadAttribute
USE MOD_2D                ,ONLY: ExpandArrayTo3D
USE MOD_IO_HDF5
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: file_exists
LOGICAL             :: doAvgDir(3,nElems)
INTEGER             :: i,j,k,iElem
INTEGER             :: N_testFilter
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

!#if FV_ENABLED
!CALL CollectiveStop(__STAMP__,"The Dynamic Smagorinsky model is not tested for FV yet!.")
!#endif

! Allocate necessary arrays
ALLOCATE(yWall(0:PP_N,0:PP_N,0:PP_NZ,0:FV_SIZE,nElems))
ALLOCATE(fd   (0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(DeltaRatio(nElems))

ALLOCATE(Cdes2(0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(IntElem(0:PP_N,0:PP_N,0:PP_N,nElems))
ALLOCATE(doFilterDir(3,nElems))

! Readin custom limits of eddy viscosity
muSGS_limits(1) = 0.

! Build filter matrix for test filter (modal cut-off filter)
N_testFilter = GETINT('N_testFilter')
IF ((N_testFilter .LT. 1) .AND. (PP_N .ge. 3) ) THEN
  N_testFilter = PP_N-2
ELSEIF ((N_testFilter .LT. 1) .AND. (PP_N .eq. 2) ) THEN
  N_testFilter = 1
ENDIF
IF (N_testFilter .LT. 1) CALL CollectiveStop(__STAMP__,&
    "Please use at least N=2 for the dynamic Smagorinsky model.")

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

! get wall distance for zonal averaging approach if needed
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
ELSE 
  SWRITE(UNIT_stdOut, *) "WARNING: No walldistance file found! DDES shielding not working!"
  CALL CollectiveStop(__STAMP__,'Please use POSTI to compute wall distance first!')
ENDIF

! Find in which x, y, z direction are the i, j ,k index pointing, and
! then decide which index to filter
IF(testcase.EQ.'channel') THEN
  ! Channel testcase, filter in wall-parallel directions
  DO iElem=1,nElems
    doFilterDir(:,iElem) = (/.TRUE.,.FALSE.,.TRUE./)
    doAvgDir(   :,iElem) = (/.TRUE.,.FALSE.,.TRUE./)
  END DO !iElem
ELSE
  DO iElem=1,nElems
    ! Default, filter in all directions
    doFilterDir(:,iElem) = (/.TRUE.,.TRUE.,.TRUE./)
    doAvgDir(   :,iElem) = (/.TRUE.,.TRUE.,.TRUE./)
  END DO !iElem
ENDIF

! Build integration weights 
DO iElem=1,nElems
  IntElem(:,:,:,iElem) = 1./sJ(:,:,:,iElem,0)
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    IntElem(i,j,k,iElem) = IntElem(i,j,k,iElem)*wGP(i)*wGP(j)*wGP(k)
  END DO; END DO; END DO
END DO

! Calculate the filter width deltaS: deltaS=( Cell volume )^(1/3) / ( PP_N+1 )
DO iElem=1,nElems
  CellVol = 0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CellVol = CellVol + wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO
  DeltaS(iElem) = CellVol**(1./3.)  / (REAL(PP_N)+1.)
END DO

! Compute DeltaRatio=(Delta/Delta) as ratio of filter widths (alpha by Germano)
DeltaRatio(:) = (REAL(PP_N+1)/REAL(N_testFilter+1))**2

DynSmagorinskyInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT Dynamic Smagorinsky DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitDynSmagorinsky

!===================================================================================================================================
!> Compute Dynamic Smagorinsky Eddy-Visosity
!===================================================================================================================================
PPURE SUBROUTINE DynSmagorinsky_Point(gradUx,gradUy,gradUz,UPrim,Cdes2,y,Delta,hmax,fd,muSGS)
! MODULES
USE MOD_Equation_Vars,  ONLY: Cmu, sqrt6
USE MOD_EddyVisc_Vars,  ONLY: CDES0
USE MOD_VISCOSITY 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx, gradUy, gradUz   !> Gradients in x,y,z directions
REAL,DIMENSION(PP_nVarPrim)   ,INTENT(IN)  :: UPrim                    !> pointwise conserved variable
REAL                          ,INTENT(IN)  :: Cdes2                    !> constant factor (CDES*deltaS)**2
REAL                          ,INTENT(IN)  :: y                        !> pointwise wall distance
REAL                          ,INTENT(IN)  :: Delta                    !> pointwise cell spacing
REAL                          ,INTENT(IN)  :: hmax                     !> local grid size
REAL                          ,INTENT(OUT) :: fd                       !> fd
REAL                          ,INTENT(OUT) :: muSGS                    !> pointwise eddyviscosity
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

! Already take the square root of 2 into account here
#if PP_dim==2
magS = SQRT( &
    2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2) + &
    (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))**2)
dUdU  = gradUx(LIFT_VEL1)**2 + gradUx(LIFT_VEL2)**2 &
      + gradUy(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 
#else
magS = SQRT( &
    2. * (gradUx(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2) + &
    (gradUy(LIFT_VEL1) + gradUx(LIFT_VEL2))**2 + &
    (gradUz(LIFT_VEL1) + gradUx(LIFT_VEL3))**2 + &
    (gradUz(LIFT_VEL2) + gradUy(LIFT_VEL3))**2)

dUdU  = gradUx(LIFT_VEL1)**2 + gradUx(LIFT_VEL2)**2 + gradUx(LIFT_VEL3)**2 &
      + gradUy(LIFT_VEL1)**2 + gradUy(LIFT_VEL2)**2 + gradUy(LIFT_VEL3)**2 &
      + gradUz(LIFT_VEL1)**2 + gradUz(LIFT_VEL2)**2 + gradUz(LIFT_VEL3)**2
#endif

kPos    = MAX( UPrim(TKE), 1.e-16 )
gPos    = MAX( UPrim(OMG), 1.e-16 )
muTOrig = Cmu * UPrim(DENS) * kPos * gPos**2

muTLim = MIN(muTOrig,  UPrim(DENS) * kPos / MAX(sqrt6 * magS, 1.e-16))

lRANS = SQRT(mutLim * sRho * Cmu * gPos**2 )

rd    = ( Cmu * kPos * gPos**2 + muS * sRho) / ((kappa * y)**2 * SQRT(MAX(1.e-16, dUdU)))
fd    = 1.0 - TANH((8.0 * rd)**3.0)

lLES = SQRT(Cdes2) * (fd * Delta + (1. - fd) * hmax)

lDDES = lRANS - fd * MAX(0., lRANS-lLES)

muSGS = UPrim(DENS) * lDDES**2 / MAX( Cmu * gPos**2, 1.e-16 )

END SUBROUTINE DynSmagorinsky_Point

!===================================================================================================================================
!> Compute Dynamic Smagorinsky Eddy-Visosity for the volume
!===================================================================================================================================
SUBROUTINE DynSmagorinsky_Volume()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY: nElems, Elem_hmx
USE MOD_EddyVisc_Vars,     ONLY: Cdes2, muSGS, yWall, DeltaS, fd
USE MOD_Lifting_Vars,      ONLY: gradUx, gradUy, gradUz
USE MOD_DG_Vars,           ONLY: UPrim, U
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
!===================================================================================================================================
! Compute Dynamic Smagorinsky Coefficients
CALL Compute_Cd(U)
CALL Apply_Clim()

DO iElem = 1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CALL DynSmagorinsky_Point(gradUx(:,i,j,k,iElem), gradUy(:,i,j,k,iElem),       gradUz(:,i,j,k,iElem), &
                              UPrim(:,i,j,k,iElem),  Cdes2(i,j,k,iElem),          yWall(i,j,k,0,iElem),  &
                              DeltaS(iElem),         Elem_hmx(iElem),             fd(i,j,k,iElem),       &
                              muSGS(1,i,j,k,iElem)                                )
  END DO; END DO; END DO ! i,j,k
END DO
END SUBROUTINE DynSmagorinsky_Volume

!===============================================================================================================================
!> Compute the lower bound for the dynamic coefficient
!===============================================================================================================================
SUBROUTINE Apply_Clim()
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,     ONLY: Cmu
USE MOD_Mesh_Vars,         ONLY: nElems, Elem_hmx
USE MOD_EddyVisc_Vars,     ONLY: CDES0, Cdes2
USE MOD_DG_Vars,           ONLY: UPrim
USE MOD_VISCOSITY
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
REAL                :: eta, ratio, diss, Clim
REAL                :: nuS
!===================================================================================================================================

! Limit Dynamic Smagorinsky Coefficients
DO iElem = 1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    nuS     = VISCOSITY_TEMPERATURE(UPrim(TEMP,i,j,k,iElem)) / UPrim(DENS,i,j,k,iElem)
    diss    = MAX( UPrim(TKE,i,j,k,iElem) / UPrim(OMG,i,j,k,iElem)**2 , 1.e-16)
    eta     = ( nuS**3 / diss )**0.25
    ratio   = Elem_hmx(iElem) / eta
    Clim    = 0.5 * CDES0 * ( MAX( MIN( (ratio-23.0)/7.0, 1.0), 0.0) + MAX( MIN( (ratio-65.0)/25.0, 1.0), 0.0) )
    Cdes2(i,j,k,iElem) = MAX( Cdes2(i,j,k,iElem), Clim**2 )
  END DO; END DO; END DO ! i,j,k
END DO

END SUBROUTINE Apply_Clim

!===============================================================================================================================
!> Compute TKE indicator needed for the modification of Smagorinskys viscosity.
!===============================================================================================================================
SUBROUTINE Compute_Cd(U_in)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Equation_Vars,  ONLY: Cmu
USE MOD_EddyVisc_Vars ,ONLY: Cdes2,DeltaRatio,DeltaS,FilterMat_testFilter
USE MOD_EddyVisc_Vars ,ONLY: doFilterDir,IntElem
USE MOD_Lifting_Vars  ,ONLY: gradUx,gradUy,gradUz
USE MOD_Mesh_Vars     ,ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems) !< Conservative volume solution
!-------------------------------------------------------------------------------------------------------------------------------
! aOCAL VARIABLES
INTEGER                                      :: i,j,k,l,m,iElem
REAL,DIMENSION(3,3,0:PP_N,0:PP_N,0:PP_N)     :: S_lm, M_lm, L_lm, gradV
REAL,DIMENSION(  3,0:PP_N,0:PP_N,0:PP_N)     :: V,V_filtered
REAL,DIMENSION(    0:PP_N,0:PP_N,0:PP_N)     :: omega, omega_filtered
REAL,DIMENSION(    0:PP_N,0:PP_N,0:PP_N)     :: MM ,ML ,divV
REAL                                         :: MMa,MLa
!===============================================================================================================================
DO iElem=1,nElems
  ! TODO: Use UPrim here!!
  V(1,:,:,:) = U_in(MOM1,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  V(2,:,:,:) = U_in(MOM2,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)
  V(3,:,:,:) = U_in(MOM3,:,:,:,iElem)/U_in(DENS,:,:,:,iElem)

  ! Store gradients in matrix for readability and filtering
  gradV(1:3,1,:,:,:) = gradUx(VELV,:,:,:,iElem)
  gradV(1:3,2,:,:,:) = gradUy(VELV,:,:,:,iElem)
  gradV(1:3,3,:,:,:) = gradUz(VELV,:,:,:,iElem)
  divV(:,:,:) = 1./3.*(gradV(1,1,:,:,:)+gradV(2,2,:,:,:)+gradV(3,3,:,:,:))

  !                             _
  ! Filter velocities to obtain u
  V_Filtered = V
  CALL Filter_Selective(3,FilterMat_testFilter,V_filtered,doFilterDir(:,iElem))

  !              _ _   __
  ! Compute L = -u*u + uu
  DO l=1,3
    DO m=1,3
      L_lm(l,m,:,:,:) = V(l,:,:,:)*V(m,:,:,:) ! uu
    END DO                                                                               ! __
    CALL Filter_Selective(3,FilterMat_testFilter,L_lm(l,1:3,:,:,:),doFilterDir(:,iElem)) ! uu
  END DO
  DO l=1,3
    DO m=1,3                                                                      !     __   _ _
      L_lm(l,m,:,:,:) = L_lm(l,m,:,:,:) - V_filtered(l,:,:,:)*V_filtered(m,:,:,:) ! L = uu - u*u
    END DO
  END DO

  !             _____            _ _   ____
  ! Compute M=-(Delta/Delta)**2* w S +  w S

  ! Compute S
  DO m=1,3; DO l=1,3
    S_lm(l,m,:,:,:) = 0.5*(gradV(l,m,:,:,:)+gradV(m,l,:,:,:))
  END DO; END DO

  ! Compute omega
  omega(:,:,:) = 1. / (Cmu * MAX(U_in(RHOG,:,:,:,iElem)/U_in(DENS,:,:,:,iElem),1.e-16)**2 )

  ! Correct for compressibility
  S_lm(1,1,:,:,:) = S_lm(1,1,:,:,:) - divV(:,:,:)
  S_lm(2,2,:,:,:) = S_lm(2,2,:,:,:) - divV(:,:,:)
  S_lm(3,3,:,:,:) = S_lm(3,3,:,:,:) - divV(:,:,:)

  ! Save first term  |S|S
  DO m=1,3
    DO l=1,3
      M_lm(l,m,:,:,:) = DeltaS(iElem)**2 * omega(:,:,:)*S_lm(l,m,:,:,:) ! |S|S
    END DO ! l                                                                             ____
    CALL Filter_Selective(3,FilterMat_testFilter,M_lm(1:3,m,:,:,:),doFilterDir(:,iElem))  ! |S|S
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

  ! filter omega
  omega_filtered(:,:,:) = omega(:,:,:)
  CALL Filter_Selective(1,FilterMat_testFilter,omega_filtered(:,:,:),doFilterDir(:,iElem))
  !          _

  ! Correct for compressibility
  S_lm(1,1,:,:,:) = S_lm(1,1,:,:,:) - divV(:,:,:)
  S_lm(2,2,:,:,:) = S_lm(2,2,:,:,:) - divV(:,:,:)
  S_lm(3,3,:,:,:) = S_lm(3,3,:,:,:) - divV(:,:,:)

  !             ____    _____              _ _
  ! Compute M = |w|S - (Delta/Delta)**2 * |w|S
  DO m=1,3; DO l=1,3
    M_lm(l,m,:,:,:) = M_lm(l,m,:,:,:) - DeltaRatio(iElem) * DeltaS(iElem)**2 * omega_filtered(:,:,:)*S_lm(l,m,:,:,:)
  END DO; END DO ! l,m

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
  MMa = 0.; MLa = 0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    MMa = MMa + MM(i,j,k)*IntElem(i,j,k,iElem)
    MLa = MLa + ML(i,j,k)*IntElem(i,j,k,iElem)
  END DO; END DO; END DO ! i,j,k
  Cdes2(:,:,:,iElem) = 0.5*MLa/MMa

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
DynSmagorinskyInitIsDone = .FALSE.

SDEALLOCATE(yWall)
SDEALLOCATE(fd)
SDEALLOCATE(DeltaRatio)

SDEALLOCATE(Cdes2)
SDEALLOCATE(IntElem)
SDEALLOCATE(doFilterDir)

END SUBROUTINE FinalizeDynSmagorinsky

END MODULE MOD_DynSmagorinsky
