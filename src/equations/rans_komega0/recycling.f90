#include "flexi.h"
#include "eos.h"

MODULE MOD_Recycling
IMPLICIT NONE
PUBLIC
CONTAINS

SUBROUTINE DefineParametersRecycling
USE MOD_Globals
USE MOD_ReadInTools

CALL prms%SetSection("Recycling")
CALL prms%CreateLogicalOption("doRecycling", "", "F")
CALL prms%CreateIntOption( "inftyRefState",  "Refstate required for freestream value.")
CALL prms%CreateIntOption("nDofsInY", "Number of DOFs in Y direction")
CALL prms%CreateIntOption("nDofsInZ", "Number of DOFs in Z direction")
CALL prms%CreateRealOption("recycling_xrec", "")
CALL prms%CreateRealOption("recycling_momentum_thickness", "")
CALL prms%CreateRealOption("yMatchingTolerance", "tolerance for finding match Y", "1.0E-6")
CALL prms%CreateRealOption("zMatchingTolerance", "tolerance for finding match Z", "1.0E-6")

END SUBROUTINE DefineParametersRecycling

SUBROUTINE InitRecycling
USE MOD_Global
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: Elem_xGP, Face_xGP, nElems, nBCSides, BoundaryType, BC
USE MOD_DG_Vars, ONLY: U, UPrim, UPrim_master
USE MOD_Equation_Vars, ONLY: RefStatePrim
USE MOD_ReadInTools
USE MOD_Recycling_Vars
IMPLICIT NONE

INTEGER :: i,j,k,p,q,nn,iELem,iSide,iFile,iosi, nLocalMatched, localInd, nGlobalMatched
REAL :: tmp, recycMinDist, minCrdGlb

doRecycling = GETLOGICAL("doRecycling")
IF (.NOT. doRecycling) THEN
  RETURN
END IF

freeStreamRefState=GETINT("inftyRefState")
nY      = GETINT("nDofsInY")
nZ      = GETINT("nDofsInZ")
xRecyc  = GETREAL("recycling_xrec")
thetaD  = GETREAL("recycling_momentum_thickness")
yMatchingTolerance = GETREAL("yMatchingTolerance")
zMatchingTolerance = GETREAL("zMatchingTolerance")
u_inf = RefStatePrim(VEL1,freeStreamRefState)
! Twall   = GETREAL("recycling_Twall")

! initialize recycling plane
! find the nearest location
tmp = 1.0E10
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( ((Elem_xGP(1,i,j,k,iElem) - xRecyc).LT.tmp) .AND. ((Elem_xGP(1,i,j,k,iElem) - xRecyc) .GT. 0.0) ) THEN
      tmp = Elem_xGP(1,i,j,k,iElem) - xRecyc
    END IF
END DO; END DO; END DO; END DO

CALL MPI_ALLREDUCE(tmp, recycMinDist, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FLEXI, iError)
recycMinDist = recycMinDist + xRecyc

ALLOCATE( indexMap(1:3, 0:PP_N, 0:PP_N, 0:PP_NZ, 1:nElems) )
indexMap = 0

nLocalMatched = 0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( ABS(Elem_xGP(1,i,j,k,iElem) - recycMinDist) .LT. 1.e-6 ) THEN
      nLocalMatched = nLocalMatched + 1
      indexMap(1,i,j,k,iELem) = 3
    ELSE
      indexMap(1,i,j,k,iELem) = -1
    END IF
END DO; END DO; END DO; END DO

CALL MPI_ALLREDUCE(nLocalMatched, nGlobalMatched, 1, MPI_INT, MPI_SUM, MPI_COMM_FLEXI, iError)
IF (nGlobalMatched.NE.(nY*nZ)) THEN
  CALL Abort(__STAMP__,'Number of sampled DOFs in the recycling plane is not equal to nY*nZ')
END IF

SWRITE(*,*) "Number of sampled DoFs in the recycling plane = ", nGlobalMatched

ALLOCATE( y(nY) )
ALLOCATE( z(nZ) )
y = -1.0E6
z = -1.0E6

! for the recycling plane
! use Y information to get IndexY mapping without involving MPIALLGATHER, reduce indexMap(1,i,j,k,iELem) = 3 to 2
nn = 0
DO localInd=1,nY
  nn = 0
  tmp = 1.0E6
  DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( indexMap(1,i,j,k,iELem) .EQ. 3 ) THEN
      tmp = MIN(tmp, Elem_xGP(2,i,j,k,iELem))
    END IF
  END DO; END DO; END DO; END DO
  y(localInd) = tmp
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, y(localInd), 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FLEXI, iError)
  DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( indexMap(1,i,j,k,iELem) .EQ. 3 ) THEN
      IF (ABS(Elem_xGP(2,i,j,k,iELem) - y(localInd)) .LE. yMatchingTolerance) THEN
        indexMap(2,i,j,k,iELem) = localInd
        indexMap(1,i,j,k,iELem) = 2
        nn = nn + 1
      END IF
    END IF
  END DO; END DO; END DO; END DO
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, nn, 1, MPI_INT, MPI_SUM, MPI_COMM_FLEXI, iError)
  SWRITE(*,*) "No. DoFs collected totally for ", localInd, "-th Y layer is ", nn
END DO
! check if all DOFs in the recycling plane are matched in Y coordinates
nn = 0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  IF ( indexMap(1,i,j,k,iELem) .EQ. 3 ) THEN
    nn = nn + 1
  END IF
END DO; END DO; END DO; END DO
IF ( nn .GT. 0) THEN
  CALL Abort(__STAMP__,'There exists DOFs in the recycling plane that are not matched in Y coordinates')
END IF
DO localInd = 1, nY
  SWRITE(*,*) "y(i) = ", localInd, y(localInd)
  IF (y(localInd) .EQ. 1.0E6) THEN
    CALL Abort(__STAMP__,'y layer not found')
  END IF
END DO

! use Z information to get IndexZ mapping without involving MPIALLGATHER, reduce indexMap(1,iELem) = 2 to 1
DO localInd=1,nZ
  nn = 0
  tmp = 1.0E6
  DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( indexMap(1,i,j,k,iELem) .EQ. 2 ) THEN
      tmp = MIN(tmp, Elem_xGP(3,i,j,k,iELem))
    END IF
  END DO; END DO; END DO; END DO
  z(localInd) = tmp
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, z(localInd),  1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FLEXI, iError)
  DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    IF ( indexMap(1,i,j,k,iELem) .EQ. 2 ) THEN
      IF ( ABS(Elem_xGP(3,i,j,k,iELem) - z(localInd)) .LE. zMatchingTolerance ) THEN
        indexMap(3,i,j,k,iELem) = localInd
        indexMap(1,i,j,k,iELem) = 1
        nn = nn + 1
      END IF
    END IF
  END DO; END DO; END DO; END DO
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, nn, 1, MPI_INT, MPI_SUM, MPI_COMM_FLEXI, iError)
  SWRITE(*,*) "No. DoFs collected totally for ", localInd, "-th Z layer is ", nn
END DO
! check if all DOFs in the recycling plane are matched in Z coordinates
nn = 0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  IF ( indexMap(1,i,j,k,iELem) .EQ. 2 ) THEN
    nn = nn + 1
  END IF
END DO; END DO; END DO; END DO
IF ( nn .GT. 0) THEN
  CALL Abort(__STAMP__,'There exists DOFs in the recycling plane that are not matched in Z coordinates')
END IF
DO localInd = 1, nZ
  SWRITE(*,*) "z(i) = ", localInd, z(localInd)
  IF (z(localInd) .EQ. 1.0E6) THEN
    CALL Abort(__STAMP__,'z layer not found')
  END IF
END DO

! initialize inflow plane
nInflow = 0
DO iSide=1,nBCSides; DO q=0,PP_NZ; DO p=0,PP_N
  IF (BoundaryType(BC(iSide),BC_TYPE).EQ.32) THEN
    nInflow = nInflow + 1
  END IF
END DO; END DO; END DO
ALLOCATE(im(3,nInflow))
nInflow = 0
DO iSide=1,nBCSides; DO q=0,PP_NZ; DO p=0,PP_N
  IF (BoundaryType(BC(iSide),BC_TYPE).EQ.32) THEN
    nInflow = nInflow + 1
    im(:,nInflow) = (/p,q,iSide/)
  END IF
END DO; END DO; END DO

ALLOCATE(rim(0:PP_N,0:PP_NZ,nBCSides))
rim = 0
DO nn=1,nInflow
  rim(im(1,nn),im(2,nn),im(3,nn)) = nn
END DO

ALLOCATE(img(2,nInflow))
img = -1
DO nn=1,nInflow
  DO j=1,nY
    IF (ABS(Face_xGP(2,im(1,nn),im(2,nn),0,im(3,nn)) - y(j)).LE.yMatchingTolerance) THEN
      img(1,nn) = j
      EXIT
    END IF
  END DO
  IF ( img(1,nn) .LT. 0) THEN
    CALL Abort(__STAMP__,'there is a inflow DOF that is not matched in Y coordinate')
  END IF
  DO k=1,nZ
    IF (ABS(Face_xGP(3,im(1,nn),im(2,nn),0,im(3,nn)) - z(k)).LE.zMatchingTolerance) THEN
      img(2,nn) = k
      EXIT
    END IF
  END DO
  IF ( img(2,nn) .LT. 0) THEN
    SWRITE(*,*) "current DOF's Z coordinate is = ", Face_xGP(3,im(1,nn),im(2,nn),0,im(3,nn))
    CALL Abort(__STAMP__,'there is a inflow DOF that is not matched in Z coordinate')
  END IF
END DO

ALLOCATE(recycl_U_global(PP_nVar,nY,nZ))
ALLOCATE(recycl_UPrim_global(PP_nVarPrim,nY,nZ))
ALLOCATE(recycl_U_mean(PP_nVar,nY))
ALLOCATE(recycl_UPrim_mean(PP_nVarPrim,nY))
ALLOCATE(recycl_UPrim_fluc(PP_nVarPrim,nY,nZ))
recycl_U_global     = 0.0
recycl_UPrim_global = 0.0
recycl_U_mean       = 0.0
recycl_UPrim_mean   = 0.0
recycl_UPrim_fluc   = 0.0

END SUBROUTINE InitRecycling

SUBROUTINE Recycling
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: Elem_xGP, Face_xGP, nElems, nBCSides
USE MOD_DG_Vars, ONLY: U, UPrim, UPrim_master
USE MOD_TimeDisc_Vars,ONLY: dt, RKb, CurrentStage
USE MOD_EOS, ONLY: ConsToPrim, PrimToCons
USE MOD_Viscosity
USE MOD_Recycling_Vars
IMPLICIT NONE

INTEGER :: i,j,k,p,q,nn,iELem,iSide,iFile,ios
REAL :: tmp,u1,u2
REAL :: tauw,rhow

! RETURN

recycl_U_global = 0.0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  IF ( indexMap(1,i,j,k,iELem) .EQ. 1 ) THEN
    recycl_U_global(:, indexMap(2,i,j,k,iELem), indexMap(3,i,j,k,iELem) ) = U(:,i,j,k,iELem)
  END IF
END DO; END DO; END DO; END DO
nn = PP_nVar * nY * nZ
CALL MPI_ALLREDUCE(MPI_IN_PLACE, recycl_U_global, nn, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FLEXI, iError)

DO k=1,nZ; DO j=1,nY
  CALL ConsToPrim(recycl_UPrim_global(:,j,k), recycl_U_global(:,j,k))
END DO; END DO

recycl_U_mean = 0.0
DO k=1,nZ; DO j=1,nY
  recycl_U_mean(:,j) = recycl_U_mean(:,j) + recycl_U_global(:,j,k) / nZ
END DO; END DO

DO j=1,nY
  CALL ConsToPrim(recycl_UPrim_mean(:,j), recycl_U_mean(:,j))
END DO

DO k=1,nZ; DO j=1,nY
  recycl_UPrim_fluc(:,j,k) = recycl_UPrim_global(:,j,k) - recycl_UPrim_mean(:,j)
END DO; END DO

tmp = 0.0
DO j=nY-5,nY
  tmp = tmp + recycl_UPrim_mean(VEL1,j)
END DO
tmp = tmp / 6.0
! tmp = u_inf

! from top to bottom, because sometimes fluctuation changes delta99 estimation if from wall
d99R = 0.0
DO k = 2, (nY - 1)
  j = nY - k
  IF ( (recycl_UPrim_mean(VEL1,j) .LE. (0.99 * tmp)) .AND. (recycl_UPrim_mean(VEL1,j+1) .GT. (0.99 * tmp))) THEN
    u1 = recycl_UPrim_mean(VEL1,j)
    u2 = recycl_UPrim_mean(VEL1,j+1)
    d99R = ((u2 - (0.99 * tmp)) * y(j) + ((0.99 * tmp) - u1) * y(j+1)) / (u2 - u1)
    EXIT
  ENDIF
END DO
! if not found or uniform initialization, just start with scaling ratio of 1
IF (d99R .LT. 1e-8) THEN
  d99R = thetaD
END IF

#if PP_NodeType != 2
! TODO(Shimushu): remove hard-coded variables
rhow = ((y(2) - 0.0) * recycl_UPrim_mean(DENS,1) + (0.0 - y(1)) * recycl_UPrim_mean(DENS,2)) / (y(2) - y(1))
Twall = ((y(2) - 0.0) * recycl_UPrim_mean(TEMP,1) + (0.0 - y(1)) * recycl_UPrim_mean(TEMP,2)) / (y(2) - y(1))
tauw = VISCOSITY_TEMPERATURE(Twall) * recycl_UPrim_mean(VEL1,1) / y(1)
#else
rhow = recycl_UPrim_mean(DENS,1)
Twall = recycl_UPrim_mean(TEMP,1)
tauw = VISCOSITY_TEMPERATURE(Twall) * recycl_UPrim_mean(VEL1,2) / y(2)
#endif
uTauR = SQRT(tauw / rhow)
! for hard start
IF (uTauR .LT. 1.0e-8) THEN
  uTauR = 1.0
ENDIF

! TODO(Shimushu): remove hard-coded variables
d99I = MIN( thetaD, d99R ) ! avoid wrong instantaneous scaling ratio for natural development
uTauI = uTauR * (d99R / d99I)**(0.1)

utau_ratio = uTauI / uTauR

! NOTE(Shimushu): assume rho_w,I = rho_w,R
! TODO(Shimushu): remove hard-coded variables
dnuR = VISCOSITY_TEMPERATURE(Twall) / (rhow * uTauR)
dnuI = VISCOSITY_TEMPERATURE(Twall) / (rhow * uTauI)

! thetaR = MIN(thetaR,thetaD)

recycling_called = .TRUE.

END SUBROUTINE Recycling

SUBROUTINE FillState(UPrim_boundary,UPrim_master,iSide)
USE MOD_Preproc
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY: Face_xGP
USE MOD_EOS_Vars, ONLY: R
USE MOD_EOS, ONLY: ConsToPrim
USE MOD_Equation_Vars  ,ONLY: RefStatePrim
USE MOD_Recycling_Vars
USE MOD_ExactFunc_Vars, ONLY: BCData, BCLength
IMPLICIT NONE

REAL,INTENT(OUT)   :: UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)    :: UPrim_master(  PP_nVarPrim,0:PP_N,0:PP_NZ)
INTEGER,INTENT(IN) :: iSide

REAL, PARAMETER :: cutEta = 1.25

INTEGER :: p,q
INTEGER :: lower,upper,i
REAL    :: UCons(PP_nVar)

REAL    :: UMeanInner(PP_nVar)
REAL    :: UFlucInner(PP_nVar)
REAL    :: UMeanOuter(PP_nVar)
REAL    :: UFlucOuter(PP_nVar)
REAL    :: UPrimMeanInner(PP_nVarPrim)
REAL    :: UPrimMeanOuter(PP_nVarPrim)
REAL    :: UPrimFlucInner(PP_nVarPrim)
REAL    :: UPrimFlucOuter(PP_nVarPrim)
REAL    :: UPrimMean(PP_nVarPrim)
REAL    :: UPrimMean2(PP_nVar)
REAL    :: UPrimFluc(PP_nVarPrim)
REAL    :: yplus,eta
REAL    :: weta
REAL    :: bdamp
REAL    :: fdamp

IF (.NOT.recycling_called) THEN
  UPrim_boundary = UPrim_master
  RETURN
END IF

DO q=0,PP_NZ; DO p=0,PP_N
  CALL BinarySearch(PP_nVarPrim, recycl_UPrim_mean(:,:),                            UPrimMeanInner, (dnuR / dnuI) * Face_xGP(2,p,q,0,iSide))
  CALL BinarySearch(PP_nVarPrim, recycl_UPrim_mean(:,:),                            UPrimMeanOuter, (d99R / d99I) * Face_xGP(2,p,q,0,iSide))
  CALL BinarySearch(PP_nVarPrim, recycl_UPrim_fluc(:,:,nZ+1-img(2,rim(p,q,iSide))), UPrimFlucInner, (dnuR / dnuI) * Face_xGP(2,p,q,0,iSide))
  CALL BinarySearch(PP_nVarPrim, recycl_UPrim_fluc(:,:,nZ+1-img(2,rim(p,q,iSide))), UPrimFlucOuter, (d99R / d99I) * Face_xGP(2,p,q,0,iSide))

  eta = Face_xGP(2,p,q,0,iSide) / d99I
  weta = 0.5 * (1.0 + TANH(4.0 * (eta - 0.2) / ((1.0 - 2.0 * 0.2) * eta + 0.2)) / TANH(4.0))
  bdamp = 0.5 * (1.0 - TANH(4.0 * (eta - 2.0)))

  IF (eta .GT. cutEta) THEN
    UPrim_boundary(:,p,q)=RefStatePrim(:,freeStreamRefState)
  ELSEIF (eta .GT. 1.0) THEN
    weta = 1.0
    UPrimFlucOuter = 0.0

    UPrimMeanInner(VEL1) = UPrimMeanInner(VEL1) * utau_ratio
    UPrimMeanInner(TKE ) = UPrimMeanInner(TKE ) * utau_ratio**2.0
    UPrimMeanInner(OMG ) = UPrimMeanInner(OMG ) * utau_ratio**2.0

    UPrimMeanOuter(VEL1) = UPrimMeanOuter(VEL1) * utau_ratio + u_inf * (1.0 - utau_ratio)

    UPrimFlucInner(VELV) = UPrimFlucInner(VELV) * utau_ratio

    UPrimFlucInner(TKE ) = UPrimFlucInner(TKE ) * utau_ratio**2.0
    UPrimFlucInner(OMG ) = UPrimFlucInner(OMG ) * utau_ratio**2.0

    UPrimMean = UPrimMeanInner * (1.0 - weta) + UPrimMeanOuter * weta
    UPrimFluc = UPrimFlucInner * (1.0 - weta) + UPrimFlucOuter * weta

    UPrim_boundary(:,p,q) = (UPrim_boundary(:,p,q) * (cutEta-eta) + ( UPrimMean + UPrimFluc ) * (eta-1.0) ) / (cutEta - 1.0 )
    UPrim_boundary(VEL3,p,q) = -UPrim_boundary(VEL3,p,q)
  ELSE
    UPrimMeanInner(VEL1) = UPrimMeanInner(VEL1) * utau_ratio
    UPrimMeanInner(TKE ) = UPrimMeanInner(TKE ) * utau_ratio**2.0
    UPrimMeanInner(OMG ) = UPrimMeanInner(OMG ) * utau_ratio**2.0

    UPrimMeanOuter(VEL1) = UPrimMeanOuter(VEL1) * utau_ratio + u_inf * (1.0 - utau_ratio)

    UPrimFlucInner(VELV) = UPrimFlucInner(VELV) * utau_ratio
    UPrimFlucOuter(VELV) = UPrimFlucOuter(VELV) * utau_ratio

    UPrimFlucInner(TKE ) = UPrimFlucInner(TKE ) * utau_ratio**2.0
    UPrimFlucInner(OMG ) = UPrimFlucInner(OMG ) * utau_ratio**2.0

    UPrimMean = UPrimMeanInner * (1.0 - weta) + UPrimMeanOuter * weta
    UPrimFluc = UPrimFlucInner * (1.0 - weta) + UPrimFlucOuter * weta

    UPrim_boundary(:,p,q) = UPrimMean + UPrimFluc
    UPrim_boundary(VEL3,p,q) = -UPrim_boundary(VEL3,p,q)
  ENDIF

  UPrim_boundary(TKE,p,q) = MAX( UPrim_boundary(TKE,p,q), 1.e-6 )
  UPrim_boundary(OMG,p,q) = MAX( UPrim_boundary(OMG,p,q), 1.e-6 )

END DO; END DO

END SUBROUTINE FillState

SUBROUTINE BinarySearch(nVars, U, UDest, yDest, isUPrim)
USE MOD_Preproc
USE MOD_Globals
USE MOD_Recycling_Vars

INTEGER,INTENT(IN) :: nVars
REAL,INTENT(IN)    :: U(nVars,nY)
REAL,INTENT(OUT)   :: UDest(nVars)
REAL,INTENT(IN)    :: yDest
LOGICAL,OPTIONAL,INTENT(IN) :: isUPrim

INTEGER :: lower, upper, i
LOGICAL :: isUPrim2

IF (.NOT. PRESENT(isUPrim)) THEN
  isUPrim2 = .TRUE.
ELSE
  isUPrim2 = isUPrim
END IF

#if PP_NodeType != 2
IF (isUPrim2) THEN
  IF (yDest.LE.y(1)) THEN
    ! avoid extrapolation at the first layer
    UDest(DENS) = U(DENS,1)
    UDest(PRES) = U(PRES,1)
    UDest(TEMP) = U(TEMP,1)
    UDest(VELV) = U(VELV,1)
    UDest(TKE ) = U(TKE ,1)
    UDest(OMG ) = U(OMG ,1)
    RETURN
  END IF
END IF
#endif

IF (yDest.GE.y(nY)) THEN
  UDest = U(:,nY)
  RETURN
END IF

lower = 1
upper = nY
DO
  IF ((upper - lower) .EQ. 1) THEN
    UDest = ( &
        (y(upper) - yDest) * U(:,lower) + &
        (yDest - y(lower)) * U(:,upper) ) / &
        (y(upper) - y(lower))
    EXIT
  END IF
  i = (upper + lower) / 2
  IF (yDest .EQ. y(i)) THEN
    UDest = U(:,i)
    EXIT
  ELSE IF (yDest .GT. y(i)) THEN
    lower = i
  ELSE
    upper = i
  END IF
END DO

END SUBROUTINE BinarySearch

END MODULE MOD_Recycling
