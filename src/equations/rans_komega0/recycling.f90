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
CALL prms%CreateRealOption("recycling_xrec", "")
CALL prms%CreateRealOption("recycling_dens", "")
CALL prms%CreateRealOption("recycling_vel1", "")
CALL prms%CreateRealOption("recycling_pres", "")
CALL prms%CreateRealOption("recycling_freq", "")
CALL prms%CreateRealOption("recycling_momentum_thickness", "")
! CALL prms%CreateRealOption("recycling_Twall", "")

END SUBROUTINE DefineParametersRecycling

SUBROUTINE InitRecycling
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: Elem_xGP, Face_xGP, nElems, nBCSides, BoundaryType, BC
USE MOD_DG_Vars, ONLY: U, UPrim, UPrim_master
USE MOD_ReadInTools
USE MOD_Recycling_Vars
IMPLICIT NONE

INTEGER :: i,j,k,p,q,nn,iELem,iSide,iFile,ios
REAL :: tmp

doRecycling = GETLOGICAL("doRecycling")
IF (.NOT. doRecycling) THEN
  RETURN
END IF

xRecyc  = GETREAL("recycling_xrec")
thetaD  = GETREAL("recycling_momentum_thickness")
rho_inf = GETREAL("recycling_dens")
u_inf   = GETREAL("recycling_vel1")
p_inf   = GETREAL("recycling_pres")
timeavgFeq = GETREAL("recycling_freq")
! Twall   = GETREAL("recycling_Twall")

OPEN(newunit=iFile, file="inflowData.dat", status="old", action="read")
READ(iFile, *, iostat=ios) nInflowData
ALLOCATE(inflowData(1+PP_nVar,nInflowData))
READ(iFile, *, iostat=ios) inflowData
CLOSE(iFile)

! read y
nY = 0
OPEN(newunit=iFile, file="y.dat", status="old", action="read")
DO
  READ(iFile, *, iostat=ios) tmp
  IF (ios.NE.0) THEN
    EXIT
  END IF
  nY = nY + 1
END DO
REWIND(iFile)
ALLOCATE(y(nY))
DO nn=1,nY
  READ(iFile, *, iostat=ios) y(nn)
END DO
CLOSE(iFile)

! read z
nZ = 0
OPEN(newunit=iFile, file="z.dat", status="old", action="read")
DO
  READ(iFile, *, iostat=ios) tmp
  IF (ios.NE.0) THEN
    EXIT
  END IF
  nZ = nZ + 1
END DO
REWIND(iFile)
ALLOCATE(z(nZ))
DO nn=1,nZ
  READ(iFile, *, iostat=ios) z(nn)
END DO
CLOSE(iFile)

! initialize recycling plane
nRecycl = 0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  IF (ABS(Elem_xGP(1,i,j,k,iElem) - xRecyc).LE.1.0E-6) THEN
    nRecycl = nRecycl + 1
  END IF
END DO; END DO; END DO; END DO
ALLOCATE(rm(4,nRecycl))
nRecycl = 0
DO iELem=1,nElems; DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
  IF (ABS(Elem_xGP(1,i,j,k,iElem) - xRecyc).LE.1.0E-6) THEN
    nRecycl = nRecycl + 1
    rm(:,nRecycl) = (/i,j,k,iElem/)
  END IF
END DO; END DO; END DO; END DO

ALLOCATE(rmg(2,nRecycl))
DO nn=1,nRecycl
  DO j=1,nY
    IF (ABS(Elem_xGP(2,rm(1,nn),rm(2,nn),rm(3,nn),rm(4,nn)) - y(j)).LE.1.0E-6) THEN
      rmg(1,nn) = j
      EXIT
    END IF
  END DO
  DO k=1,nZ
    IF (ABS(Elem_xGP(3,rm(1,nn),rm(2,nn),rm(3,nn),rm(4,nn)) - z(k)).LE.1.0E-6) THEN
      rmg(2,nn) = k
      EXIT
    END IF
  END DO
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
img = 0
DO nn=1,nInflow
  DO j=1,nY
    IF (ABS(Face_xGP(2,im(1,nn),im(2,nn),0,im(3,nn)) - y(j)).LE.1.0E-6) THEN
      img(1,nn) = j
      EXIT
    END IF
  END DO
  DO k=1,nZ
    IF (ABS(Face_xGP(3,im(1,nn),im(2,nn),0,im(3,nn)) - z(k)).LE.1.0E-6) THEN
      img(2,nn) = k
      EXIT
    END IF
  END DO
END DO

CALL MPI_ALLREDUCE(nInflow, p, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_FLEXI, iError)
SWRITE(*,*) "nInflow = ", p
CALL MPI_ALLREDUCE(nRecycl, q, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_FLEXI, iError)
SWRITE(*,*) "nRecycl = ", q
IF (p.NE.q) THEN
  CALL Abort(__STAMP__,'nInflow /= nRecycl, check your recycling_xrec')
END IF

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

INQUIRE(file="recycling_mean.dat", exist=timeavgInit)
IF (timeavgInit) THEN
  OPEN(newunit=iFile,file="recycling_mean.dat",status="old",action="read")
  READ(iFile,*,iostat=ios) recycl_U_mean
  CLOSE(iFile)
  DO k=1,nZ; DO j=1,nY
    recycl_U_global(:,j,k) = recycl_U_mean(:,j)
  END DO; END DO
END IF

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
DO nn=1,nRecycl
  recycl_U_global(:,rmg(1,nn),rmg(2,nn)) = U(:,rm(1,nn),rm(2,nn),rm(3,nn),rm(4,nn))
END DO

nn = PP_nVar * nY * nZ
CALL MPI_ALLREDUCE(MPI_IN_PLACE, recycl_U_global, nn, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FLEXI, iError)

DO k=1,nZ; DO j=1,nY
  CALL ConsToPrim(recycl_UPrim_global(:,j,k), recycl_U_global(:,j,k))
END DO; END DO

IF (.NOT.timeavgInit) THEN
  timeavgInit = .TRUE.

  recycl_U_mean = 0.0
  DO k=1,nZ; DO j=1,nY
    recycl_U_mean(:,j) = recycl_U_mean(:,j) + recycl_U_global(:,j,k)
  END DO; END DO
  recycl_U_mean = recycl_U_mean / nZ
ELSE
  ! tmp = 1.0 - EXP(-RKb(CurrentStage)*dt/timeavgFeq)
  ! tmp = (RKb(CurrentStage) * dt) / timeavgFeq
  tmp = timeavgFeq

  recycl_U_mean = (1.0 - tmp) * recycl_U_mean
  DO k=1,nZ; DO j=1,nY
    recycl_U_mean(:,j) = recycl_U_mean(:,j) + tmp * recycl_U_global(:,j,k) / nZ
  END DO; END DO
END IF

DO j=1,nY
  CALL ConsToPrim(recycl_UPrim_mean(:,j), recycl_U_mean(:,j))
END DO

DO k=1,nZ; DO j=1,nY
  recycl_UPrim_fluc(:,j,k) = recycl_UPrim_global(:,j,k) - recycl_UPrim_mean(:,j)
END DO; END DO

tmp = 0.0
DO j=nY-5,nY
  ! tmp = tmp + recycl_U_mean(MOM1,j) / recycl_U_mean(DENS,j)
  tmp = tmp + recycl_UPrim_mean(VEL1,j)
END DO
tmp = tmp / 6
! tmp = u_inf

d99R = 0.0
DO j=2,nY
  IF (recycl_UPrim_mean(VEL1,j) .GE. (0.99 * tmp)) THEN
    u1 = recycl_UPrim_mean(VEL1,j-1)
    u2 = recycl_UPrim_mean(VEL1,j-0)
    d99R = ((u2 - (0.99 * tmp)) * y(j-1) + ((0.99 * tmp) - u1) * y(j)) / (u2 - u1)
    EXIT
  END IF
END DO

IF (d99R .EQ. 0.0) THEN
  SWRITE(*,*) "u_inf = ", tmp
  SWRITE(*,*) "u = ", recycl_UPrim_mean(VEL1,:)
  CALL Abort(__STAMP__,"d99 not found")
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

! TODO(Shimushu): remove hard-coded variables
d99I = thetaD
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
USE MOD_Recycling_Vars
USE MOD_ExactFunc_Vars, ONLY: BCData, BCLength
IMPLICIT NONE

REAL,INTENT(OUT)   :: UPrim_boundary(PP_nVarPrim,0:PP_N,0:PP_NZ)
REAL,INTENT(IN)    :: UPrim_master(  PP_nVarPrim,0:PP_N,0:PP_NZ)
INTEGER,INTENT(IN) :: iSide

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

  IF (eta .GE. 1.0) THEN
    weta = 1.0
    UPrimFlucOuter = 0.0
  END IF

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

  UPrimFluc(DENS) = bdamp * UPrimFluc(DENS)
  ! UPrimFluc(TEMP) = bdamp * UPrimFluc(TEMP)
  UPrimFluc(VELV) = bdamp * UPrimFluc(VELV)

  UPrimFluc(DENS) = MAX(MIN(UPrimFluc(DENS), 4.0 * UPrimMean(DENS)), -0.75 * UPrimMean(DENS))
  ! UPrimFluc(TEMP) = MAX(MIN(UPrimFluc(TEMP), 4.0 * UPrimMean(TEMP)), -0.75 * UPrimMean(TEMP))

  ! UPrim_boundary(:,p,q) = (UPrimMeanInner + UPrimFlucInner) * (1.0 - weta) + (UPrimMeanOuter + UPrimFlucOuter) * weta
  UPrim_boundary(:,p,q) = UPrimMean + UPrimFluc

  UPrim_boundary(PRES,p,q) = p_inf
  UPrim_boundary(TEMP,p,q) = UPrim_boundary(PRES,p,q) / (R * UPrim_boundary(DENS,p,q))
  ! UPrim_boundary(PRES,p,q) = UPrim_boundary(DENS,p,q) * R * UPrim_boundary(TEMP,p,q)

  UPrim_boundary(VEL3,p,q) = -UPrim_boundary(VEL3,p,q)
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
    UDest(DENS) = U(DENS,1) + (yDest - y(1)) / (y(2) - y(1)) * (U(DENS,2) - U(DENS,1))
    UDest(PRES) = U(PRES,1) + (yDest - y(1)) / (y(2) - y(1)) * (U(PRES,2) - U(PRES,1))
    UDest(TEMP) = U(TEMP,1) + (yDest - y(1)) / (y(2) - y(1)) * (U(TEMP,2) - U(TEMP,1))
    UDest(VELV) = (yDest / y(1)) * U(VELV,1)
    UDest(TKE ) = (yDest / y(1)) * U(TKE ,1)
    UDest(OMG ) = (yDest / y(1)) * U(OMG ,1)
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
