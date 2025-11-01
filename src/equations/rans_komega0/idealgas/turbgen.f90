#include "flexi.h"
#include "eos.h"

MODULE MOD_TurbGen_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: useTurbGen = .FALSE.

REAL                :: RefState(PP_nVarPrim)

REAL,ALLOCATABLE    :: inflowData(:,:)
INTEGER             :: nRows
INTEGER,PARAMETER   :: nColumns = 7 + 4

INTEGER             :: nTurbGenSides
INTEGER,ALLOCATABLE :: sideToIdx(:)

INTEGER,PARAMETER   :: nFourierModes = 600
REAL,ALLOCATABLE    :: kappaHat(:,:)
REAL,ALLOCATABLE    :: sigmaHat(:,:)
REAL,ALLOCATABLE    :: psiHat(:)
REAL,ALLOCATABLE    :: magnitude(:,:,:,:)

! REAL,ALLOCATABLE    :: uHat(:,:,:,:)
! REAL,ALLOCATABLE    :: omegaHat(:,:,:,:)
! REAL,ALLOCATABLE    :: psiHatHat(:,:,:,:)
! REAL,ALLOCATABLE    :: kHat(:,:,:,:,:)
! REAL,ALLOCATABLE    :: sigmaHatHat(:,:,:,:,:)

! REAL
END MODULE MOD_TurbGen_Vars

MODULE MOD_TurbGen
! MODULES
IMPLICIT NONE
PUBLIC
CONTAINS

SUBROUTINE InitTurbGen()
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_TurbGen_Vars
  USE MOD_Mesh_Vars, ONLY: nBCSides,BC,BoundaryType,Face_xGP,nBCs
  IMPLICIT NONE

  INTEGER :: p, q, iSide
  INTEGER :: iMode
  REAL    :: StateExtended(nColumns)
  ! REAL    :: rho,u,v,w,p,k,o,R11,R22,R33,R12
  REAL    :: le, le_max
  REAL    :: k_min, k_max
  REAL    :: k_ener, k_curr
  REAL    :: magnitude_sum, delta_k
  INTEGER :: iBC

  nTurbGenSides = 0
  ALLOCATE(sideToIdx(nBCSides))
  DO iSide=1,nBCSides
    IF (BoundaryType(BC(iSide),BC_TYPE) .EQ. 31) THEN
      nTurbGenSides = nTurbGenSides + 1
      sideToIdx(iSide) = nTurbGenSides
    END IF
  END DO

  IF (nTurbGenSides .GT. 0) THEN
    useTurbGen = .TRUE.
  END IF
#if USE_MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, useTurbGen, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_FLEXI, iError)
#endif

  IF (.NOT. useTurbGen) RETURN

  ! RefState
  DO iBC=1,nBCs
    IF (BoundaryType(iBC,BC_TYPE) .EQ. 31) THEN
      RefState = BoundaryType(iBC,BC_STATE)
      EXIT
    END IF
  END DO

  CALL GenerateRandomNumbers()

  CALL ReadInflowData()

  ! le_max, k_min
  le_max = 0.0
  DO iSide=1,nBCSides
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL BinarySearch(StateExtended, Face_xGP(2,p,q,0,iSide))
      le = MIN( &
          2.0 * Face_xGP(2,p,q,0,iSide), &
          3.0 * (SQRT(StateExtended(6)) / (0.09 * StateExtended(7))))
      le_max = MAX(le_max, le)
    END DO; END DO
  END DO
#if USE_MPI
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, le_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FLEXI, iError)
#endif
  k_min = 0.5 * 6.28 / le_max

  ! k_max
  k_max = 6.28 / lEta(RefState)

  IF (k_max .LT. k_min) THEN
    CALL Abort(__STAMP__, "k_max < k_min")
  END IF

  DO iMode=1,nFourierModes
    kappaHat(:,iMode) = EXP(LOG(k_min) + (iMode - 1) / (nFourierModes - 1) * (LOG(k_max) - LOG(k_min))) * kappaHat(:,iMode)
  END DO

  ALLOCATE(magnitude(nFourierModes,0:PP_N,0:PP_NZ,nTurbGenSides))
  DO iSide=1,nBCSides
    IF (BoundaryType(BC(iSide),BC_TYPE) .NE. 31) CYCLE
    DO q=0,PP_NZ; DO p=0,PP_N
      CALL BinarySearch(StateExtended, Face_xGP(2,p,q,0,iSide))
      le = MIN( &
          2.0 * Face_xGP(2,p,q,0,iSide), &
          3.0 * (SQRT(StateExtended(6)) / (0.09 * StateExtended(7))))
      k_ener = 6.28 / le
      magnitude_sum = 0.0
      DO iMode=1,nFourierModes
        k_curr = EXP(LOG(k_min) + (iMode - 1) / (nFourierModes - 1) * (LOG(k_max) - LOG(k_min)))
        delta_k = &
            EXP(LOG(k_min) + (iMode - 1 + 0.5) / (nFourierModes - 1) * (LOG(k_max) - LOG(k_min))) - &
            EXP(LOG(k_min) + (iMode - 1 - 0.5) / (nFourierModes - 1) * (LOG(k_max) - LOG(k_min)))
        magnitude(iMode,p,q,sideToIdx(iSide)) = Ek(k_curr, k_ener, k_max) * delta_k
        magnitude_sum = magnitude_sum + magnitude(iMode,p,q,sideToIdx(iSide))
      END DO
      magnitude(:,p,q,sideToIdx(iSide)) = magnitude(:,p,q,sideToIdx(iSide)) / magnitude_sum
      magnitude(:,p,q,sideToIdx(iSide)) = 2.0 * SQRT(1.5) * SQRT(magnitude(:,p,q,sideToIdx(iSide)))
    END DO; END DO
  END DO
END SUBROUTINE InitTurbGen

SUBROUTINE CalcProfile(UPrim, X, t, p, q, iSide)
  USE MOD_TurbGen_Vars
  USE MOD_EOS_Vars, ONLY: sKappaM1,Kappa,KappaM1,R
  IMPLICIT NONE
  REAL,INTENT(OUT)   :: UPrim(PP_nVarPrim)
  REAL,INTENT(IN)    :: X(3)
  REAL,INTENT(IN)    :: t
  INTEGER,INTENT(IN) :: p
  INTEGER,INTENT(IN) :: q
  INTEGER,INTENT(IN) :: iSide
  REAL :: StateExtend(nColumns)
  REAL :: State(PP_nVarPrim)
  REAL :: X_new(3)
  REAL :: alpha
  INTEGER :: iMode
  REAL :: V(3)
  REAL :: R11, R22, R33, R12, AA(3,3)
  REAL :: c, Ma, pt, Tt

  CALL BinarySearch(StateExtend, X(2))
  State(1:7) = StateExtend(1:7)
  State(OMG) = 1.0 / SQRT(0.09 * State(7))
  ! State(TKE) = State(6)
  State(TKE) = 1.0e-8
  State(TEMP) = State(PRES) / (State(DENS) * R)
  UPrim = State

  c  = SQRT(kappa * UPrim(PRES) / UPrim(DENS))
  Ma = SQRT(DOT_PRODUCT(UPrim(VELV), UPrim(VELV))) / c
  pt = UPrim(PRES) * ((1.0 + 0.5 * (kappa - 1) * Ma * Ma)**(kappa * sKappaM1))
  Tt = UPrim(TEMP) * (1.0 + 0.5 * (kappa - 1) * Ma * Ma)

  X_new = X - t * State(VELV)
  V = 0.0
  DO iMode=1,nFourierModes
    alpha = DOT_PRODUCT(X_new, kappaHat(:,iMode)) + psiHat(iMode)
    V = V + magnitude(iMode,p,q,sideToIdx(iSide)) * COS(alpha) * sigmaHat(:,iMode)
  END DO
  R11 = StateExtend(PP_nVar+1)
  R22 = StateExtend(PP_nVar+2)
  R33 = StateExtend(PP_nVar+3)
  R12 = StateExtend(PP_nVar+4)
  AA = 0.0
  AA(1,1) = SQRT(MAX(R11, 1.0e-16))
  AA(2,1) = R12 / AA(1,1)
  AA(2,2) = SQRT(MAX(R22 - AA(2,1)**2, 1.0e-16))
  AA(3,3) = SQRT(MAX(R33, 1.0e-16))
  UPrim(VELV) = UPrim(VELV) + MATMUL(AA, V)

  UPrim(PRES) = pt / ((1.0 + 0.5 * (kappa - 1) * Ma * Ma)**(kappa * sKappaM1))
  UPrim(TEMP) = Tt / ( 1.0 + 0.5 * (kappa - 1) * Ma * Ma)
  UPrim(DENS) = UPrim(PRES) / (R * UPrim(TEMP))
END SUBROUTINE CalcProfile


SUBROUTINE FinalizeTurbGen()
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_TurbGen_Vars
  IMPLICIT NONE
  SDEALLOCATE(inflowData)
  SDEALLOCATE(sideToIdx)
  SDEALLOCATE(kappaHat)
  SDEALLOCATE(sigmaHat)
  SDEALLOCATE(psiHat)
  SDEALLOCATE(magnitude)
END SUBROUTINE FinalizeTurbGen


REAL FUNCTION lEta(RefState)
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_Viscosity
  REAL,INTENT(IN) :: RefState(PP_nVarPrim)
  REAL :: nu, rho, k, g
  rho = RefState(DENS)
  k   = Refstate(TKE)
  g   = Refstate(OMG)
  nu = VISCOSITY_TEMPERATURE(RefState(TEMP)) / rho
  lEta = (nu**3 / (0.09 * k * g**2))**0.25
END FUNCTION lEta


REAL FUNCTION Ek(k, kEner, kEta)
  IMPLICIT NONE
  REAL,INTENT(IN) :: k
  REAL,INTENT(IN) :: kEner
  REAL,INTENT(IN) :: kEta
  Ek = (k / kEner)**4 / (1.0 + (k / kEner)**2)**(17.0/6.0)
  Ek = Ek * EXP(-2.0 * (k / kEta)**2)
END FUNCTION Ek

SUBROUTINE GenerateRandomNumbers()
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_TurbGen_Vars, ONLY: nFourierModes, kappaHat, sigmaHat, psiHat
  IMPLICIT NONE
  REAL :: buffer(4)
  REAL :: phi, theta, alpha
  INTEGER :: iMode
  ALLOCATE(kappaHat(3,nFourierModes))
  ALLOCATE(sigmaHat(3,nFourierModes))
  ALLOCATE(psiHat(nFourierModes))
  IF (MPIRoot) THEN
    DO iMode=1,nFourierModes
      CALL RANDOM_NUMBER(buffer)
      theta      = ACOS(1.0 - 2.0 * buffer(1))
      alpha      = 6.28 * buffer(2)
      phi        = 6.28 * buffer(3)
      psiHat(iMode) = 6.28 * buffer(4)
      kappaHat(:,iMode) = (/SIN(theta) * COS(phi), SIN(theta) * SIN(phi), COS(theta)/)
      sigmaHat(:,iMode) = (/COS(alpha), SIN(alpha), 0.0/)
      buffer(1)      = sigmaHat(3,iMode)
      buffer(2)      = sigmaHat(1,iMode)
      sigmaHat(3,iMode) = buffer(1) * COS(theta) - buffer(2) * SIN(theta)
      sigmaHat(1,iMode) = buffer(1) * SIN(theta) + buffer(2) * COS(theta)
      buffer(1)      = sigmaHat(1,iMode)
      buffer(2)      = sigmaHat(2,iMode)
      sigmaHat(1,iMode) = buffer(1) * COS(phi) - buffer(2) * SIN(phi)
      sigmaHat(2,iMode) = buffer(1) * SIN(phi) + buffer(2) * COS(phi)
    END DO
  END IF
#if USE_MPI
  CALL MPI_BCAST(kappaHat, 3*nFourierModes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FLEXI, iError)
  CALL MPI_BCAST(sigmaHat, 3*nFourierModes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FLEXI, iError)
  CALL MPI_BCAST(psiHat,     nFourierModes, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FLEXI, iError)
#endif /*MPI*/
END SUBROUTINE GenerateRandomNumbers


SUBROUTINE ReadInflowData()
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_TurbGen_Vars, ONLY: inflowData, nRows, nColumns
  IMPLICIT NONE
  INTEGER :: iFile
  OPEN(newunit=iFile, file="inflowDistribution.dat", status="old", action="read")
  READ(iFile, *) nRows
  ALLOCATE(inflowData(0:nColumns,nRows))
  READ(iFile, *) inflowData
  CLOSE(iFile)
END SUBROUTINE ReadInflowData


SUBROUTINE BinarySearch(result, y)
  USE MOD_Globals
  USE MOD_PreProc
  USE MOD_TurbGen_Vars, ONLY: inflowData, nRows, nColumns
  IMPLICIT NONE
  REAL,INTENT(IN)  :: y
  REAL,INTENT(OUT) :: result(nColumns)
  INTEGER :: lower, upper, i

  lower = 1
  upper = nRows

  DO WHILE (.TRUE.)
    IF (upper - lower .EQ. 1) THEN
      result = ( &
          (inflowData(0,upper) - y) * inflowData(1:nColumns,lower) + &
          (y - inflowData(0,lower)) * inflowData(1:nColumns,upper)) / &
          (inflowData(0,upper) - inflowData(0,lower))
      EXIT
    END IF
    i = (upper + lower) / 2
    IF (y .EQ. inflowData(0,i)) THEN
      result = inflowData(1:nColumns,i)
      EXIT
    ELSE IF (y .GT. inflowData(0,i)) THEN
      lower = i
    ELSE
      upper = i
    END IF
  END DO
END SUBROUTINE BinarySearch

END MODULE MOD_TurbGen