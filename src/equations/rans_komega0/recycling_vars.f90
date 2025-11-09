#include "flexi.h"
#include "eos.h"

MODULE MOD_Recycling_Vars
IMPLICIT NONE
PUBLIC
SAVE

LOGICAL             :: doRecycling = .FALSE.

LOGICAL             :: recyclingInit = .FALSE.

INTEGER             :: nInflow ! number of DOFs on the inflow plane in the core
INTEGER,ALLOCATABLE :: rm(:,:) ! from recycle index to i,j,k,iElem (4,nRecycl)
INTEGER,ALLOCATABLE :: im(:,:) ! from inflow index to p,q,iSide (3,nInflow)

INTEGER,ALLOCATABLE :: rim(:,:,:) ! from p,q,iSide to iInflow

INTEGER,ALLOCATABLE :: img(:,:) ! from index to global plane index (2,nInflow)
REAL,ALLOCATABLE    :: recycl_UPrim_global(:,:,:) ! global variables (8,Ny,Nz)
REAL,ALLOCATABLE    :: inflow_UPrim_global(:,:,:) ! global variables (8,Ny,Nz)
REAL,ALLOCATABLE    :: recycl_U_global(:,:,:) ! global variables (8,Ny,Nz)
REAL,ALLOCATABLE    :: inflow_U_global(:,:,:) ! global variables (8,Ny,Nz)

REAL,ALLOCATABLE    :: recycl_UPrim_mean(:,:) ! global variables (8,Ny)
REAL,ALLOCATABLE    :: recycl_UPrim_fluc(:,:,:) ! global variables (8,Ny)
REAL,ALLOCATABLE    :: recycl_U_mean(:,:) ! global variables (8,Ny)
REAL,ALLOCATABLE    :: recycl_U_fluc(:,:,:) ! global variables (8,Ny)

INTEGER             :: freeStreamRefState = 1

INTEGER             :: nY
INTEGER             :: nZ
INTEGER,ALLOCATABLE :: indexMap(:,:,:,:,:)
REAL,ALLOCATABLE    :: y(:)
REAL,ALLOCATABLE    :: z(:)
REAL                :: xRecyc

INTEGER             :: nInflowData
REAL,ALLOCATABLE    :: inflowData(:,:)

REAL                :: yMatchingTolerance
REAL                :: zMatchingTolerance

REAL                :: thetaR = 1.0E15
REAL                :: thetaD = 0.1 * 4.0E-3
REAL                :: d99R, d99I
REAL                :: uTauR, uTauI
REAL                :: dnuR, dnuI

LOGICAL             :: timeavgInit = .FALSE.

REAL                :: rho_inf = 0.1
REAL                :: u_inf = 823.6
REAL                :: p_inf = 270.0

REAL                :: timeavgFeq = 1.0E-6

LOGICAL             :: recycling_called = .FALSE.

REAL                :: utau_ratio

REAL                :: Twall

END MODULE MOD_Recycling_Vars