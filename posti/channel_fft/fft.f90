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

MODULE MOD_FFT
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitFFT
  MODULE PROCEDURE InitFFT
END INTERFACE

INTERFACE FinalizeFFT
  MODULE PROCEDURE FinalizeFFT
END INTERFACE

INTERFACE PerformFFT
  MODULE PROCEDURE PerformFFT
END INTERFACE

INTERFACE FFTOutput
  MODULE PROCEDURE FFTOutput
END INTERFACE

INTERFACE PrimStateAtFFTCoords
  MODULE PROCEDURE PrimStateAtFFTCoords
END INTERFACE

PUBLIC:: InitFFT,FinalizeFFT,PerformFFT,FFTOutput,PrimStateAtFFTCoords

CONTAINS

!===================================================================================================================================
!> Initialize FFT. Read in user-defined parameters, allocate arrays for solution and auxilliary FFT variables.
!> Prepare the Vandermonde to interpolate from the state to the FFT grid. Prepare coordinates of FFT grid.
!===================================================================================================================================
SUBROUTINE InitFFT()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE FFTW3
USE MOD_FFT_Vars
USE MOD_Commandline_Arguments
USE MOD_Interpolation,           ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_ReadInTools,             ONLY: GETREAL,GETINT,GETLOGICAL,GETSTR
USE MOD_Mesh_Vars,               ONLY: nElems_IJK,nElems
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Interpolation_Vars ,     ONLY: NodeType,NodeTypeVISUInner
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
!===================================================================================================================================
! Read in user-defined parameters
OutputFormat = GETINT('OutputFormat','0')
NCalc  = GETINT('NCalc')          ! Polynomial degree to perfrom DFFT on
Re_tau = GETREAL('Re_tau')        ! Reynolds number bases on friction velocity and channel half height
ReadMean = GETLOGICAL('ReadMean') ! Read TimeAvg instead of State
OutputNodeType = GETSTR('OutputNodeType', 'VISU_INNER')
NMean = -1
IF (ReadMean) THEN
  NMean = GETINT('NMean', '-1')
END IF
IF (NMean .EQ. -1) THEN
  NMean = PP_nVar
ENDIF
Permutation = GETINT('Permutation', '0')
ALLOCATE(PIJK(3,0:NCalc,0:NCalc,0:NCalc))
DO k=0,NCalc
  DO j=0,NCalc
    DO i=0,NCalc
      IF (Permutation.EQ.0) THEN
        PIJK(:,i,j,k) = (/i,j,k/)
      ELSE IF (Permutation.EQ.1) THEN
        PIJK(:,i,j,k) = (/i,k,j/)
      ELSE IF (Permutation.EQ.2) THEN
        PIJK(:,i,j,k) = (/j,i,k/)
      ELSE IF (Permutation.EQ.3) THEN
        PIJK(:,i,j,k) = (/j,k,i/)
      ELSE IF (Permutation.EQ.4) THEN
        PIJK(:,i,j,k) = (/k,i,j/)
      ELSE IF (Permutation.EQ.5) THEN
        PIJK(:,i,j,k) = (/k,j,i/)
      END IF
    END DO
  END DO
END DO

! Allocate solution array in DG vars, used to store state
ALLOCATE(U(NMean,0:PP_N,0:PP_N,0:PP_N,nElems))

! Allocate auxilliary arrays for FFT
N_FFT(1)=(NCalc+1)*nElems_IJK(1)
N_FFT(2)=(NCalc+1)*nElems_IJK(2)
N_FFT(3)=(NCalc+1)*nElems_IJK(3)
nSamplesI = (NCalc+1)*nElems_IJK(1)
nSamples_specI  = INT((nSamplesI)/2)+1
nSamplesK = (NCalc+1)*nElems_IJK(3)
nSamples_specK  = INT((nSamplesK)/2)+1
ALLOCATE(inI(nSamplesI))
ALLOCATE(outI(nSamplesI))
ALLOCATE(inK(nSamplesK))
ALLOCATE(outK(nSamplesK))
ALLOCATE(Ex_uu(N_FFT(2),nSamples_specI))
ALLOCATE(Ex_vv(N_FFT(2),nSamples_specI))
ALLOCATE(Ex_ww(N_FFT(2),nSamples_specI))
ALLOCATE(Ex_pp(N_FFT(2),nSamples_specI))
ALLOCATE(Ez_uu(N_FFT(2),nSamples_specK))
ALLOCATE(Ez_vv(N_FFT(2),nSamples_specK))
ALLOCATE(Ez_ww(N_FFT(2),nSamples_specK))
ALLOCATE(Ez_pp(N_FFT(2),nSamples_specK))
Ex_uu=0.;Ex_vv=0.;Ex_ww=0.;Ex_pp=0.
Ez_uu=0.;Ez_vv=0.;Ez_ww=0.;Ez_pp=0.
ALLOCATE(U_FFT(1:NMean,1:N_FFT(1),1:N_FFT(2),1:N_FFT(3)))
ALLOCATE(MS_PSD(N_FFT(2),4))
ALLOCATE(MS_t(N_FFT(2),3))
ALLOCATE(M_t(N_FFT(2),NMean))
! #if EQNSYSNR == 2
! ! navierstokes
! ALLOCATE(M_t(N_FFT(2),5))
! #elif EQNSYSNR == 4
! ! rans-kg
! ALLOCATE(M_t(N_FFT(2),7))
! #endif
MS_t=0.
MS_PSD=0.
M_t=0.

! Get vandermonde from computation N on Gauss or Gauss-Lobatto points to equidistand FFT points
ALLOCATE(VdmGaussEqui(0:NCalc,0:PP_N))
CALL GetVandermonde(PP_N,NodeType,NCalc,OutputNodeType,VdmGaussEqui)
ALLOCATE(xGP(0:NCalc))
ALLOCATE(wGP(0:NCalc))
CALL GetNodesAndWeights(NCalc,OutputNodeType,xGP,wGP)
wGP = wGP * REAL(NCalc + 1) / 2.0

CALL DFFTW_PLAN_DFT_1D(planI,nSamplesI,inI,outI,FFTW_FORWARD,FFTW_ESTIMATE)
CALL DFFTW_PLAN_DFT_1D(planK,nSamplesK,inK,outK,FFTW_FORWARD,FFTW_ESTIMATE)
CALL FFTCoords()

END SUBROUTINE InitFFT

!===================================================================================================================================
!> Interpolate the coordinates of the mesh to the FFT grid.
!===================================================================================================================================
SUBROUTINE FFTCoords
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FFT_Vars
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars             ,ONLY: Elem_IJK,Elem_xGP,nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: ii,jj,kk,i,j,k,iElem
INTEGER                             :: iGlob,jGlob,kGlob
REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: X_aux
!===================================================================================================================================
ALLOCATE(X_aux(1:3,0:NCalc,0:NCalc,0:NCalc))
ALLOCATE(X_FFT(1:3,1:N_FFT(1),1:N_FFT(2),1:N_FFT(3)))
X_FFT=0.
DO iElem=1,nGlobalElems
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  CALL ChangeBasis3D(3,PP_N,NCalc,VdmGaussEqui,Elem_xGP(:,:,:,:,iElem),X_aux)
  ! DO k=0,NCalc
  !   kGlob=(kk-1)*(NCalc+1)+k+1
  !   DO j=0,NCalc
  !     jGlob=(jj-1)*(NCalc+1)+j+1
  !     DO i=0,NCalc
  !       iGlob=(ii-1)*(NCalc+1)+i+1
  !       X_FFT(:,iGlob,jGlob,kGlob)=X_aux(:,i,j,k)
  !     END DO
  !   END DO
  ! END DO
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        X_FFT(:,iGlob,jGlob,kGlob)=X_aux(:,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k))
      END DO
    END DO
  END DO
END DO

DEALLOCATE(X_aux)
END SUBROUTINE FFTCoords

!===================================================================================================================================
!> Interpolate the state to the FFT grid and convert solution to primite variables.
!===================================================================================================================================
SUBROUTINE PrimStateAtFFTCoords
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_FFT_Vars
USE MOD_ChangeBasis                 ,ONLY: ChangeBasis3D
USE MOD_Mesh_Vars                   ,ONLY: Elem_IJK,nGlobalElems
USE MOD_DG_Vars                     ,ONLY: U
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: ii,jj,kk,i,j,k,iElem
INTEGER                             :: iGlob,jGlob,kGlob
REAL                                :: kappaM1
REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: U_aux
!===================================================================================================================================
! 1:5 does not work for RANS eq
ALLOCATE(U_aux(1:NMean,0:NCalc,0:NCalc,0:NCalc))

KappaM1=0.4
U_FFT=0.
DO iElem=1,nGlobalElems
  ii=Elem_IJK(1,iElem)
  jj=Elem_IJK(2,iElem)
  kk=Elem_IJK(3,iElem)
  CALL ChangeBasis3D(NMean,PP_N,NCalc,VdmGaussEqui,U(:,:,:,:,iElem),U_aux)
  DO k=0,NCalc
    kGlob=(kk-1)*(NCalc+1)+k+1
    DO j=0,NCalc
      jGlob=(jj-1)*(NCalc+1)+j+1
      DO i=0,NCalc
        iGlob=(ii-1)*(NCalc+1)+i+1
        U_FFT(1,iGlob,jGlob,kGlob)=U_aux(1,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k))
        U_FFT(2:4,iGlob,jGlob,kGlob)= &
            U_aux(2:4,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) / &
            U_aux(1,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k))
#if EQNSYSNR == 2
        ! navierstokes
        U_FFT(5,iGlob,jGlob,kGlob)=KappaM1*( &
            U_aux(5,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) - &
            U_aux(1,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) * &
            SUM(U_FFT(2:4,iGlob,jGlob,kGlob)**2)*0.5)
#elif EQNSYSNR == 4
        ! rans-kg
        U_FFT(5,iGlob,jGlob,kGlob)=KappaM1*( &
            U_aux(5,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) - &
            U_aux(6,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) - &
            U_aux(1,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) * &
            SUM(U_FFT(2:4,iGlob,jGlob,kGlob)**2)*0.5)
        U_FFT(6:7,iGlob,jGlob,kGlob)= & 
            U_aux(6:7,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k)) / &
            U_aux(1,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k))
#endif
      IF (NMean .GT. PP_nVar) THEN
        U_FFT(PP_nVar+1:NMean,iGlob,jGlob,kGlob)= &
            U_aux(PP_nVar+1:NMean,PIJK(1,i,j,k),PIJK(2,i,j,k),PIJK(3,i,j,k))
      END IF
      END DO
    END DO
  END DO
END DO

DEALLOCATE(U_Aux)
END SUBROUTINE PrimStateAtFFTCoords

!===================================================================================================================================
!> Calls the actual FFT routines, called for each statefile.
!===================================================================================================================================
SUBROUTINE PerformFFT()
! MODULES                                                                                                                          !
USE FFTW3
USE MOD_FFT_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: j,k,ii,jj,kk
!===================================================================================================================================
DO j=1,N_FFT(2)
  DO k=1,N_FFT(3)
    ! line FFT along i-direction (x in our channel testcases)
    ! output is amplitude squares of U-spectrum, i.E. U**2(k)
    ! other variables/directions accordingly
    CALL FFT(planI,inI,outI,nSamplesI,nSamples_specI,U_FFT(2,1:nSamplesI,j,k),Ex_uu(j,1:nSamples_specI))
    CALL FFT(planI,inI,outI,nSamplesI,nSamples_specI,U_FFT(3,1:nSamplesI,j,k),Ex_vv(j,1:nSamples_specI))
    CALL FFT(planI,inI,outI,nSamplesI,nSamples_specI,U_FFT(4,1:nSamplesI,j,k),Ex_ww(j,1:nSamples_specI))
    CALL FFT(planI,inI,outI,nSamplesI,nSamples_specI,U_FFT(5,1:nSamplesI,j,k),Ex_pp(j,1:nSamples_specI))
  END DO
END DO
DO j=1,N_FFT(2)
  DO k=1,N_FFT(1)
    CALL FFT(planK,inK,outK,nSamplesK,nSamples_specK,U_FFT(2,k,j,1:nSamplesK),Ez_uu(j,1:nSamples_specK))
    CALL FFT(planK,inK,outK,nSamplesK,nSamples_specK,U_FFT(3,k,j,1:nSamplesK),Ez_vv(j,1:nSamples_specK))
    CALL FFT(planK,inK,outK,nSamplesK,nSamples_specK,U_FFT(4,k,j,1:nSamplesK),Ez_ww(j,1:nSamples_specK))
    CALL FFT(planK,inK,outK,nSamplesK,nSamples_specK,U_FFT(5,k,j,1:nSamplesK),Ez_pp(j,1:nSamples_specK))
  END DO
END DO
DO j=1,N_FFT(2)
  ! For the mean of u,v,w,p compute the sum over J-Planes
  ! (==Y-Planes in our channel testcase) and add to the
  ! sum of previous state files
  DO ii=1,N_FFT(1); DO kk=1,N_FFT(3)
    M_t(j,1)=M_t(j,1)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(2,ii,j,kk)
    M_t(j,2)=M_t(j,2)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(3,ii,j,kk)
    M_t(j,3)=M_t(j,3)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(4,ii,j,kk)
    M_t(j,4)=M_t(j,4)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(5,ii,j,kk)
    M_t(j,5)=M_t(j,5)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(1,ii,j,kk)      
  END DO; END DO
  ! M_t(j,1)=M_t(j,1)+SUM(U_FFT(2,1:N_FFT(1),j,1:N_FFT(3)))
  ! M_t(j,2)=M_t(j,2)+SUM(U_FFT(3,1:N_FFT(1),j,1:N_FFT(3)))
  ! M_t(j,3)=M_t(j,3)+SUM(U_FFT(4,1:N_FFT(1),j,1:N_FFT(3)))
  ! M_t(j,4)=M_t(j,4)+SUM(U_FFT(5,1:N_FFT(1),j,1:N_FFT(3)))
  ! M_t(j,5)=M_t(j,5)+SUM(U_FFT(1,1:N_FFT(1),j,1:N_FFT(3)))
#if EQNSYSNR == 2
  ! navierstokes
#elif EQNSYSNR == 4
  ! rans-kg
  DO ii=1,N_FFT(1); DO kk=1,N_FFT(3)
    M_t(j,6)=M_t(j,6)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(6,ii,j,kk)
    M_t(j,7)=M_t(j,7)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(7,ii,j,kk)
  END DO; END DO
  ! M_t(j,6)=M_t(j,6)+SUM(U_FFT(6,1:N_FFT(1),j,1:N_FFT(3)))
  ! M_t(j,7)=M_t(j,7)+SUM(U_FFT(7,1:N_FFT(1),j,1:N_FFT(3)))
#endif
  IF (NMean .GT. PP_nVar) THEN
    DO k=PP_nVar+1,NMean
      DO ii=1,N_FFT(1); DO kk=1,N_FFT(3)
        M_t(j,k)=M_t(j,K)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(k,ii,j,kk)
      END DO; END DO
      ! M_t(j,k)=M_t(j,k)+SUM(U_FFT(k,1:N_FFT(1),j,1:N_FFT(3)))
    END DO
  END IF
  !  For mean uv
  DO ii=1,N_FFT(1); DO kk=1,N_FFT(3)
    MS_t(j,1)=MS_t(j,1)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(2,ii,j,kk)*U_FFT(3,ii,j,kk)
    MS_t(j,2)=MS_t(j,2)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(2,ii,j,kk)*U_FFT(4,ii,j,kk)
    MS_t(j,3)=MS_t(j,3)+wGP(MOD(ii-1,NCalc+1))*wGP(MOD(kk-1,NCalc+1))*U_FFT(3,ii,j,kk)*U_FFT(4,ii,j,kk)
  END DO; END DO
  ! MS_t(j,1)=MS_t(j,1)+SUM(U_FFT(2,1:N_FFT(1),j,1:N_FFT(3))*U_FFT(3,1:N_FFT(1),j,1:N_FFT(3)))
  ! MS_t(j,2)=MS_t(j,2)+SUM(U_FFT(2,1:N_FFT(1),j,1:N_FFT(3))*U_FFT(4,1:N_FFT(1),j,1:N_FFT(3)))
  ! MS_t(j,3)=MS_t(j,3)+SUM(U_FFT(3,1:N_FFT(1),j,1:N_FFT(3))*U_FFT(4,1:N_FFT(1),j,1:N_FFT(3)))
END DO

END SUBROUTINE PerformFFT

!===================================================================================================================================
!> Low-level wrapper routine to DFFTW call.
!> Does single sided FFT spectrum of input U_in
!> Computes square of amplitude spectrum and adds it to output variable u_hat
!===================================================================================================================================
SUBROUTINE FFT(plan,in,out,nSamples,nSamples_spec,U_in,U_hat)
! MODULES
USE MOD_Globals
USE FFTW3
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)       ::plan
INTEGER,INTENT(IN)               ::nSamples
INTEGER,INTENT(IN)               ::nSamples_spec
COMPLEX,INTENT(INOUT)            ::in(nSamples)
COMPLEX,INTENT(INOUT)            ::out(nSamples)
REAL,INTENT(IN)                  ::U_in(1:nSamples)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)               ::U_hat(nSamples_spec)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
in=U_in
CALL DFFTW_EXECUTE_DFT(plan, in, out)
!for square of single sided spectrum we need "2*amplitude**2"
out(2:nSamples_spec-1)=2*(1./REAL(nSamples)*ABS(out(2:nSamples_spec-1)))**2
out(1)=ABS(out(1)/REAL(nSamples))**2 !mean value, not two sided (unique)
IF(MOD(nSamples,2).EQ.0.)THEN
  !Even number of samples,single nyquist stored in one index
  out(nSamples_spec)=ABS(out(nSamples_spec)/REAL(nSamples))**2
ELSE
  !Odd number of samples, nyquist stored two sided (just as any feq)
  out(nSamples_spec)=2.*(ABS(out(nSamples_spec))/REAL(nSamples))**2
END IF
!add squared amplitude
U_hat(2:nSamples_spec)=REAL(U_hat(2:nSamples_spec)+out(2:nSamples_spec))
!sum of mean square into first index
U_hat(1)=REAL(U_hat(1)+SUM(out(1:nSamples_spec)))

END SUBROUTINE FFT

!===================================================================================================================================
!> Average the results and write the output files
!===================================================================================================================================
SUBROUTINE FFTOutput()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_FFT_Vars
USE MOD_Commandline_Arguments
USE MOD_IO_HDF5
USE MOD_HDF5_Output
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,Fileunit_EK,offsetVar,iVar,nVal,nSamples
CHARACTER(LEN=255) :: FileName_EK
CHARACTER(LEN=255) :: ZoneTitle
CHARACTER(LEN=255) :: tmp255
LOGICAL            :: connected
REAL,ALLOCATABLE   :: PointData(:,:)
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesFFT(:)
!===================================================================================================================================
!Sum up the lower and upper channel halfes
PRINT *, "The first direction"
DO i=1,N_FFT(1)
  PRINT *, X_FFT(1,i,1,1), X_FFT(2,i,1,1), X_FFT(3,i,1,1)
END DO
PRINT *, "The second direction"
DO i=1,N_FFT(2)
  PRINT *, X_FFT(1,1,i,1), X_FFT(2,1,i,1), X_FFT(3,1,i,1)
END DO
PRINT *, "The third direction"
DO i=1,N_FFT(3)
  PRINT *, X_FFT(1,1,1,i), X_FFT(2,1,1,i), X_FFT(3,1,1,i)
END DO
DO j=N_FFT(2)/2+1,N_FFT(2)
  Ex_uu(j,:)=(Ex_uu(N_FFT(2)-j+1,:)+Ex_uu(j,:))
  Ex_vv(j,:)=(Ex_vv(N_FFT(2)-j+1,:)+Ex_vv(j,:))
  Ex_ww(j,:)=(Ex_ww(N_FFT(2)-j+1,:)+Ex_ww(j,:))
  Ex_pp(j,:)=(Ex_pp(N_FFT(2)-j+1,:)+Ex_pp(j,:))
  Ez_uu(j,:)=(Ez_uu(N_FFT(2)-j+1,:)+Ez_uu(j,:))
  Ez_vv(j,:)=(Ez_vv(N_FFT(2)-j+1,:)+Ez_vv(j,:))
  Ez_ww(j,:)=(Ez_ww(N_FFT(2)-j+1,:)+Ez_ww(j,:))
  Ez_pp(j,:)=(Ez_pp(N_FFT(2)-j+1,:)+Ez_pp(j,:))
  MS_t(j,:)=(MS_t(N_FFT(2)-j+1,:)-MS_t(j,:))
  M_t(j,1)=(M_t(N_FFT(2)-j+1,1)+M_t(j,1))
  M_t(j,2:3)=(M_t(N_FFT(2)-j+1,2:3)-M_t(j,2:3))
  M_t(j,4:5)=(M_t(N_FFT(2)-j+1,4:5)+M_t(j,4:5))
#if EQNSYSNR == 2
  ! navierstokes
#elif EQNSYSNR == 4
  ! rans-kg
  M_t(j,6:7)=(M_t(N_FFT(2)-j+1,6:7)+M_t(j,6:7))
#endif
  IF (NMean .GT. PP_nVar) THEN
    M_t(j,PP_nVar+1:NMean)=(M_t(N_FFT(2)-j+1,PP_nVar+1:NMean)+M_t(j,PP_nVar+1:NMean))
  END IF
END DO
!Mean of u,v,w,p
M_t=M_t/((nArgs-1)*N_FFT(1)*N_FFT(3)*2)
!Mean squares fluctuations in Y-planes of u,v,w,p
MS_PSD(:,1)=Ex_uu(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,1))**2
MS_PSD(:,2)=Ex_vv(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,2))**2
MS_PSD(:,3)=Ex_ww(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,3))**2
MS_PSD(:,4)=Ex_pp(:,1)/((nArgs-1)*N_FFT(3)*2)-(M_t(:,4))**2
!Mean square of uv,uw,vw in Y-planes, not computed from spectra as
!there are non: mathematically equal.
MS_t(:,1)=  MS_t(:,1)  /((nArgs-1)*N_FFT(1)*N_FFT(3)*2)-(M_t(:,1)*M_t(:,2))
MS_t(:,2)=  MS_t(:,2)  /((nArgs-1)*N_FFT(1)*N_FFT(3)*2)-(M_t(:,1)*M_t(:,3))
MS_t(:,3)=  MS_t(:,3)  /((nArgs-1)*N_FFT(1)*N_FFT(3)*2)-(M_t(:,2)*M_t(:,3))
!TWO Times meansquare summed so far over time and z/x => amplitude squares wanted
!z/x :: NFFT(3/1) ;  time :: nArgs-1 ; 1/2 for correct mean square
Ex_uu(:,2:nSamples_SpecI)=Ex_uu(:,2:nSamples_SpecI)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_vv(:,2:nSamples_SpecI)=Ex_vv(:,2:nSamples_SpecI)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_ww(:,2:nSamples_SpecI)=Ex_ww(:,2:nSamples_SpecI)/(N_FFT(3)*REAL(nArgs-1)*2)
Ex_pp(:,2:nSamples_SpecI)=Ex_pp(:,2:nSamples_SpecI)/(N_FFT(3)*REAL(nArgs-1)*2)
Ez_uu(:,2:nSamples_SpecK)=Ez_uu(:,2:nSamples_SpecK)/(N_FFT(1)*REAL(nArgs-1)*2)
Ez_vv(:,2:nSamples_SpecK)=Ez_vv(:,2:nSamples_SpecK)/(N_FFT(1)*REAL(nArgs-1)*2)
Ez_ww(:,2:nSamples_SpecK)=Ez_ww(:,2:nSamples_SpecK)/(N_FFT(1)*REAL(nArgs-1)*2)
Ez_pp(:,2:nSamples_SpecK)=Ez_pp(:,2:nSamples_SpecK)/(N_FFT(1)*REAL(nArgs-1)*2)

SELECT CASE(OutputFormat)
  CASE(0) !Write mean square fluctuations data and mean velocity to Tecplot format
    IF (Nmean .NE. PP_nVar) THEN
      CALL CollectiveStop(__STAMP__,'ERROR - Tecplot output is not supported when NMean != PP_nVar.')
    ENDIF
    FileUnit_EK=155
    INQUIRE(UNIT=FileUnit_EK, OPENED=connected)
    IF (Connected) THEN
      DO
        FileUnit_EK=Fileunit_EK+1
        INQUIRE(UNIT=FileUnit_EK, OPENED=connected)
        IF (.NOT. connected) EXIT
      END DO
    END IF
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_MS',time)
    FileName_EK=TRIM(Filename_EK)//'.dat'
    OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
    WRITE(FileUnit_EK,*)'TITLE     = "MeanSquares "'
    WRITE(FileUnit_EK,'(a)')'VARIABLES ="yPlus"'
    WRITE(FileUnit_EK,'(a)')'"uu"'
    WRITE(FileUnit_EK,'(a)')'"vv"'
    WRITE(FileUnit_EK,'(a)')'"ww"'
    WRITE(FileUnit_EK,'(a)')'"uv"'
    WRITE(FileUnit_EK,'(a)')'"uw"'
    WRITE(FileUnit_EK,'(a)')'"vw"'
    WRITE(FileUnit_EK,'(a)')'"u_mean"'
    WRITE(FileUnit_EK,'(a)')'"v_mean"'
    WRITE(FileUnit_EK,'(a)')'"w_mean"'
    WRITE(FileUnit_EK,'(a)')'"p_mean"'
    WRITE(FileUnit_EK,'(a)')'"rho_mean"'
#if EQNSYSNR == 2
    ! navierstokes
#elif EQNSYSNR == 4
    ! rans-kg
    WRITE(FileUnit_EK,'(a)')'"k_mean"'
    WRITE(FileUnit_EK,'(a)')'"g_mean"'
#endif
    WRITE(FileUnit_EK,*)'ZONE T="',ProjectName,'"'
    WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
    WRITE(FileUnit_EK,*)' I=',N_FFT(2)/2,', J=1, K=1, ZONETYPE=Ordered'
    WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
#if EQNSYSNR == 2
    ! navierstokes
    WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
#elif EQNSYSNR == 4
    ! rans-kg
    WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
#endif

    DO j=N_FFT(2)/2+1,N_FFT(2)
#if EQNSYSNR == 2
      ! navierstokes
      WRITE(FileUnit_EK,'(12(E20.12,X))')(1-(X_FFT(2,1,j,1)))*Re_tau,MS_PSD(j,1),MS_PSD(j,2),MS_PSD(j,3),MS_t(j,1),MS_t(j,2),&
        MS_t(j,3),M_t(j,1),M_t(j,2),M_t(j,3),M_t(j,4),M_t(j,5)
#elif EQNSYSNR == 4
      ! rans-kg
      WRITE(FileUnit_EK,'(12(E20.12,X))')(1-(X_FFT(2,1,j,1)))*Re_tau,MS_PSD(j,1),MS_PSD(j,2),MS_PSD(j,3),MS_t(j,1),MS_t(j,2),&
        MS_t(j,3),M_t(j,1),M_t(j,2),M_t(j,3),M_t(j,4),M_t(j,5),M_t(j,6),M_t(j,7)
#endif

    END DO
    CLOSE(FILEUnit_EK)
    !-------------------------------------------------
    !write energy spectra in x-direction
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_x',time)
    FileName_EK=TRIM(Filename_EK)//'.dat'
    OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
    WRITE(FileUnit_EK,*)'TITLE     = "Energy Spectra_x "'
    WRITE(FileUnit_EK,'(a)')'VARIABLES ="k"'
    WRITE(FileUnit_EK,'(a)')'"E_uu_x"'
    WRITE(FileUnit_EK,'(a)')'"E_vv_x"'
    WRITE(FileUnit_EK,'(a)')'"E_ww_x"'
    WRITE(FileUnit_EK,'(a)')'"E_pp_x"'
    DO j=N_FFT(2)/2+1,N_FFT(2)
      WRITE(FileUnit_EK,'(A8,A,A7,I4,A1)')'ZONE T="',TRIM(ProjectName),',yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau),'"'
      WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
      WRITE(FileUnit_EK,*)' I=',nSamples_specI,', J=1, K=1, ZONETYPE=Ordered'
      WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
      WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
      DO k=1,nSamples_specI
        WRITE(FileUnit_EK,'((I5.5,X),4(E20.12,X))')k-1,Ex_uu(j,k),Ex_vv(j,k),Ex_ww(j,k),Ex_pp(j,k)
      END DO
    END DO
    CLOSE(FILEUnit_EK)
    !-------------------------------------------------
    !write energy spectra in z-direction
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_z',time)
    FileName_EK=TRIM(Filename_EK)//'.dat'
    OPEN(FileUnit_Ek,FILE=Filename_EK,STATUS="REPLACE")
    WRITE(FileUnit_EK,*)'TITLE     = "Energy Spectra_z "'
    WRITE(FileUnit_EK,'(a)')'VARIABLES ="k"'
    WRITE(FileUnit_EK,'(a)')'"E_uu_z"'
    WRITE(FileUnit_EK,'(a)')'"E_vv_z"'
    WRITE(FileUnit_EK,'(a)')'"E_ww_z"'
    WRITE(FileUnit_EK,'(a)')'"E_pp_z"'
    DO j=N_FFT(2)/2+1,N_FFT(2)
      WRITE(FileUnit_EK,'(A8,A,A7,I4,A1)')'ZONE T="',TRIM(ProjectName),',yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau),'"'
      WRITE(FileUnit_EK,'(a)')' STRANDID=0, SOLUTIONTIME=0'
      WRITE(FileUnit_EK,*)' I=',nSamples_specK,', J=1, K=1, ZONETYPE=Ordered'
      WRITE(FileUnit_EK,'(a)')' DATAPACKING=POINT'
      WRITE(FileUnit_EK,'(a)')' DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
      DO k=1,nSamples_specK
        WRITE(FileUnit_EK,'((I5.5,X),4(E20.12,X))')2*(k-1),Ez_uu(j,k),Ez_vv(j,k),Ez_ww(j,k),Ez_pp(j,k)
      END DO
    END DO
    CLOSE(FILEUnit_EK)
  CASE(2) !Write mean square fluctuations data and mean velocity to HDF5 format
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_MS',time)
    FileName_EK=TRIM(Filename_EK)//'.h5'

    nVal = 7 + Nmean
! #if EQNSYSNR == 2
!     ! navierstokes
!     nVal=12
! #elif EQNSYSNR == 4
!     ! rans-kg
!     nVal=14
! #endif
    nSamples=N_FFT(2)/2
    ALLOCATE(VarNamesFFT(nVal))
    ALLOCATE(PointData(nVal,nSamples))
    VarNamesFFT(1) ='yPlus'
    VarNamesFFT(2) ='uu'
    VarNamesFFT(3) ='vv'
    VarNamesFFT(4) ='ww'
    VarNamesFFT(5) ='uv'
    VarNamesFFT(6) ='uw'
    VarNamesFFT(7) ='vw'
    VarNamesFFT(8) ='u_mean'
    VarNamesFFT(9) ='v_mean'
    VarNamesFFT(10)='w_mean'
    VarNamesFFT(11)='p_mean'
    VarNamesFFT(12)='rho_mean'
#if EQNSYSNR == 2
    ! navierstokes
#elif EQNSYSNR == 4
    ! rans-kg
    VarNamesFFT(13)='k_mean'
    VarNamesFFT(14)='g_mean'
#endif
    IF (Nmean.GT.PP_nVar) THEN
      DO i=PP_nVar+1,Nmean
        VarNamesFFT(i)='MEAN_VARIABLE'
      END DO
    END IF
    i=0
    DO j=N_FFT(2)/2+1,N_FFT(2)
      k=N_FFT(2)-i
      i=i+1
      PointData(1,i)    = (X_FFT(2,1,i,1)+1)*Re_tau
      PointData(2:4,i)  = MS_PSD(k,1:3)
      PointData(5:7,i)  = MS_t(k,1:3)
      PointData(8:7+NMean,i) = M_t(k,1:NMean)
! #if EQNSYSNR == 2
!       ! navierstokes
!       PointData(8:12,i) = M_t(k,1:5)
! #elif EQNSYSNR == 4
!       ! rans-kg
!       PointData(8:14,i) = M_t(k,1:7)
! #endif
    END DO
    CALL OpenDataFile(Filename_EK,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    ! Write dataset attributes
    CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/CHARACTER(LEN=255)::'MeanSquares'/))
    tmp255=TRIM(ProjectName)
    CALL WriteAttribute(File_ID,'ProjectName' ,1,StrScalar=(/tmp255/))
    CALL WriteAttribute(File_ID,'Time'        ,1,RealScalar=time)
    CALL WriteAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNamesFFT)
    WRITE(ZoneTitle,'(A)')'MeanSquares'
    offsetVar=0
    DO iVar=1,nVal
      CALL WriteArray(TRIM(ZoneTitle),2,(/nVal,nSamples/),(/1,nSamples/),(/offsetVar,0/),&
                      .FALSE.,RealArray=PointData(iVar,:))
      offsetVar=offsetVar+1
    END DO       ! iVar
    SDEALLOCATE(VarNamesFFT)
    SDEALLOCATE(PointData)
    CALL CloseDataFile()
    !-------------------------------------------------
    !write energy spectra in x-direction
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_x',time)
    FileName_EK=TRIM(Filename_EK)//'.h5'

    nVal=5
    nSamples=nSamples_specK
    ALLOCATE(VarNamesFFT(nVal))
    ALLOCATE(PointData(nVal,nSamples))
    VarNamesFFT(1) ='k'
    VarNamesFFT(2) ='E_uu_x'
    VarNamesFFT(3) ='E_vv_x'
    VarNamesFFT(4) ='E_ww_x'
    VarNamesFFT(5) ='E_pp_x'
    CALL OpenDataFile(Filename_EK,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/CHARACTER(LEN=255)::'MeanSquares'/))
    tmp255=TRIM(ProjectName)
    CALL WriteAttribute(File_ID,'ProjectName' ,1,StrScalar=(/tmp255/))
    CALL WriteAttribute(File_ID,'Time'        ,1,RealScalar=time)
    CALL WriteAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNamesFFT)
    DO j=N_FFT(2)/2+1,N_FFT(2)
      offsetVar=0
      DO k=1,nSamples_specK
        PointData(1,k) = 2*(k-1)
        PointData(2,k) = Ex_uu(j,k)
        PointData(3,k) = Ex_vv(j,k)
        PointData(4,k) = Ex_ww(j,k)
        PointData(5,k) = Ex_pp(j,k)
      END DO
      ! Write dataset attributes
      offsetVar=0
      WRITE(ZoneTitle,'(A,I4)')'yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau)
      DO iVar=1,nVal
        CALL WriteArray(TRIM(ZoneTitle),2,(/nVal,nSamples/),(/1,nSamples/),(/offsetVar,0/),&
                        .FALSE.,RealArray=PointData(iVar,:))
        offsetVar=offsetVar+1
      END DO
    END DO
    SDEALLOCATE(VarNamesFFT)
    SDEALLOCATE(PointData)
    !-------------------------------------------------
    !write energy spectra in z-direction
    FileName_EK=TIMESTAMP(TRIM(ProjectName)//'_EnergySpectra_z',time)
    FileName_EK=TRIM(Filename_EK)//'.h5'

    nVal=5
    nSamples=nSamples_specK
    ALLOCATE(VarNamesFFT(nVal))
    ALLOCATE(PointData(nVal,nSamples))
    VarNamesFFT(1) ='k'
    VarNamesFFT(2) ='E_uu_z'
    VarNamesFFT(3) ='E_vv_z'
    VarNamesFFT(4) ='E_ww_z'
    VarNamesFFT(5) ='E_pp_z'
    CALL OpenDataFile(Filename_EK,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)
    CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/CHARACTER(LEN=255)::'MeanSquares'/))
    tmp255=TRIM(ProjectName)
    CALL WriteAttribute(File_ID,'ProjectName' ,1,StrScalar=(/tmp255/))
    CALL WriteAttribute(File_ID,'Time'        ,1,RealScalar=time)
    CALL WriteAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNamesFFT)
    DO j=N_FFT(2)/2+1,N_FFT(2)
      offsetVar=0
      DO k=1,nSamples_specK
        PointData(1,k) = 2*(k-1)
        PointData(2,k) = Ez_uu(j,k)
        PointData(3,k) = Ez_vv(j,k)
        PointData(4,k) = Ez_ww(j,k)
        PointData(5,k) = Ez_pp(j,k)
      END DO
      ! Write dataset attributes
      offsetVar=0
      WRITE(ZoneTitle,'(A,I4)')'yPlus=',INT((1-ABS(X_FFT(2,1,j,1)))*Re_tau)
      DO iVar=1,nVal
        CALL WriteArray(TRIM(ZoneTitle),2,(/nVal,nSamples/),(/1,nSamples/),(/offsetVar,0/),&
                        .FALSE.,RealArray=PointData(iVar,:))
        offsetVar=offsetVar+1
      END DO
    END DO
    SDEALLOCATE(VarNamesFFT)
    SDEALLOCATE(PointData)

END SELECT
END SUBROUTINE FFTOutput


!===================================================================================================================================
!> Finalize FFT
!===================================================================================================================================
SUBROUTINE FinalizeFFT()
! MODULES                                                                                                                          !
USE MOD_FFT_Vars
USE MOD_DG_Vars,     ONLY: U
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(U)
SDEALLOCATE(U_FFT)
SDEALLOCATE(X_FFT)
SDEALLOCATE(Ex_uu)
SDEALLOCATE(Ex_vv)
SDEALLOCATE(Ex_pp)
SDEALLOCATE(Ez_uu)
SDEALLOCATE(Ez_vv)
SDEALLOCATE(Ez_pp)
SDEALLOCATE(inI)
SDEALLOCATE(outI)
SDEALLOCATE(inK)
SDEALLOCATE(outK)
SDEALLOCATE(VdmGaussEqui)
SDEALLOCATE(MS_t)
SDEALLOCATE(MS_PSD)
SDEALLOCATE(M_t)
SDEALLOCATE(PIJK)
SDEALLOCATE(xGP)
SDEALLOCATE(wGP)
END SUBROUTINE FinalizeFFT

END MODULE MOD_FFT
