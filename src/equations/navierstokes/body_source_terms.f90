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

!==================================================================================================================================
!> Subroutines providing the source terms of all the transport equations
!==================================================================================================================================
Module MOD_BodySourceTerms
!MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersBodySourceTerms
  MODULE PROCEDURE DefineParametersBodySourceTerms
END INTERFACE

INTERFACE InitBodySourceTerms
  MODULE PROCEDURE InitBodySourceTerms
END INTERFACE

INTERFACE AddBodySourceTerms
  MODULE PROCEDURE AddBodySourceTerms
END INTERFACE

PUBLIC::DefineParametersBodySourceTerms
PUBLIC::InitBodySourceTerms
PUBLIC::AddBodySourceTerms
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters of exact functions
!==================================================================================================================================
SUBROUTINE DefineParametersBodySourceTerms()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("BodyForce")
CALL prms%CreateIntFromStringOption('BodyForceType', "Body force to be added to momentum equation.", "none")
CALL addStrListEntry('BodyForceType','none',-1)
CALL addStrListEntry('BodyForceType','fixedPressureGradient',1)

CALL prms%CreateRealArrayOption( 'bodyForceVector',  "pressure gradient or force vector (f1, f2, f3) required for CASE(1)")

END SUBROUTINE DefineParametersBodySourceTerms

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitBodySourceTerms()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_ExactFunc_Vars
USE MOD_Mesh_Vars          ,ONLY: nElems
USE MOD_Equation_Vars      ,ONLY: IniBodyForce

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE

!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT BODY FORCE...'

IniBodyForce = GETINTFROMSTR('BodyForceType')
SELECT CASE (IniBodyForce)
CASE(1) ! fixed pressure gradient
  BodyForceVector = GETREALARRAY('bodyForceVector',3)
CASE DEFAULT 
END SELECT

#if PP_dim==2 
SELECT CASE (IniBodyForce)
CASE(1) ! fixed pressure gradient
  IF(BodyForceVector(3).NE.0.)THEN
    CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set BodyForceVector(3) = 0!')
  ENDIF
END SELECT
#endif 

SWRITE(UNIT_stdOut,'(A)')' INIT BODY FORCE DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE InitBodySourceTerms

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!> dt is only needed to compute the time dependent boundary values for the RK scheme
!> for each function resu and the first and second time derivative resu_t and resu_tt have to be defined (is trivial for constants)
!==================================================================================================================================
SUBROUTINE AddBodySourceTerms(Ut,t)
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY: U
USE MOD_Mesh_Vars          ,ONLY: nElems,Elem_xGP,sJ
USE MOD_Equation_Vars      ,ONLY: IniBodyForce
USE MOD_Exactfunc_Vars
#if FV_ENABLED
USE MOD_ChangeBasisByDim   ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars            ,ONLY: FV_Vdm,FV_Elems
#endif
  
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

!----------------------------------------------------------------------------------------------------------------------------------
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: Ut_src2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
INTEGER             :: i,j,k,iElem
INTEGER             :: FV_Elem
!==================================================================================================================================

SELECT CASE (IniBodyForce)
CASE(1) ! fixed pressure gradient
DO iElem=1,nElems

#if FV_ENABLED
  FV_Elem = FV_Elems(iElem)
#else
  FV_Elem = 0
#endif

#if FV_ENABLED
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    Ut_src(MOM1:MOM3,i,j,k) = BodyForceVector
  END DO ; END DO; END DO ! i,j,k
  IF (FV_Elems(iElem).GT.0) THEN ! FV elem
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(MOM1:MOM3,i,j,k,iElem) = Ut(MOM1:MOM3,i,j,k,iElem)+Ut_src2(MOM1:MOM3,i,j,k)/sJ(i,j,k,iElem,FV_Elem)
      Ut(ENER,i,j,k,iElem) = Ut(ENER,i,j,k,iElem) + &
          dot_product( &
              Ut_src2(MOM1:MOM3, i, j, k) / Ut_src2(DENS, i, j, k), &
              BodyForceVector &
          ) / sJ(i, j, k, iElem, FV_Elem)
    END DO; END DO; END DO ! i,j,k
  ELSE
#endif 
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(MOM1:MOM3,i,j,k,iElem) = Ut(MOM1:MOM3,i,j,k,iElem)+BodyForceVector/sJ(i,j,k,iElem,FV_Elem)
      Ut(ENER,i,j,k,iElem) = Ut(ENER,i,j,k,iElem) + &
          dot_product( &
              U(MOM1:MOM3, i, j, k, iElem) / U(DENS, i, j, k, iElem), &
              BodyForceVector &
          ) / sJ(i, j, k, iElem, FV_Elem)
    END DO; END DO; END DO ! i,j,k
#if FV_ENABLED
  END IF
#endif
END DO ! element loop
END SELECT

END SUBROUTINE AddBodySourceTerms 


END MODULE MOD_BodySourceTerms