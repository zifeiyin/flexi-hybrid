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

! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE

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
USE MOD_Mesh_Vars          ,ONLY: nElems

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES



END SUBROUTINE AddBodySourceTerms 


END MODULE MOD_BodySourceTerms










