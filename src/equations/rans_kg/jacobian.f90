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

!===================================================================================================================================
!> Contains the routines for the calculation of the analytical flux jacobians of the NS-kg equations
!===================================================================================================================================
MODULE MOD_Jacobian
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE dConsdPrimTemp
  MODULE PROCEDURE dConsdPrimTemp
END INTERFACE

INTERFACE dPrimTempdCons
  MODULE PROCEDURE dPrimTempdCons
END INTERFACE

INTERFACE EvalAdvFluxJacobianPoint
  MODULE PROCEDURE EvalAdvFluxJacobianPoint
END INTERFACE

PUBLIC::EvalAdvFluxJacobian,EvalAdvFluxJacobianPoint
#if PARABOLIC
PUBLIC::EvalDiffFluxJacobian
PUBLIC::EvalFluxGradJacobian
#endif
PUBLIC::dConsdPrimTemp,dPrimTempdCons
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Navier-Stokes-Problem:
!> The Jacobian of the advective Flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalAdvFluxJacobian(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars           ,ONLY:nDOFElem
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar,nDOFElem),INTENT(IN)              :: U               !< Conservative solution
REAL,DIMENSION(PP_nVarPrim,nDOFElem),INTENT(IN)          :: UPrim           !< Primitive solution
REAL,DIMENSION(PP_nVar,PP_nVar,nDOFElem),INTENT(OUT)     :: fJac,gJac,hJac  !< Derivative of the physical fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!===================================================================================================================================
! TODO(Shimushu): implement this
STOP '! jacobian.f90 is not ready for use.'

DO i=1,nDOFElem
  CALL EvalAdvFluxJacobianPoint(U(:,i),UPrim(:,i),fJac(:,:,i),gJac(:,:,i),hJac(:,:,i))
END DO !i
END SUBROUTINE EvalAdvFluxJacobian

!===================================================================================================================================
!> Navier-Stokes-Problem:
!> The Jacobian of the advective Flux with respect to the conservative variables U
!===================================================================================================================================
PPURE SUBROUTINE EvalAdvFluxJacobianPoint(U,UPrim,fJac,gJac,hJac)
! MODULES
USE MOD_EOS_Vars          ,ONLY:Kappa,KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar),INTENT(IN)              :: U               !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)          :: UPrim           !< Primitive solution
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: fJac,gJac,hJac  !< Derivative of the physical fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE EvalAdvFluxJacobianPoint

#if PARABOLIC
!===================================================================================================================================
!> The Jacobian of the diffusion flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalDiffFluxJacobian(nDOF_loc,U,UPrim,gradUx,gradUy,gradUz,fJac,gJac,hJac &
#if EDDYVISCOSITY
                                ,muSGS &
#endif
                                )
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars     ,ONLY:s23,s43
USE MOD_Viscosity
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN)                                   :: nDOF_loc             !< number of degrees of freedom
REAL,DIMENSION(PP_nVar        ,nDOF_loc),INTENT(IN)  :: U                    !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim    ,nDOF_loc),INTENT(IN)  :: UPrim                !< solution in primitive variables
REAL,DIMENSION(PP_nVarLifting ,nDOF_loc),INTENT(IN)  :: gradUx,gradUy,gradUz !< primitive gradients
REAL,DIMENSION(PP_nVar,PP_nVar,nDOF_loc),INTENT(OUT) :: fJac,gJac,hJac       !< Derivative of the physical diffusive fluxes
#if EDDYVISCOSITY
REAL,DIMENSION(1              ,nDOF_loc),INTENT(IN)  :: muSGS                !< eddy viscosity
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE EvalDiffFluxJacobian

!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(nDOF_loc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz &
#if EDDYVISCOSITY
                               ,muSGS &
#endif
                               )
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Equation_Vars,ONLY:s43,s23
USE MOD_EOS_Vars,     ONLY:cp,Pr
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,ONLY: PrSGS
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                       :: nDOF_loc !< number of degrees of freedom
REAL,DIMENSION(PP_nVar    ,nDOF_loc),INTENT(IN)          :: U        !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim,nDOF_loc),INTENT(IN)          :: UPrim    !< solution in primitive variables
#if EDDYVISCOSITY
REAL,DIMENSION(1          ,nDOF_loc),INTENT(IN)          :: muSGS    !< eddy viscosity
#endif
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVarPrim,nDOF_loc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz !<
                                                            !> Gradient of the diffusive Cartesian fluxes (iVar,i,j,k) w.r.t. grad
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE EvalFluxGradJacobian
#endif /*PARABOLIC*/

!===================================================================================================================================
!> The Jacobian of the transformation from primitive to conservative variables
!===================================================================================================================================
SUBROUTINE dConsdPrim(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY: KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE dConsdPrim

!===================================================================================================================================
!> The Jacobian of the transformation from primitive (including temperature) to conservative variables
!===================================================================================================================================
SUBROUTINE dConsdPrimTemp(UPrim,Jac)
! MODULES
USE MOD_PreProc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVarPrim),INTENT(OUT) :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE dConsdPrimTemp

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dPrimdCons(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY:KappaM1
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE dPrimdCons

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive (including temperature) variables
!===================================================================================================================================
SUBROUTINE dPrimTempdCons(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY:KappaM1,R
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,PP_nVar),INTENT(OUT) :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

! TODO(Shimushu): implement this

END SUBROUTINE dPrimTempdCons

END MODULE MOD_Jacobian
