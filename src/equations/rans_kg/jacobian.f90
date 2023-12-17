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
!> Contains the routines for the calculation of the analytical flux jacobians of the Navier-Stokes equations
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
USE MOD_Equation_Vars     ,ONLY:Cmu
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar),INTENT(IN)              :: U               !< Conservative solution
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)          :: UPrim           !< Primitive solution
REAL,DIMENSION(PP_nVar,PP_nVar),INTENT(OUT)     :: fJac,gJac,hJac  !< Derivative of the physical fluxes (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: KappaM2
REAL    :: uv,uu,vv,absu,v1,v2,srho,k,g
REAL    :: a1,phi
#if PP_dim==3
REAL    :: uw,vw,ww,v3
#endif
!===================================================================================================================================
KappaM2   = Kappa-2.

srho = 1./UPrim(DENS)
v1   = UPrim(VEL1)
v2   = UPrim(VEL2)
k    = UPrim(TKE)
g    = UPrim(OMG)
uv   = UPrim(VEL1)*UPrim(VEL2)
uu   = UPrim(VEL1)*UPrim(VEL1)
vv   = UPrim(VEL2)*UPrim(VEL2)
#if PP_dim==3
v3   = UPrim(VEL3)
uw   = UPrim(VEL1)*UPrim(VEL3)
vw   = UPrim(VEL2)*UPrim(VEL3)
ww   = UPrim(VEL3)*UPrim(VEL3)
absu = uu+vv+ww
phi  = kappaM1*0.5*absu
a1   = kappa * U(ENER)*sRho - phi


!                       drho         drhov1       drhov2        drhov3     drhoH    drhok   drhog
fJac(1,1:7)= (/          0.,             1.,          0.,           0.,        0.,     0.,   0. /)
fJac(2,1:7)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,  -kappaM1*v3,   kappaM1,     0.,   0. /)
fJac(3,1:7)= (/         -uv,             v2,          v1,           0.,        0.,     0.,   0. /)
fJac(4,1:7)= (/         -uw,             v3,          0.,           v1,        0.,     0.,   0. /)
fJac(5,1:7)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,  -kappaM1*uw,  kappa*v1,     0.,   0. /)
fJac(6,1:7)= (/       -v1*k,              k,          0.,           0.,        0.,     v1,   0. /)
fJac(7,1:7)= (/       -v1*g,              g,          0.,           0.,        0.,     0.,   v1 /)

gJac(1,1:7)= (/          0.,           0.,              1.,          0.,       0.,    0.,   0. /)
gJac(2,1:7)= (/         -uv,           v2,              v1,          0.,       0.,    0.,   0. /)
gJac(3,1:7)= (/      phi-vv,  -kappaM1*v1, (1.-kappaM2)*v2, -kappaM1*v3,  kappaM1,    0.,   0. /)
gJac(4,1:7)= (/         -vw,           0.,              v3,          v2,       0.,    0.,   0. /)
gJac(5,1:7)= (/ v2*(phi-a1),  -kappaM1*uv,   a1-kappaM1*vv, -kappaM1*vw, kappa*v2,    0.,   0. /)
gJac(6,1:7)= (/       -v2*k,           0.,               k,          0.,       0.,    v2,   0. /)
gJac(7,1:7)= (/       -v2*g,           0.,               g,          0.,       0.,    0.,   v2 /)


hJac(1,1:7)= (/          0.,          0.,           0.,              1.,       0.,    0.,   0. /)
hJac(2,1:7)= (/         -uw,          v3,           0.,              v1,       0.,    0.,   0. /)
hJac(3,1:7)= (/         -vw,          0.,           v3,              v2,       0.,    0.,   0. /)
hJac(4,1:7)= (/      phi-ww, -kappaM1*v1,  -kappaM1*v2, (1.-kappaM2)*v3,  kappaM1,    0.,   0. /)
hJac(5,1:7)= (/ v3*(phi-a1), -kappaM1*uw,  -kappaM1*vw,   a1-kappaM1*ww, kappa*v3,    0.,   0. /)
gJac(6,1:7)= (/       -v3*k,           0.,          0.,               k,       0.,    v3,   0. /)
gJac(7,1:7)= (/       -v3*g,           0.,          0.,               g,       0.,    0.,   v3 /)
#else
absu=uu+vv
phi  = kappaM1*0.5*absu
a1   = kappa * U(ENER)*sRho - phi


fJac(1,1:7)= (/          0.,             1.,          0.,           0.,        0.,    0.,   0. /)
fJac(2,1:7)= (/      phi-uu, (1-kappaM2)*v1, -kappaM1*v2,           0.,   kappaM1,    0.,   0. /)
fJac(3,1:7)= (/         -uv,             v2,          v1,           0.,        0.,    0.,   0. /)
fJac(4,1:7)= (/          0.,             0.,          0.,           0.,        0.,    0.,   0. /)
fJac(5,1:7)= (/ v1*(phi-a1),  a1-kappaM1*uu, -kappaM1*uv,           0.,  kappa*v1,    0.,   0. /)
fJac(6,1:7)= (/       -v1*k,              k,          0.,           0.,        0.,    v1,   0. /)
fJac(7,1:7)= (/       -v1*g,              g,          0.,           0.,        0.,    0.,   v1 /)


gJac(1,1:7)= (/          0.,          0.,              1.,          0.,        0.,    0.,   0. /)
gJac(2,1:7)= (/         -uv,          v2,              v1,          0.,        0.,    0.,   0. /)
gJac(3,1:7)= (/      phi-vv, -kappaM1*v1, (1.-kappaM2)*v2,          0.,   kappaM1,    0.,   0. /)
gJac(4,1:7)= (/          0.,          0.,              0.,          0.,        0.,    0.,   0. /)
gJac(5,1:7)= (/ v2*(phi-a1), -kappaM1*uv,   a1-kappaM1*vv,          0.,  kappa*v2,    0.,   0. /)
gJac(6,1:7)= (/       -v2*k,          0.,               k,          0.,        0.,    v2,   0. /)
gJac(7,1:7)= (/       -v2*g,          0.,               g,          0.,        0.,    0.,   v2 /)

hJac(:,:)=0.
#endif
END SUBROUTINE EvalAdvFluxJacobianPoint

#if PARABOLIC
!===================================================================================================================================
!> The Jacobian of the diffusion flux with respect to the conservative variables U
!===================================================================================================================================
SUBROUTINE EvalDiffFluxJacobian(nDOF_loc,U,UPrim,gradUx,gradUy,gradUz,fJac,gJac,hJac )
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars     ,ONLY:s23,s43
USE MOD_Equation_Vars, ONLY: PrTurb,Comega1,Comega2,Cmu,sigmaK,sigmaG
USE MOD_EOS_Vars,      ONLY: cp
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
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: dir,i
REAL                :: muS
#if PP_dim==3
REAL                :: tau(3,3)
#else
REAL                :: tau(2,2)
#endif
REAL                :: invRho
REAL                :: muTurb, muEff, lambda, kDiffEff, gDiffEff
REAL                :: dmuTurb_drhok, dmuTurb_drhog
!===================================================================================================================================
fJac = 0.
gJac = 0.
hJac = 0.

DO i=1,nDOF_loc
  muS = VISCOSITY_TEMPERATURE(UPrim(TEMP,i))

  ! Add turbulent viscosity and diffusivity
  invRho = 1.0 / U(DENS,i)
  muTurb = Cmu * U(RHOK,i) * U(RHOG,i) * U(RHOG,i) * invRho * invRho  
  muEff  = MAX(muS, muS + muTurb)
  lambda = MAX(lambda,lambda+muTurb*cp/PrTurb)
  !diffusivity of turbulence variables
  kDiffEff = MAX(muS, muS + muTurb * sigmaK)
  gDiffEff = MAX(muS, muS + muTurb * sigmaG)

  ! compute dmuTurb_d
  dmuTurb_drhok = Cmu * U(RHOG,i) * U(RHOG,i) * invRho * invRho 
  dmuTurb_drhog = 2.0 * U(RHOK,i) * U(RHOG,i) * invRho * invRho
#if PP_dim==3
    tau(1,1) = ( s43 * gradUx(LIFT_VEL1,i) - s23 * gradUy(LIFT_VEL2,i) - s23 * gradUz(LIFT_VEL3,i)) ! 4/3*u_x-2/3*v_y -2/3*w*z
    tau(2,2) = (-s23 * gradUx(LIFT_VEL1,i) + s43 * gradUy(LIFT_VEL2,i) - s23 * gradUz(LIFT_VEL3,i)) !-2/3*u_x+4/3*v_y -2/3*w*z
    tau(3,3) = (-s23 * gradUx(LIFT_VEL1,i) - s23 * gradUy(LIFT_VEL2,i) + s43 * gradUz(LIFT_VEL3,i)) !-2/3*u_x-2/3*v_y +4/3*w*z
    tau(1,2) = (gradUy(LIFT_VEL1,i) + gradUx(LIFT_VEL2,i))               !(u_y+v_x)
    tau(2,1) = tau(1,2)
    tau(1,3) = (gradUz(LIFT_VEL1,i) + gradUx(LIFT_VEL3,i))               !(u_z+w_x)
    tau(3,1) = tau(1,3)
    tau(2,3) = (gradUz(LIFT_VEL2,i) + gradUy(LIFT_VEL3,i))               !(y_z+w_y)
    tau(3,2) = tau(2,3)
#else
    tau(1,1) = ( s43 * gradUx(LIFT_VEL1,i) - s23 * gradUy(LIFT_VEL2,i))  ! 4/3*u_x-2/3*v_y -2/3*w*z
    tau(2,2) = (-s23 * gradUx(LIFT_VEL1,i) + s43 * gradUy(LIFT_VEL2,i))  !-2/3*u_x+4/3*v_y -2/3*w*z
    tau(1,2) = (gradUy(LIFT_VEL1,i) + gradUx(LIFT_VEL2,i))               !(u_y+v_x)
    tau(2,1) = tau(1,2)
#endif

  ! For RANS-k-g, the momentum equations depend on the k-g working variable via the effective viscosity
  ! e.g. F_{rho*u} = -tau_xx * (mu+muTurb)
  fJac(2,6:7,i) = -1. * (/ dmuTurb_drhok, dmuTurb_drhog /) * tau(1,1)
  fJac(3,6:7,i) = -1. * (/ dmuTurb_drhok, dmuTurb_drhog /) * tau(1,2)
#if PP_dim==3
  fJac(4,6:7,i) = -1. * (/ dmuTurb_drhok, dmuTurb_drhog /) * tau(1,3)
#endif

  ! dF^d(5)/dU(1) = (u*tau_(1,1)*muEff + v*tau_(1,2)*muEff + w*tau_(1,3)*muEff)/rho
  ! dF^d(5)/dU(2) = - tau_(1,1)*muEff/rho
  ! dF^d(5)/dU(3) = - tau_(1,2)*muEff/rho
  ! dF^d(5)/dU(4) = - tau_(1,3)*muEff/rho
  ! dF^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    fJac(5,1,i) = fJac(5,1,i) + tau(1,dir)*UPrim(1+dir,i) * muEff
  END DO
  fJac(5,1,i) = fJac(5,1,i)/ U(1,i) 
  fJac(5,2,i) = -tau(1,1)  / U(1,i) * muEff
  fJac(5,3,i) = -tau(1,2)  / U(1,i) * muEff
#if PP_dim==3
  fJac(5,4,i) = -tau(1,3)  / U(1,i) * muEff
  ! The energy equation depends on the k-g working variable through both the effective viscosity and the effective thermal
  ! conductivity
  fJac(5,6:7,i) = -1. * (/ dmuTurb_drhok, dmuTurb_drhog /) &
                * ( (tau(1,1)*UPrim(2,i) + tau(1,2)*UPrim(3,i) + tau(1,3)*UPrim(4,i)) + cp/PrTurb * gradUx(LIFT_TEMP,i) )
#elif 
  fJac(5,6:7,i) = -1.* (/dmuTurb_drhok, dmuTurb_drhog /) & 
                * ( (tau(1,1)*UPrim(2,i) + tau(1,2)*UPrim(3,i)) + cp/PrTurb * gradUx(LIFT_TEMP,i) )
#endif
  ! k-g
  fJac(6,6,i) = -sigmaK * dmuTurb_drhok * gradUx(LIFT_TKE,i)
  fJac(6,7,i) = -sigmaK * dmuTurb_drhog * gradUx(LIFT_TKE,i)
  fJac(7,6,i) = -sigmaG * dmuTurb_drhok * gradUx(LIFT_OMG,i)
  fJac(7,7,i) = -sigmaG * dmuTurb_drhog * gradUx(LIFT_OMG,i)


  gJac(2,6:7,i) = -tau(2,1) * (/ dmuTurb_drhok, dmuTurb_drhog /) 
  gJac(3,6:7,i) = -tau(2,2) * (/ dmuTurb_drhok, dmuTurb_drhog /)
#if PP_dim==3
  gJac(4,6:7,i) = -tau(2,3) * (/ dmuTurb_drhok, dmuTurb_drhog /)
#endif
  ! dG^d(1:4)/dU(:) = 0, only the energy equation is directly depending on the conservative solution!
  ! dG^d(5)/dU(1) = (u*tau_(2,1) + v*tau_(2,2) + w*tau_(2,3))/rho
  ! dG^d(5)/dU(2) = - tau_(2,1)/rho
  ! dG^d(5)/dU(3) = - tau_(2,2)/rho
  ! dG^d(5)/dU(4) = - tau_(2,3)/rho
  ! dG^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    gJac(5,1,i) = gJac(5,1,i) + tau(2,dir)*UPrim(1+dir,i) * muEff
  END DO
  gJac(5,1,i) = gJac(5,1,i)/ U(1,i)
  gJac(5,2,i) = -tau(2,1)  / U(1,i) * muEff
  gJac(5,3,i) = -tau(2,2)  / U(1,i) * muEff
#if PP_dim==3
  gJac(5,4,i) = -tau(2,3)  / U(1,i) * muEff
#endif
  gJac(6,6,i) = -sigmaK * dmuTurb_drhok * gradUy(LIFT_TKE,i)
  gJac(6,7,i) = -sigmaK * dmuTurb_drhog * gradUy(LIFT_TKE,i)
  gJac(7,6,i) = -sigmaG * dmuTurb_drhok * gradUy(LIFT_OMG,i)
  gJac(7,7,i) = -sigmaG * dmuTurb_drhog * gradUy(LIFT_OMG,i)

#if PP_dim==3
  hJac(2,6:7,i) = -tau(3,1) * (/ dmuTurb_drhok, dmuTurb_drhog /)
  hJac(3,6:7,i) = -tau(3,2) * (/ dmuTurb_drhok, dmuTurb_drhog /)
  hJac(4,6:7,i) = -tau(3,3) * (/ dmuTurb_drhok, dmuTurb_drhog /)
  ! dH^d(1:4)/dU(:) = 0, only the energy equation is directly depending on the conservative solution!
  ! dH^d(5)/dU(1) = (u*tau_(3,1) + v*tau_(3,2) + w*tau_(3,3))/rho
  ! dH^d(5)/dU(2) = - tau_(3,1)/rho
  ! dH^d(5)/dU(3) = - tau_(3,2)/rho
  ! dH^d(5)/dU(4) = - tau_(3,3)/rho
  ! dH^d(5)/dU(5) = 0.
  DO dir=1,PP_dim
    hJac(5,1,i) = hJac(5,1,i) + tau(3,dir)*UPrim(1+dir,i) * muEff
  END DO
  hJac(5,1,i) = hJac(5,1,i)/ U(1,i)
  hJac(5,2,i) = -tau(3,1)  / U(1,i) * muEff
  hJac(5,3,i) = -tau(3,2)  / U(1,i) * muEff
  hJac(5,4,i) = -tau(3,3)  / U(1,i) * muEff
#endif
END DO

END SUBROUTINE EvalDiffFluxJacobian

!===================================================================================================================================
!> Computes the volume derivative of the analytical diffusive flux with respect to the gradient of U: d(F^v)/dQ, Q=grad U
!===================================================================================================================================
SUBROUTINE EvalFluxGradJacobian(nDOF_loc,U,UPrim,fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz )
! MODULES
USE MOD_PreProc
USE MOD_Viscosity
USE MOD_Equation_Vars,ONLY:s43,s23
Use MOD_Equation_Vars,ONLY: PrTurb, Comega1, Comega2, Cmu, sigmaK, sigmaG
USE MOD_EOS_Vars,     ONLY:cp,Pr

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                       :: nDOF_loc !< number of degrees of freedom
REAL,DIMENSION(PP_nVar    ,nDOF_loc),INTENT(IN)          :: U        !< solution in conservative variables
REAL,DIMENSION(PP_nVarPrim,nDOF_loc),INTENT(IN)          :: UPrim    !< solution in primitive variables
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,PP_nVarPrim,nDOF_loc),INTENT(OUT) :: fJacQx,fJacQy,fJacQz,gJacQx,gJacQy,gJacQz,hJacQx,hJacQy,hJacQz !<
                                                            !> Gradient of the diffusive Cartesian fluxes (iVar,i,j,k) w.r.t. grad
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
REAL                :: muS,lambda
REAL                :: muTurb
REAL                :: kDiffEff, gDiffEff
REAL                :: dmuTurb_dk, dmuTurb_dg
REAL                :: v1, v2, v3
!===================================================================================================================================
DO i=1,nDOF_loc
  muS    = VISCOSITY_TEMPERATURE(UPrim(TEMP,i))
  lambda = THERMAL_CONDUCTIVITY_H(muS)
  !Add turbulent sub grid scale viscosity to mu

  v1 = UPrim(VEL1,i) 
  v2 = UPrim(VEL2,i) 
  v3 = UPrim(VEL3,i) 

  ! add turbulence part
  muTurb = Cmu * UPrim(DENS,i) * UPrim(TKE,i) * UPrim(OMG,i) * UPrim(OMG,i)
  muS    = MAX( muS, muS + muTurb )
  lambda = MAX( lambda, lambda + muTurb * cp / PrTurb )

  kDiffEff = MAX(muS, muS + muTurb * sigmaK)
  gDiffEff = MAX(muS, muS + muTurb * sigmaG)

  dmuTurb_dk = Cmu * UPrim(DENS,i) * UPrim(OMG,i) * UPrim(OMG,i)
  dmuTurb_dg = 2.0 * Cmu * UPrim(DENS,i) * UPrim(TKE,i) * UPrim(OMG,i)

#if PP_dim==3
  !              grad(rho)         grad(v1)        grad(v2)      grad(v3) grad(p)     grad(T)    grad(k)   grad(g)  grad(nut)
  ! derivatives of diffusive flux in x-direction
  fJacQx(1,1:9,i) = 0.
  fJacQx(2,1:9,i) = (/ 0.,        -muS*s43,             0.,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQx(3,1:9,i) = (/ 0.,              0.,           -muS,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQx(4,1:9,i) = (/ 0.,              0.,             0.,         -muS,     0.,         0.,        0.,         0.,   0.  /)
  fJacQx(5,1:9,i) = (/ 0.,     -muS*s43*v1,        -muS*v2,      -muS*v3,     0.,    -lambda,        0.,         0.,   0.  /)
  fJacQx(6,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0., -kDiffEff,         0.,   0.  /)
  fJacQx(7,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0.,        0.,  -gDiffEff,   0.  /)

  fJacQy(1,1:9,i) = 0.
  fJacQy(2,1:9,i) = (/ 0.,              0.,        muS*s23,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQy(3,1:9,i) = (/ 0.,            -muS,             0.,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQy(4,1:9,i) = 0.
  fJacQy(5,1:9,i) = (/ 0.,         -muS*v2,     muS*s23*v1,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQy(6,1:9,i) = 0.
  fJacQy(7,1:9,i) = 0.

  fJacQz(1,1:9,i) = 0.
  fJacQz(2,1:9,i) = (/ 0.,              0.,             0.,      muS*s23,     0.,         0.,        0.,         0.,   0.  /)
  fJacQz(3,1:9,i) = 0.
  fJacQz(4,1:9,i) = (/ 0.,            -muS,             0.,           0.,     0.,         0.,        0.,         0.,   0.  /)
  fJacQz(5,1:9,i) = (/ 0.,         -muS*v3,             0.,   muS*s23*v1,     0.,         0.,        0.,         0.,   0.  /)
  fJacQz(6,1:9,i) = 0.
  fJacQz(7,1:9,i) = 0.


  ! derivatives of diffusive flux in y-direction
  gJacQx(1,1:9,i) = 0.
  gJacQx(2,1:9,i) = (/ 0.,              0.,           -muS,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQx(3,1:9,i) = (/ 0.,         muS*s23,             0.,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQx(4,1:9,i) = 0.
  gJacQx(5,1:9,i) = (/ 0.,      muS*s23*v2,        -muS*v1,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQx(6,1:9,i) = 0.
  gJacQx(7,1:9,i) = 0.

  gJacQy(1,1:9,i) = 0.
  gJacQy(2,1:9,i) = (/ 0.,            -muS,             0.,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQy(3,1:9,i) = (/ 0.,              0.,       -muS*s43,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQy(4,1:9,i) = (/ 0.,              0.,             0.,         -muS,     0.,         0.,        0.,         0.,   0.  /)
  gJacQy(5,1:9,i) = (/ 0.,         -muS*v1,    -muS*s43*v2,      -muS*v3,     0.,    -lambda,        0.,         0.,   0.  /)
  gJacQy(6,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0., -kDiffEff,         0.,   0.  /)
  gJacQy(7,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0.,        0.,  -gDiffEff,   0.  /)

  gJacQz(1,1:9,i) = 0.
  gJacQz(2,1:9,i) = 0.
  gJacQz(3,1:9,i) = (/ 0.,              0.,             0.,      muS*s23,     0.,         0.,        0.,         0.,   0.  /)
  gJacQz(4,1:9,i) = (/ 0.,              0.,           -muS,           0.,     0.,         0.,        0.,         0.,   0.  /)
  gJacQz(5,1:9,i) = (/ 0.,              0.,        -muS*v3,   muS*s23*v2,     0.,         0.,        0.,         0.,   0.  /)
  gJacQz(6,1:9,i) = 0.
  gJacQz(7,1:9,i) = 0.

  ! derivatives of diffusive flux in z-direction
  !              grad(rho)         grad(v1)        grad(v2)      grad(v3) grad(p)     grad(T)    grad(k)   grad(g)  grad(nut)
  hJacQx(1,1:9,i) = 0.
  hJacQx(2,1:9,i) = (/ 0.,              0.,             0.,         -muS,     0.,         0.,         0.,        0.,   0.  /)
  hJacQx(3,1:9,i) = 0.
  hJacQx(4,1:9,i) = (/ 0.,         muS*s23,             0.,           0.,     0.,         0.,         0.,        0.,   0.  /)
  hJacQx(5,1:9,i) = (/ 0.,      muS*s23*v3,             0.,      -muS*v1,     0.,         0.,         0.,        0.,   0.  /)
  hJacQx(6,1:9,i) = 0.
  hJacQx(7,1:9,i) = 0.

  hJacQy(1,1:9,i) = 0.
  hJacQy(2,1:9,i) = 0.
  hJacQy(3,1:9,i) = (/ 0.,              0.,             0.,         -muS,     0.,         0.,         0.,         0.,  0.  /)
  hJacQy(4,1:9,i) = (/ 0.,              0.,        muS*s23,           0.,     0.,         0.,         0.,         0.,  0.  /)
  hJacQy(5,1:9,i) = (/ 0.,              0.,     muS*s23*v3,      -muS*v2,     0.,         0.,         0.,         0.,  0.  /)
  hJacQy(6,1:9,i) = 0.
  hJacQy(7,1:9,i) = 0.

  hJacQz(1,1:9,i) = 0.
  hJacQz(2,1:9,i) = (/ 0.,            -muS,             0.,           0.,     0.,         0.,         0.,         0.,  0.  /)
  hJacQz(3,1:9,i) = (/ 0.,              0.,           -muS,           0.,     0.,         0.,         0.,         0.,  0.  /)
  hJacQz(4,1:9,i) = (/ 0.,              0.,             0.,     -muS*s43,     0.,         0.,         0.,         0.,  0.  /)
  hJacQz(5,1:9,i) = (/ 0.,         -muS*v1,        -muS*v2,  -muS*s43*v3,     0.,    -lambda,         0.,         0.,  0.  /)
  hJacQz(6,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0.,  -kDiffEff,         0.,  0.  /)
  hJacQz(7,1:9,i) = (/ 0.,              0.,             0.,           0.,     0.,         0.,         0.,  -gDiffEff,  0.  /)
#else
  ! derivatives of diffusive flux in x-direction
  fJacQx(1,1:9,i) = 0.
  fJacQx(2,1:9,i) = (/ 0.,        -muS*s43,             0.,           0.,     0.,         0.,         0.,         0.,  0.  /)
  fJacQx(3,1:9,i) = (/ 0.,              0.,           -muS,           0.,     0.,         0.,         0.,         0.,  0.  /)
  fJacQx(4,1:9,i) = 0.
  fJacQx(5,1:9,i) = (/ 0.,     -muS*s43*v1,        -muS*v2,           0.,     0.,    -lambda,         0.,         0.,  0.  /)
  fJacQx(6,1:9,i) = (/ 0.,              0.,                           0.,     0.,         0.,  -kDiffEff,         0.,  0.  /)
  fJacQx(7,1:9,i) = (/ 0.,              0.,                           0.,     0.,         0.,         0.,  -gDiffEff,  0.  /)

  fJacQy(1,1:9,i) = 0.
  fJacQy(2,1:9,i) = (/ 0.,              0.,        muS*s23,           0.,     0.,         0.,         0.,         0.,  0.  /)
  fJacQy(3,1:9,i) = (/ 0.,            -muS,             0.,           0.,     0.,         0.,         0.,         0.,  0.  /)
  fJacQy(4,1:9,i) = 0.
  fJacQy(5,1:9,i) = (/ 0.,         -muS*v2,     muS*s23*v1,           0.,     0.,         0.,         0.,         0.,  0.  /)
  fJacQy(6,1:9,i) = 0.
  fJacQy(7,1:9,i) = 0.

  fJacQz(:,:,i) = 0.

  ! derivatives of diffusive flux in y-direction
  gJacQx(1,1:9,i) = 0.
  gJacQx(2,1:9,i) = (/ 0.,              0.,           -muS,        0.,        0.,         0.,         0.,         0.,  0.  /)
  gJacQx(3,1:9,i) = (/ 0.,         muS*s23,             0.,        0.,        0.,         0.,         0.,         0.,  0.  /)
  gJacQx(4,1:9,i) = 0.
  gJacQx(5,1:9,i) = (/ 0.,      muS*s23*v2,        -muS*v1,        0.,        0.,         0.,         0.,         0.,  0.  /)
  gJacQx(6,1:9,i) = 0.
  gJacQx(7,1:9,i) = 0.

  gJacQy(1,1:9,i) = 0.
  gJacQy(2,1:9,i) = (/ 0.,            -muS,             0.,        0.,        0.,         0.,         0.,         0.,  0.   /)
  gJacQy(3,1:9,i) = (/ 0.,              0.,       -muS*s43,        0.,        0.,         0.,         0.,         0.,  0.   /)
  gJacQy(4,1:9,i) = 0.
  gJacQy(5,1:9,i) = (/ 0.,         -muS*v1,    -muS*s43*v2,        0.,        0.,    -lambda,         0.,         0.,  0.   /)
  gJacQy(6,1:9,i) = (/ 0.,              0.,             0.,        0.,        0.,         0.,  -kDiffEff,         0.,  0.   /)
  gJacQy(7,1:9,i) = (/ 0.,              0.,             0.,        0.,        0.,         0.,         0.,  -gDiffEff,  0.   /)

  gJacQz(:,:,i) = 0.

  ! derivatives of diffusive flux in z-direction
  hJacQx(:,:,i) = 0.
  hJacQy(:,:,i) = 0.
  hJacQz(:,:,i) = 0.
#endif
END DO
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
REAL,DIMENSION(PP_nVar,PP_nVarPrim),INTENT(OUT)     :: Jac      !< cons to prim Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                            :: UE(PP_2Var),dpdrho,dedrho
REAL                                            :: rho, v1, v2, v3, tk, tg
!===================================================================================================================================
UE(EXT_PRIM) = UPrim
UE(EXT_DENS) = UPrim(1)
UE(EXT_SRHO) = 1./UE(EXT_DENS)

rho = UPrim(DENS)
v1  = UPrim(VEL1)
v2  = UPrim(VEL2) 
v3  = UPrim(VEL3)  
tk  = UPrim(TKE)
tg  = UPrim(OMG)     

#if PP_dim == 3
dedrho = 0.5 * ( v1*v1 + v2*v2 + v3*v3 )
#else
dedrho = 0.5 * ( v1*v1 + v2*v2 ) 
#endif

dpdrho = KappaM1*0.5*( v1*v1 + v2*v2 + v3*v3 )
#if PP_dim == 3
dedrho = ( v1*v1 + v2*v2 + v3*v3 ) - dpdrho / KappaM1
#else
dedrho = ( v1*v1 + v2*v2 ) - dpdrho / KappaM1
#endif

! column 6:9 are temperature, tke, g, nut
!                  rho,            v1                 v2                 v3         p      T       k .     g     nut  
Jac(1,1:9) = (/     1.,            0.,                0.,                0.,         0.,    0.,     0.,     0.,    0. /)
Jac(2,1:9) = (/     v1,           rho,                0.,                0.,         0.,    0.,     0.,     0.,    0. /)
Jac(3,1:9) = (/     v2,            0.,               rho,                0.,         0.,    0.,     0.,     0.,    0. /)
#if PP_dim == 3
Jac(4,1:9) = (/     v3,            0.,                0.,               rho,         0.,    0.,     0.,     0.,    0. /)
Jac(5,1:9) = (/ dedrho,        rho*v1,            rho*v2,            rho*v3, 1./KappaM1,    0.,    rho,     0.,    0. /)
#else
Jac(4,1:9) = 0.
Jac(5,1:9) = (/ dedrho,        rho*v1,            rho*v2,                0., 1./KappaM1,    0.,    rho,     0.,    0. /)
#endif
! dependency of the rho*k and rho*g (=> conservative variable)
Jac(6,1:9) = (/     tk,            0.,                0.,                0.,          0., 0.,      rho,     0.,    0. /)
Jac(7,1:9) = (/     tg,            0.,                0.,                0.,          0., 0.,      rho,     0.,    0. /)

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
! fill jacobian without temperature
CALL dConsdPrim(UPrim,Jac(1:7,1:9))
Jac(1:7,6) = 0.
END SUBROUTINE dConsdPrimTemp

!===================================================================================================================================
!> The Jacobian of the transformation from conservative to primitive variables
!===================================================================================================================================
SUBROUTINE dPrimdCons(UPrim,Jac)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars                 ,ONLY:KappaM1
USE MOD_Equation_Vars            ,ONLY:Cmu
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim)        ,INTENT(IN)  :: UPrim    !< primitive state vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,PP_nVar),INTENT(OUT)     :: Jac      !< prim to cons Jacobian
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                            :: invRho 
REAL                                            :: rho, v1, v2, v3, tk, tg
REAL                                            :: dHdRho
!===================================================================================================================================

invRho = 1./UPrim(DENS)
rho = UPrim(DENS)
v1  = UPrim(VEL1)
v2  = UPrim(VEL2) 
v3  = UPrim(VEL3)  
tk  = UPrim(TKE)
tg  = UPrim(OMG) 

#if PP_dim == 3
dHdRho = KappaM1*0.5*(v1*v1+v2*v2+v3*v3) 
#else
dHdRho = KappaM1*0.5*(v1*v1+v2*v2)
#endif

!                      rho          rho*v1        rho*v2        rho*v3    rho*h,    rho*k,      rho*g
Jac(1,1:7)= (/          1.,             0.,           0.,           0.,      0.,        0.,        0.   /)
Jac(2,1:7)= (/  -v1*invRho,         invRho,           0.,           0.,      0.,        0.,        0.   /)
Jac(3,1:7)= (/  -v2*invRho,             0.,       invRho,           0.,      0.,        0.,        0.   /)
#if PP_dim == 3
Jac(4,1:7)= (/  -v3*invRho,             0.,           0.,       invRho,      0.,        0.,        0.   /)
Jac(5,1:7)= (/      dHdRho,    -v1*KappaM1,  -v2*KappaM1,  -v3*KappaM1, KappaM1,  -kappaM1,        0.   /)
#else
Jac(4,1:7)= 0.
Jac(5,1:7)= (/      dHdRho,    -v1*KappaM1,  -v2*KappaM1,           0., KappaM1,  -KappaM1,        0.   /)
#endif

! turbulence kinetic energy
!                      rho          rho*v1        rho*v2        rho*v3    rho*h,    rho*k,      rho*g
Jac(7,1:7) = (/ -tk * invRho,           0.,           0.,           0.,      0.,   invRho,         0.  /) 

! turbulence g
Jac(8,1:7) = (/ -tg * invRho,           0.,           0.,           0.,      0.,       0.,     invRho  /)

! turbulence nut
Jac(9,1)   = -3. * Cmu * tk * tg * tg * invRho
Jac(9,2:5) = 0.
Jac(9,6)   = Cmu * tg * tg * invRho
Jac(9,7)   = 2. * Cmu * tk * tg * invRho


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
REAL                                            :: sRhoR,dpdU(6)
!===================================================================================================================================
! fill jacobian without temperature
CALL dPrimdCons(UPrim,Jac(1:9,1:7))

! fill jacobian of transformation to temperature
#if PP_dim==3
dpdU(1)   =  KappaM1*(0.5*(UPrim(VEL1)**2+UPrim(VEL2)**2+UPrim(VEL3)**2) + UPrim(TKE) )
#else
dpdU(1)   =  KappaM1*(0.5*(UPrim(VEL1)**2+UPrim(VEL2)**2) + UPrim(TKE) )
#endif
dpdU(2:4) = -KappaM1*UPrim(VELV)
dpdU(5)   =  KappaM1
dpdU(6)   = -KappaM1
sRhoR     =  1./(R*UPrim(DENS))

Jac(6,1)   = dpdU(1  )*sRhoR-UPrim(PRES)*sRhoR/UPrim(DENS)
Jac(6,2:4) = dpdU(2:4)*sRhoR
Jac(6,5)   = dpdU(5  )*sRhoR
Jac(6,6)   = dpdU(6  )*sRhoR
Jac(6,7)   = 0.

END SUBROUTINE dPrimTempdCons

END MODULE MOD_Jacobian
