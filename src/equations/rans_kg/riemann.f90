!=================================================================================================================================
! Copyright (c) 2010-2024 Prof. Claus-Dieter Munz
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PPURE SUBROUTINE RiemannInt(F_L,F_R,U_LL,U_RR,F)
    REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(RiemannInt),POINTER :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver

INTEGER,PARAMETER      :: PRM_RIEMANN_SAME          = -1
INTEGER,PARAMETER      :: PRM_RIEMANN_LF            = 1
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLC          = 2
INTEGER,PARAMETER      :: PRM_RIEMANN_ROE           = 3
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEL2         = 32
INTEGER,PARAMETER      :: PRM_RIEMANN_ROEENTROPYFIX = 33
INTEGER,PARAMETER      :: PRM_RIEMANN_HLL           = 4
INTEGER,PARAMETER      :: PRM_RIEMANN_HLLE          = 5
INTEGER,PARAMETER      :: PRM_RIEMANN_SLAU          = 6
INTEGER,PARAMETER      :: PRM_RIEMANN_SLAU2         = 62
#ifdef SPLIT_DG
INTEGER,PARAMETER      :: PRM_RIEMANN_CH            = 7
INTEGER,PARAMETER      :: PRM_RIEMANN_Average       = 0
#endif

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann_Point
  MODULE PROCEDURE Riemann_Side
END INTERFACE

#if PARABOLIC
INTERFACE ViscousFlux
  MODULE PROCEDURE ViscousFlux_Point
  MODULE PROCEDURE ViscousFlux_Side
END INTERFACE
#endif

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE

PUBLIC::DefineParametersRiemann
PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
#if PARABOLIC
PUBLIC::ViscousFlux
#endif
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
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
CALL prms%SetSection("Riemann")
CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: LF, HLLC, Roe, RoeEntropyFix, HLL, HLLE, HLLEM", &
                                                 "RoeEntropyFix")
CALL addStrListEntry('Riemann','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('Riemann','hllc',         PRM_RIEMANN_HLLC)
CALL addStrListEntry('Riemann','RoeL2',        PRM_RIEMANN_ROEL2)
CALL addStrListEntry('Riemann','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('Riemann','hll',          PRM_RIEMANN_HLL)
CALL addStrListEntry('Riemann','hlle',         PRM_RIEMANN_HLLE)
CALL addStrListEntry('Riemann','slau',         PRM_RIEMANN_SLAU)
CALL addStrListEntry('Riemann','slau2',        PRM_RIEMANN_SLAU2)
#ifdef SPLIT_DG
CALL addStrListEntry('Riemann','ch',           PRM_RIEMANN_CH)
CALL addStrListEntry('Riemann','avg',          PRM_RIEMANN_Average)
#endif
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: Same, LF, Roe, RoeEntropyFix, "//&
                                                 "HLL, HLLE, HLLEM",&
                                                 "Same")
CALL addStrListEntry('RiemannBC','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','hllc',         PRM_RIEMANN_HLLC)
CALL addStrListEntry('RiemannBC','RoeL2',        PRM_RIEMANN_ROEL2)
CALL addStrListEntry('RiemannBC','roeentropyfix',PRM_RIEMANN_ROEENTROPYFIX)
CALL addStrListEntry('RiemannBC','hll',          PRM_RIEMANN_HLL)
CALL addStrListEntry('RiemannBC','hlle',         PRM_RIEMANN_HLLE)
CALL addStrListEntry('RiemannBC','slau',         PRM_RIEMANN_SLAU)
CALL addStrListEntry('RiemannBC','slau2',        PRM_RIEMANN_SLAU2)
#ifdef SPLIT_DG
CALL addStrListEntry('RiemannBC','ch',           PRM_RIEMANN_CH)
CALL addStrListEntry('RiemannBC','avg',          PRM_RIEMANN_Average)
#endif
CALL addStrListEntry('RiemannBC','same',         PRM_RIEMANN_SAME)
END SUBROUTINE DefineParametersRiemann


!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
!==================================================================================================================================
#ifndef SPLIT_DG
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLC)
  Riemann_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_HLL)
  Riemann_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLE)
  Riemann_pointer => Riemann_HLLE
CASE(PRM_RIEMANN_SLAU)
  Riemann_pointer => Riemann_SLAU
CASE(PRM_RIEMANN_SLAU2)
  Riemann_pointer => Riemann_SLAU2
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE(PRM_RIEMANN_HLLC)
  RiemannBC_pointer => Riemann_HLLC
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_HLL)
  RiemannBC_pointer => Riemann_HLL
CASE(PRM_RIEMANN_HLLE)
  RiemannBC_pointer => Riemann_HLLE
CASE(PRM_RIEMANN_SLAU)
  Riemann_pointer => Riemann_SLAU
CASE(PRM_RIEMANN_SLAU2)
  Riemann_pointer => Riemann_SLAU2
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT

#else
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  Riemann_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_SLAU)
  Riemann_pointer => Riemann_SLAU
CASE(PRM_RIEMANN_SLAU2)
  Riemann_pointer => Riemann_SLAU2
CASE(PRM_RIEMANN_CH)
  Riemann_pointer => Riemann_CH
CASE(PRM_RIEMANN_Average)
  Riemann_pointer => Riemann_FluxAverage
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE(PRM_RIEMANN_ROEL2)
  Riemann_pointer => Riemann_RoeL2
CASE(PRM_RIEMANN_ROEENTROPYFIX)
  RiemannBC_pointer => Riemann_RoeEntropyFix
CASE(PRM_RIEMANN_SLAU)
  Riemann_pointer => Riemann_SLAU
CASE(PRM_RIEMANN_SLAU2)
  Riemann_pointer => Riemann_SLAU2
CASE(PRM_RIEMANN_CH)
  Riemann_pointer => Riemann_CH
CASE(PRM_RIEMANN_Average)
  RiemannBC_pointer => Riemann_FluxAverage
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT
#endif /*SPLIT_DG*/
END SUBROUTINE InitRiemann


!==================================================================================================================================
!> Computes the numerical flux for a side calling the flux calculation pointwise.
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann_Side(Nloc,FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                          :: Nloc      !< local polynomial degree
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_L       !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: U_R       !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L   !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_R   !< primitive solution at right side of the interface
REAL,DIMENSION(3          ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv,t1,t2  !> normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: FOut      !< advective flux
LOGICAL,INTENT(IN)                                          :: doBC      !< marker whether side is a BC side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
PROCEDURE(RiemannInt),POINTER :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

DO j=0,ZDIM(Nloc); DO i=0,Nloc
  ! Momentum has to be rotated using the normal system individual for each
  ! left state: U_L
  U_LL(EXT_DENS)=U_L(DENS,i,j)
  U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
  U_LL(EXT_ENER)=U_L(ENER,i,j)
  U_LL(EXT_PRES)=UPrim_L(PRES,i,j)
  U_LL(EXT_RHOK)=U_L(RHOK,i,j)
  U_LL(EXT_RHOG)=U_L(RHOG,i,j)
  U_LL(EXT_TKE )=UPrim_L(TKE,i,j)
  U_LL(EXT_OMG )=UPrim_L(OMG,i,j)

  ! rotate velocity in normal and tangential direction
  U_LL(EXT_VEL1)=DOT_PRODUCT(UPrim_L(VELV,i,j),nv(:,i,j))
  U_LL(EXT_VEL2)=DOT_PRODUCT(UPrim_L(VELV,i,j),t1(:,i,j))
  U_LL(EXT_MOM1)=U_LL(EXT_DENS)*U_LL(EXT_VEL1)
  U_LL(EXT_MOM2)=U_LL(EXT_DENS)*U_LL(EXT_VEL2)
#if PP_dim==3
  U_LL(EXT_VEL3)=DOT_PRODUCT(UPrim_L(VELV,i,j),t2(:,i,j))
  U_LL(EXT_MOM3)=U_LL(EXT_DENS)*U_LL(EXT_VEL3)
#else
  U_LL(EXT_VEL3)=0.
  U_LL(EXT_MOM3)=0.
#endif
  ! right state: U_R
  U_RR(EXT_DENS)=U_R(DENS,i,j)
  U_RR(EXT_SRHO)=1./U_RR(EXT_DENS)
  U_RR(EXT_ENER)=U_R(ENER,i,j)
  U_RR(EXT_PRES)=UPrim_R(PRES,i,j)
  U_RR(EXT_RHOK)=U_R(RHOK,i,j)
  U_RR(EXT_RHOG)=U_R(RHOG,i,j)
  U_RR(EXT_TKE )=UPrim_R(TKE,i,j)
  U_RR(EXT_OMG )=UPrim_R(OMG,i,j)
  ! rotate momentum in normal and tangential direction
  U_RR(EXT_VEL1)=DOT_PRODUCT(UPRIM_R(VELV,i,j),nv(:,i,j))
  U_RR(EXT_VEL2)=DOT_PRODUCT(UPRIM_R(VELV,i,j),t1(:,i,j))
  U_RR(EXT_MOM1)=U_RR(EXT_DENS)*U_RR(EXT_VEL1)
  U_RR(EXT_MOM2)=U_RR(EXT_DENS)*U_RR(EXT_VEL2)
#if PP_dim==3
  U_RR(EXT_VEL3)=DOT_PRODUCT(UPRIM_R(VELV,i,j),t2(:,i,j))
  U_RR(EXT_MOM3)=U_RR(EXT_DENS)*U_RR(EXT_VEL3)
#else
  U_RR(EXT_VEL3)=0.
  U_RR(EXT_MOM3)=0.
#endif

#ifndef SPLIT_DG
  CALL EvalEulerFlux1D_fast(U_LL,F_L)
  CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

  CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j)=F(DENS)
  Fout(MOMV,i,j)=nv(:,i,j)*F(MOM1)  &
               + t1(:,i,j)*F(MOM2)  &
#if PP_dim==3
               + t2(:,i,j)*F(MOM3)
#else
               + 0.
#endif
  Fout(ENER,i,j)=F(ENER)
  Fout(RHOK,i,j)=F(RHOK)
  Fout(RHOG,i,j)=F(RHOG)
END DO; END DO

END SUBROUTINE Riemann_Side


!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotated: don't use it after this routine!!
!> Attention 2: numerical flux is backrotated at the end of the routine!!
!==================================================================================================================================
SUBROUTINE Riemann_Point(FOut,U_L,U_R,UPrim_L,UPrim_R,nv,t1,t2,doBC)
! MODULES
USE MOD_Flux         ,ONLY:EvalEulerFlux1D_fast
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_L        !< conservative solution at left side of the interface
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U_R        !< conservative solution at right side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_L    !< primitive solution at left side of the interface
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim_R    !< primitive solution at right side of the interface
REAL,DIMENSION(3          ),INTENT(IN)  :: nv,t1,t2   !< normal vector and tangential vectors at side
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: FOut       !< advective flux
LOGICAL,INTENT(IN)                      :: doBC       !< marker whether side is a BC side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar) :: F_L,F_R,F
REAL,DIMENSION(PP_2Var) :: U_LL,U_RR
PROCEDURE(RiemannInt),POINTER :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

! Momentum has to be rotated using the normal system individual for each
! left state: U_L
U_LL(EXT_DENS)=U_L(DENS)
U_LL(EXT_SRHO)=1./U_LL(EXT_DENS)
U_LL(EXT_ENER)=U_L(ENER)
U_LL(EXT_PRES)=UPrim_L(PRES)
U_LL(EXT_RHOK)=U_L(RHOK)
U_LL(EXT_RHOG)=U_L(RHOG)
U_LL(EXT_TKE )=UPrim_L(TKE)
U_LL(EXT_OMG )=UPrim_L(OMG)

! rotate velocity in normal and tangential direction
U_LL(EXT_VEL1)=DOT_PRODUCT(UPrim_L(VELV),nv(:))
U_LL(EXT_VEL2)=DOT_PRODUCT(UPrim_L(VELV),t1(:))
U_LL(EXT_MOM1)=U_LL(EXT_DENS)*U_LL(EXT_VEL1)
U_LL(EXT_MOM2)=U_LL(EXT_DENS)*U_LL(EXT_VEL2)
#if PP_dim==3
U_LL(EXT_VEL3)=DOT_PRODUCT(UPrim_L(VELV),t2(:))
U_LL(EXT_MOM3)=U_LL(EXT_DENS)*U_LL(EXT_VEL3)
#else
U_LL(EXT_VEL3)=0.
U_LL(EXT_MOM3)=0.
#endif
! right state: U_R
U_RR(EXT_DENS)=U_R(DENS)
U_RR(EXT_SRHO)=1./U_RR(EXT_DENS)
U_RR(EXT_ENER)=U_R(ENER)
U_RR(EXT_PRES)=UPrim_R(PRES)
U_RR(EXT_RHOK)=U_R(RHOK)
U_RR(EXT_RHOG)=U_R(RHOG)
U_RR(EXT_TKE )=UPrim_R(TKE)
U_RR(EXT_OMG )=UPrim_R(OMG)
! rotate momentum in normal and tangential direction
U_RR(EXT_VEL1)=DOT_PRODUCT(UPRIM_R(VELV),nv(:))
U_RR(EXT_VEL2)=DOT_PRODUCT(UPRIM_R(VELV),t1(:))
U_RR(EXT_MOM1)=U_RR(EXT_DENS)*U_RR(EXT_VEL1)
U_RR(EXT_MOM2)=U_RR(EXT_DENS)*U_RR(EXT_VEL2)
#if PP_dim==3
U_RR(EXT_VEL3)=DOT_PRODUCT(UPRIM_R(VELV),t2(:))
U_RR(EXT_MOM3)=U_RR(EXT_DENS)*U_RR(EXT_VEL3)
#else
U_RR(EXT_VEL3)=0.
U_RR(EXT_MOM3)=0.
#endif

# ifndef SPLIT_DG
CALL EvalEulerFlux1D_fast(U_LL,F_L)
CALL EvalEulerFlux1D_fast(U_RR,F_R)
#endif /*SPLIT_DG*/

CALL Riemann_loc(F_L,F_R,U_LL,U_RR,F)

! Back rotate the normal flux into Cartesian direction
Fout(DENS)=F(DENS)
Fout(MOMV)=nv(:)*F(MOM1)  &
          +t1(:)*F(MOM2)  &
#if PP_dim==3
          +t2(:)*F(MOM3)
#else
          + 0.
#endif
Fout(ENER)=F(ENER)
Fout(RHOK)=F(RHOK)
Fout(RHOG)=F(RHOG)
END SUBROUTINE Riemann_Point


#if PARABOLIC
!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Side(Nloc,F,UPrim_L,UPrim_R, &
                            gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                           ,muSGS_L,muSGS_R &
#endif
                           )
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
USE MOD_Lifting_Vars ,ONLY: diffFluxX_L,diffFluxY_L,diffFluxZ_L
USE MOD_Lifting_Vars ,ONLY: diffFluxX_R,diffFluxY_R,diffFluxZ_R
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                                             :: Nloc     !< local polynomial degree
                                                               !> solution in primitive variables at left/right side of interface
REAL,DIMENSION(PP_nVarPrim   ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: UPrim_L,UPrim_R
                                                               !> solution gradients in x/y/z-direction left/right of interface
REAL,DIMENSION(PP_nVarLifting,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ,0:Nloc,0:ZDIM(Nloc)),INTENT(OUT) :: F   !< viscous flux
#if EDDYVISCOSITY
REAL,DIMENSION(1             ,0:Nloc,0:ZDIM(Nloc)),INTENT(IN)  :: muSGS_L,muSGS_R   !> eddy viscosity left/right of the interface
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                       :: p,q
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux3D(Nloc,UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                                ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  &
#if EDDYVISCOSITY
                   ,muSGS_L &
#endif
      )
CALL EvalDiffFlux3D(Nloc,UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                                ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  &
#if EDDYVISCOSITY
                   ,muSGS_R&
#endif
      )
! Arithmetic mean of the fluxes
DO q=0,ZDIM(Nloc); DO p=0,Nloc
  F(:,p,q)=0.5*(nv(1,p,q)*(diffFluxX_L(:,p,q)+diffFluxX_R(:,p,q)) &
               +nv(2,p,q)*(diffFluxY_L(:,p,q)+diffFluxY_R(:,p,q)) &
               +nv(3,p,q)*(diffFluxZ_L(:,p,q)+diffFluxZ_R(:,p,q)))
END DO; END DO
END SUBROUTINE ViscousFlux_Side

!==================================================================================================================================
!> Computes the viscous NSE diffusion fluxes in all directions to approximate the numerical flux
!> Actually not a Riemann solver, only here for coding reasons
!==================================================================================================================================
SUBROUTINE ViscousFlux_Point(F,UPrim_L,UPrim_R, &
                             gradUx_L,gradUy_L,gradUz_L,gradUx_R,gradUy_R,gradUz_R,nv &
#if EDDYVISCOSITY
                            ,muSGS_L,muSGS_R &
#endif
                            )
! MODULES
USE MOD_Flux         ,ONLY: EvalDiffFlux3D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> solution in primitive variables at left/right side of the interface
REAL,DIMENSION(PP_nVarPrim   ),INTENT(IN)  :: UPrim_L,UPrim_R
                                           !> solution gradients in x/y/z-direction left/right of the interface
REAL,DIMENSION(PP_nVarLifting),INTENT(IN)  :: gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R
REAL,DIMENSION(3             ),INTENT(IN)  :: nv  !< normal vector
REAL,DIMENSION(PP_nVar       ),INTENT(OUT) :: F   !< viscous flux
#if EDDYVISCOSITY
REAL,INTENT(IN)                            :: muSGS_L,muSGS_R    !> eddy viscosity left/right of the interface
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)  :: diffFluxX_L,diffFluxY_L,diffFluxZ_L
REAL,DIMENSION(PP_nVar)  :: diffFluxX_R,diffFluxY_R,diffFluxZ_R
!==================================================================================================================================
! Don't forget the diffusion contribution, my young padawan
! Compute NSE Diffusion flux
CALL EvalDiffFlux3D(UPrim_L,   gradUx_L,   gradUy_L,   gradUz_L  &
                           ,diffFluxX_L,diffFluxY_L,diffFluxZ_L  &
#if EDDYVISCOSITY
                   ,muSGS_L &
#endif
      )
CALL EvalDiffFlux3D(UPrim_R,   gradUx_R,   gradUy_R,   gradUz_R  &
                           ,diffFluxX_R,diffFluxY_R,diffFluxZ_R  &
#if EDDYVISCOSITY
                   ,muSGS_R&
#endif
      )
! Arithmetic mean of the fluxes
F(:)=0.5*(nv(1)*(diffFluxX_L(:)+diffFluxX_R(:)) &
         +nv(2)*(diffFluxY_L(:)+diffFluxY_R(:)) &
         +nv(3)*(diffFluxZ_L(:)+diffFluxZ_R(:)))
END SUBROUTINE ViscousFlux_Point
#endif /* PARABOLIC */


!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
PPURE SUBROUTINE Riemann_LF(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )
#ifndef SPLIT_DG
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(CONS) - U_LL(CONS)))
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! compute surface flux
F = F - 0.5*LambdaMax*(U_RR(CONS) - U_LL(CONS))
#endif /*SPLIT_DG*/
END SUBROUTINE Riemann_LF

!=================================================================================================================================
!> Harten-Lax-Van-Leer Riemann solver resolving contact discontinuity
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLLC(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: KappaM1!,kappa
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                           !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                           !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F    !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL    :: RoeVel(3),RoeH,Roec,absVel
REAL    :: Ssl,Ssr,SStar
REAL    :: U_Star(PP_nVar),EStar
REAL    :: sMu_L,sMu_R
#if DECOUPLE==0
REAL    :: RoeK
#endif
!REAL    :: c_L,c_R
!=================================================================================================================================
! HLLC flux

! Version A: Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - SPEEDOFSOUND_HE(U_LL)
!Ssr = U_RR(EXT_VEL1) + SPEEDOFSOUND_HE(U_RR)

! Version B: Basic Davis estimate for wave speed
!c_L = SPEEDOFSOUND_HE(U_LL)
!c_R = SPEEDOFSOUND_HE(U_RR)
!Ssl = MIN(U_LL(EXT_VEL1) - c_L,U_RR(EXT_VEL1) - c_R)
!Ssr = MAX(U_LL(EXT_VEL1) + c_L,U_RR(EXT_VEL1) + c_R)

! Version C: Better Roe estimate for wave speeds Davis, Einfeldt
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L       )     * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
#if DECOUPLE==0
RoeK      = (SqrtRho_R*U_RR(EXT_TKE)  + SqrtRho_L*U_LL(EXT_TKE) ) * sSqrtRho
Roec      = SQRT(KappaM1*(RoeH-0.5*absVel-RoeK))
#else
Roec      = SQRT(KappaM1*(RoeH-0.5*absVel))
#endif
Ssl       = RoeVel(1) - Roec
Ssr       = RoeVel(1) + Roec

! TODO(Shimushu): Fix HLLC, considering k into account.

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  sMu_L = Ssl - U_LL(EXT_VEL1)
  sMu_R = Ssr - U_RR(EXT_VEL1)
  SStar = (U_RR(EXT_PRES) - U_LL(EXT_PRES) + U_LL(EXT_MOM1)*sMu_L - U_RR(EXT_MOM1)*sMu_R) / (U_LL(EXT_DENS)*sMu_L - U_RR(EXT_DENS)*sMu_R)
  IF ((Ssl .LE. 0.).AND.(SStar .GE. 0.)) THEN
    EStar  = TOTALENERGY_HE(U_LL) + (SStar-U_LL(EXT_VEL1))*(SStar + U_LL(EXT_PRES)*U_LL(EXT_SRHO)/sMu_L)
    U_Star = U_LL(EXT_DENS) * sMu_L/(Ssl-SStar) * (/ 1., SStar, U_LL(EXT_VEL2:EXT_VEL3), EStar, U_LL(EXT_TKE:EXT_OMG) /)
    F=F_L+Ssl*(U_Star-U_LL(CONS))
  ELSE
    EStar  = TOTALENERGY_HE(U_RR) + (SStar-U_RR(EXT_VEL1))*(SStar + U_RR(EXT_PRES)*U_RR(EXT_SRHO)/sMu_R)
    U_Star = U_RR(EXT_DENS) * sMu_R/(Ssr-SStar) * (/ 1., SStar, U_RR(EXT_VEL2:EXT_VEL3), EStar, U_RR(EXT_TKE:EXT_OMG) /)
    F=F_R+Ssr*(U_Star-U_RR(CONS))
  END IF
END IF ! subsonic case
END SUBROUTINE Riemann_HLLC

!=================================================================================================================================
!> Roe's approximate Riemann solver using the Harten and Hymen II entropy fix, see
!> Pelanti, Marica & Quartapelle, Luigi & Vigevano, L & Vigevano, Luigi. (2018):
!>  A review of entropy fixes as applied to Roe's linearization.
!=================================================================================================================================
PPURE SUBROUTINE Riemann_RoeEntropyFix(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
#ifdef SPLIT_DG
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iVar
REAL                    :: c_L,c_R
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL                    :: RoeVel(3),RoeH,Roec,RoeDens
REAL,DIMENSION(5)       :: r1,r2,r3,r4,r5,a,al,ar,Delta_U,Alpha  ! Roe eigenvectors
REAL                    :: tmp,da
REAL                    :: Ma_loc,cl,cr,lambdaMax
!=================================================================================================================================
c_L       = SPEEDOFSOUND_HE(U_LL)
c_R       = SPEEDOFSOUND_HE(U_RR)
H_L       = TOTALENTHALPY_NS_HE(U_LL)
H_R       = TOTALENTHALPY_NS_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)
RoeDens   = SQRT(U_LL(EXT_DENS)*U_RR(EXT_DENS))
! Roe+Pike version of Roe Riemann solver

! calculate jump
Delta_U(DELTA_U1)   = U_RR(EXT_DENS) - U_LL(EXT_DENS)
Delta_U(DELTA_UV)   = U_RR(EXT_VELV) - U_LL(EXT_VELV)
Delta_U(DELTA_U5)   = U_RR(EXT_PRES) - U_LL(EXT_PRES)

! low Mach-Number fix
Ma_loc = MIN(SQRT(absVel)/(Roec*SQRT(kappa)),1.0)
Delta_U(DELTA_UV) = Delta_U(DELTA_UV) * Ma_loc

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate wave strenghts
tmp      = 0.5/(Roec*Roec)
Alpha(1) = tmp*(Delta_U(DELTA_U5)-RoeDens*Roec*Delta_U(DELTA_U2))
Alpha(2) = Delta_U(DELTA_U1) - Delta_U(DELTA_U5)*2.*tmp
Alpha(3) = RoeDens*Delta_U(DELTA_U3)
Alpha(4) = RoeDens*Delta_U(DELTA_U4)
Alpha(5) = tmp*(Delta_U(DELTA_U5)+RoeDens*Roec*Delta_U(DELTA_U2))

! Harten+Hyman entropy fix (apply only for acoustic waves, don't fix r)

al(1) = U_LL(EXT_VEL1) - c_L
al(2) = U_LL(EXT_VEL1)
al(3) = U_LL(EXT_VEL1)
al(4) = U_LL(EXT_VEL1)
al(5) = U_LL(EXT_VEL1) + c_L
ar(1) = U_RR(EXT_VEL1) - c_R
ar(2) = U_RR(EXT_VEL1)
ar(3) = U_RR(EXT_VEL1)
ar(4) = U_RR(EXT_VEL1)
ar(5) = U_RR(EXT_VEL1) + c_R
! HH1
!IF(ABS(a(1)).LT.da1) a(1)=da1
!IF(ABS(a(5)).LT.da5) a(5)=da5
! HH2
DO iVar=1,5
  da = MAX(0.,a(iVar)-al(iVar),ar(iVar)-a(iVar))

  IF(ABS(a(iVar)).LT.da) THEN
    a(iVar)=0.5*(a(iVar)*a(iVar)/da+da)
  ELSE
    a(iVar) = ABS(a(iVar))
  END IF
END DO

cl = SQRT(kappa * U_LL(EXT_PRES) / U_LL(EXT_DENS))
cr = SQRT(kappa * U_RR(EXT_PRES) / U_RR(EXT_DENS))
lambdaMax = MAX(ABS(U_LL(EXT_VEL1)) + cl, ABS(U_RR(EXT_VEL1)) + cr)

#ifndef SPLIT_DG
F = 0.5 * (F_L + F_R)
#else
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
#endif
! assemble Roe flux
F(1:5) = F(1:5) - 0.5 * ( &
    Alpha(1) * a(1) * r1 + &
    Alpha(2) * a(2) * r2 + &
    Alpha(3) * a(3) * r3 + &
    Alpha(4) * a(4) * r4 + &
    Alpha(5) * a(5) * r5)
F(ENER) = F(ENER) - 0.5 * lambdaMax * (U_RR(EXT_RHOK) - U_LL(EXT_RHOK))
F(RHOK) = F(RHOK) - 0.5 * lambdaMax * (U_RR(EXT_RHOK) - U_LL(EXT_RHOK))
F(RHOG) = F(RHOG) - 0.5 * lambdaMax * (U_RR(EXT_RHOG) - U_LL(EXT_RHOG))
END SUBROUTINE Riemann_RoeEntropyFix


!=================================================================================================================================
!> low mach number Roe's approximate Riemann solver according to Oßwald(2015)
!=================================================================================================================================
PPURE SUBROUTINE Riemann_RoeL2(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars  ,ONLY: kappaM1,kappa
#ifdef SPLIT_DG
USE MOD_SplitFlux ,ONLY: SplitDGSurface_pointer
#endif /*SPLIT_DG*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: H_L,H_R
REAL                    :: SqrtRho_L,SqrtRho_R,sSqrtRho
REAL                    :: RoeVel(3),RoeH,Roec,absVel
REAL                    :: Ma_loc ! local Mach-Number
REAL,DIMENSION(5)       :: a,r1,r2,r3,r4,r5  ! Roe eigenvectors
REAL                    :: Alpha1,Alpha2,Alpha3,Alpha4,Alpha5,Delta_U(5+1)
REAL                    :: cl,cr,lambdaMax
!=================================================================================================================================
! Roe flux
H_L       = TOTALENTHALPY_NS_HE(U_LL)
H_R       = TOTALENTHALPY_NS_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))

sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
RoeH      = (SqrtRho_R*H_R+SqrtRho_L*H_L) * sSqrtRho
Roec      = ROEC_RIEMANN_H(RoeH,RoeVel)

! mean eigenvalues and eigenvectors
a  = (/ RoeVel(1)-Roec, RoeVel(1), RoeVel(1), RoeVel(1), RoeVel(1)+Roec      /)
r1 = (/ 1.,             a(1),      RoeVel(2), RoeVel(3), RoeH-RoeVel(1)*Roec /)
r2 = (/ 1.,             RoeVel(1), RoeVel(2), RoeVel(3), 0.5*absVel          /)
r3 = (/ 0.,             0.,        1.,        0.,        RoeVel(2)           /)
r4 = (/ 0.,             0.,        0.,        1.,        RoeVel(3)           /)
r5 = (/ 1.,             a(5),      RoeVel(2), RoeVel(3), RoeH+RoeVel(1)*Roec /)

! calculate differences
Delta_U(1:5) = U_RR(1:5) - U_LL(1:5)
Delta_U(5)    = Delta_U(5) - (U_RR(EXT_RHOK) - U_LL(EXT_RHOK))
Delta_U(DELTA_U6)   = Delta_U(DELTA_U5)-(Delta_U(DELTA_U3)-RoeVel(2)*Delta_U(DELTA_U1))*RoeVel(2) - &
                      (Delta_U(DELTA_U4)-RoeVel(3)*Delta_U(DELTA_U1))*RoeVel(3)

! low Mach-Number fix
Ma_loc = MIN(SQRT(absVel)/(Roec*SQRT(kappa)),1.0)
Delta_U(DELTA_UV) = Delta_U(DELTA_UV) * Ma_loc

! calculate factors
Alpha3 = Delta_U(DELTA_U3) - RoeVel(2)*Delta_U(DELTA_U1)
Alpha4 = Delta_U(DELTA_U4) - RoeVel(3)*Delta_U(DELTA_U1)
Alpha2 = ALPHA2_RIEMANN_H(RoeH,RoeVel,Roec,Delta_U)
Alpha1 = 0.5/Roec * (Delta_U(DELTA_U1)*(RoeVel(1)+Roec) - Delta_U(DELTA_U2) - Roec*Alpha2)
Alpha5 = Delta_U(DELTA_U1) - Alpha1 - Alpha2

cl = SQRT(kappa * U_LL(EXT_PRES) / U_LL(EXT_DENS))
cr = SQRT(kappa * U_RR(EXT_PRES) / U_RR(EXT_DENS))
lambdaMax = MAX(ABS(U_LL(EXT_VEL1)) + cl, ABS(U_RR(EXT_VEL1)) + cr)

#ifndef SPLIT_DG
! assemble Roe flux
F = 0.5 * (F_L + F_R)
F(1:5)=F(1:5) - 0.5 * ( &
       Alpha1*ABS(a(1))*r1 + &
       Alpha2*ABS(a(2))*r2 + &
       Alpha3*ABS(a(3))*r3 + &
       Alpha4*ABS(a(4))*r4 + &
       Alpha5*ABS(a(5))*r5)
#else
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
! assemble Roe flux
F(1:5) = F(1:5) - 0.5 * (Alpha1*ABS(a(1))*r1 + &
                         Alpha2*ABS(a(2))*r2 + &
                         Alpha3*ABS(a(3))*r3 + &
                         Alpha4*ABS(a(4))*r4 + &
                         Alpha5*ABS(a(5))*r5)
#endif /*SPLIT_DG*/
F(ENER) = F(ENER) - 0.5 * lambdaMax * (U_RR(EXT_RHOK) - U_LL(EXT_RHOK))
F(RHOK) = F(RHOK) - 0.5 * lambdaMax * (U_RR(EXT_RHOK) - U_LL(EXT_RHOK))
F(RHOG) = F(RHOG) - 0.5 * lambdaMax * (U_RR(EXT_RHOG) - U_LL(EXT_RHOG))
END SUBROUTINE Riemann_RoeL2

!=================================================================================================================================
!> Standard Harten-Lax-Van-Leer Riemann solver without contact discontinuity
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLL(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars, ONLY: KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roek,Roec
REAL    :: Ssl,Ssr
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
Roek      = (SqrtRho_R*U_RR(EXT_TKE)  + SqrtRho_L*U_LL(EXT_TKE))  * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH-Roek,RoeVel)
! HLL flux
! Basic Davis estimate for wave speed
!Ssl = U_LL(EXT_VEL1) - c_L
!Ssr = U_RR(EXT_VEL1) + c_R
! Better Roe estimate for wave speeds Davis, Einfeldt
Ssl = RoeVel(1) - Roec
Ssr = RoeVel(1) + Roec
! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLL

!=================================================================================================================================
!> Harten-Lax-Van-Leer-Einfeldt Riemann solver
!=================================================================================================================================
PPURE SUBROUTINE Riemann_HLLE(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: SqrtRho_L,SqrtRho_R,sSqrtRho,absVel
REAL    :: RoeVel(3),RoeH,Roek,Roec
REAL    :: Ssl,Ssr,beta
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
SqrtRho_L = SQRT(U_LL(EXT_DENS))
SqrtRho_R = SQRT(U_RR(EXT_DENS))
sSqrtRho  = 1./(SqrtRho_L+SqrtRho_R)
! Roe mean values
RoeVel    = (SqrtRho_R*U_RR(EXT_VELV) + SqrtRho_L*U_LL(EXT_VELV)) * sSqrtRho
RoeH      = (SqrtRho_R*H_R            + SqrtRho_L*H_L)            * sSqrtRho
Roek      = (SqrtRho_R*U_RR(EXT_TKE)  + SqrtRho_L*U_LL(EXT_TKE))  * sSqrtRho
absVel    = DOT_PRODUCT(RoeVel,RoeVel)
Roec      = ROEC_RIEMANN_H(RoeH-Roek,RoeVel)
! HLLE flux (positively conservative)
beta=BETA_RIEMANN_H()
SsL=MIN(RoeVel(1)-Roec,U_LL(EXT_VEL1) - beta*SPEEDOFSOUND_HE(U_LL), 0.)
SsR=MAX(RoeVel(1)+Roec,U_RR(EXT_VEL1) + beta*SPEEDOFSOUND_HE(U_RR), 0.)

! positive supersonic speed
IF(Ssl .GE. 0.)THEN
  F=F_L
! negative supersonic speed
ELSEIF(Ssr .LE. 0.)THEN
  F=F_R
! subsonic case
ELSE
  F=(Ssr*F_L-Ssl*F_R+Ssl*Ssr*(U_RR(EXT_CONS)-U_LL(EXT_CONS)))/(Ssr-Ssl)
END IF ! subsonic case
END SUBROUTINE Riemann_HLLE

!=================================================================================================================================
!> SLAU
!=================================================================================================================================
PPURE SUBROUTINE Riemann_SLAU(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: cL, cR, c0_5, VL2, VR2, VAvg, MM, chi, MaL, MaR, fL, fR
REAL    :: pTilde, g, Vn, massFlux
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
cL        = SPEEDOFSOUND_HE(U_LL)
cR        = SPEEDOFSOUND_HE(U_RR)
c0_5      = 0.5 * (cL + cR) ! Eq. 2.3h
VL2       = DOT_PRODUCT(U_LL(EXT_VELV), U_LL(EXT_VELV))
VR2       = DOT_PRODUCT(U_RR(EXT_VELV), U_RR(EXT_VELV))
VAvg      = SQRT(0.5 * (VL2 + VR2))
MM        = MIN(1.0, VAvg / c0_5) ! Eq. 2.3e
chi       = (1.0 - MM)**2 ! Eq. 2.3d
MaL       = U_LL(EXT_VEL1) / c0_5 ! Eq. 2.3g
MaR       = U_RR(EXT_VEL1) / c0_5 ! Eq. 2.3g
IF (ABS(MaL) .GE. 1.0) THEN
  fL = 0.5 * (1.0 + SIGN(1.0, MaL)) ! Eq. 2.3f
ELSE
  fL = 0.25 * (MaL + 1.0)**2 * (2.0 - MaL) ! Eq. 2.3f
END IF
IF (ABS(MaR) .GE. 1.0) THEN
  fR = 0.5 * (1.0 - SIGN(1.0, MaR)) ! Eq. 2.3f
ELSE
  fR = 0.25 * (MaR - 1.0)**2 * (2.0 + MaR) ! Eq. 2.3f
END IF
! Eq. 2.3c
pTilde    = 0.5 * (U_LL(EXT_PRES) + U_RR(EXT_PRES)) + &
            0.5 * (fL - fR) * (U_LL(EXT_PRES) - U_RR(EXT_PRES)) + &
            (1.0 - chi) * (fL + fR - 1.0) * 0.5 * (U_LL(EXT_PRES) + U_RR(EXT_PRES))
! Eq. 2.3l
g         = -MAX(MIN(MaL, 0.0), -1.0) * MIN(MAX(MaR, 0.0), 1.0)
Vn        = (ABS(U_LL(EXT_MOM1)) + ABS(U_RR(EXT_MOM1))) / (U_LL(EXT_DENS) + U_RR(EXT_DENS)) ! Eq. 2.3k
! Eq. 2.3i
massFlux  = 0.5 * (( &
              U_LL(EXT_MOM1) + U_RR(EXT_MOM1) - &
              Vn * (U_RR(EXT_DENS) - U_LL(EXT_DENS))) * (1.0 - g) - &
              chi / c0_5 * (U_RR(EXT_PRES) - U_LL(EXT_PRES)))
IF (massFlux .GE. 0.0) THEN
  F = massFlux * (/1.0, U_LL(EXT_VELV), H_L, U_LL(EXT_TKE:EXT_OMG)/) ! Eq. 2.3a
ELSE
  F = massFlux * (/1.0, U_RR(EXT_VELV), H_R, U_RR(EXT_TKE:EXT_OMG)/) ! Eq. 2.3a
END IF
F(MOM1) = F(MOM1) + pTilde ! Eq. 2.3a
END SUBROUTINE Riemann_SLAU

!=================================================================================================================================
!> SLAU2
!=================================================================================================================================
PPURE SUBROUTINE Riemann_SLAU2(F_L,F_R,U_LL,U_RR,F)
!=================================================================================================================================
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,KappaM1
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                               !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                               !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F        !< resulting Riemann flux
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: H_L,H_R
REAL    :: cL, cR, c0_5, VL2, VR2, VAvg, MM, chi, MaL, MaR, fL, fR
REAL    :: pTilde, g, Vn, VnPlus, VnMinus, massFlux
!=================================================================================================================================
H_L       = TOTALENTHALPY_HE(U_LL)
H_R       = TOTALENTHALPY_HE(U_RR)
cL        = SPEEDOFSOUND_HE(U_LL)
cR        = SPEEDOFSOUND_HE(U_RR)
c0_5      = 0.5 * (cL + cR) ! Eq. 2.3h
! p. 1694, Shima and Kitamura, 2011
! c0_5      = SQRT(kappa*(U_LL(EXT_PRES)+U_RR(EXT_PRES))/(U_LL(EXT_DENS)+U_RR(EXT_DENS)))
VL2       = DOT_PRODUCT(U_LL(EXT_VELV), U_LL(EXT_VELV))
VR2       = DOT_PRODUCT(U_RR(EXT_VELV), U_RR(EXT_VELV))
VAvg      = SQRT(0.5 * (VL2 + VR2))
MM        = MIN(1.0, VAvg / c0_5) ! Eq. 2.3e
chi       = (1.0 - MM)**2 ! Eq. 2.3d
MaL       = U_LL(EXT_VEL1) / c0_5 ! Eq. 2.3g
MaR       = U_RR(EXT_VEL1) / c0_5 ! Eq. 2.3g
! Adapted from line 766,
! https://github.com/su2code/SU2/blob/master/SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp
IF (ABS(MaL) .LE. 1.0) THEN
  fL = 0.25 * (MaL + 1.0)**2 * (2.0 - MaL) ! Eq. 2.3f
ELSE IF (MaL .GE. 0.0) THEN
  fL = 1.0 ! Eq. 2.3f
ELSE
  fL = 0.0 ! Eq. 2.3f
END IF
IF (ABS(MaR) .LE. 1.0) THEN
  fR = 0.25 * (MaR - 1.0)**2 * (2.0 + MaR) ! Eq. 2.3f
ELSE IF (MaR .GE. 0.0) THEN
  fR = 0.0 ! Eq. 2.3f
ELSE
  fR = 1.0 ! Eq. 2.3f
END IF
! IF (ABS(MaL) .GE. 1.0) THEN
!   fL = 0.5 * (1.0 + SIGN(1.0, MaL)) ! Eq. 2.3f
! ELSE
!   fL = 0.25 * (MaL + 1.0)**2 * (2.0 - MaL) ! Eq. 2.3f
! END IF
! IF (ABS(MaR) .GE. 1.0) THEN
!   fR = 0.5 * (1.0 - SIGN(1.0, MaR)) ! Eq. 2.3f
! ELSE
!   fR = 0.25 * (MaR - 1.0)**2 * (2.0 + MaR) ! Eq. 2.3f
! END IF
! Eq. 2.3c
pTilde    = 0.5 * (U_LL(EXT_PRES) + U_RR(EXT_PRES)) + &
            0.5 * (fL - fR) * (U_LL(EXT_PRES) - U_RR(EXT_PRES)) + &
            VAvg * (fL + fR - 1.0) * 0.5 * (U_LL(EXT_DENS) + U_RR(EXT_DENS)) * c0_5
! Eq. 2.3l
g         = -MAX(MIN(MaL, 0.0), -1.0) * MIN(MAX(MaR, 0.0), 1.0)
Vn        = (ABS(U_LL(EXT_MOM1)) + ABS(U_RR(EXT_MOM1))) / (U_LL(EXT_DENS) + U_RR(EXT_DENS)) ! Eq. 2.3k
VnPlus    = (1.0 - g) * Vn + g * ABS(U_LL(EXT_VEL1)) ! Eq. 2.3j
VnMinus   = (1.0 - g) * Vn + g * ABS(U_RR(EXT_VEL1)) ! Eq. 2.3j
! Eq. 2.3i
massFlux  = 0.5 * ( &
              U_LL(EXT_DENS) * (U_LL(EXT_VEL1) + VnPlus) + &
              U_RR(EXT_DENS) * (U_RR(EXT_VEL1) - VnMinus) - &
              chi / c0_5 * (U_RR(EXT_PRES) - U_LL(EXT_PRES)))
! Eq. 29, Shima and Kitamura, 2011
! massFlux  = 0.5 * (( &
!               U_LL(EXT_MOM1) + U_RR(EXT_MOM1) - &
!               Vn * (U_RR(EXT_DENS) - U_LL(EXT_DENS))) * (1.0 - g) - &
!               chi / c0_5 * (U_RR(EXT_PRES) - U_LL(EXT_PRES)))
IF (massFlux .GE. 0.0) THEN
  F = massFlux * (/1.0, U_LL(EXT_VELV), H_L, U_LL(EXT_TKE:EXT_OMG)/) ! Eq. 2.3a
ELSE
  F = massFlux * (/1.0, U_RR(EXT_VELV), H_R, U_RR(EXT_TKE:EXT_OMG)/) ! Eq. 2.3a
END IF
F(MOM1) = F(MOM1) + pTilde ! Eq. 2.3a
END SUBROUTINE Riemann_SLAU2

#ifdef SPLIT_DG
!==================================================================================================================================
!> Riemann solver using purely the average fluxes
!==================================================================================================================================
PPURE SUBROUTINE Riemann_FluxAverage(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)
END SUBROUTINE Riemann_FluxAverage

!==================================================================================================================================
!> kinetic energy preserving and entropy consistent flux according to Chandrashekar (2012)
!==================================================================================================================================
PPURE SUBROUTINE Riemann_CH(F_L,F_R,U_LL,U_RR,F)
! MODULES
USE MOD_EOS_Vars      ,ONLY: Kappa,sKappaM1
USE MOD_SplitFlux     ,ONLY: SplitDGSurface_pointer,GetLogMean
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
                                                !> extended solution vector on the left/right side of the interface
REAL,DIMENSION(PP_2Var),INTENT(IN) :: U_LL,U_RR
                                                !> advection fluxes on the left/right side of the interface
REAL,DIMENSION(PP_nVar),INTENT(IN) :: F_L,F_R
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F         !< resulting Riemann flux
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                               :: LambdaMax
REAL                               :: beta_LL,beta_RR   ! auxiliary variables for the inverse Temperature
REAL                               :: rhoMean           ! auxiliary variable for the mean density
REAL                               :: uMean,vMean,wMean ! auxiliary variable for the average velocities
REAL                               :: betaLogMean       ! auxiliary variable for the logarithmic mean inverse temperature
!==================================================================================================================================
! Lax-Friedrichs
LambdaMax = MAX( ABS(U_RR(EXT_VEL1)),ABS(U_LL(EXT_VEL1)) ) + MAX( SPEEDOFSOUND_HE(U_LL),SPEEDOFSOUND_HE(U_RR) )

! average quantities
rhoMean = 0.5*(U_LL(EXT_DENS) + U_RR(EXT_DENS))
uMean   = 0.5*(U_LL(EXT_VEL1) + U_RR(EXT_VEL1))
vMean   = 0.5*(U_LL(EXT_VEL2) + U_RR(EXT_VEL2))
wMean   = 0.5*(U_LL(EXT_VEL3) + U_RR(EXT_VEL3))

! inverse temperature
beta_LL = 0.5*U_LL(EXT_DENS)/U_LL(EXT_PRES)
beta_RR = 0.5*U_RR(EXT_DENS)/U_RR(EXT_PRES)

! logarithmic mean
CALL GetLogMean(beta_LL,beta_RR,betaLogMean)

! get split flux
CALL SplitDGSurface_pointer(U_LL,U_RR,F)

!compute flux
F(DENS:MOM3) = F(DENS:MOM3) - 0.5*LambdaMax*(U_RR(EXT_DENS:EXT_MOM3)-U_LL(EXT_DENS:EXT_MOM3))
F(ENER)      = F(ENER)      - 0.5*LambdaMax*( &
         (U_RR(EXT_DENS)-U_LL(EXT_DENS))*(0.5*sKappaM1/betaLogMean +0.5*(U_RR(EXT_VEL1)*U_LL(EXT_VEL1)+U_RR(EXT_VEL2)*U_LL(EXT_VEL2)+U_RR(EXT_VEL3)*U_LL(EXT_VEL3))) &
         +rhoMean*uMean*(U_RR(EXT_VEL1)-U_LL(EXT_VEL1)) + rhoMean*vMean*(U_RR(EXT_VEL2)-U_LL(EXT_VEL2)) + rhoMean*wMean*(U_RR(EXT_VEL3)-U_LL(EXT_VEL3)) &
         +0.5*rhoMean*sKappaM1*(1./beta_RR - 1./beta_LL))
F(RHOK:RHOG) = F(RHOK:RHOG) - 0.5*LambdaMax*(U_RR(EXT_RHOK:EXT_RHOG)-U_LL(EXT_RHOK:EXT_RHOG))

END SUBROUTINE Riemann_CH
#endif /*SPLIT_DG*/


!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann

END MODULE MOD_Riemann
