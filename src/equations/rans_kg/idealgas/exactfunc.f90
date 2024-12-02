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

!==================================================================================================================================
!> Subroutines providing exactly evaluated functions used in initialization or boundary conditions.
!==================================================================================================================================
MODULE MOD_Exactfunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersExactFunc
  MODULE PROCEDURE DefineParametersExactFunc
END INTERFACE

INTERFACE InitExactFunc
  MODULE PROCEDURE InitExactFunc
END INTERFACE

INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE

INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE


PUBLIC::DefineParametersExactFunc
PUBLIC::InitExactFunc
PUBLIC::ExactFunc
PUBLIC::CalcSource
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters of exact functions
!==================================================================================================================================
SUBROUTINE DefineParametersExactFunc()
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
CALL prms%SetSection("Exactfunc")
CALL prms%CreateIntFromStringOption('IniExactFunc', "Exact function to be used for computing initial solution.")
CALL addStrListEntry('IniExactFunc','testcase'          ,-1)
CALL addStrListEntry('IniExactFunc','testcase'          ,0)
CALL addStrListEntry('IniExactFunc','refstate'          ,1)
CALL addStrListEntry('IniExactFunc','sinedens'          ,2)
CALL addStrListEntry('IniExactFunc','sinedensx'         ,21)
CALL addStrListEntry('IniExactFunc','lindens'           ,3)
CALL addStrListEntry('IniExactFunc','sinevel'           ,4)
CALL addStrListEntry('IniExactFunc','sinevelx'          ,41)
CALL addStrListEntry('IniExactFunc','sinevely'          ,42)
CALL addStrListEntry('IniExactFunc','sinevelz'          ,43)
CALL addStrListEntry('IniExactFunc','roundjet'          ,5)
CALL addStrListEntry('IniExactFunc','cylinder'          ,6)
CALL addStrListEntry('IniExactFunc','shuvortex'         ,7)
CALL addStrListEntry('IniExactFunc','couette'           ,8)
CALL addStrListEntry('IniExactFunc','cavity'            ,9)
CALL addStrListEntry('IniExactFunc','shock'             ,10)
CALL addStrListEntry('IniExactFunc','sod'               ,11)
CALL addStrListEntry('IniExactFunc','dmr'               ,13)
CALL addStrListEntry('IniExactFunc','harmonicgausspulse',14)
CALL addStrListEntry('IniExactFunc','blast_shock'       ,111)
CALL addStrListEntry('IniExactFunc','turbulentchannel'  ,517)
CALL addStrListEntry('IniExactFunc','readfromfile'      ,301)
#if PARABOLIC
CALL addStrListEntry('IniExactFunc','blasius'  ,1338)
#endif
CALL prms%CreateRealArrayOption(    'AdvVel',       "Advection velocity (v1,v2,v3) required for exactfunction CASE(2,21,4,8)")
CALL prms%CreateRealOption(         'IniAmplitude', "Amplitude for synthetic test case")
CALL prms%CreateRealOption(         'IniFrequency', "Frequency for synthetic test case")
CALL prms%CreateRealOption(         'MachShock',    "Parameter required for CASE(10)", '1.5')
CALL prms%CreateRealOption(         'PreShockDens', "Parameter required for CASE(10)", '1.0')
CALL prms%CreateRealArrayOption(    'IniCenter',    "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealArrayOption(    'IniAxis',      "Shu Vortex CASE(7) (x,y,z)")
CALL prms%CreateRealOption(         'IniHalfwidth', "Shu Vortex CASE(7)", '0.2')
CALL prms%CreateRealOption(         'JetRadius',    "Roundjet CASE(5/33)", '1.0')
CALL prms%CreateRealOption(         'JetEnd',       "Roundjet CASE(5/33)", '10.0')
CALL prms%CreateRealOption(         'Ramping',      "Subsonic mass inflow CASE(28)"  , '1.0')
CALL prms%CreateRealOption(         'P_Parameter',  "Couette-Poiseuille flow CASE(8)", '0.0')
CALL prms%CreateRealOption(         'U_Parameter',  "Couette-Poiseuille flow CASE(8)", '0.01')
CALL prms%CreateRealOption(         'AmplitudeFactor',         "Harmonic Gauss Pulse CASE(14)", '0.1')
CALL prms%CreateRealOption(         'HarmonicFrequency',       "Harmonic Gauss Pulse CASE(14)", '400')
CALL prms%CreateRealOption(         'SigmaSqr',                "Harmonic Gauss Pulse CASE(14)", '0.1')
CALL prms%CreateRealOption(         'Re_tau',                  "Re_tau CASE(517)")
CALL prms%CreateStringOption(       'BCfile',                  "BC file CASE(301)")
#if PARABOLIC
CALL prms%CreateRealOption(         'delta99_in',              "Blasius boundary layer CASE(1338)")
CALL prms%CreateRealArrayOption(    'x_in',                    "Blasius boundary layer CASE(1338)")
#endif

CALL prms%CreateLogicalOption(      'useBCFile',               "Use BC file CASE(301)", 'F')

! Extra source term
CALL prms%SetSection("SourceTerm")
CALL prms%CreateIntFromStringOption('IniSourceTerm', "Source term function to be used for computing source term.", "None")
CALL addStrListEntry('IniSourceTerm','None'             ,0)
CALL addStrListEntry('IniSourceTerm','ConstantBodyForce',1)

CALL prms%CreateRealArrayOption('ConstantBodyForce', "Constant body force to be added, IniSourceTerm==ConstantBodyForce")
CALL prms%CreateRealOption('ConstantHeatSource',     "Constant heat source")
CALL prms%CreateRealOption('Fluctuation',            "Fluctation of the constant body force (-1 < f < 1)")

END SUBROUTINE DefineParametersExactFunc

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitExactFunc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_ExactFunc_Vars
USE MOD_Equation_Vars      ,ONLY: IniExactFunc,IniRefState,IniSourceTerm,ConstantBodyForce,ConstantHeatSource,Fluctuation
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

IniExactFunc = GETINTFROMSTR('IniExactFunc')
IniRefState  = GETINT('IniRefState', "-1")
! Read in boundary parameters
SELECT CASE (IniExactFunc)
  CASE(2,21) ! sinus
    AdvVel            = GETREALARRAY('AdvVel',3)
    IniFrequency      = GETREAL('IniFrequency','0.5')
    IniAmplitude      = GETREAL('IniAmplitude','0.3')
  CASE(3) ! synthetic test cases
    AdvVel            = GETREALARRAY('AdvVel',3)
  CASE(4,41,42,43) ! synthetic test cases
    AdvVel            = GETREALARRAY('AdvVel',3)
    IniFrequency      = GETREAL('IniFrequency','1.0')
    IniAmplitude      = GETREAL('IniAmplitude','0.1')
  CASE(7) ! Shu Vortex
    IniCenter         = GETREALARRAY('IniCenter',3,'(/0.,0.,0./)')
    IniAxis           = GETREALARRAY('IniAxis',3,'(/0.,0.,1./)')
    IniAmplitude      = GETREAL('IniAmplitude','0.2')
    IniHalfwidth      = GETREAL('IniHalfwidth','0.2')
  CASE(8) ! couette-poiseuille flow
    P_Parameter       = GETREAL('P_Parameter')
    U_Parameter       = GETREAL('U_Parameter')
  CASE(10) ! shock
    MachShock         = GETREAL('MachShock')
    PreShockDens      = GETREAL('PreShockDens')
  CASE(14)
    HarmonicFrequency = GETREAL('HarmonicFrequency')
    AmplitudeFactor   = GETREAL('AmplitudeFactor')
    SiqmaSqr          = GETREAL('SigmaSqr')
  CASE(517) ! turbulent channel
    Re_tau            = GETREAL('Re_tau')
  ! CASE(301) ! readfromfile
  !   BCfile            = GETSTR('BCFile')
  !   open(newunit=BCFileID, file=BCfile, status="old", action="read")
  !   read(BCFileID, *) BCLength
  !   ALLOCATE(BCData(1+PP_nVar,BCLength))
  !   read(BCFileID, *) BCData
  !   ! print * , BCData
  !   close(BCFileID)
#if PARABOLIC
  CASE(1338) ! Blasius boundary layer solution
    delta99_in      = GETREAL('delta99_in')
    x_in            = GETREALARRAY('x_in',2,'(/0.,0./)')
    BlasiusInitDone = .TRUE. ! Mark Blasius init as done so we don't read the parameters again in BC init
#endif
  CASE DEFAULT
    ! Everything defined, do nothing
END SELECT ! IniExactFunc

#if PP_dim==2
SELECT CASE (IniExactFunc)
CASE(43) ! synthetic test cases
  CALL CollectiveStop(__STAMP__,'The selected exact function is not available in 2D!')
CASE(2,3,4,41,42) ! synthetic test cases
  IF(AdvVel(3).NE.0.) THEN
    CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set AdvVel(3) = 0!')
  END IF
END SELECT
#endif

IF(GETLOGICAL('useBCFile'))THEN
  BCfile            = GETSTR('BCFile')
  open(newunit=BCFileID, file=BCfile, status="old", action="read")
  read(BCFileID, *) BCLength
  ALLOCATE(BCData(1+PP_nVar,BCLength))
  read(BCFileID, *) BCData
  close(BCFileID)
ENDIF

IniSourceTerm = GETINTFROMSTR('IniSourceTerm')
SELECT CASE (IniSourceTerm)
CASE(0) ! None
CASE(1) ! ConstantBodyForce
  ConstantBodyForce = GETREALARRAY('ConstantBodyForce', 3)
  ConstantHeatSource = GETREAL("ConstantHeatSource", "0.0")
  print*, "ConstantHeatSource = ",ConstantHeatSource
  Fluctuation = GETREAL('Fluctuation', '0.0')
#if PP_dim==2
  IF(ConstantBodyForce(3).NE.0.) THEN
    CALL CollectiveStop(__STAMP__,'You are computing in 2D! Please set ConstantBodyForce(3) = 0!')
  END IF
#endif
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,'Unknown IniSourceTerm!')
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_stdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!> t is the actual time
!> dt is only needed to compute the time dependent boundary values for the RK scheme
!> for each function resu and the first and second time derivative resu_t and resu_tt have to be defined (is trivial for constants)
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu,RefStateOpt)
! MODULES
USE MOD_Preproc        ,ONLY: PP_PI
USE MOD_Globals        ,ONLY: Abort
USE MOD_Mathtools      ,ONLY: CROSS
USE MOD_Eos_Vars       ,ONLY: Kappa,sKappaM1,KappaM1,KappaP1,R
USE MOD_Exactfunc_Vars ,ONLY: IniCenter,IniHalfwidth,IniAmplitude,IniFrequency,IniAxis,AdvVel
USE MOD_Exactfunc_Vars ,ONLY: MachShock,PreShockDens
USE MOD_Exactfunc_Vars ,ONLY: P_Parameter,U_Parameter
USE MOD_Exactfunc_Vars ,ONLY: JetRadius,JetEnd
USE MOD_Exactfunc_Vars ,ONLY: Re_tau
USE MOD_Exactfunc_Vars ,ONLY: BCLength,BCData
USE MOD_Equation_Vars  ,ONLY: IniRefState,RefStateCons,RefStatePrim
USE MOD_Timedisc_Vars  ,ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
USE MOD_TestCase       ,ONLY: ExactFuncTestcase
USE MOD_EOS            ,ONLY: PrimToCons,ConsToPrim
#if PARABOLIC
USE MOD_Eos_Vars       ,ONLY: mu0
USE MOD_Exactfunc_Vars ,ONLY: delta99_in,x_in
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction          !< determines the exact function
REAL,INTENT(IN)                 :: x(3)                   !< physical coordinates
REAL,INTENT(IN)                 :: tIn                    !< solution time (Runge-Kutta stage)
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< state in conservative variables
INTEGER,INTENT(IN),OPTIONAL     :: RefStateOpt            !< refstate to be used for exact func
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: RefState
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar),ov ! state in conservative variables
REAL                            :: Frequency,Amplitude
REAL                            :: Omega
REAL                            :: Vel(3),Cent(3),a
REAL                            :: Prim(PP_nVarPrim)
REAL                            :: r_len
REAL                            :: Ms,xs
REAL                            :: Resul(PP_nVar),Resur(PP_nVar)
REAL                            :: random
REAL                            :: du, dTemp, RT, r2       ! aux var for SHU VORTEX,isentropic vortex case 12
REAL                            :: pi_loc,phi,radius       ! needed for cylinder potential flow
REAL                            :: h,sRT,pexit,pentry   ! needed for Couette-Poiseuille
REAL                            :: y_plus
INTEGER                         :: flag, upper, lower
#if PARABOLIC
! needed for blasius BL
INTEGER                         :: nSteps,i
REAL                            :: eta,deta,deta2,f,fp,fpp,fppp,fbar,fpbar,fppbar,fpppbar
REAL                            :: x_eff(3),x_offset(3)
#endif
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration
IF (PRESENT(RefStateOpt)) THEN
  RefState = RefStateOpt
ELSE
  RefState = IniRefState
END IF

Resu   =0.
Resu_t =0.
Resu_tt=0.

! Determine the value, the first and the second time derivative
SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(0)
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)
CASE(1) ! constant
  Resu = RefStateCons(:,RefState)
CASE(2) ! sinus
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(DENS)   = 1.
  prim(VELV) = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)*(1.+Amplitude*SIN(Omega*SUM(cent(1:3))))
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(DENS)=-Amplitude*cos(Omega*SUM(cent(1:3)))*ov
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
    ! g''(t)
    Resu_tt(DENS)=-Amplitude*sin(Omega*SUM(cent(1:3)))*ov**2.
    Resu_tt(MOMV)=Resu_tt(DENS)*prim(VELV)
    Resu_tt(ENER)=0.5*SUM(Resu_tt(MOMV)*prim(VELV))
  END IF
CASE(21) ! sinus x
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=2.*PP_Pi*Frequency
  ! base flow
  prim(DENS)   = 1.
  prim(VELV) = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)*(1.+Amplitude*SIN(Omega*cent(1)))
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e

  IF(fullBoundaryOrder)THEN
    ov=Omega*SUM(vel)
    ! g'(t)
    Resu_t(DENS)=-Amplitude*cos(Omega*cent(1))*ov
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
    ! g''(t)
    Resu_tt(DENS)=-Amplitude*sin(Omega*cent(1))*ov**2.
    Resu_tt(MOMV)=Resu_tt(DENS)*prim(VELV)
    Resu_tt(ENER)=0.5*SUM(Resu_tt(MOMV)*prim(VELV))
  END IF
CASE(3) ! linear in rho
  ! base flow
  prim(DENS)   = 100.
  prim(VELV)   = AdvVel
  prim(PRES)   = 1.
  Vel=prim(VELV)
  cent=x-Vel*tEval
  prim(DENS)=prim(DENS)+SUM(AdvVel*cent)
  ! g(t)
  Resu(DENS)=prim(DENS) ! rho
  Resu(MOMV)=prim(DENS)*prim(VELV) ! rho*vel
  Resu(ENER)=prim(PRES)*sKappaM1+0.5*SUM(Resu(MOMV)*prim(VELV)) ! rho*e
  IF(fullBoundaryOrder)THEN
    ! g'(t)
    Resu_t(DENS)=-SUM(Vel)
    Resu_t(MOMV)=Resu_t(DENS)*prim(VELV) ! rho*vel
    Resu_t(ENER)=0.5*SUM(Resu_t(MOMV)*prim(VELV))
  END IF
CASE(4) ! oblique sine wave (in x,y,z for 3D calculations, and x,y for 2D)
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi

  ! g(t)
#if (PP_dim == 3)
  Resu(DENS:MOM3) = 2.+ Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
  Resu(DENS:MOM2) = 2.+ Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
  Resu(MOM3)   = 0.
#endif
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    ! g'(t)
#if (PP_dim == 3)
    Resu_t(DENS:MOM3)=(-a)*Amplitude*cos(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
    Resu_t(DENS:MOM2)=(-a)*Amplitude*cos(Omega*SUM(x(1:PP_dim)) - a*tEval)
    Resu_t(MOM3)=0.
#endif
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
#if (PP_dim == 3)
    Resu_tt(DENS:MOM3)=-a*a*Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
#else
    Resu_tt(DENS:MOM2)=-a*a*Amplitude*sin(Omega*SUM(x(1:PP_dim)) - a*tEval)
    Resu_tt(MOM3)=0.
#endif
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
CASE(41) ! SINUS in x
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(1)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS:MOM1)=2.+ Amplitude*sin(Omega*x(1) - a*tEval)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_t = 0.
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS:MOM1)=(-a)*Amplitude*cos(Omega*x(1) - a*tEval)
    Resu_t(ENER)=2.*Resu(1)*Resu_t(1)
    ! g''(t)
    Resu_tt(DENS:MOM1)=-a*a*Amplitude*sin(Omega*x(1) - a*tEval)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
CASE(42) ! SINUS in y
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(2)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS)=2.+ Amplitude*sin(Omega*x(2) - a*tEval)
  Resu(MOM2)=Resu(DENS)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS)=(-a)*Amplitude*cos(Omega*x(2) - a*tEval)
    Resu_t(MOM2)=Resu_t(DENS)
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
    Resu_tt(DENS)=-a*a*Amplitude*sin(Omega*x(2) - a*tEval)
    Resu_tt(MOM2)=Resu_tt(DENS)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
#if PP_dim==3
CASE(43) ! SINUS in z
  Frequency=IniFrequency
  Amplitude=IniAmplitude
  Omega=PP_Pi*Frequency
  a=AdvVel(3)*2.*PP_Pi
  ! g(t)
  Resu = 0.
  Resu(DENS)=2.+ Amplitude*sin(Omega*x(3) - a*tEval)
  Resu(MOM3)=Resu(DENS)
  Resu(ENER)=Resu(DENS)*Resu(DENS)
  IF(fullBoundaryOrder)THEN
    Resu_tt = 0.
    ! g'(t)
    Resu_t(DENS)=(-a)*Amplitude*cos(Omega*x(3) - a*tEval)
    Resu_t(MOM3)=Resu_t(DENS)
    Resu_t(ENER)=2.*Resu(DENS)*Resu_t(DENS)
    ! g''(t)
    Resu_tt(DENS)=-a*a*Amplitude*sin(Omega*x(3) - a*tEval)
    Resu_tt(MOM3)=Resu_tt(DENS)
    Resu_tt(ENER)=2.*(Resu_t(DENS)*Resu_t(DENS) + Resu(DENS)*Resu_tt(DENS))
  END IF
#endif
CASE(5) !Roundjet Bogey Bailly 2002, Re=65000, x-axis is jet axis
  prim(DENS)  =1.
  prim(VELV)  =0.
  prim(PRES)  =1./Kappa
  prim(TEMP)  = prim(PRES)/(prim(DENS)*R)
  ! Jet inflow (from x=0, diameter 2.0)
  ! Initial jet radius: rj=1.
  ! Momentum thickness: delta_theta0=0.05=1/20
  ! Re=65000
  ! Uco=0.
  ! Uj=0.9
  r_len=SQRT((x(2)*x(2)+x(3)*x(3)))
  prim(VEL1)=0.9*0.5*(1.+TANH((JetRadius-r_len)/JetRadius*10.))
  CALL RANDOM_NUMBER(random)
  ! Random disturbance +-5%; uniform distribution between -1,1
  random=0.05*2.*(random-0.5)
  prim(VEL1)=prim(VEL1)+random*prim(VEL1)
  prim(VEL2)=x(2)/r_len*0.5*random*prim(VEL1)
  prim(VEL3)=x(3)/r_len*0.5*random*prim(VEL1)
  CALL PrimToCons(prim,ResuL)
  prim(VELV)  =0.
  CALL PrimToCons(prim,ResuR)
  ! after x/r0=10 blend to ResuR
  Resu=ResuL+(ResuR-ResuL)*0.5*(1.+tanh(x(1)/JetRadius-JetEnd))
CASE(6)  ! Cylinder flow
  IF(tEval .EQ. 0.)THEN   ! Initialize potential flow
    prim(DENS)=RefStatePrim(DENS,RefState)  ! Density
    prim(VEL3)=0.                           ! VelocityZ=0. (2D flow)
    ! Calculate cylinder coordinates (0<phi<Pi/2)
    pi_loc=ASIN(1.)*2.
    IF(x(1) .LT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc+phi
      ELSE
        phi=pi_loc-phi
      END IF
    ELSEIF(x(1) .GT. 0.)THEN
      phi=ATAN(ABS(x(2))/ABS(x(1)))
      IF(x(2) .LT. 0.) phi=2.*pi_loc-phi
    ELSE
      IF(x(2) .LT. 0.)THEN
        phi=pi_loc*1.5
      ELSE
        phi=pi_loc*0.5
      END IF
    END IF
    ! Calculate radius**2
    radius=x(1)*x(1)+x(2)*x(2)
    ! Calculate velocities, radius of cylinder=0.5
    prim(VEL1)=RefStatePrim(VEL1,RefState)*(COS(phi)**2*(1.-0.25/radius)+SIN(phi)**2*(1.+0.25/radius))
    prim(VEL2)=RefStatePrim(VEL1,RefState)*(-2.)*SIN(phi)*COS(phi)*0.25/radius
    ! Calculate pressure, RefState(2)=u_infinity
    prim(PRES)=RefStatePrim(PRES,RefState) + &
            0.5*prim(DENS)*(RefStatePrim(VEL1,RefState)*RefStatePrim(VEL1,RefState)-prim(VEL1)*prim(VEL1)-prim(VEL2)*prim(VEL2))
    prim(TEMP) = prim(PRES)/(prim(DENS)*R)
  ELSE  ! Use RefState as BC
    prim=RefStatePrim(:,RefState)
  END IF  ! t=0
  CALL PrimToCons(prim,resu)
CASE(7) ! SHU VORTEX,isentropic vortex
  ! base flow
  prim=RefStatePrim(:,RefState)  ! Density
  ! ini-Parameter of the Example
  vel=prim(VELV)
  RT=prim(PRES)/prim(DENS) !ideal gas
  cent=(iniCenter+vel*tEval)    !centerpoint time dependant
  cent=x-cent                   ! distance to centerpoint
  cent=CROSS(iniAxis,cent)      !distance to axis, tangent vector, length r
  cent=cent/iniHalfWidth        !Halfwidth is dimension 1
  r2=SUM(cent*cent) !
  du = IniAmplitude/(2.*PP_Pi)*exp(0.5*(1.-r2))   ! vel. perturbation
  dTemp = -kappaM1/(2.*kappa*RT)*du**2            ! adiabatic
  prim(DENS)=prim(DENS)*(1.+dTemp)**(1.*skappaM1) !rho
  prim(VELV)=prim(VELV)+du*cent(:)                !v
#if PP_dim == 2
  prim(VEL3)=0.
#endif
  prim(PRES) = prim(PRES)*(1.+dTemp)**(kappa/kappaM1) !p
  prim(TEMP) = prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,resu)
CASE(8) !Couette-Poiseuille flow between plates: exact steady lamiar solution with height=1 !
        !(Gao, hesthaven, Warburton)
  RT=1. ! Hesthaven: Absorbing layers for weakly compressible flows
  sRT=1./RT
  ! size of domain must be [-0.5, 0.5]^2 -> (x*y)
  h=0.5
  prim(VEL1)     = U_parameter*(0.5*(1.+x(2)/h) + P_parameter*(1-(x(2)/h)**2))
  prim(VEL2:VEL3)= 0.
  pexit=0.9996
  pentry=1.0004
  prim(PRES)= ( ( x(1) - (-0.5) )*( pexit - pentry) / ( 0.5 - (-0.5)) ) + pentry
  prim(DENS)=prim(PRES)*sRT
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,Resu)
CASE(9) !lid driven cavity flow from Gao, Hesthaven, Warburton
        !"Absorbing layers for weakly compressible flows", to appear, JSC, 2016
        ! Special "regularized" driven cavity BC to prevent singularities at corners
        ! top BC assumed to be in x-direction from 0..1
  Prim = RefStatePrim(:,RefState)
  IF (x(1).LT.0.2) THEN
    prim(VEL1)=1000*4.9333*x(1)**4-1.4267*1000*x(1)**3+0.1297*1000*x(1)**2-0.0033*1000*x(1)
  ELSEIF (x(1).LE.0.8) THEN
    prim(VEL1)=1.0
  ELSE
    prim(VEL1)=1000*4.9333*x(1)**4-1.8307*10000*x(1)**3+2.5450*10000*x(1)**2-1.5709*10000*x(1)+10000*0.3633
  ENDIF
  CALL PrimToCons(prim,Resu)
CASE(10) ! shock
  prim=0.

  ! pre-shock
  prim(DENS) = PreShockDens
  Ms         = MachShock

  prim(PRES)=prim(DENS)/Kappa
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,Resur)

  ! post-shock
  prim(VEL2)=prim(DENS) ! temporal storage of pre-shock density
  prim(DENS)=prim(DENS)*((KappaP1)*Ms*Ms)/(KappaM1*Ms*Ms+2.)
  prim(PRES)=prim(PRES)*(2.*Kappa*Ms*Ms-KappaM1)/(KappaP1)
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  IF (prim(VEL1) .EQ. 0.0) THEN
    prim(VEL1)=Ms*(1.-prim(VEL2)/prim(DENS))
  ELSE
    prim(VEL1)=prim(VEL1)*prim(VEL2)/prim(DENS)
  END IF
  prim(3)=0. ! reset temporal storage
  CALL PrimToCons(prim,Resul)
  xs=5.+Ms*tEval ! 5. bei 10x10x10 Rechengebiet
  ! Tanh boundary
  Resu=-0.5*(Resul-Resur)*TANH(5.0*(x(1)-xs))+Resur+0.5*(Resul-Resur)
CASE(11) ! Sod Shock tube
  xs = 0.5
  IF (X(1).LE.xs) THEN
    Resu = RefStateCons(:,1)
  ELSE
    Resu = RefStateCons(:,2)
  END IF
CASE(111) ! Sedov blast wave
  xs = 0.5
  Ms = SQRT(SUM((X(1:3)-1.5)**2))
  IF ((Ms.LE.xs).AND.(Ms.NE.0)) THEN
    prim(DENS)      = 1.3416
    prim(VEL1:VEL3) = 0.3615*(X(1:3)-1.5)/Ms
    prim(PRES)      = 1.5133
  ELSE
    prim(DENS)      = 1.
    prim(VELV)      = 0.
    prim(PRES)      = 1.
  END IF
  prim(TEMP)=prim(PRES)/(prim(DENS)*R)
  CALL PrimToCons(prim,resu)
CASE(12) ! Shu Osher density fluctuations shock wave interaction
  IF (x(1).LT.-4.0) THEN
    prim(DENS)      = 3.857143
    prim(VEL1)      = 2.629369
    prim(VEL2:VEL3) = 0.
    prim(PRES)      = 10.33333
  ELSE
    prim(DENS)      = 1.+0.2*SIN(5.*x(1))
    prim(VELV)      = 0.
    prim(PRES)      = 1.
  END IF
  CALL PrimToCons(prim,resu)
CASE(13) ! DoubleMachReflection (see e.g. http://www.astro.princeton.edu/~jstone/Athena/tests/dmr/dmr.html )
  IF (x(1).EQ.0.) THEN
    prim = RefStatePrim(:,1)
  ELSE IF (x(1).EQ.4.0) THEN
    prim = RefStatePrim(:,2)
  ELSE
    IF (x(1).LT.1./6.+(x(2)+20.*t)*1./3.**0.5) THEN
      prim = RefStatePrim(:,1)
    ELSE
      prim = RefStatePrim(:,2)
    END IF
  END IF
  CALL PrimToCons(prim,resu)
CASE(14) ! harmonic gauss pulse
  Resu = RefStateCons(:,RefState)
CASE(517) ! turbulent channel
  IF (x(2).LT.0.) THEN
    y_plus = Re_tau * (x(2) + 1.)
  ELSE
    y_plus = Re_tau * (1. - x(2))
  END IF
  Prim = RefStatePrim(:,RefState)
  Prim(VELV) = 0.
  Prim(VEL1) = 1. / 0.41 * LOG(1. + 0.41 * y_plus) + 7.8 * (1. - EXP(-y_plus / 11.) - y_plus / 11. * EXP(-y_plus / 3.))
  Prim(VEL1) = Prim(VEL1) * RefStatePrim(VEL1, RefState)
  Amplitude = 0.1 * Prim(VEL1)

  ! copied from src/testcase/channel/testcase.f90
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(20.0*PP_PI*(x(2)/(2.0)))*SIN(20.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(30.0*PP_PI*(x(2)/(2.0)))*SIN(30.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(35.0*PP_PI*(x(2)/(2.0)))*SIN(35.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(40.0*PP_PI*(x(2)/(2.0)))*SIN(40.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(45.0*PP_PI*(x(2)/(2.0)))*SIN(45.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL1)=Prim(VEL1)+Amplitude*SIN(50.0*PP_PI*(x(2)/(2.0)))*SIN(50.0*PP_PI*(x(3)/(2*PP_PI)))

  Prim(VEL2)=Prim(VEL2)+Amplitude*SIN(30.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(30.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL2)=Prim(VEL2)+Amplitude*SIN(35.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(35.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL2)=Prim(VEL2)+Amplitude*SIN(40.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(40.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL2)=Prim(VEL2)+Amplitude*SIN(45.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(45.0*PP_PI*(x(3)/(2*PP_PI)))
  Prim(VEL2)=Prim(VEL2)+Amplitude*SIN(50.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(50.0*PP_PI*(x(3)/(2*PP_PI)))

  Prim(VEL3)=Prim(VEL3)+Amplitude*SIN(30.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(30.0*PP_PI*(x(2)/(2.0)))
  Prim(VEL3)=Prim(VEL3)+Amplitude*SIN(35.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(35.0*PP_PI*(x(2)/(2.0)))
  Prim(VEL3)=Prim(VEL3)+Amplitude*SIN(40.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(40.0*PP_PI*(x(2)/(2.0)))
  Prim(VEL3)=Prim(VEL3)+Amplitude*SIN(45.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(45.0*PP_PI*(x(2)/(2.0)))
  Prim(VEL3)=Prim(VEL3)+Amplitude*SIN(50.0*PP_PI*(x(1)/(4*PP_PI)))*SIN(50.0*PP_PI*(x(2)/(2.0)))

  CALL PrimToCons(Prim, Resu)
CASE(301) ! readfromfile
  lower = 1
  upper = BCLength
  DO WHILE (.TRUE.)
    IF (upper - lower .EQ. 1) THEN
      Prim(1:PP_nVar) = ( &
          (BCData(1,upper) - x(2)) * BCData(2:(1+PP_nVar),lower) + &
          (x(2) - BCData(1,lower)) * BCData(2:(1+PP_nVar),upper) ) / &
          (BCData(1,upper) - BCData(1,lower))
      EXIT
    END IF
    i = (upper + lower) / 2
    IF (x(2) .EQ. BCData(1,i)) THEN
      Prim(1:PP_nVar) = BCData(2:(1+PP_nVar),i)
      EXIT
    ELSE IF (x(2) .GT. BCData(1,i)) THEN
      lower = i
    ELSE
      upper = i
    END IF
  END DO
      Prim(OMG) = 1. / SQRT( 0.09 * Prim(7) )
      Prim(TKE) = Prim(6)
      Prim(TEMP) = 0.
      CALL PrimToCons(Prim, Resu)
#if PARABOLIC
CASE(1338) ! blasius
  prim=RefStatePrim(:,RefState)
  ! calculate equivalent x for Blasius flat plate to have delta99_in at x_in
  x_offset(1)=(delta99_in/5)**2*prim(DENS)*prim(VEL1)/mu0-x_in(1)
  x_offset(2)=-x_in(2)
  x_offset(3)=0.
  x_eff=x+x_offset
  IF(x_eff(2).GT.0 .AND. x_eff(1).GT.0) THEN
    ! scale bl position in physical space to reference space, eta=5 is ~99% bl thickness
    eta=x_eff(2)*(prim(DENS)*prim(VEL1)/(mu0*x_eff(1)))**0.5

    deta=0.02 ! step size
    nSteps=CEILING(eta/deta)
    deta =eta/nSteps
    deta2=0.5*deta

    f=0.
    fp=0.
    fpp=0.332 ! default literature value, don't change if you don't know what you're doing
    fppp=0.
    !Blasius boundary layer
    DO i=1,nSteps
      ! predictor
      fbar    = f   + deta * fp
      fpbar   = fp  + deta * fpp
      fppbar  = fpp + deta * fppp
      fpppbar = -0.5*fbar*fppbar
      ! corrector
      f       = f   + deta2 * (fp   + fpbar)
      fp      = fp  + deta2 * (fpp  + fppbar)
      fpp     = fpp + deta2 * (fppp + fpppbar)
      fppp    = -0.5*f*fpp
    END DO
    prim(VEL2)=0.5*(mu0*prim(VEL1)/prim(DENS)/x_eff(1))**0.5*(fp*eta-f)
    prim(VEL1)=RefStatePrim(VEL1,RefState)*fp
  ELSE
    IF(x_eff(2).LE.0) THEN
      prim(VEL1)=0.
    END IF
  END IF
  CALL PrimToCons(prim,resu)
#endif
END SELECT ! ExactFunction
#if PP_dim==2
Resu(MOM3)=0.
#endif

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
#if PP_dim==2
  Resu_t=0.
  Resu_tt=0.
#endif
  SELECT CASE(CurrentStage)
  CASE(1)
    ! resu = g(t)
  CASE(2)
    ! resu = g(t) + dt/3*g'(t)
    Resu=Resu + dt*RKc(2)*Resu_t
  CASE(3)
    ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
    Resu=Resu + RKc(3)*dt*Resu_t + RKc(2)*RKb(2)*dt*dt*Resu_tt
  CASE DEFAULT
    ! Stop, works only for 3 Stage O3 LS RK
    CALL Abort(__STAMP__,&
               'Exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT
END IF


END SUBROUTINE ExactFunc

!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: sJ,nElems
USE MOD_EOS_Vars         ,ONLY: mu0
USE MOD_Equation_Vars    ,ONLY: s43,s23,Cmu,Comega1,Comega2,invSigmaG,invSigmaK
USE MOD_DG_Vars          ,ONLY: U
USE MOD_Lifting_Vars     ,ONLY: gradUx,gradUy,gradUz
USE MOD_EOS              ,ONLY: ConsToPrim
USE MOD_Equation_Vars    ,ONLY: IniSourceTerm,ConstantBodyForce,ConstantHeatSource,Fluctuation
#if FV_ENABLED
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisVolume
USE MOD_FV_Vars          ,ONLY: FV_Vdm,FV_Elems
#endif
#if PP_VISC == 1
USE MOD_Viscosity        ,ONLY: muSuth
#endif
USE MOD_EddyVisc_Vars    ,ONLY: muSGS,prodK,dissK,prodG,dissG,crossG,SijUij,dGidGi
! USE MOD_EddyVisc_Vars    ,ONLY: muSGS
USE MOD_EOS_Vars         ,ONLY: cp,Pr,kappa
USE MOD_EddyVisc_Vars    ,ONLY: PrSGS,fd
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< current solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< DG time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,iElem
REAL                :: Ut_src(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
REAL                :: UPrim(PP_nVarPrim)
REAL                :: muS,muT,kPos,gPos,muTOrig
REAL                :: muEffG,invR,invG
#if FV_ENABLED
REAL                :: Ut_src2(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ)
#endif
REAL                :: Sxx,Sxy,Syy,SijGradU
#if PP_dim == 3
REAL                :: Sxz, Syz, Szz
#endif
REAL                :: dGdG
REAL                :: bodyForce(3)
REAL                :: lambda
REAL                :: comp_f, comp_c, comp_t, comp_qx, comp_qy, comp_qz, comp_q, comp_u, comp_M, comp_utau, comp_Mtau, comp_Bq
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    CALL ConsToPrim(UPrim, U(:,i,j,k,iElem))

    muS = VISCOSITY_TEMPERATURE(UPrim(TEMP))
    muT = muSGS(1,i,j,k,iElem)
    kPos    = MAX( UPrim(TKE), 1.e-16 )
    gPos    = MAX( UPrim(OMG), 1.e-16 )

    muTOrig = MIN( Cmu * UPrim(DENS) * kPos * gPos**2, 10000.0 * muS )

    invR   = 1. / MAX( 0.01 * muS, muTOrig )
    invG   = SQRT(Cmu * UPrim(DENS) * kPos * invR)

    muEffG = (muS + invSigmaG * muTOrig)

    lambda = muS*cp/Pr + muT*cp/PrSGS

#if PP_dim==2
    ASSOCIATE( &
        ux => gradUx(LIFT_VEL1,i,j,k,iElem), uy => gradUy(LIFT_VEL1,i,j,k,iElem), &
        vx => gradUx(LIFT_VEL2,i,j,k,iElem), vy => gradUy(LIFT_VEL2,i,j,k,iElem), &
        kx => gradUx(LIFT_TKE ,i,j,k,iElem), ky => gradUy(LIFT_TKE ,i,j,k,iElem), &
        gx => gradUx(LIFT_OMG ,i,j,k,iElem), gy => gradUy(LIFT_OMG ,i,j,k,iElem))

    Sxx = 0.5 * (s43 * ux - s23 * vy)
    Syy = 0.5 * (s43 * vy - s23 * ux)
    Sxy = 0.5 * (uy + vx)

    SijGradU = Sxx * ux + Sxy * uy + Sxy * vx + Syy * vy

    dGdG = gx * gx + gy * gy

    comp_qx = lambda * gradUx(LIFT_TEMP,i,j,k,iElem) + 2.0 * (muS + muT) * (UPrim(VEL1) * Sxx + UPrim(VEL2) * Sxy)
    comp_qy = lambda * gradUy(LIFT_TEMP,i,j,k,iElem) + 2.0 * (muS + muT) * (UPrim(VEL1) * Sxy + UPrim(VEL2) * Syy)
    ! comp_q = lambda * SQRT(gradUx(LIFT_TEMP,i,j,k,iElem)**2+gradUy(LIFT_TEMP,i,j,k,iElem)**2)
    comp_q = SQRT(comp_qx**2 + comp_qy**2)
    comp_t = (muS + muT) * SQRT(2.0 * (Sxx**2 + Syy**2) + 4.0 * Sxy**2)
    comp_u = SQRT(UPrim(VEL1)**2 + UPrim(VEL2)**2)

    END ASSOCIATE
#else
    ASSOCIATE( &
        ux => gradUx(LIFT_VEL1,i,j,k,iElem), uy => gradUy(LIFT_VEL1,i,j,k,iElem), uz => gradUz(LIFT_VEL1,i,j,k,iElem), &
        vx => gradUx(LIFT_VEL2,i,j,k,iElem), vy => gradUy(LIFT_VEL2,i,j,k,iElem), vz => gradUz(LIFT_VEL2,i,j,k,iElem), &
        wx => gradUx(LIFT_VEL3,i,j,k,iElem), wy => gradUy(LIFT_VEL3,i,j,k,iElem), wz => gradUz(LIFT_VEL3,i,j,k,iElem), &
        kx => gradUx(LIFT_TKE ,i,j,k,iElem), ky => gradUy(LIFT_TKE ,i,j,k,iElem), kz => gradUz(LIFT_TKE ,i,j,k,iElem), &
        gx => gradUx(LIFT_OMG ,i,j,k,iElem), gy => gradUy(LIFT_OMG ,i,j,k,iElem), gz => gradUz(LIFT_OMG ,i,j,k,iElem))

    Sxx = 0.5 * (s43 * ux - s23 * (vy + wz))
    Syy = 0.5 * (s43 * vy - s23 * (ux + wz))
    Szz = 0.5 * (s43 * wz - s23 * (ux + vy))
    Sxy = 0.5 * (uy + vx)
    Sxz = 0.5 * (uz + wx)
    Syz = 0.5 * (vz + wy)

    SijGradU = &
        Sxx * ux + Sxy * uy + Sxz * uz + &
        Sxy * vx + Syy * vy + Syz * vz + &
        Sxz * wx + Syz * wy + Szz * wz

    dGdG = gx * gx + gy * gy + gz * gz

    comp_qx = lambda * gradUx(LIFT_TEMP,i,j,k,iElem) + 2.0 * (muS + muT) * (UPrim(VEL1) * Sxx + UPrim(VEL2) * Sxy + UPrim(VEL3) * Sxz)
    comp_qy = lambda * gradUy(LIFT_TEMP,i,j,k,iElem) + 2.0 * (muS + muT) * (UPrim(VEL1) * Sxy + UPrim(VEL2) * Syy + UPrim(VEL3) * Syz)
    comp_qz = lambda * gradUz(LIFT_TEMP,i,j,k,iElem) + 2.0 * (muS + muT) * (UPrim(VEL1) * Sxz + UPrim(VEL2) * Syz + UPrim(VEL3) * Szz)
    ! comp_q = lambda * SQRT(gradUx(LIFT_TEMP,i,j,k,iElem)**2+gradUy(LIFT_TEMP,i,j,k,iElem)**2+gradUz(LIFT_TEMP,i,j,k,iElem)**2)
    comp_q = SQRT(comp_qx**2 + comp_qy**2 + comp_qz**2)
    comp_t = (muS + muT) * SQRT(2.0 * (Sxx**2 + Syy**2 + Szz**2) + 4.0 * (Sxy**2 + Syz**2 + Sxz**2))
    comp_u = SQRT(UPrim(VEL1)**2 + UPrim(VEL2)**2 + UPrim(VEL3)**2)

    END ASSOCIATE
#endif

    comp_c    = SPEEDOFSOUND_H(UPrim(PRES),(1.0/UPrim(DENS)))
    comp_M    = comp_u / comp_c
    comp_utau = SQRT(comp_t / UPrim(DENS))
    comp_Mtau = comp_utau / comp_c
    comp_Bq   = comp_q / (UPrim(DENS) * Cp * UPrim(TEMP) * comp_utau)
    CALL CalcCompf(comp_f, comp_M, comp_Mtau, comp_Bq)

    IF (ALLOCATED(fd)) THEN
      comp_f = fd(i,j,k,iElem) + (1.0 - fd(i,j,k,iElem)) * comp_f
    END IF

    SijUij(i,j,k,iElem) = SijGradU
    dGidGi(i,j,k,iElem) = dGdG

    prodK (1,i,j,k,iElem) = 2. * muT * SijGradU
    IF (UPrim(TKE).GE.0.0) THEN
      dissK (1,i,j,k,iElem) = +(Cmu * (UPrim(DENS) * UPrim(TKE))**2 * invR)
    ELSE
      dissK (1,i,j,k,iElem) = -(Cmu * (UPrim(DENS) * UPrim(TKE))**2 * invR)
    END IF

    ! prodG (1,i,j,k,iElem) = comp_f * Comega2 * UPrim(DENS)**2 * kPos * gPos * 0.5 * invR
    prodG (1,i,j,k,iElem) = comp_f * Comega2 * UPrim(DENS) / (2 * Cmu) * invG
    dissG (1,i,j,k,iElem) = comp_f * Comega1 * Cmu * UPrim(DENS) * gPos**3 * SijGradU
    ! crossG(1,i,j,k,iElem) = 3.0 * muEffG * Cmu * UPrim(DENS) * kPos * gPos * invR * dGdG
    crossG(1,i,j,k,iElem) = muEffG * 3.0 * invG * dGdG

    Ut_src(RHOK,i,j,k) = prodK(1,i,j,k,iElem) - dissK(1,i,j,k,iElem)
    Ut_src(RHOG,i,j,k) = prodG(1,i,j,k,iElem) - dissG(1,i,j,k,iElem) - crossG(1,i,j,k,iElem)
  END DO; END DO; END DO ! i,j,k

#if FV_ENABLED
  IF (FV_Elems(iElem).GT.0) THEN ! FV elem
    CALL ChangeBasisVolume(2,PP_N,PP_N,FV_Vdm,Ut_src(RHOK:RHOG,:,:,:),Ut_src2(RHOK:RHOG,:,:,:))
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(RHOK:RHOG,i,j,k,iElem) = Ut(RHOK:RHOG,i,j,k,iElem)+Ut_src2(RHOK:RHOG,i,j,k)/sJ(i,j,k,iElem,1)
    END DO; END DO; END DO ! i,j,k
  ELSE
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(RHOK:RHOG,i,j,k,iElem) = Ut(RHOK:RHOG,i,j,k,iElem)+Ut_src(RHOK:RHOG,i,j,k)/sJ(i,j,k,iElem,0)
    END DO; END DO; END DO ! i,j,k
  END IF
#else
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    Ut(RHOK:RHOG,i,j,k,iElem) = Ut(RHOK:RHOG,i,j,k,iElem)+Ut_src(RHOK:RHOG,i,j,k)/sJ(i,j,k,iElem,0)
  END DO; END DO; END DO ! i,j,k
#endif

END DO

SELECT CASE (IniSourceTerm)
CASE(1) ! ConstantBodyForce
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut_src(:,i,j,k) = 0.
      IF (Fluctuation .NE. 0.0) THEN
        CALL RANDOM_NUMBER(bodyForce)
        bodyForce = ConstantBodyForce * (1.0 + Fluctuation * (2.0 * bodyForce - 1.0))
      ELSE
        bodyForce = ConstantBodyForce
      ENDIF
      Ut_src(MOMV,i,j,k) = bodyForce
      Ut_src(ENER,i,j,k) = dot_product(U(MOMV,i,j,k,iElem), bodyForce) / U(DENS,i,j,k,iElem) + ConstantHeatSource
    END DO; END DO; END DO ! i,j,k

#if FV_ENABLED
    IF (FV_Elems(iElem).GT.0) THEN ! FV elem
    CALL ChangeBasisVolume(PP_nVar,PP_N,PP_N,FV_Vdm,Ut_src(:,:,:,:),Ut_src2(:,:,:,:))
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src2(:,i,j,k)/sJ(i,j,k,iElem,1)
      END DO; END DO; END DO ! i,j,k
    ELSE
      DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
        Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
      END DO; END DO; END DO ! i,j,k
    END IF
#else
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem)+Ut_src(:,i,j,k)/sJ(i,j,k,iElem,0)
    END DO; END DO; END DO ! i,j,k
#endif

  END DO
END SELECT

END SUBROUTINE CalcSource

SUBROUTINE CalcCompf(f,M,Mtau,Bq)
IMPLICIT NONE
REAL,INTENT(OUT)    :: f
REAL,INTENT(IN)     :: M
REAL,INTENT(IN)     :: Mtau
REAL,INTENT(IN)     :: Bq

REAL                :: c1,c2,c3,c4,phi

c1 = 60.0
c2 = 5.0
c3 = 15.0/4.0 * EXP(-M/10.0)
c4 = -21.0/4.0 * EXP(-M/5.0)
phi = TANH(3.0/4.0 * M)

f = (1.0 - phi) * EXP(-c1 * Mtau**2 - c2 * Bq) + phi * EXP(c3 * Mtau + c4 * Bq)

END SUBROUTINE CalcCompf

END MODULE MOD_Exactfunc
