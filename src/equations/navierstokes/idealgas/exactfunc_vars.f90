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

!==================================================================================================================================
!> Contains the parameters needed for the different initializations of the Navier Stokes equations
!==================================================================================================================================
MODULE MOD_Exactfunc_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
REAL              :: MachShock         !< Shock Mach speed for ExactFunction = 6 (shock)
REAL              :: PreShockDens      !< Pre-shock density for ExactFunction = 6 (shock)
REAL              :: AdvVel(3)         !< Advection Velocity for the test cases
REAL              :: IniFrequency      !< parameter used for synthetic test cases
REAL              :: IniAmplitude      !< parameter used for synthetic test cases
REAL              :: IniCenter(3)      !< parameter used for Shu vortex
REAL              :: IniAxis(3)        !< parameter used for Shu vortex
REAL              :: IniHalfwidth      !< parameter used for Shu vortex
REAL              :: P_Parameter       !< parameter for Couette-Poiseuille flow
REAL              :: U_Parameter       !< parameter for Couette-Poiseuille flow
REAL              :: HarmonicFrequency !< frequency of the harmonic gauss pulse
REAL              :: AmplitudeFactor   !< amplitudefactor of the harmoinic gauss pulse
REAL              :: SiqmaSqr          !< Sigma squared of the gaussian distribution
REAL              :: JetRadius         !< parameter for Roundjet
REAL              :: JetEnd            !< parameter for Roundjet
LOGICAL           :: RoundJetInitDone=.FALSE.
REAL              :: Re_tau            !< parameter for turbulent channel
#if PARABOLIC
REAL              :: delta99_in        !< boundary layer thickness for Blasius solution
REAL              :: x_in(2)           !< inflow position for Blasius solution
LOGICAL           :: BlasiusInitDone = .FALSE. !< Flag indicating that the parameters for Blasius have been read (they can be read
                                               !< both in exact func init and in BC init)
#endif
!==================================================================================================================================

END MODULE MOD_Exactfunc_Vars
