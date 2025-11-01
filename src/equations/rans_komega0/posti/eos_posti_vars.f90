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
!> This module contains the dependency tables for the quantities to be visualized by the posti routines
!==================================================================================================================================
MODULE MOD_EOS_Posti_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

#if PARABOLIC
#define CUT(x)
INTEGER,PARAMETER :: nVarDepEOS=41+4
#else
#define CUT(x) x!
INTEGER,PARAMETER :: nVarDepEOS=24
#endif
! ATTENTION: The first     7(previously 5) variables must be the conservative ones
!            The following 8(previously 5) variables must be the primitive ones
!           E
!           n
!           e                                                                         W
!           r                                                                         a
!           g                                                                         l
!           y                              E                        V N               l
!           S                      V       n       P                o o               F
!           t                      e     E t   T   r                r r               r W
!           a                      l     n h   o   e                t m               i a             
! W         g                      o     e a   t   s                i a               c l             t
! i         n               k   t  c V   r l   a T s                c l         W W W t l             u
! t         a               i   u  i e   g p   l o u                i i         a a a i H             r
! h         t             T n t r  t l   y y   T t r                t z         l l l o e             b
! D         i             e e u b  y o   S S   e a e          V V V y e   D Q   l l l n a             o
! G   M M M o     V V V   m t r u  M c   t t   m l T          o o o M d   i C S F F F M t             m
! O   o o o n D D e e e P p i b l  a i   a a   p P i          r r r a H   l r c r r r a T             e
! p D m m m D e e l l l r e c u e  g t   g g E e r m          t t t g e L a i h i i i g r             g
! e e e e e e n n o o o e r e l n  n y   n n n r e e          i i i n l a t t l c c c n a             a
! r n n n n n s s c c c s a n e c  i S   a a t a s D          c c c i i m a e i t t t i n       m   c _
! a s t t t s i i i i i s t e n e  t o M t t r t s e          i i i t c b t r e i i i t s       u   d r
! t i u u u i t t t t t u u r c m  u u a i i o u u r          t t t u i d i i r o o o u f       S   e e
! o t m m m t y y y y y r r g e u  d n c o o p r r i          y y y d t a o o e n n n d e x y z G f s a
! r y X Y Z y k g X Y Z e e y g t  e d h n n y e e v          X Y Z e y 2 n n n X Y Z e r + + + S d 2 l
INTEGER,DIMENSION(1:nVarDepEOS,0:nVarDepEOS),PARAMETER :: DepTableEOS = TRANSPOSE(RESHAPE(&
(/&
  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1  Density
  0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2  MomentumX
  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3  MomentumY
  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4  MomentumZ
  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !5  EnergyStagnationDensity
  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !6  DensityK
  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !7  DensityG
  0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !8  VelocityX
  0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !9  VelocityY
  0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !10 VelocityZ
  0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !11 Pressure
  0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !12 Temperature
  0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !13 turbK
  0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !14 turbOmega
  0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !15 turbMut
  0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !16 VelocityMagnitude
  0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !17 VelocitySound
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !18 Mach
  0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !19 EnergyStagnation
  0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 0,0,0,1,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !20 EnthalpyStagnation
  0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !21 Entropy
  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, 1,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !22 TotalTemperature
  0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 1,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !23 TotalPressure
  1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !24 PressureTimeDeriv
#if PARABOLIC
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !25 VorticityX
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !26 VorticityY
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !27 VorticityZ
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !28 VorticityMagnitude
  1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0, 1,0,0,0,0,0,0,0,0         ,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !29 NormalizedHelicity
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !30 Lambda2
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !31 Dilatation
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !32 QCriterion
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !33 Schlieren
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !34 WallFrictionX
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !35 WallFrictionY
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !36 WallFrictionZ
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0 ,& !37 WallFrictionMagnitude
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !38 WallHeatTransfer
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0 ,& !39 x+
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0 ,& !40 y+
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0 ,& !41 z+
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !42 muSGS
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !43 fd
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !44 CDes2
  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  & !45 turbomega_real
#endif
/),(/nVarDepEOS+1,nVarDepEOS/)))

! Mark all quantities that can be calculated exclusively on the surface
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepSurfaceOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0 &
/)

! Mark all quantities that can be calculated exclusively in the volume and must be prolonged to the surface from the volume
! NOTE(Shimushu): may have issues when dealing with muSGS etc., need to set to 1
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepVolumeOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1 &
/)

#if FV_ENABLED && FV_RECONSTRUCT
!           E
!           n
!           e
!           r
!           g
!           y
!           S
!           t
!           a
! W         g
! i         n               k   t
! t         a               i   u
! h         t             T n t r
! D         i             e e u b
! G   M M M o     V V V   m t r u
! O   o o o n D D e e e P p i b l
! p D m m m D e e l l l r e c u e
! e e e e e e n n o o o e r e l n
! r n n n n n s s c c c s a n e c
! a s t t t s i i i i i s t e n e
! t i u u u i t t t t t u u r c n
! o t m m m t y y y y y r r g e u
! r y X Y Z y k g X Y Z e e y g t
INTEGER,DIMENSION(PP_nVar,0:nVarDepEOS),PARAMETER :: DepTablePrimToCons =TRANSPOSE(RESHAPE(&
(/&
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1 Density
  0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2 MomentumX
  0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3 MomentumY
  0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4 MomentumZ
  0,1,0,0,0,0,0,0,1,1,1,1,0,1,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !5 EnergyStagnationDensity
  0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !6 RHOK
  0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  & !7 RHOG
/),(/nVarDepEOS+1,7/)))
#endif
#undef CUT

CHARACTER(LEN=255),DIMENSION(nVarDepEOS),PARAMETER :: DepNames = &
(/ CHARACTER(LEN=255) ::    &
"Density"                  ,& !1
"MomentumX"                ,& !2
"MomentumY"                ,& !3
"MomentumZ"                ,& !4
"EnergyStagnationDensity"  ,& !5
"DensityK"                 ,& !6
"DensityG"                 ,& !7
"VelocityX"                ,& !8
"VelocityY"                ,& !9
"VelocityZ"                ,& !10
"Pressure"                 ,& !11
"Temperature"              ,& !12
"turbK"                    ,& !13
"turbOmega"                ,& !14
"turbMut"                  ,& !15
"VelocityMagnitude"        ,& !16
"VelocitySound"            ,& !17
"Mach"                     ,& !18
"EnergyStagnation"         ,& !19
"EnthalpyStagnation"       ,& !20
"Entropy"                  ,& !21
"TotalTemperature"         ,& !22
"TotalPressure"            ,& !23
"PressureTimeDeriv"         & !24
#if PARABOLIC
,"VorticityX"              ,& !25
"VorticityY"               ,& !26
"VorticityZ"               ,& !27
"VorticityMagnitude"       ,& !28
"NormalizedHelicity"       ,& !29
"Lambda2"                  ,& !30
"Dilatation"               ,& !31
"QCriterion"               ,& !32
"Schlieren"                ,& !33
"WallFrictionX"            ,& !34
"WallFrictionY"            ,& !35
"WallFrictionZ"            ,& !36
"WallFrictionMagnitude"    ,& !37
"WallHeatTransfer"         ,& !38
"x+"                       ,& !39
"y+"                       ,& !40
"z+"                       ,& !41
"muSGS"                    ,& !42
"fd"                       ,& !43
"CDes2"                    ,& !44
"TurbOmega_real"            & !45
#endif /*PARABOLIC*/
/)

END MODULE MOD_EOS_Posti_Vars
