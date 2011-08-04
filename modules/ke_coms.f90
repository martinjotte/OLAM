!===============================================================================
! OLAM version 4.0

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
Module ke_coms

!STC---------------------------------------------------------------------------
!STC---------------------------------------------------------------------------
!       Setting of the empirical constants and parameters
!       to be used for the E-l and E-eps turbulence closures
! (S. Trini Castelli)
!                       TABLE OF THE VALUES
!
! C_MI, C_EPS           :   empirical constants for E-l closure 
! C_MI, C1_EPS, C2_EPS  :   empirical constants for E-eps closure 
! ALF_TKE, ALF_EPS      :   coefficients for the proportionality of 
!                           diffusivities Km/Ke and Km/Keps 
! ALF_THT               :   Prandtl-Schimdt turb. num. for temperature, 
!                           coefficient for the proportionality of i
!                           diffusivity Km/Kth
! C_EPS_Y,C_ust_tk      :   Constants for the bottom boundary condition of 
!                           TKE and EPS in E-l and E-eps closures
! IOPZL                 :   option flag for the aymptotic mixing length 
!                          = 1  constant value (i.e., deriv. from Ying, 1992)
!                          = 2  Mellor-Yamada 1974 formulation
!                          = 3  Zilitinkevich-Laitkman 1965 formulation
! AL0_CONST             :   Constant value for the asymptotic mixing length
!                          (IOPZ = 1 ) (i.e., Ying's value=0.336 m, 1992)
!STC---------------------------------------------------------------------------
real, parameter ::                    &
        C_EPS    = 0.08               &   ! RUSVAL value
    ,   C_MI     = 0.42               &   ! RUSVAL value
!
!_alternative    ,   C_EPS =0.17               &   ! Ying's value
!_alternative    ,   C_MI  =0.55               &   ! Ying's value
!_alternative    ,   C_MI  =0.55               &   ! DET-ETL (1984) LAB
!_alternative    ,   C_MI  =0.40               &   ! DET-ETL (1984) ABL
!_alternative    ,   C_MI  =0.43               &   ! DUYNKERKE (1988)  
!_alternative    ,   C_MI  =0.41               &   ! TAYLOR (1997)     
!
    ,   C1_EPS   = 1.22               &   !  RUSVAL value
    ,   C2_EPS   = 1.92               &   !  used in RUSVAL 
!
!_alternative    ,   C1_EPS   = 1.44               &   !  DET - ETL (1984) LAB
!_alternative    ,   C2_EPS   = 1.92               &   !  DET - ETL (1984) LAB
!_alternative    ,   C1_EPS   = 1.13               &   !  DET - ETL (1984) ABL
!_alternative    ,   C2_EPS   = 1.90               &   !  DET - ETL (1984) ABL
!_alternative    ,   C1_EPS   = 1.46               &   !  DUYNKERKE (1988)    
!_alternative    ,   C2_EPS   = 1.83               &   !  DUYNKERKE (1988)    
!_alternative    ,   C1_EPS   = 1.44               &   !  TAYLOR (1997)       
!_alternative    ,   C2_EPS   = 1.92               &   !  TAYLOR (1997)       
!_alternative    ,   C1_EPS   = 1.22               &   ! RUSVAL value
!_alternative    ,   C2_EPS   = 1.92               &   ! used in RUSVAL 
!
    ,   ALF_TKE  = 1.                 &   ! used in RUSVAL (DET-ETL 1984, LAB)
    ,   ALF_EPS  = 0.77                  &   ! used in RUSVAL (DET-ETL 1984, LAB)
!
!_alternative    ,   ALF_TKE  = 1.35               &   ! DET - ETL (1984) ABL
!_alternative    ,   ALF_EPS  = 0.77               &   ! DET - ETL (1984) ABL
!_alternative    ,   ALF_TKE  = 1.                 &   ! DUYNKERKE (1988)
!_alternative    ,   ALF_EPS  = 0.42               &   ! DUYNKERKE (1988)
!_alternative    ,   ALF_TKE  = 1.                 &   ! TAYLOR (1997) 
!_alternative    ,   ALF_EPS  = 0.51               &   ! TAYLOR (1997) 
!
    ,   ALF_THT = 1.11 
!
real, parameter ::  C_ust_tk  = 0.16      ! Panofsky 1977 and DET-ETL (1984) ABL
!alternative   ,    C_ust_tk  = 0.3                    ! DET-ETL (1984) LAB
!
integer, parameter :: IOPZL  = 1
!!!real, parameter :: AL0_CONST = 201.6    ! NBL = 600 m
real, parameter :: AL0_CONST = 100.8    ! NBL = 300 m

real, parameter :: coef_km     = c_mi*c_eps

End Module
