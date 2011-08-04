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
Module consts_coms

! Portable single (real*4) and double (real*8) precision kinds:
integer, parameter :: r4 = selected_real_kind(6,37)
integer, parameter :: r8 = selected_real_kind(13,300)

real, parameter :: rdry     = 287.               ! dry air gas constant [J/(kg K)]
real, parameter :: cp       = 1004.              ! dry air spec heat at const P [J/(kg K)]
real, parameter :: cv       = 717.               ! dry air spec heat at const vol [J/(kg K)]
real, parameter :: rvap     = 461.               ! water vapor gas constant [J/(kg K)]
real, parameter :: p00      = 1.e5               ! reference pressure [Pa]
real, parameter :: t00      = 273.15             ! 0 Celsius temperature [K]
real, parameter :: grav     = 9.81               ! acceleration of gravity [m/s^2]
real, parameter :: pio180   = 3.1415927 / 180.   ! radians per degree
real, parameter :: piu180   = 180. / 3.1415927   ! degrees per radian
real, parameter :: pi1      = 3.1415927           
real, parameter :: pi2      = 3.1415927 * 2.      
real, parameter :: pi4      = 3.1415927 * 4.      
real, parameter :: pio2     = 3.1415927 * .5      
real, parameter :: pio4     = 3.1415927 * .25     
real, parameter :: erad     = 6367.e3            ! Earth radius [m]
real, parameter :: erad2    = erad * 2.            
real, parameter :: erad4    = erad * 4.           
real, parameter :: eradsq   = erad * erad         
real, parameter :: erador5  = erad * .4472135955 ! second factor is sqrt(.2)
real, parameter :: eradi    = 1. / erad           
real, parameter :: vonk     = 0.40               ! von Karman constant
real, parameter :: tkmin    = 5.e-4              ! minimum allowed TKE
real, parameter :: alvl     = 2.50e6             ! latent heat of evaporation [J/kg]
real, parameter :: alvi     = 2.834e6            ! latent heat of sublimation [J/kg]
real, parameter :: alli     = 0.334e6            ! latent heat of fusion [J/kg]
real, parameter :: alli1000 = 0.334e9            ! latent heat of fusion [J/(1000 kg)]
real, parameter :: alvl2    = 6.25e12            ! alvl * alvl
real, parameter :: alvi2    = 8.032e12           ! alvi * alvi
real, parameter :: cliq     = 4186.              ! spec heat of liquid [J/(kg K)]
real, parameter :: cice     = 2106.              ! spec heat of ice [J/(kg K)]
real, parameter :: cliq1000 = cliq * 1000.       ! spec heat of liquid [J/(1000 kg K)]
real, parameter :: cice1000 = cice * 1000.       ! spec heat of ice [J/(1000 kg K)]
real, parameter :: solar    = 1.3533e3           ! solar constant [W/m^2]
real, parameter :: stefan   = 5.6696e-8          ! Stefan-Boltzmann constant
real, parameter :: dlat     = 111120.            ! Earth 1-latitude-degree dist [m]
real, parameter :: omega    = 7.29212e-5           ! Earth rotational angular velocity
real, parameter :: omega2   = 7.29212e-5 * 2.      ! omega * 2
real, parameter :: p00i     = 1. / p00            
real, parameter :: rocp     = rdry / cp           
real, parameter :: cpor     = cp / rdry           
real, parameter :: rocv     = rdry / cv           
real, parameter :: cvor     = cv / rdry           
real, parameter :: cpocv    = cp / cv             
real, parameter :: cvocp    = cv / cp             
real, parameter :: cpi      = 1. / cp             
real, parameter :: cpi4     = 4. * cpi            
real, parameter :: cp253i   = cpi / 253.           
real, parameter :: cliqi    = 1. / cliq             
real, parameter :: cicei    = 1. / cice             
real, parameter :: allii    = 1. / alli           
real, parameter :: alvlor   = alvl / rdry         
real, parameter :: alvlocp  = alvl / cp           
real, parameter :: alviocp  = alvi / cp           
real, parameter :: grav2    = 2. * grav           
real, parameter :: gravo2   = .5 * grav           
real, parameter :: gravi    = 1.0 / grav           
real, parameter :: cpog     = cp / grav
real, parameter :: gocp     = grav / cp
real, parameter :: gordry   = grav / rdry
real, parameter :: rdryog   = rdry / grav
real, parameter :: eps_vap  = rdry / rvap          
real, parameter :: eps_virt = (rvap - rdry) / rdry 
real, parameter :: p00k     = 26.870941          ! = p00 ** rocp  
real, parameter :: p00ki    = 1. / p00k           
real, parameter :: pc1      = .9967937e-2        ! pc1 = (1./p00) ** rocv
                                      ! old form:  pc1 = (1./p00) ** rocp

End Module
