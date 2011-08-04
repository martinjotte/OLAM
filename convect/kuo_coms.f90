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
Module kuo_coms

integer, parameter :: nkp=100

integer :: kcon, klcl, klfc, ketl, kct, igo, kmt, icprtfl, icpltfl
               
real    :: zmid, cdzmin, dzlow, dzhigh, plcl, tlcl, dzlcl, zlcl, garea
real    :: wconmin, contim, preff, envshr, supply, cprecip

real, dimension(nkp) :: uepe,vepe,wepe,wpe,ze,te,the,pe,rte,pke,rhoe,thve,zc, &
                        rve,thee,qvct1,qvct2,qvct3,qvct4,  &
                        vheat,vmois,vmdry, &
                        theu,rsu,thu,tu,thd,wtd, ftcone, frcone, thcone

End Module kuo_coms
