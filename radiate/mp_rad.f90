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
SUBROUTINE SHRADP(NZP,RVR,DN0R,DZR,PIRD,COSZ,ALBEDO  &
                 ,SOLAR,FTHR,RSHORT)
implicit none                 
!+----------------------------------------------------------------------
!     Shortwave radiation parameterization described in Mahrer and
!     Pielke(1977).
!
!       Arguments:
!       ----------
!
!       Input:  NZP    - number of vertical levels
!               RVR    - water vapor at each level
!               DN0R   - air density at each level
!               DZR    - inverse of delta z = 1./(Z(K)-Z(K-1))
!                        where Z is at the vapor levels
!               PIRD   - Exner function (p/p00)**(R/Cp) at each level
!               SC     - scratch array at least 2*NZP long
!               COSZ   - cosine of the zenith angle
!               ALBEDO - albedo of the ground surface
!               SOLAR  - the solar constant
!
!       Output: FTHR   - radiation tendency on potential temperture.
!               RSHORT - downward shortwave flux on a flat surface at
!                        the ground
!
!+----------------------------------------------------------------------
integer :: nzp
real :: COSZ,ALBEDO,SOLAR,RSHORT
real :: RVR(NZP),DN0R(NZP),PIRD(NZP),FTHR(NZP),DZR(NZP)

real :: SC1(NZP),sc2(nzp)  ! automatic arrays

real, parameter :: CP=1.004E7

integer :: nz,k
real :: raysct,rdcon1,vabs

NZ=NZP-1

!     Rayleigh scattering (numerator in SQRT should be
!        SQRT((.000949*P+.051)/COSZ), but is ignored. (P in mb)

RAYSCT=1.021-.0824*SQRT(1./COSZ)

!     Vapor path length

DO K=1,NZ
  SC1(K)=(RVR(K)*DN0R(K)+RVR(K+1)*DN0R(K+1))*.5/DZR(K)
  SC1(K)=MAX(SC1(K),1E-10)
ENDDO

DO K=1,NZ
  sc2(k) = sum(sc1(k:nz)) + 1.e-10
ENDDO
SC2(NZP)=0.

!     Shortwave heating by vapor absorbtion

RDCON1=.0231*SOLAR/CP
DO K=2,NZ
  FTHR(K)=(RDCON1*(SC2(K)/COSZ)**(-.7)*COSZ*RVR(K))/PIRD(K)
ENDDO
VABS=.077*(SC2(1)/COSZ)**.3

!     Shortwave on a flat surface

RSHORT=MAX(SOLAR*COSZ*(1.-ALBEDO)*(RAYSCT-VABS),0.)

RETURN
END

!===============================================================================

SUBROUTINE LWRADP(NZP,TEMPRD,RVR,DN0R,DZZR,PIRD,FTHR,RLONG)
implicit none
!-----------------------------------------------------------------------
!     Longwave radiation parameterization described in Mahrer and
!     Pielke (1977).  Does not include any cloud effects.
!
!       Arguments:
!       ----------
!
!       Input:  NZP    - number of vertical levels
!               TEMPRD - temperature in Kelvin at each level
!               RVR    - water vapor at each level
!               DN0R   - air density at each level
!               DZZR   - inverse of delta z = 1./(ZZ(K)-ZZ(K-1))
!                        where ZZ is staggered midway between T levels
!               PIRD   - Exner function (p/p00)**(R/Cp) at each level
!               SC     - scratch array at least 20*NZP long
!
!       Output: FTHR   - radiation tendency on potential temperture.
!               RLONG  - downward longwave flux at the ground
!
!-----------------------------------------------------------------------
integer :: nzp
real :: RLONG
real :: RVR(NZP),DN0R(NZP),TEMPRD(NZP),FTHR(NZP),DZZR(NZP),PIRD(NZP)
real sc(nzp,18)  ! automatic array

real, parameter :: G=980.,CP=1.004E7,STEFAN=5.6696E-5,R=.287E7,P00=1E6
integer, parameter :: IV1=1,IV2=2,IV3=3,IV4=4,IV5=5,IV6=6,IV7=7,IV8=8  &
                     ,IV9=9,IV10=10,IV11=11,IV12=12,IV13=13,IV14=14  &
                     ,IV15=15,IV16=16,IV17=17,IV18=18
integer :: nz,nz1,k
real :: c1,c2                     

sc(1:nzp,1:18) = 0.

NZ=NZP-1
NZ1=NZ-1

!                   COMPUTE UPWARD AND DOWNWARD VAPOR PATH

DO K=2,NZ
  SC(K,IV1)=RVR(K)*DN0R(K)/DZZR(K)
  SC(K,IV1)=MAX(SC(K,IV1),1E-10)
ENDDO
SC(NZP,IV1)=SC(NZ,IV1)

SC(1,IV2)=0.
DO K=2,NZ
  SC(K,IV2)=SC(K-1,IV2)+SC(K,IV1)
ENDDO
SC(NZ,IV3)=SC(NZP,IV1)
DO K=NZ1,1,-1
  SC(K,IV3)=SC(K+1,IV3)+SC(K+1,IV1)
ENDDO

!                          WATER VAPOR EMISSIVITY CALCULATION

DO K=1,NZ
  SC(K,IV4)=LOG10(SC(K,IV2)+1E-30)
  SC(K,IV5)=LOG10(SC(K,IV3)+1E-30)
ENDDO

DO K=1,NZ
  IF(SC(K,IV4).LE.-4.)  &
    SC(K,IV6)=.1129*LOG10(1.+12.63*SC(K,IV2))
  IF(SC(K,IV4).LE.-3.0.AND.SC(K,IV4).GT.-4.0)  &
    SC(K,IV6)=.104*SC(K,IV4)+.440
  IF(SC(K,IV4).LE.-1.5.AND.SC(K,IV4).GT.-3.0)  &
    SC(K,IV6)=.121*SC(K,IV4)+.491
  IF(SC(K,IV4).LE.-1.0.AND.SC(K,IV4).GT.-1.5)  &
    SC(K,IV6)=.146*SC(K,IV4)+.527
  IF(SC(K,IV4).LE. 0.0.AND.SC(K,IV4).GT.-1.0)  &
    SC(K,IV6)=.161*SC(K,IV4)+.542
  IF(SC(K,IV4).GT. 0.)  &
    SC(K,IV6)=.136*SC(K,IV4)+.542

  IF(SC(K,IV5).LE.-4.)  &
    SC(K,IV7)=.1129*LOG10(1.+12.63*SC(K,IV3))
  IF(SC(K,IV5).LE.-3.0.AND.SC(K,IV5).GT.-4.0)  &
    SC(K,IV7)=.104*SC(K,IV5)+.440
  IF(SC(K,IV5).LE.-1.5.AND.SC(K,IV5).GT.-3.0)  &
    SC(K,IV7)=.121*SC(K,IV5)+.491
  IF(SC(K,IV5).LE.-1.0.AND.SC(K,IV5).GT.-1.5)  &
    SC(K,IV7)=.146*SC(K,IV5)+.527
  IF(SC(K,IV5).LE. 0.0.AND.SC(K,IV5).GT.-1.0)  &
    SC(K,IV7)=.161*SC(K,IV5)+.542
  IF(SC(K,IV5).GT. 0.)  &
    SC(K,IV7)=.136*SC(K,IV5)+.542
ENDDO

!                           CO2 path lengths and emissivities

C1=.0004148239
C2=C1*G
DO K=2,NZ
  SC(K,IV11)=C2*DN0R(K)/DZZR(K)
ENDDO
SC(NZP,IV11)=C1*PIRD(NZP)**(CP/R)*P00

SC(1,IV12)=0.
DO K=2,NZ
  SC(K,IV12)=SC(K-1,IV12)+SC(K,IV11)
ENDDO
SC(NZ,IV13)=SC(NZP,IV11)
DO K=NZ1,1,-1
  SC(K,IV13)=SC(K+1,IV13)+SC(K+1,IV11)
ENDDO

DO K=1,NZ
  SC(K,IV8)=.185*(1.-EXP(-.3919*SC(K,IV12)**.4))
  SC(K,IV9)=.185*(1.-EXP(-.3919*SC(K,IV13)**.4))
ENDDO

!                        Add CO2 and H2O emissivities, find SIG(T**4)

DO K=1,NZP
  SC(K,IV14)=SC(K,IV8)+SC(K,IV6)
  SC(K,IV15)=SC(K,IV9)+SC(K,IV7)
  SC(K,IV16)=STEFAN*TEMPRD(K)**4
ENDDO
SC(NZP,IV16)=SC(NZP,IV16)*SC(NZ,IV15)

!                       Calculate upward and downward divergences

DO K=2,NZ
  SC(K,IV17)=(SC(K,IV16)-SC(1,IV16))*(SC(K,IV14)-SC(K-1,IV14))
  SC(K,IV18)=(SC(NZP,IV16)-SC(K,IV16))*(SC(K,IV15)-SC(K-1,IV15))
ENDDO
SC(1,IV18)=SC(NZP,IV16)*(1.-SC(1,IV15))+SC(2,IV16)*SC(2,IV15)

DO K=2,NZ
  FTHR(K)=-(SC(K,IV17)+SC(K,IV18))*DZZR(K)/(CP*DN0R(K)*PIRD(K))
ENDDO
RLONG=SC(1,IV18)

!------------------------------------------------------------------
!      PRINT 6667,(K,TEMPRD(K),SC(K,IV6),SC(K,IV7),SC(K,IV8),SC(K,IV9)
!     +  ,SC(K,IV13),FTHR(K)*24.*3600.,K=NZP,1,-1)
6667 FORMAT(' LONGWAVE-T,EMISSUR,EMISSDR,EMISSUC,EMISSDC,PATHD,D/DAY',  &
      /,(I3,7E10.3))
!      PRINT 6668,(K,TEMPRD(K),SC(K,IV2),SC(K,IV3),SC(K,IV12),SC(K,IV13)
!     +  ,K=NZP,1,-1)
6668 FORMAT(' LONGWAVE-T,UP VAP,DN VAP,UP CO2, DN CO2',  &
      /,(I3,5E10.3))
!      PRINT*,'  LONGWAVE DOWN ',-SC(1,IV18)*1E-3

RETURN
END

