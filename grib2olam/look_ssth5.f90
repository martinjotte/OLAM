subroutine plotfield(nxll,nyll,array)
  integer, intent(inout) :: nxll, nyll
  real,    intent(inout) :: array(nxll,nyll)

  write(*,*) maxval(array)
  write(*,*) minval(array)

  call opngks
  call gks_colors(1)

   call plotback
   call tileplot(nxll,nyll,array)
   call frame

   call clsgks
 end subroutine plotfield

!****************************************************************************

subroutine tileplot(nii,njj,slab)

integer :: nii,njj
real, dimension(nii,njj) :: slab

dimension rif(6),rjf(6),ric(7),rjc(7),dst(6),ind(8)
dimension iasf(18)
!real, dimension(365) :: xmplot,ymplot
data iasf / 18*1 /
character*8 number,numbr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section adapted from OLAM colortables

integer                :: nvals  ! Number of palette values in each table
integer, dimension(50) :: ipal   ! Set of color palette indices
real, dimension(50)    :: vals   ! Set of delimiting contour values

data nvals/44/  ! Color Table (designed for SST)
data ipal(1:44)/  &
   144,143,142,100,101,102,103,104,105,106,  &
   107,108,109,110,111,112,113,114,115,116,  &
   117,118,119,120,121,122,123,124,125,126,  &
   127,128,129,130,131,132,133,134,135,136,  &
   137,138,139,140                           /
data vals(1:44)/  &
   -5.,-4.,-3.,-2.,-1., 0., 1., 2., 3., 4.,  &
    5., 6., 7., 8., 9.,10.,11.,12.,13.,14.,  &
   15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,  &
   25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,  &
   35.,36.,37.,38.                           /
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! turn off the clipping indicator.

call gsclip (0)

! set all the gks aspect source flags to "individual".

call gsasf (iasf)

! force solid fill.

call gsfais (1)

!call set(.01,.87,.01,.99,0.,361.,0.,361.,1)
call set(.01,.87,.01,.99,0.,real(nii+1),0.,real(nii+1),1)

!!do i=1,365
!!   xmplot(i)=float(i-1)
!!   ymplot(i)=float(i-1)
!!enddo

call sfseti ('type of fill',0)

do j = njj,1,-1
   do i = 1,nii
      rif(1) = real(i-1)
      rif(2) = real(i)
      rjf(1) = real(j-1)
      rjf(3) = real(j)

      rif(3) = rif(2)
      rif(4) = rif(1)
      rif(5) = rif(1)
      rjf(2) = rjf(1)
      rjf(4) = rjf(3)
      rjf(5) = rjf(1)

      ival = 1
      do while (slab(i,j) > vals(ival) .and. ival < nvals)
         ival = ival + 1
      enddo
      icolor = ipal(ival)
      call sfsgfa (rif,rjf,4,dst,6,ind,8,icolor)
   enddo
enddo

!  Draw grid lines

call sflush
call gsplci(10)
call sflush
call gstxci(10)
call sflush

!!do i = 2,362,10
!!   call frstpt(xmplot(i),ymplot(1))
!!   call vector(xmplot(i),ymplot(183))
!!enddo
!!
!!do j = 2,182,10
!!   call frstpt(xmplot(1),ymplot(j))
!!   call vector(xmplot(363),ymplot(j))
!!enddo

call sflush

! draw a color bar for the plot.

call set (0.,1.,0.,1.,0.,1.,0.,1.,1)

call gsplci(8)
call gstxci(8)

rif(1) = .88
rif(2) = .91
rif(3) = .91
rif(4) = .88
rif(5) = rif(1)

bsize = -.63

yinc = .98 / float(nvals)

do ibox = 1,nvals

   rjf(1) = .01 + float(ibox-1) * yinc
   rjf(2) = rjf(1)
   rjf(3) = rjf(2) + yinc
   rjf(4) = rjf(3)
   rjf(5) = rjf(1)

   write (number,'(f7.0)') vals(ibox)

   call deblank(number,numbr,nnn)
   if(numbr(len_trim(numbr):len_trim(numbr)) == '.') then
      ln = len_trim(numbr) - 1
   else
      ln = len_trim(numbr)
   endif

   if (ibox < nvals) then
      call plchlq (rif(2)+.01,rjf(3),numbr(1:ln),bsize,0.,-1.)
   endif

   call sfsgfa (rif,rjf,4,dst,6,ind,8,ipal(ibox))
!   call fillpolyg (4,rif,rjf,ipal(ibox))

enddo

return
end

!*****************************************************************************

subroutine gks_colors(iwk)

implicit none
integer :: iwk

call gscr(iwk,1,1.,0.,0.)            ! redll gscr(iwk,2,1.,1.,0.)            ! yellow
call gscr(iwk,3,0.,1.,0.)            ! green
call gscr(iwk,4,0.,1.,1.)            ! cyan
call gscr(iwk,5,.5,1.,0.)            ! yellow-green
call gscr(iwk,6,0.,0.,1.)            ! blue
call gscr(iwk,7,1.,1.,1.)            ! white
call gscr(iwk,8,0.00 , 0.00 , 0.50)  ! dark blue
call gscr(iwk,9,0.00 , 0.50 , 0.00)  ! dark green
call gscr(iwk,10,0.00 , 0.00 , 0.00) ! black
call gscr(iwk,11,1.,0.,1.)           ! purple
call gscr(iwk,12,.5,.1,.1)           ! dark red for roads
call gscr(iwk,13,.6,.6,.6)           ! gray for roads

! set of 50 colors

call hls00 (iwk, 97,145.,48., 90.)
call hls00 (iwk, 98,145.,61., 90.)  ! red oranges (for below sea level)
call hls00 (iwk, 99,145.,74., 90.)

call hls00 (iwk,100,340.,92., 75.)  ! special sea level value
call hls00 (iwk,101,340.,83., 75.)  ! blues
call hls00 (iwk,102,340.,71., 75.)
call hls00 (iwk,103,340.,60., 75.)

call hls00 (iwk,104,302.,77., 70.)
call hls00 (iwk,105,302.,55., 70.)  ! blue greens
call hls00 (iwk,106,302.,41., 70.)
call hls00 (iwk,107,302.,33., 70.)

call hls00 (iwk,108,270.,77., 65.)
call hls00 (iwk,109,270.,55., 65.)  ! greens
call hls00 (iwk,110,270.,41., 65.)
call hls00 (iwk,111,270.,33., 65.)

call hls00 (iwk,112,194.,70., 75.)
call hls00 (iwk,113,194.,50., 75.)  ! green yellows
call hls00 (iwk,114,194.,41., 75.)
call hls00 (iwk,115,194.,33., 75.)

call hls00 (iwk,116,180.,62., 90.)  ! yellow

call hls00 (iwk,117,165.,74.,100.)
call hls00 (iwk,118,165.,61.,100.)  ! oranges
call hls00 (iwk,119,165.,48.,100.)
call hls00 (iwk,120,165.,42.,100.)

call hls00 (iwk,121,150.,68., 50.)
call hls00 (iwk,122,150.,57., 50.)  ! browns
call hls00 (iwk,123,150.,44., 50.)
call hls00 (iwk,124,150.,37., 50.)

call hls00 (iwk,125,135.,80., 70.)
call hls00 (iwk,126,135.,68., 70.)  ! reds
call hls00 (iwk,127,135.,49., 70.)
call hls00 (iwk,128,135.,40., 70.)

call hls00 (iwk,129, 85.,80., 70.)
call hls00 (iwk,130, 85.,68., 70.)  ! magentas
call hls00 (iwk,131, 85.,49., 70.)
call hls00 (iwk,132, 85.,40., 70.)

call hls00 (iwk,133, 50.,80., 70.)
call hls00 (iwk,134, 50.,68., 70.)  ! violets
call hls00 (iwk,135, 50.,49., 70.)
call hls00 (iwk,136, 50.,40., 70.)

call hls00 (iwk,137, 20.,80., 70.)
call hls00 (iwk,138, 20.,70., 70.)  ! violet blues
call hls00 (iwk,139, 20.,56., 70.)
call hls00 (iwk,140, 20.,40., 70.)

call hls00 (iwk,141, 00., 90.,  0.)
call hls00 (iwk,142, 00., 77.,  0.)
call hls00 (iwk,143, 00., 60.,  0.)
call hls00 (iwk,144, 00., 43.,  0.)  ! gray shades
call hls00 (iwk,145, 00., 34.,  0.)

call hls00 (iwk,146, 00.,100.,  0.)  ! white

return
end

!***************************************************************************

subroutine hls00 (iwk,ic,h,l,s)

implicit none

integer :: iwk,ic
real :: h,l,s

real :: r,g,b

call hlsrgb (h,l,s,r,g,b)
call gscr (iwk,ic,r,g,b)

return
end

!******************************************************************

subroutine plotback

implicit none

real, dimension(5) :: rif,rjf
real, dimension(8) :: dst,ind
integer, dimension(18) :: iasf
data iasf / 18*1 /

!  This subroutine plots a white background for the entire plot.

call gsclip (0)

! Set all the GKS aspect source flags to "individual".

call gsasf (iasf)

! Force solid fill.

call gsfais (1)

call set (0.,1.,0.,1.,0.,1.,0.,1.,1)

call sfseti ('TYPE OF FILL',0)

rif(1) = 0.
rif(2) = 1.
rif(3) = 1.
rif(4) = 0.
rif(5) = 0.

rjf(1) = 0.
rjf(2) = 0.
rjf(3) = 1.
rjf(4) = 1.
rjf(5) = 0.

! The last parameter in the following call sets the color.

call sfsgfa (rif,rjf,4,dst,5,ind,8,7)

return
end

!***************************************************************************

subroutine deblank(str1,str2,nch)
implicit none
character(len=*) :: str1,str2
integer :: n,ln,nch

! strips blanks from a string and returns number of chars

str2=' '
ln=len(str1)
nch=0
do n=1,ln
   if(str1(n:n).ne.' ') then
      nch=nch+1
      str2(nch:nch)=str1(n:n)
   endif
enddo

return
end


