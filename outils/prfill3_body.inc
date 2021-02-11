integer :: i,j,k,nprxo2,iin

nprxo2 = nprx / 2

! Copy pressure-level or surface data from xx to dn array, shifting by ipoffset
! in x direction and by 2 in y direction

!$omp parallel do private(j,i,iin)
do k = 1,nprz
   do j = 1,npry
      do i = 1,nprx+4
         iin = mod(i+nprx-ipoffset-1,nprx)+1
         dn(i,j+2,k) = xx(iin,j,k)
      enddo
   enddo

! Fill in N/S boundary rows, based on whether xswlat = -90. or -90. + gdatdy/2.

   if (inproj==1 .and. abs(xswlat + 90.) < .1) then

      do i = 1,nprxo2
         dn(i,1,k)      = dn(i+nprxo2,5,k)
         dn(i,2,k)      = dn(i+nprxo2,4,k)
         dn(i,npry+3,k) = dn(i+nprxo2,npry+1,k)
         dn(i,npry+4,k) = dn(i+nprxo2,npry,k)
      enddo

      do i = nprxo2+1,nprx+4
         dn(i,1,k)      = dn(i-nprxo2,5,k)
         dn(i,2,k)      = dn(i-nprxo2,4,k)
         dn(i,npry+3,k) = dn(i-nprxo2,npry+1,k)
         dn(i,npry+4,k) = dn(i-nprxo2,npry,k)
      enddo

   elseif (inproj==2 .or. (inproj==1 .and. abs(xswlat - .5 * gdatdy + 90.) < .1)) then

      do i = 1,nprxo2
         dn(i,1,k)      = dn(i+nprxo2,4,k)
         dn(i,2,k)      = dn(i+nprxo2,3,k)
         dn(i,npry+3,k) = dn(i+nprxo2,npry+2,k)
         dn(i,npry+4,k) = dn(i+nprxo2,npry+1,k)
      enddo

      do i = nprxo2+1,nprx+4
         dn(i,1,k)      = dn(i-nprxo2,4,k)
         dn(i,2,k)      = dn(i-nprxo2,3,k)
         dn(i,npry+3,k) = dn(i-nprxo2,npry+2,k)
         dn(i,npry+4,k) = dn(i-nprxo2,npry+1,k)
      enddo

   else
      print*, 'xswlat & gdatdy combo not found ',xswlat,gdatdy
   endif
   
enddo
!$omp end parallel do
