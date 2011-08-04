/*!===============================================================================
! OLAM version 3.0_beta2

! Copyright (C) 2002-2006; All Rights Reserved; 
! Duke University, Durham, North Carolina, USA 

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

! It is requested that any scientific publications based on application of OLAM
! include the following acknowledgment:  "OLAM was developed at the 
! Edmund T. Pratt Jr. School of Engineering, Duke University."

! For additional information, including published references, please contact
! the software authors, Robert L. Walko (robert.walko@duke.edu)
! or Roni Avissar (avissar@duke.edu).
!===============================================================================*/

#if defined(CRAY) || defined(PC_NT1)

#define fh5f_open          FH5F_OPEN
#define fh5f_create        FH5F_CREATE
#define fh5f_close         FH5F_CLOSE
#define fh5d_open          FH5D_OPEN
#define fh5d_close         FH5D_CLOSE
#define fh5s_get_ndims     FH5S_GET_NDIMS
#define fh5s_get_dims      FH5S_GET_DIMS
#define fh5_prepare_read   FH5_PREPARE_READ
#define fh5d_read          FH5D_READ
#define fh5_close_read     FH5_CLOSE_READ
#define fh5_prepare_write  FH5_PREPARE_WRITE
#define fh5_write          FH5_WRITE
#define fh5_close_write    FH5_CLOSE_WRITE
#define cdc_stw            CDC_STW

#elif defined(AIX) || defined(_AIX)

#define fh5f_open_         fh5f_open
#define fh5f_create_       fh5f_create
#define fh5f_close_        fh5f_close
#define fh5d_open_         fh5d_open
#define fh5d_close_        fh5d_close
#define fh5s_get_ndims_    fh5s_get_ndims
#define fh5s_get_dims_     fh5s_get_dims
#define fh5_prepare_read_  fh5_prepare_read
#define fh5d_read_         fh5d_read
#define fh5_close_read_    fh5_close_read
#define fh5_prepare_write_ fh5_prepare_write
#define fh5_write_         fh5_write
#define fh5_close_write_   fh5_close_write
#define cdc_stw            cdc_stw

#else

#define fh5f_open          fh5f_open_
#define fh5f_create        fh5f_create_
#define fh5f_close         fh5f_close_
#define fh5d_open          fh5d_open_
#define fh5d_close         fh5d_close_
#define fh5s_get_ndims     fh5s_get_ndims_
#define fh5s_get_dims      fh5s_get_dims_
#define fh5_prepare_read   fh5_prepare_read_
#define fh5d_read          fh5d_read_
#define fh5_close_read     fh5_close_read_
#define fh5_prepare_write  fh5_prepare_write_
#define fh5_write          fh5_write_
#define fh5_close_write    fh5_close_write_
#define cdc_stw            cdc_stw_

#endif
