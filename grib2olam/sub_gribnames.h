/*!
! Copyright (C) 1991-2003  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
*/
#if defined(IBM) || defined(HP)

#define ir_popen ir_popen
#define grib_get_rec grib_get_rec

#elif defined(CRAY)

#define ir_popen IR_POPEN
#define grib_get_rec GRIB_GET_REC

#else 

#define ir_popen ir_popen_
#define grib_get_rec grib_get_rec_

#endif
