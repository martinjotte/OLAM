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
!==================================================================================*/

#include <stdlib.h>
#include <stdio.h>

static FILE *pipe;

void ir_popen (int *nbuff, char *buff, char *cmd, int *lenout) {

  int nch;

/*printf(" popen cmd: %s \n",cmd);*/

  *lenout=-1;

  if ((pipe = popen(cmd, "r")) != NULL) {

    nch =fread(buff, 1, *nbuff, pipe);

  /*printf(" popen nch: %d \n",nch);*/

    if ( nch > *nbuff ) {
      fprintf(stderr,"Maximum ir_popen buffer length exceeded: %d\n",nch);
      fprintf(stderr,"popen command: %s\n",cmd);
      fprintf(stderr,"Increase buffer size in ir_popen call\n");
      exit(1);
    }

    *lenout=nch;

    pclose(pipe);
  }
}
