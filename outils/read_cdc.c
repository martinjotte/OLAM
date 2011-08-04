/*  SUMMARY: 
      1. This programs reads the global NCEP reanalysis files (in netcdf 
      format) that you have downloaded.  
      2. It then prints out the data for a subregion over a specified time 
      period in ASCII. 
*/


#define  _GNU_SOURCE		/* for NaN */
#include "utils_sub_names.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <netcdf.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <hdf5.h>
#define NR_END 1
#define LATMIN -90
#define LATMAX 90
#define LONMIN 0
#define LONMAX 360
#define TIMESTEP 6

#define check_ncstat(stat, msg) _check_ncstat(stat, __FILE__, __LINE__, msg);
void check_cdc_files(char *infilename, int iyear);
void _check_ncstat(int stat, const char *filename, unsigned line,
		   const char *msg);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
/**********************************************************/

struct ncvar {
  char *filename;	/* file data came from */
  int ncid;			/* file id */
  int varid;  			/* variable id */
  float minlat, maxlat;		/* latitude range */
  float latstep;
  float minlon, maxlon;		/* longitude range */
  size_t nlat, nlon;
  float lonstep;
  long int mintime, maxtime;
  size_t N;			/* number of timesteps */
  float offset, scale;		/* for decoding the variable */
  short nodata;			/* code for missing data */
};
struct ncvar * new_ncvar(const char *filename, const char *varname);
void _check_ncstat(int stat, const char *filename, unsigned line,const char *msg );
float *** ncgetvar(const struct ncvar *v);
void free_ncvar(struct ncvar *var);
void ncnearest(struct ncvar *v, float lat, float lon,
	       size_t *lati, size_t *loni);
void index2latlon(struct ncvar *v, size_t lati, size_t loni,
		  float *lat, float *lon);
#define FREE_ARG char*
void free_f3tensor(float ***t, const struct ncvar *v);

/****************************************************************/

void cdc_stw(int *iyear, int *imonth, int *idate, int *itime, char *locfn, float *soilt1, float *soilt2, float *soilw1, float *soilw2, float *snow){
  FILE*outfile;
  struct ncvar *soilw010cm_var = NULL;
  struct ncvar *soilw10200cm_var = NULL;
  struct ncvar *tmp010cm_var = NULL;
  struct ncvar *tmp10200cm_var = NULL;
  struct ncvar *weasd_var = NULL;
  struct ncvar *tmpvar = NULL;
  float ***soilw010cm = NULL;
  float ***soilw10200cm = NULL;
  float ***tmp010cm = NULL;
  float ***tmp10200cm = NULL;
  float ***weasd = NULL;
  char infilename[256];
  int lati,loni;
  size_t latmini,lonmini,latmaxi,lonmaxi,rlati,rloni;
  float lat,lon;
  double time;
  unsigned int timei;
  float u;
  int beg_time[13];
  int nlat,nlon;
  float lat_min,lon_min,lat_max,lon_max,rv;
  float stime;
  int offset_from_utc,id,ib,ih;
  int ntime;
  char tmpvar_name[256];
  char mname[3];
  hid_t fileid,dsetid,mspcid,propid;
  hid_t fileid2,dsetid2,mspcid2,propid2;
  char syscmd[256];
  float var_out[1464];
  int ndims;
  int idims[1];
  hsize_t dimsc[7] = {1,1,1,1,1,1,1};
  hsize_t maxdims[7] = {1,1,1,1,1,1,1};
  hsize_t chunk_size[7] = {0,0,0,0,0,0,0};
  int i;
  char var_name[20];
  int hdferr;
  int timeo;

  /* Read in each variable */
  printf("trying tmp.0-10cm\n");
  sprintf(infilename,"%stmp.0-10cm.gauss.%d.nc",locfn,*iyear);
  check_cdc_files(infilename, *iyear);
  tmp010cm_var = new_ncvar(infilename,"tmp");
  tmp010cm = ncgetvar(tmp010cm_var);
  
  printf("trying tmp.10-200cm\n");
  sprintf(infilename,"%stmp.10-200cm.gauss.%d.nc",locfn,*iyear);
  check_cdc_files(infilename, *iyear);
  tmp10200cm_var = new_ncvar(infilename,"tmp");
  tmp10200cm = ncgetvar(tmp10200cm_var);
  
  printf("trying soilw.0-10cm\n");
  sprintf(infilename,"%ssoilw.0-10cm.gauss.%d.nc",locfn,*iyear);
  check_cdc_files(infilename, *iyear);
  soilw010cm_var = new_ncvar(infilename,"soilw");
  soilw010cm = ncgetvar(soilw010cm_var);
  
  printf("trying soilw.10-200cm\n");
  sprintf(infilename,"%ssoilw.10-200cm.gauss.%d.nc",locfn,*iyear);
  check_cdc_files(infilename, *iyear);
  soilw10200cm_var = new_ncvar(infilename,"soilw");
  soilw10200cm = ncgetvar(soilw10200cm_var);
  
  printf("trying weasd\n");
  sprintf(infilename,"%sweasd.sfc.gauss.%d.nc",locfn,*iyear);
  check_cdc_files(infilename, *iyear);
  weasd_var = new_ncvar(infilename,"weasd");
  weasd = ncgetvar(weasd_var);

  /*   The variable beg_time is used to determine what month you are in.
   */
  beg_time[0] = 0;
  beg_time[1] = 31*4;
  if(*iyear%4 == 0){
      beg_time[2] = beg_time[1] + 29 * 4;
  }else{
      beg_time[2] = beg_time[1] + 28 * 4;
  }
  beg_time[3] = beg_time[2] + 31 * 4;
  beg_time[4] = beg_time[3] + 30 * 4;
  beg_time[5] = beg_time[4] + 31 * 4;
  beg_time[6] = beg_time[5] + 30 * 4;
  beg_time[7] = beg_time[6] + 31 * 4;
  beg_time[8] = beg_time[7] + 31 * 4;
  beg_time[9] = beg_time[8] + 30 * 4;
  beg_time[10] = beg_time[9] + 31 * 4;
  beg_time[11] = beg_time[10] + 30 * 4;
  beg_time[12] = beg_time[11] + 31 * 4;
  beg_time[*imonth-1] += 4 * (*idate - 1) + *itime/6; 

  /* loop over sites */
  for(lati = 1; lati <= 73; lati++) {
      lat = 90.0 - (lati - 1) * 2.5;
      // Goes from east to west
      for(loni = 1; loni <= 144; loni++) {
	  lon = 1.25 + (loni - 1) * 2.5;
	  ncnearest(tmp010cm_var, lat, lon, &rlati, &rloni);
	  
	  timei=beg_time[*imonth-1];
	  soilt1[144*(lati-1)+loni-1] = tmp010cm[timei][rlati][rloni];
	  soilt2[144*(lati-1)+loni-1] = tmp10200cm[timei][rlati][rloni];
	  soilw1[144*(lati-1)+loni-1] = soilw010cm[timei][rlati][rloni]/0.434;
	  soilw2[144*(lati-1)+loni-1] = soilw10200cm[timei][rlati][rloni]/0.434;
	  snow[144*(lati-1)+loni-1] = weasd[timei/4][rlati][rloni];
	  
      }
  }

  /* Free memory */

  free_f3tensor(tmp010cm, tmp010cm_var);
  free_f3tensor(tmp10200cm, tmp10200cm_var);
  free_f3tensor(soilw010cm, soilw010cm_var);
  free_f3tensor(soilw10200cm, soilw10200cm_var);
  free_f3tensor(weasd, weasd_var);

  free_ncvar(tmp010cm_var);
  free_ncvar(tmp10200cm_var);
  free_ncvar(soilw010cm_var);
  free_ncvar(soilw10200cm_var);
  free_ncvar(weasd_var);
  
  return;
}

/*******************************************************************/


struct ncvar * new_ncvar(const char *filename, const char *varname) {
  struct ncvar *v = NULL;
  int varid;			/* variable id */
  int dimid;			/* dimension id */
  float range[2];		/* min and max */
  int ncstat;
  
  v = calloc(1, sizeof(struct ncvar));
  
  /* open file */
  ncstat = nc_open(filename, NC_NOWRITE, &v->ncid);
  check_ncstat(ncstat, filename);

  /* get lat/lon/time sizes */
  ncstat = nc_inq_dimid(v->ncid, "lat", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->nlat); 
  check_ncstat(ncstat, NULL);

  ncstat = nc_inq_dimid(v->ncid, "lon", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->nlon); 
  check_ncstat(ncstat, NULL);

  ncstat = nc_inq_dimid(v->ncid, "time", &dimid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_inq_dimlen(v->ncid, dimid, &v->N); 
  check_ncstat(ncstat, NULL);

  /* get mins and maxes, calculate stepsizes */
  ncstat = nc_inq_varid (v->ncid, "lat", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->minlat = range[0]; v->maxlat = range[1];
  v->latstep = (v->maxlat - v->minlat) / (v->nlat -1.0);
  
  ncstat = nc_inq_varid (v->ncid, "lon", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->minlon = range[0]; v->maxlon = range[1];
  v->lonstep = (v->maxlon - v->minlon) / (v->nlon -1.0);
  
  ncstat = nc_inq_varid (v->ncid, "time", &varid);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, varid, "actual_range", range);
  check_ncstat(ncstat, NULL);
  v->mintime = range[0]; v->maxtime = range[1];

  /* get varid */
  ncstat = nc_inq_varid (v->ncid, varname, &v->varid);
  check_ncstat(ncstat, NULL);

  /* get var unpacking info */
  ncstat = nc_get_att_float(v->ncid, v->varid, "scale_factor", &v->scale);
  check_ncstat(ncstat, NULL);
  ncstat = nc_get_att_float(v->ncid, v->varid, "add_offset", &v->offset); 
  check_ncstat(ncstat, NULL);

  /* get nodata val */
  ncstat = nc_get_att_short(v->ncid, v->varid, "missing_value", &v->nodata); 
  check_ncstat(ncstat, NULL);

  /* record filename for debugging */
  v->filename = strdup(filename);
  
  return v;
} /* new_ncvar */

void _check_ncstat(int stat, const char *filename, unsigned line,
		   const char *msg )
{
  if(stat == NC_NOERR) return; /* no error */
  fprintf(stderr, "%s:%d netCDF error: %s",
	  filename, line, nc_strerror(stat));
  if(msg != NULL) {
    fprintf(stderr, " %s", msg);
  }
  putchar('\n');
  exit(1);
}

float *** ncgetvar(const struct ncvar *v) {
  int ncstat;
  int i;
  size_t len;
  float ***arr = NULL;
  float *ptr = NULL;

  arr = f3tensor(0,  v->N-1, 0, v->nlat-1, 0, v->nlon-1);
  ptr = **arr;
  len = v->N * v->nlat * v->nlon;

  ncstat = nc_get_var_float(v->ncid, v->varid, ptr);
  check_ncstat(ncstat, NULL);
  
  /* scale the values */
  for(i = 0; i<len; i++) {
    if (ptr[i] == v->nodata) {	/* handle nodata */
      ptr[i] = NAN;
    } else {
      ptr[i] *= v->scale;
      ptr[i] += v->offset;
    }
  }
  return arr;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh){
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  t+= NR_END;
  t -= nrl;

  t[nrl] = (float **)malloc ((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  t[nrl] += NR_END;
  t[nrl] -= ncl; 

  t[nrl][ncl] = (float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  for(j=ncl+1;j<=nch;j++)t[nrl][j] = t[nrl][j-1] + ndep;
  for(i=nrl+1;i<=nrh;i++){
    t[i] = t[i-1]+ncol;
    t[i][ncl] = t[i-1][ncl] + ncol * ndep;
    for(j=ncl+1;j<=nch;j++)t[i][j] = t[i][j-1] + ndep;
  }
  return t;
}

void free_ncvar(struct ncvar *var) {
  if(var == NULL) return;
  nc_close(var->ncid);
  free(var->filename);
  free(var);
  var = NULL;
} /* free_ncvar */

void ncnearest(struct ncvar *v, float lat, float lon,
	       size_t *lati, size_t *loni)
{
  if(0)printf("lat %.1f lon %.1f > \n", lat, lon);  
  if (lon < 0) lon += 360;			/* convert to 0-360 range */
  *lati = rint( (lat - v->minlat) / v->latstep );
  *loni = rint( (lon - v->minlon) / v->lonstep );
  if( *lati >= v->nlat) *lati = v->nlat - 1;
  if( *loni >= v->nlon) *loni = v->nlon - 1;
  if(0) {
    /* debug */
    printf("lati %ld loni %ld\n", *lati, *loni);
  }
}

void index2latlon(struct ncvar *v, size_t lati, size_t loni,
		  float *lat, float *lon)
{
  *lat = v->minlat + (v->latstep * lati);
  *lon = v->minlon + (v->lonstep * loni);
  *lon += 0.5 * v->lonstep;

  if(*lon > 180) *lon -= 360;
  /* debug */
  if(0)printf("lat %.1f lon %.1f < lati %ld loni %ld\n", *lat, *lon, lati, loni);
}

void check_cdc_files(char *infilename, int iyear){
    FILE*outfile;
    outfile = fopen(infilename,"r");
    if(outfile == NULL){
	printf("\n");
	printf("==========================================================\n");
	printf("ERROR; STOPPING RUN.\n");
	printf("==========================================================\n");
	printf("\n");
	printf("Bad filename in soil temperature/water and snow initialization.\n");
	printf("File: %s does not exist.\n",infilename);
	printf("Are you sure that you downloaded it?\n");
	printf("\n");
	printf("Note that for soil temperature/water and snow initialization, \n");
	printf("you need to download the following files:\n");
	sprintf(infilename,"tmp.0-10cm.gauss.%d.nc",iyear);
	printf("%s \n",infilename);
	sprintf(infilename,"tmp.10-200cm.gauss.%d.nc",iyear);
	printf("%s \n",infilename);
	sprintf(infilename,"soilw.0-10cm.gauss.%d.nc",iyear);
	printf("%s \n",infilename);
	sprintf(infilename,"soilw.10-200cm.gauss.%d.nc",iyear);
	printf("%s \n",infilename);
	sprintf(infilename,"weasd.sfc.gauss.%d.nc",iyear);
	printf("%s \n",infilename);
	printf("\n");
	printf("These files can be downloaded from:\n");
	printf("http://www.cdc.noaa.gov/cdc/reanalysis/reanalysis.shtml \n");
	printf("\n");
	printf("Note that the tmp* and soilw* files are the 6-hourly, but the\n");
	printf("weasd files are the daily.\n");
	exit(8);
    }
}
void free_f3tensor(float ***t, const struct ncvar *v){
/* free a float f3tensor allocated by f3tensor() */

  long nrl, nrh, ncl, nch, ndl, ndh;

  nrl = 0;
  nrh = v->N-1;
  ncl = 0;
  nch = v->nlat-1;
  ndl = 0;
  ndh = v->nlon-1;

  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
