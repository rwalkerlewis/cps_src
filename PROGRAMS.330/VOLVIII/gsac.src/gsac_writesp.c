/* CHANGES:
	17 MAR 2016 - added ability to output phase, real and imaginary parts of spectra
		      to previous amplitude spectrum
	20 NOV 2019 - DEPMAX DEPMIN DEPMEN were not correctly set. The three lines were moved after
		      the call getmxmn
        17 MAR 2022 - in computing the amplitude spectra, change "float tr, ti" to "double tr, ti"
	19 DEC 2023 - Document ability to write AMplitude, PHase, ReaL and IMaginary parts. Also permit
                      phase to be determined with respect to the O marker in addition to the B marker.
        02 MAR 2024 - use IAMPH and IRLIM for the IF_TYPE
                      
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include        <libgen.h>

extern struct sacfile_ *sacdata;
extern int *sortptr;

#define    NOUT   4

static  int write_sp_ctl[NOUT];

#define	WRITESP_DFLT	0
#define	WRITESP_AM	1
#define	WRITESP_PH	2
#define	WRITESP_RL	3
#define	WRITESP_IM	4
#define	WRITESP_RLIM	5
#define	WRITESP_AMPH	6
#define	WRITESP_B	7
#define	WRITESP_O	8

#define W_AM 0
#define W_PH 1
#define W_RL 2
#define W_IM 3


struct arghdr writesparg[] = {
	{WRITESP_DFLT  , "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{WRITESP_AM    , "AM"  , RHDR, NO, 0, NO, "", -1},
	{WRITESP_PH    , "PH"  , RHDR, NO, 0, NO, "", -1},
	{WRITESP_RL    , "RL"  , RHDR, NO, 0, NO, "", -1},
	{WRITESP_IM    , "IM"  , RHDR, NO, 0, NO, "", -1},
	{WRITESP_RLIM  , "RLIM"  , RHDR, NO, 0, NO, "", 4},
	{WRITESP_AMPH  , "AMPH"  , RHDR, NO, 0, NO, "", 4},
	{WRITESP_B     , "B"  , RHDR, NO, 0, NO, "", 1},
	{WRITESP_O     , "O"  , RHDR, NO, 0, NO, "", 1},
	{0             , ""   , IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float writesp_real[10];
int   writesp_int [10];
int   writesp_yn;
int   writesp_num;
int   writesp_otime;

static float  fhdr_default = -12345.;

static float *temp = (float *)NULL;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_writesp(int ncmd, char **cmdstr)
{
	int i, j;
	/* initial debug */
        writesp_otime = NO;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* initialize */
        for ( i= 0 ; i < NOUT ; i++)
		write_sp_ctl[i] = NO;
 	/* if no option on the command line, continue but force amplitude spectra output */
	if(ncmd == 1){
		write_sp_ctl[W_AM] = YES;
		return;
	}
	if(testarg(ncmd, cmdstr, writesparg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; writesparg[i].key[0] != '\0' ; i++){
		if(writesparg[i].used > 0){
			if(writesparg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, writesp_real);
			} else if(writesparg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, writesp_int );
			} else if(writesparg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, &writesp_yn );
			} else if(writesparg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, &writesp_num );
			}
			switch(writesparg[i].id){
				case WRITESP_AM:
					write_sp_ctl[W_AM] = YES;
					break;
				case WRITESP_PH:
					write_sp_ctl[W_PH] = YES;
					break;
				case WRITESP_RL:
					write_sp_ctl[W_RL] = YES;
					break;
				case WRITESP_IM:
					write_sp_ctl[W_IM] = YES;
					break;
				case WRITESP_RLIM:
					write_sp_ctl[W_IM] = YES;
					write_sp_ctl[W_RL] = YES;
					break;
				case WRITESP_AMPH:
					write_sp_ctl[W_AM] = YES;
					write_sp_ctl[W_PH] = YES;
					break;
				case WRITESP_B:
					writesp_otime = NO ;
					break;
				case WRITESP_O:
					writesp_otime = YES ;
					break;
			}
		}
	}
	/* if nothing has been set force output of amplitude spectrum */
        for(i=0, j=0 ; i < NOUT ; i++){
		j+= write_sp_ctl[i];
	}
	if(j == 0 )
		write_sp_ctl[W_AM] = YES;
}

void gsac_exec_writesp(void)
{
	int i, j, k, ntrc, n2, n21, ncmd;
	float depmax, depmin, depmen,timeb,timeo,tshift;
        double tr, ti;
        double ttr, tti;
	double cphase, sphase;
        float freq;
	int indmax, indmin;
	char tt[1000];
	char *t;
	struct sacfile_ tsac;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
        if(gsac_control.fft == NO){
                printf("Execute FFT first before trying to write spectra\n");
                return;
        }
	/* first get the spectra even if it has already been - this is not
		too time consuming and also ensures that the spectra is of the
		current trace */
	gsac_exec_fft();
	/* now create the create the output and write the spectra - 
		allocate a structure for sac
		fill trace with spectra
		bwsac
	*/
	for(ncmd = 0 ; ncmd < NOUT ;ncmd ++){
		if(write_sp_ctl[ncmd] == YES){
			/* ensure that the sac_data file is nulled - we will free it later */
			/* get a temporary sacfile_ structure */
			for ( k=0 ; k < ntrc ; k ++){
				strcpy(tt, sacdata[k].sac_ifile_name);
				t = basename(tt);
				switch (ncmd){
					case W_AM:  
						strcat(t,".am");
						break;
					case W_PH:  ;;
						strcat(t,".ph");
						break;
					case W_RL:  ;;
						strcat(t,".rl");
						break;
					case W_IM:  ;;
						strcat(t,".im");
						break;
				}
				n2 = sacdata[k].npow2/2 ;
				n21 = n2 + 1 ;
				if(temp == NULL){
					if((temp= (float *)calloc(n21,sizeof(float)) ) == NULL){
						fprintf(stderr,"could not allocate memory for spectra file\n");
						return;
					}
				} else {
					if((temp=(float *)realloc(temp,n21*sizeof(float)))==NULL){
						fprintf(stderr,"could not allocate memory for spectra file\n");
						return;
					}
				}
                                /* implement the pahse shift to O if required */
				tshift = 0.0 ;
				if (writesp_otime == YES){
					/* check to see if the O is set */
					timeb = sacdata[k].sachdr.rhdr[H_B];
					timeo = sacdata[k].sachdr.rhdr[H_O];
					if(timeo != fhdr_default){
						tshift  = timeb - timeo ;
					} else {
						tshift = 0.0 ;
					}
				}
fprintf(stderr,"timeb %f timeo %f tshift %f df=%f n2=%d\n",timeb,timeo,tshift,sacdata[k].df,n2);
				for(i=0, j=0; i<=n2 ; i++){
					freq = i* sacdata[k].df;
					ttr = sacdata[k].sac_spectra[j++];
					tti = sacdata[k].sac_spectra[j++];
					cphase =  cos(6.2831853*freq*tshift);
					sphase = -sin(6.2831853*freq*tshift);
					tr = ttr*cphase - tti*sphase;
					ti = ttr*sphase + tti*cphase;
					switch (ncmd){
						case W_AM:  
							temp[i] = (float)sqrt(tr*tr + ti*ti);
							break;
						case W_PH:  ;;
							temp[i] = atan2(ti,tr);;
							break;
						case W_RL:  ;;
							temp[i] = tr;
							break;
						case W_IM:  ;;
							temp[i] = ti;
							break;
					}
				}
				/* copy the header from sacdata[k] */
				for(i=0 ; i < 70 ; i++)
					tsac.sachdr.rhdr[i] = sacdata[k].sachdr.rhdr[i];
				for(i=0 ; i < 40 ; i++)
					tsac.sachdr.ihdr[i] = sacdata[k].sachdr.ihdr[i];
				for(i=0 ; i < 24 ; i++)
					strncpy(tsac.sachdr.chdr[i],sacdata[k].sachdr.chdr[i],9);
				tsac.sachdr.rhdr[H_B] =  0.0 ;
				tsac.sachdr.ihdr[H_NPTS] =  n21 ;
				tsac.sachdr.rhdr[H_E] =  0.5/sacdata[k].sachdr.rhdr[H_DELTA] ;
				tsac.sachdr.rhdr[H_DELTA] = sacdata[k].df ;
				/* unset time values */
				tsac.sachdr.rhdr[H_O] = -12345. ;
				tsac.sachdr.rhdr[H_A] = -12345. ;
				tsac.sachdr.rhdr[H_T0] = -12345. ;
				tsac.sachdr.rhdr[H_T1] = -12345. ;
				tsac.sachdr.rhdr[H_T2] = -12345. ;
				tsac.sachdr.rhdr[H_T3] = -12345. ;
				tsac.sachdr.rhdr[H_T4] = -12345. ;
				tsac.sachdr.rhdr[H_T5] = -12345. ;
				tsac.sachdr.rhdr[H_T6] = -12345. ;
				tsac.sachdr.rhdr[H_T7] = -12345. ;
				tsac.sachdr.rhdr[H_T8] = -12345. ;
				tsac.sachdr.rhdr[H_T9] = -12345. ;
				switch (ncmd){
					case W_AM:  
						tsac.sachdr.ihdr[H_IFTYPE] =  ENUM_IAMPH ;
						break;
					case W_PH:  ;;
						tsac.sachdr.ihdr[H_IFTYPE] =  ENUM_IAMPH ;
						break;
					case W_RL:  ;;
						tsac.sachdr.ihdr[H_IFTYPE] =  ENUM_IRLIM ;
						break;
					case W_IM:  ;;
						tsac.sachdr.ihdr[H_IFTYPE] =  ENUM_IRLIM ;
						break;
				}
				/* define LEVEN as true */
				tsac.sachdr.ihdr[35] =  1 ;
				tsac.sac_data = temp;
				getmxmn(tsac.sac_data, n21,&depmax, &depmin, &depmen,&indmax,&indmin);
				tsac.sachdr.rhdr[H_DEPMAX] = depmax;
				tsac.sachdr.rhdr[H_DEPMIN] = depmin;
				tsac.sachdr.rhdr[H_DEPMEN] = depmen;
				tsac.sachdr.rhdr[H_TIMMAX] = tsac.sachdr.rhdr[H_B]  + ( indmax)*tsac.sachdr.rhdr[H_DELTA] ;
				tsac.sachdr.rhdr[H_TIMMIN] = tsac.sachdr.rhdr[H_B]  + ( indmin)*tsac.sachdr.rhdr[H_DELTA] ;
				bwsac(t,tsac.sachdr,tsac.sac_data);
		
			}
		}
	}
/*
	free(temp);
*/

}
