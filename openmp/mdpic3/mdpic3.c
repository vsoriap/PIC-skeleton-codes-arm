/*---------------------------------------------------------------------*/
/* Skeleton 3D Darwin OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "mdpush3.h"
#include "omplib.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy/indz = exponent which determines grid points in x/y/z: */
/* direction: nx = 2**indx, ny = 2**indy, nz = 2**indz */
   int indx =   7, indy =   7, indz =   7;
/* npx/npy/npz = number of electrons distributed in x/y/z direction */
   int npx =  384, npy =   384, npz =   384;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 3;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx/vty/vtz = thermal velocity of electrons in x/y/z direction */
   float vtx = 1.0, vty = 1.0, vtz = 1.0;
/* vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction */
   float vx0 = 0.0, vy0 = 0.0, vz0 = 0.0;
/* ax/ay/az = smoothed particle size in x/y/z direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ay = .912871, az = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 6 */
/* ipbc = particle boundary condition: 1 = periodic */
   int idimp = 6, ipbc = 1;
/* omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z */
   float omx = 0.4, omy = 0.0, omz = 0.0;
/* ndc = number of corrections in darwin iteration */
   int ndc = 1;
/* wke/we = particle kinetic/electrostatic field energy             */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
   float zero = 0.0;
/* mx/my/mz = number of grids in x/y/z in sorting tiles */
   int mx = 8, my = 8, mz = 8;
/* xtras = fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* declare scalars for standard code */
   int j, k;
   int np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh;
   int mdim, nxyzh, nxhyz, mx1, my1, mz1, mxyz1;
   int ntime, nloop, isign;
   float qbme, affp, q2m0, wpm, wpmax, wpmin;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, ntmax, npbmx, irc;
   int nvp;

/* declare arrays for standard code: */
/* part = original particle array */
   float *part = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* cue = electron current density with guard cells */
/* dcu = acceleration density with guard cells */
/* cus = transverse electric field with guard cells */
/* amu = momentum flux with guard cells */
   float *cue = NULL, *dcu = NULL, *cus = NULL, *amu = NULL;
/* exyze = smoothed total electric field with guard cells */
/* fxyze = smoothed longitudinal electric field with guard cells */
/* bxyze = smoothed magnetic field with guard cells */
   float *fxyze = NULL, *exyze = NULL, *bxyze = NULL;
/* ffc, ffe = form factor arrays for poisson solvers */
/* sct = sine/cosine table for FFT */
/* ss = scratch array for cwfft2rn */
   float complex *ffc = NULL, *ffe = NULL, *sct = NULL, *ss = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;

/* declare arrays for OpenMP (tiled) code: */
/* ppart = tiled particle array */
/* ppbuff = buffer array for reordering tiled particle array */
   float *ppart = NULL, *ppbuff = NULL;
/* kpic = number of particles in each tile */
/* ncl = number of particles departing tile in each direction */
/* ihole = location/destination of each particle departing tile */
   int *kpic = NULL, *ncl = NULL, *ihole = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0;
   double dtime;

   irc = 0;
/* nvp = number of shared memory nodes  (0=default) */
   nvp = 0;
/* printf("enter number of nodes:\n"); */
/* scanf("%i",&nvp);                   */
/* initialize for shared memory parallel processing */
   cinit_omp(nvp);

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx/ny/nz = number of grid points in x/y direction */
   np = npx*npy*npz; nx = 1L<<indx; ny = 1L<<indy; nz = 1L<<indz;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2; nzh = 1 > nz/2 ? 1 : nz/2;
   nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2;
   nxyzh = (nx > ny ? nx : ny); nxyzh = (nxyzh > nz ? nxyzh : nz)/2;
   nxhyz = nxh > ny ? nxh : ny; nxhyz = nxhyz > nz ? nxhyz : nz;
/* mx1/my1/mz1 = number of tiles in x/y/z direction */
   mx1 = (nx - 1)/mx + 1; my1 = (ny - 1)/my + 1;
   mz1 = (nz - 1)/mz + 1; mxyz1 = mx1*my1*mz1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
/* mdim = dimension of amu array */
   mdim = 2*ndim;
   qbme = qme;
   affp = ((float) nx)*((float) ny)*((float) nz)/(float ) np;

/* allocate data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   qe = (float *) malloc(nxe*nye*nze*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   dcu = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   cus = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   amu = (float *) malloc(mdim*nxe*nye*nze*sizeof(float));
   exyze = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nye*nze*sizeof(float));
   ffc = (float complex *) malloc(nxh*nyh*nzh*sizeof(float complex));
   ffe = (float complex *) malloc(nxh*nyh*nzh*sizeof(float complex));
   mixup = (int *) malloc(nxhyz*sizeof(int));
   sct = (float complex *) malloc(nxyzh*sizeof(float complex));
   ss = (float complex *) malloc(mdim*nxeh*nze*sizeof(float complex));
   kpic = (int *) malloc(mxyz1*sizeof(int));

/* prepare fft tables */
   cwfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh);
/* calculate form factor: ffc */
   isign = 0;
   cmpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
            az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* initialize electrons */
   cdistr3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,ny,nz,
           ipbc);

/* find number of particles in each of mx, my, mz tiles: */
/* updates kpic, nppmx */
   cdblkp3l(part,kpic,&nppmx,idimp,np,mx,my,mz,mx1,my1,mxyz1,&irc);
   if (irc != 0) { 
      printf("cdblkp3l error, irc=%d\n",irc);
      exit(1);
   }
/* allocate vector particle data */
   nppmx0 = (1.0 + xtras)*nppmx;
   ntmax = xtras*nppmx;
   npbmx = xtras*nppmx;
   ppart = (float *) malloc(idimp*nppmx0*mxyz1*sizeof(float));
   ppbuff = (float *) malloc(idimp*npbmx*mxyz1*sizeof(float));
   ncl = (int *) malloc(26*mxyz1*sizeof(int));
   ihole = (int *) malloc(2*(ntmax+1)*mxyz1*sizeof(int));
/* copy ordered particle data for OpenMP: updates ppart and kpic */
   cppmovin3l(part,ppart,kpic,nppmx0,idimp,np,mx,my,mz,mx1,my1,mxyz1,
              &irc);
   if (irc != 0) { 
      printf("cppmovin3l overflow error, irc=%d\n",irc);
      exit(1);
   }
/* sanity check */
   cppcheck3l(ppart,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,my1,mz1,
              &irc);
   if (irc != 0) {
      printf("%d,cppcheck3l error: irc=%d\n",ntime,irc);
      exit(1);
   }

/* find maximum and minimum initial electron density */
   for (j = 0; j < nxe*nye*nze; j++) {
      qe[j] = 0.0;
   }
   cgppost3l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,my1,
             mxyz1);
   caguard3l(qe,nx,ny,nz,nxe,nye,nze);
   cfwpminmx3(qe,qbme,&wpmax,&wpmin,nx,ny,nz,nxe,nye,nze);
   wpm = 0.5*(wpmax + wpmin)*affp;
/* accelerate convergence: update wpm */
   if (wpm <= 10.0)
      wpm = 0.75*wpm;
   printf("wpm=%f\n",wpm);
   q2m0 = wpm/affp;
/* calculate form factor: ffe */
   isign = 0;
   cmepois33((float complex *)dcu,(float complex *)cus,isign,ffe,ax,ay,
             az,affp,wpm,ci,&wf,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);

/* initialize transverse electric field */
   for (j = 0; j < ndim*nxe*nye*nze; j++) {
      cus[j] = 0.0;
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */

/* deposit current with OpenMP: updates cue */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nye*nze; j++) {
         cue[j] = 0.0;
      }
      cgjppost3l(ppart,cue,kpic,qme,zero,nppmx0,idimp,nx,ny,nz,mx,my,mz,
                 nxe,nye,nze,mx1,my1,mxyz1,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nye*nze; j++) {
         qe[j] = 0.0;
      }
      cgppost3l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,
                my1,mxyz1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with OpenMP: updates qe, cue */
      dtimer(&dtime,&itime,-1);
      caguard3l(qe,nx,ny,nz,nxe,nye,nze);
      cacguard3l(cue,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwfft3rmx((float complex *)qe,isign,mixup,sct,indx,indy,indz,nxeh,
                nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate longitudinal force/charge in fourier space with OpenMP: */
/* updates fxyze, we                                                 */ 
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cmpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,
               ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform longitudinal electric force to real space with OpenMP: */
/* updates fxyze                                                    */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwfft3rm3((float complex *)fxyze,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform current to fourier space with OpenMP: update cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwfft3rm3((float complex *)cue,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of current with OpenMP: updates cue */
      dtimer(&dtime,&itime,-1);
      cmcuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate magnetic field in fourier space with OpenMP: */
/* updates bxyze, wm                                      */
      dtimer(&dtime,&itime,-1);
      cmbbpois33((float complex *)cue,(float complex *)bxyze,ffc,ci,&wm,
                 nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform magnetic force to real space with OpenMP: updates bxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwfft3rm3((float complex *)bxyze,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* add constant to magnetic field with OpenMP: updates bxyze */
      dtimer(&dtime,&itime,-1);
      cbaddext3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* copy guard cells with OpenMP: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      ccguard3l(fxyze,nx,ny,nz,nxe,nye,nze);
      ccguard3l(bxyze,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and old transverse electric fields with OpenMP: */
/* updates exyze                                                    */
      dtimer(&dtime,&itime,-1);
      caddvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* deposit electron acceleration density and momentum flux with OpenMP: */
/* updates dcu, amu                                                     */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nye*nze; j++) {
         dcu[j] = 0.0;
      }
      for (j = 0; j < mdim*nxe*nye*nze; j++) {
         amu[j] = 0.0;
      }
      cgdjppost3l(ppart,exyze,bxyze,kpic,dcu,amu,qme,qbme,dt,idimp,
                  nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1);
/* add old scaled electric field with OpenMP: updates dcu */
      cascfguard3l(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdcjpost += time;

/* add guard cells with OpenMP: updates dcu, amu */
      dtimer(&dtime,&itime,-1);
      cacguard3l(dcu,nx,ny,nz,nxe,nye,nze);
      camcguard3l(amu,nx,ny,nz,nxe,nye,nze,mdim);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform acceleration density and momentum flux to fourier space */
/* with OpenMP: updates dcu, amu                                     */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwfft3rm3((float complex *)dcu,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
      cwfft3rmn((float complex *)amu,ss,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,mdim,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of time derivative of current with OpenMP: */
/* updates dcu                                                     */
      dtimer(&dtime,&itime,-1);
      cmadcuperp3((float complex *)dcu,(float complex *)amu,nx,ny,nz,
                  nxeh,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate transverse electric field with OpenMP: updates cus, wf */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cmepois33((float complex *)dcu,(float complex *)cus,isign,ffe,ax,
                ay,az,affp,wpm,ci,&wf,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform transverse electric field to real space with OpenMP: */
/* updates cus                                                    */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwfft3rm3((float complex *)cus,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with OpenMP: updates cus */
      dtimer(&dtime,&itime,-1);
      ccguard3l(cus,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and transverse electric fields with OpenMP: */
/* exyze = cus + fxyze, updates exyze                           */
/* cus needs to be retained for next time step                  */
      dtimer(&dtime,&itime,-1);
      caddvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* inner iteration loop */
      for (k = 0; k < ndc; k++) {

/* deposit electron current and acceleration density and momentum flux */
/* with OpenMP: updates cue, dcu, amu                                  */
         dtimer(&dtime,&itime,-1);
         for (j = 0; j < ndim*nxe*nye*nze; j++) {
            cue[j] = 0.0;
            dcu[j] = 0.0;
         }
         for (j = 0; j < mdim*nxe*nye*nze; j++) {
            amu[j] = 0.0;
         }
         cgdcjppost3l(ppart,exyze,bxyze,kpic,cue,dcu,amu,qme,qbme,dt,
                      idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,
                      my1,mxyz1);
/* add scaled electric field with OpenMP: updates dcu */
         cascfguard3l(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tdcjpost += time;

/* add guard cells for current, acceleration density, and momentum flux */
/* with OpenMP: updates cue, dcu, amu                                   */
         dtimer(&dtime,&itime,-1);
         cacguard3l(cue,nx,ny,nz,nxe,nye,nze);
         cacguard3l(dcu,nx,ny,nz,nxe,nye,nze);
         camcguard3l(amu,nx,ny,nz,nxe,nye,nze,mdim);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;

/* transform current to fourier space with OpenMP: update cue */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwfft3rm3((float complex *)cue,isign,mixup,sct,indx,indy,indz,
                   nxeh,nye,nze,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;

/* take transverse part of current with OpenMP: updates cue */
         dtimer(&dtime,&itime,-1);
         cmcuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* calculate magnetic field in fourier space with OpenMP: */
/* updates bxyze, wm                                      */
         dtimer(&dtime,&itime,-1);
         cmbbpois33((float complex *)cue,(float complex *)bxyze,ffc,ci,
                    &wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform magnetic force to real space with OpenMP: updates bxyze */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwfft3rm3((float complex *)bxyze,isign,mixup,sct,indx,indy,indz,
                   nxeh,nye,nze,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;

/* add constant to magnetic field with OpenMP: updates bxyze */
         dtimer(&dtime,&itime,-1);
         cbaddext3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform acceleration density and momentum flux to fourier space */
/* with OpenMP: updates dcu and amu                                  */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwfft3rm3((float complex *)dcu,isign,mixup,sct,indx,indy,indz,
                   nxeh,nye,nze,nxhyz,nxyzh);
         cwfft3rmn((float complex *)amu,ss,isign,mixup,sct,indx,indy,
                   indz,nxeh,nye,nze,mdim,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;
 
/* take transverse part of time derivative of current with OpenMP: */
/* updates dcu                                                     */
         dtimer(&dtime,&itime,-1);
         cmadcuperp3((float complex *)dcu,(float complex *)amu,nx,ny,nz,
                     nxeh,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* calculate convective part of transverse electric field with OpenMP: */
/* updates cus, wf                                                     */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cmepois33((float complex *)dcu,(float complex *)cus,isign,ffe,
                   ax,ay,az,affp,wpm,ci,&wf,nx,ny,nz,nxeh,nye,nze,nxh,
                   nyh,nzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* transform transverse electric field to real space with OpenMP: */
/* updates cus                                                    */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwfft3rm3((float complex *)cus,isign,mixup,sct,indx,indy,indz,
                   nxeh,nye,nze,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;
 
/* copy guard cells with OpenMP: updates bxyze, cus */
         dtimer(&dtime,&itime,-1);
         ccguard3l(bxyze,nx,ny,nz,nxe,nye,nze);
         ccguard3l(cus,nx,ny,nz,nxe,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;
 
/* add longitudinal and transverse electric fields with OpenMP: */
/* exyze = cus + fxyze, updates exyze                           */
/* cus needs to be retained for next time step                  */
         dtimer(&dtime,&itime,-1);
         caddvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

      }

/* push particles with OpenMP: */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
/* updates ppart, wke */
/*    cgbppush3l(ppart,exyze,bxyze,kpic,qbme,dt,dt,&wke,idimp,nppmx0, */
/*               nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc);   */
/* updates ppart, ncl, ihole, wke, irc */
      cgbppushf3l(ppart,exyze,bxyze,kpic,ncl,ihole,qbme,dt,dt,&wke,
                  idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,
                  mxyz1,ntmax,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) {
         printf("cgbppushf3l error: irc=%d\n",irc);
         exit(1);
      }

/* reorder particles by tile with OpenMP: */
      dtimer(&dtime,&itime,-1);
/* updates ppart, ppbuff, kpic, ncl, ihole, and irc */
/*    cpporder3l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,nz, */
/*               mx,my,mz,mx1,my1,mz1,npbmx,ntmax,&irc);            */
/* updates ppart, ppbuff, kpic, ncl, and irc */
      cpporderf3l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,my1,
                  mz1,npbmx,ntmax,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("cpporderf3l error: ntmax, irc=%d,%d\n",ntmax,irc);
         exit(1);
      }

      if (ntime==0) {
         wt = we + wm;
         printf("Initial Total Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",wt,wke,wke+wt);
         printf("Initial Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
         printf("%e %e %e\n",we,wf,wm);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime, ndc = %i,%i\n",ntime,ndc);
   wt = we + wm;
   printf("Final Total Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",wt,wke,wke+wt);
   printf("Final Electrostatic, Transverse Electric and Magnetic Field \
Energies:\n");
   printf("%e %e %e\n",we,wf,wm);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
   printf("current deposit time = %f\n",tdjpost);
   printf("current derivative deposit time = %f\n",tdcjpost);
   tdpost += tdjpost + tdcjpost;
   printf("total deposit time = %f\n",tdpost);
   printf("guard time = %f\n",tguard);
   printf("solver time = %f\n",tfield);
   printf("fft time = %f\n",tfft);
   printf("push time = %f\n",tpush);
   printf("sort time = %f\n",tsort);
   tfield += tguard + tfft;
   printf("total solver time = %f\n",tfield);
   time = tdpost + tpush + tsort;
   printf("total particle time = %f\n",time);
   wt = time + tfield;
   printf("total time = %f\n",wt);
   printf("\n");

   wt = 1.0e+09/(((float) nloop)*((float) np));
   printf("Push Time (nsec) = %f\n",tpush*wt);
   printf("Deposit Time (nsec) = %f\n",tdpost*wt);
   printf("Sort Time (nsec) = %f\n",tsort*wt);
   printf("Total Particle Time (nsec) = %f\n",time*wt);
   printf("\n");

   return 0;
}
