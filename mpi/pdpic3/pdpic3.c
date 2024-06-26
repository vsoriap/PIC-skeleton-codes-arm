/*---------------------------------------------------------------------*/
/* Skeleton 3D Darwin MPI PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "pdpush3.h"
#include "pplib3.h"

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
/* sortime = number of time steps between standard electron sorting */
   int idimp = 6, ipbc = 1, sortime = 20;
/* omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z */
   float omx = 0.4, omy = 0.0, omz = 0.0;
/* ndc = number of corrections in darwin iteration */
   int ndc = 1;
/* idps = number of partition boundaries = 4 */
/* idds = dimensionality of domain decomposition = 2 */
   int idps = 4, idds =    2;
/* wke/we = particle kinetic/electrostatic field energy             */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
   float zero = 0.0;
/* declare scalars for standard code */
   int j, k;
   int nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh, nnxe;
   int mdim, nxyzh, nxhyz, ntime, nloop, isign, ierr;
   float qbme, affp, q2m0, wpm, wpmax, wpmin;
   double np;

/* declare scalars for MPI code */
   int ntpose = 1;
   int nvpy, nvpz, nvp, idproc, kstrt, npmax, kyp, kzp;
   int kxyp, kyzp, kzyp, nypmx, nzpmx, nypmn, nzpmn, npp, nps;
   int nyzpm1, nbmax, ntmax;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
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
/* qt, qs = scalar charge density field arrays in fourier space */
   float complex *qt = NULL, *qs = NULL;
/* cut = vector current density field arrays in fourier space */
/* dcut = vector acceleration density in fourier space */
/* exyzt = vector transverse electric field in fourier space */
/* amut = tensor momentum flux in fourier space */
   float complex *cut = NULL, *dcut = NULL, *exyzt = NULL, *amut = NULL; 
/* fxyzt = vector longitudinal electric field in fourier space */
/* bxyzt = vector magnetic field array in fourier space */
   float complex *fxyzt = NULL, *bxyzt = NULL;
/* fxyzs = vector field array in fourier space */
/* amus = tensor field array in fourier space */
   float complex *fxyzs = NULL, *amus = NULL;
/* ffc, ffe = form factor arrays for poisson solvers */
/* sct = sine/cosine table for FFT */
/* ss = scratch array for cwfft3rn */
   float complex *ffc = NULL, *ffe = NULL, *sct = NULL, *ss = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* ihole = location of hole left in particle arrays */
   int *ihole = NULL;
/* npic = scratch array for reordering particles */
   int *npic = NULL;
   double wtot[7], work[7];
   int info[7];

/* declare arrays for MPI code: */
/* bs/br = complex send/receive buffers for data transpose */
   float complex *bs = NULL, *br = NULL;
/* sbufl/sbufr = particle buffers sent to nearby processors */
/* rbufl/rbufr = particle buffers received from nearby processors */
   float *sbufl = NULL, *sbufr = NULL, *rbufl = NULL, *rbufr = NULL;
/* edges[0:1] = lower:upper y boundaries of particle partition */
/* edges[2:3] = back:front z boundaries of particle partition */
   float *edges = NULL;
/* nyzp[0:1] = number of primary (complete) gridpoints in y/z */
/* noff[0:1] = lowermost global gridpoint in y/z */
   int *nyzp = NULL, *noff = NULL;
/* scr/scs = guard cell buffers received/sent from nearby processors */
   float *scr = NULL, *scs = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0;
   float tmov = 0.0;
   float tfft[2] = {0.0,0.0};
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
   np = ((double) npx)*((double) npy)*((double) npz);
/* nx/ny/nz = number of grid points in x/y direction */
   nx = 1L<<indx; ny = 1L<<indy; nz = 1L<<indz;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2; nzh = 1 > nz/2 ? 1 : nz/2;
   nxe = nx + 2; nye = ny + 2; nze = nz + 2;
   nxeh = nxe/2;  nnxe = ndim*nxe;
   nxyzh = (nx > ny ? nx : ny); nxyzh = (nxyzh > nz ? nxyzh : nz)/2;
   nxhyz = nxh > ny ? nxh : ny; nxhyz = nxhyz > nz ? nxhyz : nz;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
/* mdim = dimension of amu array */
   mdim = 2*ndim;
   qbme = qme;
   affp = ((float) nx)*((float) ny)*((float) nz)/(float ) np;

/* nvp = number of MPI ranks */
/* initialize for distributed memory parallel processing */
   cppinit2(&idproc,&nvp,argc,argv);
   kstrt = idproc + 1;
/* obtain 2D partition (nvpy,nvpz) from nvp: */
/* nvpy/nvpz = number of processors in y/z */
   cfcomp32(nvp,nx,ny,nz,&nvpy,&nvpz,&ierr);
   if (ierr != 0) {
      if (kstrt==1) {
         printf("cfcomp32 error: nvp,nvpy,nvpz=%d,%d,%d\n",nvp,nvpy,nvpz);
      }
      goto L3000;
   }

/* initialize data for MPI code */
   edges = (float *) malloc(idps*sizeof(float));
   nyzp = (int *) malloc(idds*sizeof(int));
   noff = (int *) malloc(idds*sizeof(int));
/* calculate partition variables: */
/* edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn                          */
/* edges[0:1] = lower:upper boundary of particle partition in y           */
/* edges[2:3] = back:front boundary of particle partition in z            */
/* nyzp[0:1] = number of primary (complete) gridpoints in y/z             */
/* noff[0:1] = lowermost global gridpoint in y/z in particle partition    */
/* nypmx = maximum size of particle partition in y, including guard cells */
/* nzpmx = maximum size of particle partition in z, including guard cells */
/* nypmn = minimum value of nyzp[0]                                       */
/* nzpmn = minimum value of nyzp[1]                                       */
   cpdicomp32l(edges,nyzp,noff,&nypmx,&nzpmx,&nypmn,&nzpmn,ny,nz,kstrt,
               nvpy,nvpz,idps,idds);
   if (kstrt==1) {
      if (nypmn < 1) {
                     printf("combination not supported nvpy,ny= %d,%d\n",
                            nvpy,ny);
      }
      if (nzpmn < 1) {
                     printf("combination not supported nvpz,nz= %d,%d\n",
                            nvpz,nz);
      }
   }
   if ((nypmn < 1) || (nzpmn < 1)) goto L3000;
/* initialize additional scalars for MPI code */
/* kyp = number of complex grids in each field partition in y direction */
   kyp = (ny - 1)/nvpy + 1;
/* kzp = number of complex grids in each field partition in z direction */
   kzp = (nz - 1)/nvpz + 1;
/* kxyp = number of complex grids in each field partition in x direction */
/* in transposed data */
   kxyp = (nxh - 1)/nvpy + 1;
/* kyzp = number of complex grids in each field partition in y direction, */
/* in transposed data */
   kyzp = (ny - 1)/nvpz + 1; kzyp = kyzp > kyp ? kyzp : kyp;
/* dimension for scratch array for reordering particles */
   nyzpm1 = (kyp + 1)*(kzp + 1);
/* npmax = maximum number of electrons in each partition */
   npmax = (np/nvp)*1.25;
/* nbmax = size of buffer for passing particles between processors */
   nbmax = 0.1*npmax;
/* ntmax = size of ihole buffer for particles leaving processor */
   ntmax = 2*nbmax;

/* allocate data for standard code */
   part = (float *) malloc(idimp*npmax*sizeof(float));
   if (sortime > 0)
      part2 = (float *) malloc(idimp*npmax*sizeof(float));
   qe = (float *) malloc(nxe*nypmx*nzpmx*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   dcu = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   cus = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   amu = (float *) malloc(mdim*nxe*nypmx*nzpmx*sizeof(float));
   exyze = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   qt = (float complex *) malloc(nze*kxyp*kyzp*sizeof(float complex));
   qs = (float complex *) malloc(nye*kxyp*nzpmx*sizeof(float complex));
   exyzt = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                    *sizeof(float complex));
   fxyzt = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                    *sizeof(float complex));
   bxyzt = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                    *sizeof(float complex));
   cut = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                  *sizeof(float complex));
   dcut = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                   *sizeof(float complex));
   amut = (float complex *) malloc(mdim*nze*kxyp*kyzp
                                   *sizeof(float complex));
   fxyzs = (float complex *) malloc(ndim*nye*kxyp*nzpmx
                                    *sizeof(float complex));
   amus = (float complex *) malloc(mdim*nye*kxyp*nzpmx
                                   *sizeof(float complex));
   ffc = (float complex *) malloc(nzh*kxyp*kyzp*sizeof(float complex));
   ffe = (float complex *) malloc(nzh*kxyp*kyzp*sizeof(float complex));
   mixup = (int *) malloc(nxhyz*sizeof(int));
   sct = (float complex *) malloc(nxyzh*sizeof(float complex));
   ihole = (int *) malloc((ntmax+1)*2*sizeof(int));
   npic = (int *) malloc(nyzpm1*sizeof(int));
   ss = (float complex *) malloc(mdim*nxeh*sizeof(float complex));

/* allocate data for MPI code */
   bs = (float complex *) malloc(mdim*kxyp*kzyp*kzp*sizeof(float complex));
   br = (float complex *) malloc(mdim*kxyp*kzyp*kzp*sizeof(float complex));
   sbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   sbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   scr = (float *) malloc(mdim*nxe*nypmx*sizeof(float));
   scs = (float *) malloc(mdim*nxe*2*nzpmx*sizeof(float));

/* prepare fft tables */
   cwpfft32rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh);
/* calculate form factor: ffc */
   isign = 0;
   cppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,&we,nx,ny,nz,kstrt,nvpy,
             nvpz,nze,kxyp,kyzp,nzh);

/* initialize electrons */
   nps = 1;
   npp = 0;
   cpdistr32(part,edges,&npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,nx,
             ny,nz,idimp,npmax,idps,ipbc,&ierr);
/* check for particle initialization error */
   if (ierr != 0) {
      if (kstrt==1) {
         printf("particle initialization error: ierr=%d\n",ierr);
      }
      goto L3000;
   }

/* find maximum and minimum initial electron density */
   for (j = 0; j < nxe*nypmx*nzpmx; j++) {
      qe[j] = 0.0;
   }
   cppgpost32l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx,nzpmx,idds);
   cppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds);
   cppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,
                 idds);

   cppfwpminmx32(qe,nyzp,qbme,&wpmax,&wpmin,nx,nxe,nypmx,nzpmx,idds);
   wtot[0] = wpmax;
   wtot[1] = -wpmin;
   cppdmax(wtot,work,2);
   wpmax = wtot[0];
   wpmin = -wtot[1];
   wpm = 0.5*(wpmax + wpmin)*affp;
/* accelerate convergence: update wpm */
   if (wpm <= 10.0)
      wpm = 0.75*wpm;
   if (kstrt==1)
      printf("wpm=%f\n",wpm);
   q2m0 = wpm/affp;
/* calculate form factor: ffe */
   isign = 0;
   cppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,&wf,nx,ny,nz,
                kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh);

/* initialize transverse electric field */
   for (j = 0; j < ndim*nxe*nypmx*nzpmx; j++) {
      cus[j] = 0.0;
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    if (kstrt==1) printf("ntime = %i\n",ntime); */

/* deposit current with standard procedure: updates cue and ihole */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx*nzpmx; j++) {
         cue[j] = 0.0;
      }
      cppgjpost32l(part,cue,edges,npp,noff,ihole,qme,zero,nx,ny,nz,
                   idimp,npmax,nxe,nypmx,nzpmx,idps,idds,ntmax,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nypmx*nzpmx; j++) {
         qe[j] = 0.0;
      }
      cppgpost32l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates qe, cue */
      dtimer(&dtime,&itime,-1);
      cppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds);
      cppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,
                    idds);
      cppacguard32xl(cue,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      cppnacguard32l(cue,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                     nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: updates qt */
/* modifies qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft32r((float complex *)qe,qs,qt,bs,br,isign,ntpose,mixup,sct,
                 &ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,
                 kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* calculate longitudinal force/charge in fourier space with standard */
/* procedure: updates fxyzt, we                                       */ 
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,&we,nx,ny,nz,kstrt,nvpy,
                nvpz,nze,kxyp,kyzp,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform longitudinal electric force to real space with standard */
/* procedure: updates fxyze, modifies fxyzt                          */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft32r3((float complex *)fxyze,fxyzs,fxyzt,bs,br,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* transform current to fourier space with standard procedure: update cut */
/* modifies cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft32r3((float complex *)cue,fxyzs,cut,bs,br,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* take transverse part of current with standard procedure: updates cut */
      dtimer(&dtime,&itime,-1);
      cppcuperp32(cut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates bxyzt, wm                                                  */
      dtimer(&dtime,&itime,-1);
      cppbbpoisp332(cut,bxyzt,ffc,ci,&wm,nx,ny,nz,kstrt,nvpy,nvpz,nze,
                    kxyp,kyzp,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze, modifies bxyzt                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft32r3((float complex *)bxyze,fxyzs,bxyzt,bs,br,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* add constant to magnetic field with standard procedure: updates bxyze */
      dtimer(&dtime,&itime,-1);
      cppbaddext32(bxyze,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* copy guard cells with standard procedure: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      cppncguard32l(fxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,
                    idds);
      cppcguard32xl(fxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      cppncguard32l(bxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,
                    idds);
      cppcguard32xl(bxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and old transverse electric fields with standard */
/* procedure: updates exyze                                          */
      dtimer(&dtime,&itime,-1);
      cppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* deposit electron acceleration density and momentum flux with */
/* standard procedure: updates dcu, amu                         */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx*nzpmx; j++) {
         dcu[j] = 0.0;
      }
      for (j = 0; j < mdim*nxe*nypmx*nzpmx; j++) {
         amu[j] = 0.0;
      }
      cppgdjpost32l(part,exyze,bxyze,npp,noff,dcu,amu,qme,qbme,dt,idimp,
      npmax,nxe,nypmx,nzpmx,idds);
/* add old scaled electric field with standard procedure: updates dcu */
      cppascfguard32l(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdcjpost += time;

/* add guard cells with standard procedure: updates dcu, amu */
      dtimer(&dtime,&itime,-1);
      cppacguard32xl(dcu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      cppnacguard32l(dcu,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                     nzpmx,idds);
      cppacguard32xl(amu,nyzp,nx,mdim,nxe,nypmx,nzpmx,idds);
      cppnacguard32l(amu,scs,scr,nyzp,mdim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                     nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcut, amut, modifies dcu, amu    */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft32r3((float complex *)dcu,fxyzs,dcut,bs,br,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      tfft[1] += ttp;
      cwppfft32rn((float complex *)amu,amus,amut,bs,br,ss,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,mdim,
                  nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* take transverse part of time derivative of current with standard */
/* procedure: updates dcut                                          */
      dtimer(&dtime,&itime,-1);
      cppadcuperp32(dcut,amut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate transverse electric field with standard procedure: */
/* updates exyzt, wf                                              */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,&wf,nx,ny,
                   nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform transverse electric field to real space with standard */
/* procedure: updates cus, modifies exyzt                          */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft32r3((float complex *)cus,fxyzs,exyzt,bs,br,isign,ntpose,
                  mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,
                  nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* copy guard cells with standard procedure: updates cus */
      dtimer(&dtime,&itime,-1);
      cppncguard32l(cus,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds);
      cppcguard32xl(cus,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxyze, updates exyze                 */
/* cus needs to be retained for next time step                   */
      dtimer(&dtime,&itime,-1);
      cppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* inner iteration loop */
      for (k = 0; k < ndc; k++) {

/* deposit electron current and acceleration density and momentum flux */
/* with standard procedure: updates cue, dcu, amu                      */
         dtimer(&dtime,&itime,-1);
         for (j = 0; j < ndim*nxe*nypmx*nzpmx; j++) {
            cue[j] = 0.0;
            dcu[j] = 0.0;
         }
         for (j = 0; j < mdim*nxe*nypmx*nzpmx; j++) {
            amu[j] = 0.0;
         }
         cppgdcjpost32l(part,exyze,bxyze,npp,noff,cue,dcu,amu,qme,qbme,
                        dt,idimp,npmax,nxe,nypmx,nzpmx,idds);
/* add scaled electric field with standard procedure: updates dcu */
         cppascfguard32l(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,idds);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tdcjpost += time;

/* add guard cells for current, acceleration density, and momentum flux */
/* with standard procedure: updates cue, dcu, amu                       */
         dtimer(&dtime,&itime,-1);
         cppacguard32xl(cue,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
         cppnacguard32l(cue,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,
                        nypmx,nzpmx,idds);
         cppacguard32xl(dcu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
         cppnacguard32l(dcu,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,
                        nypmx,nzpmx,idds);
         cppacguard32xl(amu,nyzp,nx,mdim,nxe,nypmx,nzpmx,idds);
         cppnacguard32l(amu,scs,scr,nyzp,mdim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                     nzpmx,idds);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;

/* transform current to fourier space with standard procedure: update cut */
/* modifies cue */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwppfft32r3((float complex *)cue,fxyzs,cut,bs,br,isign,ntpose,
                     mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,
                     nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,
                     kzyp,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;

/* take transverse part of current with standard procedure: updates cut */
         dtimer(&dtime,&itime,-1);
         cppcuperp32(cut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates bxyzt, wm                                                  */
         dtimer(&dtime,&itime,-1);
         cppbbpoisp332(cut,bxyzt,ffc,ci,&wm,nx,ny,nz,kstrt,nvpy,nvpz,nze,
                       kxyp,kyzp,nzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze, modifies bxyzt                                   */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwppfft32r3((float complex *)bxyze,fxyzs,bxyzt,bs,br,isign,
                     ntpose,mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,
                     nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,
                     nzpmx,kzyp,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;

/* add constant to magnetic field with standard procedure: updates bxyze */
         dtimer(&dtime,&itime,-1);
         cppbaddext32(bxyze,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx,idds);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcut, amut, modifies dcu, amu    */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwppfft32r3((float complex *)dcu,fxyzs,dcut,bs,br,isign,ntpose,
                     mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,
                     nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,
                     kzyp,nxhyz,nxyzh);
         tfft[1] += ttp;
         cwppfft32rn((float complex *)amu,amus,amut,bs,br,ss,isign,
                     ntpose,mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,
                     nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,
                     nzpmx,kzyp,mdim,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;
 
/* take transverse part of time derivative of current with standard */
/* procedure: updates dcut                                          */
         dtimer(&dtime,&itime,-1);
         cppadcuperp32(dcut,amut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* calculate convective part of transverse electric field with standard */
/* procedure: updates exyzt, wf                                         */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,&wf,nx,
                      ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* transform transverse electric field to real space with standard */
/* procedure: updates cus, modifies exyzt                          */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwppfft32r3((float complex *)cus,fxyzs,exyzt,bs,br,isign,ntpose,
                     mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,
                     nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,
                     kzyp,nxhyz,nxyzh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;
 
/* copy guard cells with standard procedure: updates bxyze, cus */
         dtimer(&dtime,&itime,-1);
         cppncguard32l(bxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,
                       idds);
         cppcguard32xl(bxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
         cppncguard32l(cus,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,
                       idds);
         cppcguard32xl(cus,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;
 
/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxyze, updates exyze                 */
/* cus needs to be retained for next time step                   */
         dtimer(&dtime,&itime,-1);
         cppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

      }

/* push particles with standard procedure: updates part, wke and ihole */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      cppgbpush32l(part,exyze,bxyze,edges,npp,noff,ihole,qbme,dt,dt,&wke,
                   nx,ny,nz,idimp,npmax,nxe,nypmx,nzpmx,idps,idds,ntmax,
                   ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
/* check for ihole overflow error */
      if (ihole[0] < 0) {
         ierr = -ihole[0];
         printf("ihole overflow error: ntmax,ih=%d,%d\n",ntmax,ierr);
         cppabort();
         goto L3000;
      }

/* move electrons into appropriate spatial regions: updates part, npp */
      dtimer(&dtime,&itime,-1);
      cppmove32(part,edges,&npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,nz,
                kstrt,nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* check for particle manager error */
      if (info[0] != 0) {
         ierr = info[0];
         if (kstrt==1) {
            printf("particle manager error: ierr=%d\n",ierr);
         }
         goto L3000;
      }

/* sort particles by cell for standard procedure */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            cppdsortp32yzl(part,part2,npic,npp,noff,nyzp,idimp,npmax,
                           nyzpm1,idds);
/* exchange pointers */
            tpart = part;
            part = part2;
            part2 = tpart;
            dtimer(&dtime,&itime,1);
            time = (float) dtime;
            tsort += time;
         }
      }

/* energy diagnostic */
      wt = we + wm;
      wtot[0] = wt;
      wtot[1] = wke;
      wtot[2] = 0.0;
      wtot[3] = wke + wt;
      wtot[4] = we;
      wtot[5] = wf;
      wtot[6] = wm;
      cppdsum(wtot,work,7);
      wke = wtot[1];
      we = wtot[4];
      wf = wtot[5];
      wm = wtot[6];
      if ((ntime==0) && (kstrt==1)) {
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

   if (kstrt==1) {
      printf("ntime, ndc = %i,%i\n",ntime,ndc);
      printf("MPI nodes nvpy, nvpz = %i,%i\n",nvpy,nvpz);
      wt = we + wm;
      printf("Final Total Field, Kinetic and Total Energies:\n");
      printf("%e %e %e\n",wt,wke,wke+wt);
      printf("Final Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
      printf("%e %e %e\n",we,wf,wm);

      printf("\n");
      printf("deposit time = %f\n",tdpost);
      printf("current deposit time = %f\n",tdjpost);
      printf("current derivative deposit time = %f\n",tdcjpost);
      tdpost += tdjpost + tdcjpost;
      printf("total deposit time = %f\n",tdpost);
      printf("guard time = %f\n",tguard);
      printf("solver time = %f\n",tfield);
      printf("fft and transpose time = %f,%f\n",tfft[0],tfft[1]);
      printf("push time = %f\n",tpush);
      printf("particle move time = %f\n",tmov);
      printf("sort time = %f\n",tsort);
      tfield += tguard + tfft[0];
      printf("total solver time = %f\n",tfield);
      tsort += tmov;
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
   }

L3000:
   cppexit();
   return 0;
}
