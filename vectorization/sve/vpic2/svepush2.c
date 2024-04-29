/* SVE C Library for Skeleton 2D Electrostatic Vector PIC Code */
/* written by Víctor Soria, BSC */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "svepush2.h"
#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif /* __ARM_FEATURE_SVE */

/*--------------------------------------------------------------------*/
void csvexiscan2(int *isdata, int nths) {
}

/*--------------------------------------------------------------------*/
void csvegpush2lt(float part[], float fxy[], float qbm, float dt,
                   float *ek, int idimp, int nop, int npe, int nx,
                   int ny, int nxv, int nyv, int ipbc) {
/* for 2d code, this subroutine updates particle co-ordinates and
   velocities using leap-frog scheme in time and first-order linear
   interpolation in space, with various boundary conditions.
   vector version using guard cells
   44 flops/particle, 12 loads, 4 stores
   input: all, output: part, ek
   equations used are:
   vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
   vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
   where q/m is charge/mass, and
   x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
   fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
   the nearest grid points:
   fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
      + dx*fx(n+1,m+1))
   fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
      + dx*fy(n+1,m+1))
   where n,m = leftmost grid points and dx = x-n, dy = y-m
   part[0][n] = position x of particle n
   part[1][n] = position y of particle n
   part[2][n] = velocity vx of particle n
   part[3][n] = velocity vy of particle n
   fxy[k][j][0] = x component of force/charge at grid (j,k)
   fxy[k][j][1] = y component of force/charge at grid (j,k)
   that is, convolution of electric field over particle shape
   qbm = particle charge/mass
   dt = time interval between successive calculations
   kinetic energy/mass at time t is also calculated, using
   ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
   idimp = size of phase space = 4
   nop = number of particles
   npe = first dimension of particle array
   nx/ny = system length in x/y direction
   nxv = second dimension of field arrays, must be >= nx+1
   nyv = third dimension of field arrays, must be >= ny+1
   ipbc = particle boundary condition = (0,1,2,3) =
   (none,2d periodic,2d reflecting,mixed reflecting/periodic)
   requires SVE, part needs to be 16 byte aligned
   npe needs to be a multiple of 4
local data                                                            */
   uint32_t j, nps, nn, mm;
   float qtm, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy;
   float x, y, dx, dy, vx, vy;
   double sum1;
   svbool_t predicate;

   svint32_t v_nn, v_mm, v_it;
   svfloat32_t v_qtm, v_dt;
   svfloat32_t v_dxp, v_dyp, v_amx, v_amy, v_at;
   svfloat32_t v_x, v_y, v_dx, v_dy, v_vx, v_vy;
   svfloat32_t a, b, c, d;
   svfloat64_t v_sum1, v_d;
   uint64_t numInc = svlen_f32(v_x);
   qtm = qbm*dt;
   nps = numInc*(nop/numInc);
/* set boundary values */
   edgelx = 0.0f;
   edgely = 0.0f;
   edgerx = (float) nx;
   edgery = (float) ny;
   if (ipbc==2) {
      edgelx = 1.0f;
      edgely = 1.0f;
      edgerx = (float) (nx-1);
      edgery = (float) (ny-1);
   }
   else if (ipbc==3) {
      edgelx = 1.0f;
      edgerx = (float) (nx-1);
   }
   v_qtm = svdup_n_f32(qtm);
   v_dt = svdup_n_f32(dt);
   v_sum1 = svdup_f64(0.0);

/* vector loop over particles in blocks of 4 */
   for (j = 0; j < nps; j+=numInc) {
      predicate = svptrue_b32();
/* find interpolation weights */
/*    x = part[j];     */
/*    y = part[j+npe]; */
      v_x = svld1_f32(predicate, &part[j]);
      v_y = svld1_f32(predicate, &part[j+npe]);
/*    nn = x; */
/*    mm = y; */
      v_nn = svcvt_s32_f32_x(predicate,v_x);
      v_mm = svcvt_s32_f32_x(predicate,v_y);
/*    dxp = x - (float) nn; */
      v_dxp = svsub_f32_x(predicate,v_x,svcvt_f32_s32_x(predicate,v_nn));
/*    dyp = y - (float) mm; */
      v_dyp = svsub_f32_x(predicate,v_y,svcvt_f32_s32_x(predicate,v_mm));
/*    nn = 2*(nn + nxv*mm); */
      v_mm = svmla_n_s32_x(predicate,v_nn,v_mm,nxv);
      v_nn = svlsl_n_s32_x(predicate,v_mm,1);
/*    amx = 1.0f - dxp; */
/*    amy = 1.0f - dyp; */
      v_amx = svsubr_n_f32_x(predicate,v_dxp,1.0f);
      v_amy = svsubr_n_f32_x(predicate,v_dyp,1.0f);
/* find acceleration */
/* load fields, for lower left/right */
      v_it = svlsl_n_s32_x(predicate,v_nn,2);
      a = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      b = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      c = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      d = svld1_gather_s32offset_f32(predicate,fxy,v_it);
/*    dx = amx*fxy[nn];   */
/*    dy = amx*fxy[nn+1]; */
      v_dx = svmul_f32_x(predicate,v_amx,a);
      v_dy = svmul_f32_x(predicate,v_amx,b);
/*    dx = amy*(dxp*fxy[nn+2] + dx); */
/*    dy = amy*(dxp*fxy[nn+3] + dy); */
      v_dx = svmul_f32_x(predicate,v_amy,svadd_f32_x(predicate, svmul_f32_x(predicate,v_dxp,c),v_dx));
      v_dy = svmul_f32_x(predicate,v_amy,svadd_f32_x(predicate,svmul_f32_x(predicate,v_dxp,d),v_dy));
/*    nn += 2*nxv; */
/* load fields, for upper left/right */
      v_nn = svadd_n_s32_x(predicate,v_nn,2*nxv);
      v_it = svlsl_n_s32_x(predicate,v_nn,2);
      a = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      b = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      c = svld1_gather_s32offset_f32(predicate,fxy,v_it);
      v_it = svadd_n_s32_x(predicate,v_it,4);
      d = svld1_gather_s32offset_f32(predicate,fxy,v_it);
/*    vx = amx*fxy[nn];   */
/*    vy = amx*fxy[nn+1]; */
      a = svmul_f32_x(predicate,v_amx,a);
      b = svmul_f32_x(predicate,v_amx,b);
/*    dx += dyp*(dxp*fxy[nn+2] + vx); */
/*    dy += dyp*(dxp*fxy[nn+3] + vy); */
      a = svmul_f32_x(predicate,v_dyp,svadd_f32_x(predicate,svmul_f32_x(predicate,v_dxp,c),a));
      b = svmul_f32_x(predicate,v_dyp,svadd_f32_x(predicate,svmul_f32_x(predicate,v_dxp,d),b));
      v_dx = svadd_f32_x(predicate,v_dx,a);
      v_dy = svadd_f32_x(predicate,v_dy,b);
/* new velocity */
/*    dxp = part[j+2*npe]; */
/*    dyp = part[j+3*npe]; */
      v_dxp = svld1_f32(predicate,&part[j+2*npe]);
      v_dyp = svld1_f32(predicate,&part[j+3*npe]);
/*    vx = dxp + qtm*dx; */
/*    vy = dyp + qtm*dy; */
      v_vx = svadd_f32_x(predicate,v_dxp,svmul_f32_x(predicate,v_qtm,v_dx));
      v_vy = svadd_f32_x(predicate,v_dyp,svmul_f32_x(predicate,v_qtm,v_dy));
/* average kinetic energy */
/*    dxp += vx; */
/*    dyp += vy; */
      v_dxp = svadd_f32_x(predicate,v_dxp,v_vx);
      v_dyp = svadd_f32_x(predicate,v_dyp,v_vy);
/*    sum1 += dxp*dxp + dyp*dyp; */
      v_at = svmul_f32_x(predicate,v_dxp,v_dxp);
      v_at = svadd_f32_x(predicate,v_at,svmul_f32_x(predicate,v_dyp,v_dyp));
/* convert to double precision before accumulating */
      v_d = svcvt_f64_f32_x(predicate,v_at);
      v_sum1 = svadd_f64_x(predicate,v_sum1,v_d);
      v_at = svtrn2_f32(v_at,v_at);
      v_d = svcvt_f64_f32_x(predicate,v_at);
      v_sum1 = svadd_f64_x(predicate,v_sum1,v_d);
/* new position */
/*    dx = x + vx*dt; */
/*    dy = y + vy*dt; */
      v_dx = svadd_f32_x(predicate,v_x,svmul_f32_x(predicate,v_vx,v_dt));
      v_dy = svadd_f32_x(predicate,v_y,svmul_f32_x(predicate,v_vy,v_dt));
/* periodic boundary conditions */
      if (ipbc==1) {
/*       if (dx < edgelx) dx += edgerx; */
         svbool_t lt = svcmplt_n_f32(predicate,v_dx,edgelx);
         v_dx = svadd_n_f32_m(lt,v_dx,edgerx);
/*       if (dx >= edgerx) dx -= edgerx; */
         svbool_t geq = svcmpge_n_f32(predicate,v_dx,edgerx);
         v_dx = svsub_n_f32_m(geq,v_dx,edgerx);
/*       if (dy < edgely) dy += edgery; */
         svbool_t lt2 = svcmplt_n_f32(predicate,v_dy,edgely);
         v_dy = svadd_n_f32_m(lt2,v_dy,edgery);
/*       if (dy >= edgery) dy -= edgery; */
         svbool_t geq2 = svcmpge_n_f32(predicate,v_dy,edgery);
         v_dy = svsub_n_f32_m(geq2,v_dy,edgery);
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         svbool_t lt = svcmplt_n_f32(predicate,v_dx,edgelx);
         svbool_t beq = svcmpge_n_f32(predicate,v_dx,edgerx);
         svbool_t cond = svorr_b_z(predicate,lt,beq);
         v_dx = svsel_f32(cond,v_x,v_dx);
         v_vx = svsel_f32(cond,svneg_f32_x(predicate,v_vx), v_vx);
/*       if ((dy < edgely) || (dy >= edgery)) { */
/*          dy = y;                             */
/*          vy = -vy;                           */
/*       }                                      */
         svbool_t lt2 = svcmplt_n_f32(predicate,v_dy,edgely);
         svbool_t beq2 = svcmpge_n_f32(predicate,v_dy,edgery);
         svbool_t cond2 = svorr_b_z(predicate,lt2,beq2);
         v_dy = svsel_f32(cond2,v_y,v_dy);
         v_vy = svsel_f32(cond2,svneg_f32_x(predicate,v_vy), v_vy);
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
/*       if ((dx < edgelx) || (dx >= edgerx)) { */
/*          dx = x;                             */
/*          vx = -vx;                           */
/*       }                                      */
         svbool_t lt = svcmplt_n_f32(predicate,v_dx,edgelx);
         svbool_t beq = svcmpge_n_f32(predicate,v_dx,edgerx);
         svbool_t cond = svorr_b_z(predicate,lt,beq);
         v_dx = svsel_f32(cond,v_x,v_dx);
         v_vx = svsel_f32(cond,svneg_f32_x(predicate,v_vx), v_vx);
/*       if (dy < edgely) dy += edgery; */
         svbool_t lt2 = svcmplt_n_f32(predicate,v_dy,edgely);
         v_dy = svadd_n_f32_m(lt2,v_dy,edgery);
/*       if (dy >= edgery) dy -= edgery; */
         svbool_t beq2 = svcmpge_n_f32(predicate,v_dy,edgery);
         v_dy = svsub_n_f32_m(beq2,v_dy,edgery);
      }
/* set new position */
/*    part[j] = dx;     */
/*    part[j+npe] = dy; */
      svst1_f32(predicate,&part[j],v_dx);
      svst1_f32(predicate,&part[j+npe],v_dy);
/* set new velocity */
/*    part[j+2*npe] = vx; */
/*    part[j+3*npe] = vy; */
      svst1_f32(predicate,&part[j+2*npe],v_vx);
      svst1_f32(predicate,&part[j+3*npe],v_vy);
   }
   predicate = svptrue_b32();
   sum1 = svaddv_f64(predicate,v_sum1);
/* loop over remaining particles */
   for (j = nps; j < nop; j++) {
/* find interpolation weights */
      x = part[j];
      y = part[j+npe];
      nn = x;
      mm = y;
      dxp = x - (float) nn;
      dyp = y - (float) mm;
      nn = 2*(nn + nxv*mm);
      amx = 1.0f - dxp;
      amy = 1.0f - dyp;
/* find acceleration */
      dx = amx*fxy[nn];
      dy = amx*fxy[nn+1];
      dx = amy*(dxp*fxy[nn+2] + dx);
      dy = amy*(dxp*fxy[nn+3] + dy);
      nn += 2*nxv;
      vx = amx*fxy[nn];
      vy = amx*fxy[nn+1];
      dx += dyp*(dxp*fxy[nn+2] + vx);
      dy += dyp*(dxp*fxy[nn+3] + vy);
/* new velocity */
      dxp = part[j+2*npe];
      dyp = part[j+3*npe];
      vx = dxp + qtm*dx;
      vy = dyp + qtm*dy;
/* average kinetic energy */
      dxp += vx;
      dyp += vy;
      sum1 += dxp*dxp + dyp*dyp;
/* new position */
      dx = x + vx*dt;
      dy = y + vy*dt;
/* periodic boundary conditions */
      if (ipbc==1) {
         if (dx < edgelx) dx += edgerx;
         if (dx >= edgerx) dx -= edgerx;
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* reflecting boundary conditions */
      else if (ipbc==2) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if ((dy < edgely) || (dy >= edgery)) {
            dy = y;
            vy = -vy;
         }
      }
/* mixed reflecting/periodic boundary conditions */
      else if (ipbc==3) {
         if ((dx < edgelx) || (dx >= edgerx)) {
            dx = x;
            vx = -vx;
         }
         if (dy < edgely) dy += edgery;
         if (dy >= edgery) dy -= edgery;
      }
/* set new position */
      part[j] = dx;
      part[j+npe] = dy;
/* set new velocity */
      part[j+2*npe] = vx;
      part[j+3*npe] = vy;
   }
/* normalize kinetic energy */
/* *ek += 0.125f*sum1; */
   *ek += 0.125f*(sum1);
   return;
}

/*--------------------------------------------------------------------*/
void csvegpost2lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv) {

}

/*--------------------------------------------------------------------*/
void csvedsortp2ylt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1) {

}

/*--------------------------------------------------------------------*/
void csvecguard2l(float fxy[], int nx, int ny, int nxe, int nye) {

}

/*--------------------------------------------------------------------*/
void csveaguard2l(float q[], int nx, int ny, int nxe, int nye) {

}

/*--------------------------------------------------------------------*/
void csvepois22(float complex q[], float complex fxy[], int isign,
                 float complex ffc[], float ax, float ay, float affp,
                 float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
                 int nyhd) {

}

/*--------------------------------------------------------------------*/
void csvefft2rxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {

}

/*--------------------------------------------------------------------*/
void csvefft2rxy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int nxi, 
                 int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {

}

/*--------------------------------------------------------------------*/
void csvefft2r2x(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {

}

/*--------------------------------------------------------------------*/
void csvefft2r2y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {

}

/*--------------------------------------------------------------------*/
void csvewfft2rx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd) {

}

/*--------------------------------------------------------------------*/
void csvewfft2r2(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd) {

}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void csvexiscan2_(int *isdata, int *nths) {
   csvexiscan2(isdata,*nths);
   return;
}

/*--------------------------------------------------------------------*/
void csvegpush2lt_(float *part, float *fxy, float *qbm, float *dt,
                    float *ek, int *idimp, int *nop, int *npe, int *nx,
                    int *ny, int *nxv, int *nyv, int *ipbc) {
   csvegpush2lt(part,fxy,*qbm,*dt,ek,*idimp,*nop,*npe,*nx,*ny,*nxv,
                 *nyv,*ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void csvegpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
                    int *idimp, int *nxv, int *nyv) {
   csvegpost2lt(part,q,*qm,*nop,*npe,*idimp,*nxv,*nyv);
   return;
}

/*--------------------------------------------------------------------*/
void csvedsortp2ylt_(float *parta, float *partb, int *npic, int *idimp,
                      int *nop, int *npe, int *ny1) {
   csvedsortp2ylt(parta,partb,npic,*idimp,*nop,*npe,*ny1);
   return;
}

/*--------------------------------------------------------------------*/
void csvecguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye) {
   csvecguard2l(fxy,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csveaguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye) {
   csveaguard2l(q,*nx,*ny,*nxe,*nye);
   return;
}

/*--------------------------------------------------------------------*/
void csvepois22_(float complex *q, float complex *fxy, int *isign,
                  float complex *ffc, float *ax, float *ay, float *affp,
                  float *we, int *nx, int *ny, int *nxvh, int *nyv,
                  int *nxhd, int *nyhd) {
   csvepois22(q,fxy,*isign,ffc,*ax,*ay,*affp,we,*nx,*ny,*nxvh,*nyv,*nxhd,
               *nyhd);
   return;
}

void csvewfft2rx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csvewfft2rx(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void csvewfft2r2_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *nxhd,
                   int *nyd, int *nxhyd, int *nxyhd) {
   csvewfft2r2(f,*isign,mixup,sct,*indx,*indy,*nxhd,*nyd,*nxhyd,
                *nxyhd);
   return;
}
