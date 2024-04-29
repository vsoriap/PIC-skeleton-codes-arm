/* header file for ssepush2.c */

void csvexiscan2(int *isdata, int nths);

void csvegpush2lt(float part[], float fxy[], float qbm, float dt,
                   float *ek, int idimp, int nop, int npe, int nx,
                   int ny, int nxv, int nyv, int ipbc);

void csvegpost2lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv);

void csvedsortp2ylt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1);

void csvecguard2l(float fxy[], int nx, int ny, int nxe, int nye);

void csveaguard2l(float q[], int nx, int ny, int nxe, int nye);

void csvepois22(float complex q[], float complex fxy[], int isign,
                 float complex ffc[], float ax, float ay, float affp,
                 float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
                 int nyhd);

void csvefft2rxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csvefft2rxy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int nxi, 
                 int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csvefft2r2x(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csvefft2r2y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csvewfft2rx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);

void csvewfft2r2(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);
