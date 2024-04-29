/* header file for svelib.c */

void sve_fallocate(float **s_f, int nsize, int *irc);

void sve_callocate(float complex **s_c, int nsize, int *irc);

void sve_iallocate(int **s_i, int nsize, int *irc);

void sve_deallocate(void *s_d);

void csveiscan2(int *isdata, int nths);
