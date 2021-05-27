#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/timeb.h>
#include <iostream>

#include "Config.hpp"
#include "fftw3.h"


double read_wtime() 
{
  struct timeb T;
  (void) ftime( &T );
  
  return ((double)T.time) + (1e-3) * ((double)T.millitm);
}


/*
Build fftw wisdom to be used in Blendenpik-like variants
level 0 = lowest level of fftw calibration, takes about 20 seconds to run
level 3 = highest level of fftw calibration, could take hours
After the wisdom is built, put the wisdom path (argv[2]) into include/Config.hpp FFTW_WISDOM_FILE
If a level other than 0 is used, modify FFTW_WISDOM_FLAG as well 
*/
int main(int argc, char **argv /*
                                      argv[1], level (0-3)
                                      argv[2],  wisdom file name */) 
{
  int n, m, i;
  // mxArray *XX;
  double *X;
  fftw_plan plan;
  char type[4];
  double t0, tt0, elp, t_all0;
  char wisdom_file_name[255];
  FILE *wisdom_file;
  unsigned level_flag;
  int level;
  int tm;

  int dct_sizes[FFTW_TIMES];
  double dct_times[FFTW_TIMES];
  int dht_sizes[FFTW_TIMES];
  double dht_times[FFTW_TIMES];

  t_all0 = read_wtime();

  level = (int)(atoi(argv[1]));
  switch(level) {
  case 0:
    level_flag = FFTW_ESTIMATE;
    break;

  case 1:
    level_flag = FFTW_MEASURE;
    break;

  case 2:
    level_flag = FFTW_PATIENT;
    break;
			
  case 3:
    level_flag = FFTW_EXHAUSTIVE;
    break;
  }

  n = 10;


  for(m = FFTW_TIMES * FFTW_QUANT; m >= FFTW_QUANT; m -= FFTW_QUANT) {
    tm = m / FFTW_QUANT - 1;
    t0 = read_wtime();
    // XX = mxCreateDoubleMatrix(m, n, mxREAL);
    X = (double*) malloc(m*n*sizeof(double));
    // mxGetPr(XX);
    memset(X, 1, m * n * sizeof(double));
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_DHT, FFTW_UNALIGNED | level_flag);
    fftw_execute_r2r(plan, X, X); /* execute one time for estimate */
    fftw_destroy_plan(plan);
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_DHT, FFTW_UNALIGNED | level_flag); /* replan for estimate */
    tt0 = read_wtime();
    for(i = 0; i < n; i++)
      fftw_execute_r2r(plan, X + i * m, X + i * m);
    elp = read_wtime() - tt0;
    dht_times[tm] = elp;
    dht_sizes[tm] = m;
    fftw_destroy_plan(plan);
    free(X);

    if (tm < FFTW_TIMES-1 && dht_times[tm] > dht_times[tm + 1]) {
      dht_times[tm] = dht_times[tm + 1];
      dht_sizes[tm] = dht_sizes[tm + 1];
    } 
    printf("DHT size %d->%d run time %.2e planning time is %.2e sec\n", m, dht_sizes[tm], elp, read_wtime() - t0);
  }

  for(m = FFTW_TIMES * FFTW_QUANT; m >= FFTW_QUANT; m -= FFTW_QUANT) {
    tm = m / FFTW_QUANT - 1;
    t0 = read_wtime();
    // XX = mxCreateDoubleMatrix(m, n, mxREAL);
    // X = mxGetPr(XX);
    X = (double*) malloc(m*n*sizeof(double));
    memset(X, 1, m * n * sizeof(double));
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_REDFT10, FFTW_UNALIGNED | level_flag);
    fftw_execute_r2r(plan, X, X); /* execute one time for estimate */
    fftw_destroy_plan(plan);
    plan = fftw_plan_r2r_1d(m, X, X, FFTW_REDFT10, FFTW_UNALIGNED | level_flag);
    tt0 = read_wtime();
    for(i = 0; i < n; i++)
      fftw_execute_r2r(plan, X + i * m, X + i * m);
    elp = read_wtime() - tt0;
    dct_times[tm] = elp;
    dct_sizes[tm] = m;
    fftw_destroy_plan(plan);
    free(X);

    if (tm < FFTW_TIMES-1 && dct_times[tm] > dct_times[tm + 1]) {
      dct_times[tm] = dct_times[tm + 1];
      dct_sizes[tm] = dct_sizes[tm + 1];
    } 
    printf("DCT size %d->%d run time %.2e planning time is %.2e sec\n", m, dct_sizes[tm], elp, read_wtime() - t0);
  }

  /* Write wisdom to file. FFTW's function doesn't work for some reason... */
  // mxGetString(argv[2], wisdom_file_name, 255);	
  wisdom_file = fopen(argv[2], "w+");
  if (wisdom_file == NULL)
    std::cout << "Cannot open wisdom file" << std::endl; 
  fwrite(dht_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fwrite(dct_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fputs(fftw_export_wisdom_to_string(), wisdom_file);
  fclose(wisdom_file);

  printf("TOTAL time is %.2e sec\n", read_wtime() - t_all0); 
}