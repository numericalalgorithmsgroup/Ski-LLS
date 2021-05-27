/*
Copyright (c) 2009-2016, Haim Avron and Sivan Toledo.
Extended by Zhen Shao. 
Adapted from Avron's implementation in the Blendenpik package
Perform discrete cosin transform or discrete hartley transform from fftw 
*/

#include "blendenpik.hpp"
#include "Config.hpp"
#include "Mat.hpp"
#include "fftw_r2r.hpp"
#include <ostream>
#include <string.h>

int dct_sizes[FFTW_TIMES];
int dht_sizes[FFTW_TIMES];

void load_fftw_wisdom()
{
  FILE *wisdom_file;
  char *wisdom_string;
  long size, start;

  wisdom_file = fopen(
    FFTW_WISDOM_FILE, /* path to the FFTW wisdom file, defined in include/Config.hpp*/
    "r");
  if (wisdom_file == NULL)
  	std::cout << "Could not find wisdom file" << std::endl;
  fread(dht_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  fread(dct_sizes, sizeof(int), FFTW_TIMES, wisdom_file);
  start = ftell(wisdom_file);
  fseek(wisdom_file, 0, SEEK_END);
  size = ftell(wisdom_file) - start;
  fseek(wisdom_file, start, SEEK_SET);
  wisdom_string = (char*)malloc(sizeof(char) * size);
  fread(wisdom_string, 1, size, wisdom_file);
  fclose(wisdom_file);

  fftw_import_wisdom_from_string(wisdom_string);
}

Mat_d fast_unitary_transform(const Mat_d &A /* matrix m_0 by d */, 
	double* D /* array of size m_0 */, int wisdom){
	
	int debug=0;
	long mm, m0, i, j;
  int m,n;
	double* y;
	double scale;
	unsigned kind= FFTW_R2R_DHT;
  // unsigned kind= FFTW_R2R_DCT;

	int wis_level = FFTW_ESTIMATE;

	m0 = A.m();
	n = A.n();
	// x = A.data();

  if (kind == FFTW_R2R_DCT) {
    if (wisdom==1) {
      load_fftw_wisdom();
      wis_level = FFTW_WISDOM_FLAG;
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dct_sizes[mm / FFTW_QUANT - 1];
    } else
      m = m0;
    scale = 1 / sqrt(2 * m);

  }

  if (kind == FFTW_R2R_DHT) {
    if (wisdom==1) {
      load_fftw_wisdom();
      wis_level = FFTW_WISDOM_FLAG;
      mm = (m0 % FFTW_QUANT == 0) ? m0 : m0 + (FFTW_QUANT - (m0 % FFTW_QUANT));
      m = dht_sizes[mm / FFTW_QUANT - 1];
    } else
      m = m0;
    scale = 1 / sqrt(m);
  }

    y = (double*)malloc(m*n*sizeof(double));

    if (debug==1){
    	std::cout << "m:" << m << ", n:"<<n << std::endl;
    }

	for(i = 0; i < n; i++)/* go through columns*/ {
		for(j = 0; j < m0; j++)
		{
			y[m * i + j] = scale * D[j]* A(j,i); /* remember to add sign randomization here*/
    	  	if (debug==1){
    	  		std::cout << "i:"<< i << ", j:"<<j <<", y(i,j)"<<y[m*i+j] <<", A(j,i)" << A(j,i)
    	  			<< ", scale:" << scale << ", D[j]"<< D[j]<<std::endl;
    	  	}
		} /* go through rows*/

		if (m0 != m) /* padding with zeros */
		  memset(y + m * i + m0, 0, (m - m0) * sizeof(double));
	}

	// if (debug==1){
	// 	for (i =0; i < m*n; i++){
	// 		std::cout << y[i] <<std::endl;
	// 	}
	// }

	int r;
    r = fftw_r2r(y, y, m, n, kind, wis_level); /* in place transformation */
    if (!r){
      std::cout << "FFTW wisdom not loaded" << std::endl;
    }

  	if (kind == FFTW_R2R_DCT){
	    for (i = 0; i < n; i++)
	      y[i * m] /= sqrt(2);
  	}

  	Mat_d Y(m, n, (double*)y);
    return Y;
}