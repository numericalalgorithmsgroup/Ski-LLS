// Author: Zhen Shao
#include "cs.h"
#include "lsex_all_in_c.h"
#include "Mat.hpp"
#include "Random.hpp"
#include "flapack.h"

typedef ptrdiff_t idx_t;
Mat_d reducedQR(Mat_d &A, double rcond, long* &E){
	// input A
	// output R, truncated R factor 
	// (input) output E, s.t. A(:, E) = QR, 
	// on input, E needs to have enough workspace
	column_pivoted_qr(A, E);
	long rank =0;
	for (idx_t i=0; i< A.n(); i++){
    	if (abs( A(i,i) ) > rcond)
    	{
    		rank++;
    	}
    	else
    	{
    		break;
    	}
    }
    Mat_d R;
    R = A.submat(0, rank, 0, rank);
    return R;
}