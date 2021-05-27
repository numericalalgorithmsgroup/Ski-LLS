#include <iostream>
#include "ski-lls.h" 

int main(int argc, char **argv)
{
	// create data
	long m = 4;
	long n = 2;
	double data_A[8] = {0.306051,-1.53078,1.64493,-1.61322,
	-0.2829,0.474476,-0.586278,-0.610202}; 

	double data_b[4] = {0.0649338,0.845946,-0.0164085,0.247119};

	Mat_d A(m,n,data_A);
	Vec_d b(m,data_b);

	// parameters
	long k = NNZ_PER_COLUMN;
	double gamma = OVER_SAMPLING_RATIO;
	long max_it = MAX_IT;
	double it_tol = IT_TOL;
	double rcond = 1e-10;
	int wisdom = 0;
	int debug = 0;
	double abs_tol = ABS_TOL;

	// storage

	long it;
	int flag;
	long rank;
	double residual;   
	Vec_d x(A.n());

	// solve
	ls_dense_hashing_blendenpik(
	A, b, x, rank, flag, it, gamma, k, abs_tol, 
	rcond, it_tol, max_it, debug, wisdom);

	std::cout << "A is: "<< A << std::endl;
	std::cout << "b is: "<< b << std::endl;
	std::cout << "solution is: "<< x << std::endl;

	// compute residual
	A.mv('n', 1, x, -1, b); 
	std::cout << "Residual of ls_blendenpik_hashing is: "<<  nrm2(b) << std::endl;
	std::cout << "Iteration of ls_blendenpik_hashing is: " << it << std::endl;
	// The correct value is x = (-0.2122, 0.720)

}