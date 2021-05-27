#!/bin/bash
# Test all dense solvers on full rank dense problems
# 	   all sparse solvers on full rank sparse problems
# 	   selected solvers on rank-deficient sparse problems

matDir=./data
rankDefMatDir=./data/rank_deficient_sparse_paper
solverDir=./EXE


# dense solvers and their inputs
declare -a denseSolvers=("ls_dense_lsrn.exe" "ls_blendenpik.exe" \
	"ls_blendenpik_hashing.exe" "ls_lapack_qr.exe")
declare -a denseMatrices=("50by15DenseFullRank.mtx" "50by20DenseFullRank.mtx" "100by15DenseFullRank.mtx" \
	"100by20DenseFullRank.mtx" "120by15DenseFullRank.mtx" "120by20DenseFullRank.mtx" )
# declare -a denseSolns= //To do, get the exact solution using SVD (or can we just take LAPACK QR?)

# sparse solvers and their inputs
declare -a sparseSolvers=("ls_qr.exe" "ls_lsrn.exe" "ls_spqr_solve.exe"\
 "ls_blendenpik_hashing_sparse_input.exe" "ls_blendenpik_sparse_input.exe")
declare -a sparseMatrices=("50by15SparseFullRank.mtx" "50by20SparseFullRank.mtx" "100by15SparseFullRank.mtx" \
	"100by20SparseFullRank.mtx" "120by15SparseFullRank.mtx" "120by20SparseFullRank.mtx" )
# declare -a sparseSolns= // To do, get te exact solution using SVD (or can we just take LAPACK QR?)

# rank-deficient problem solvers and their inputs
declare -a rankDefSolvers=("ls_qr.exe" "ls_blendenpik_hashing_sparse_input.exe" \
	 "ls_lsrn.exe" "ls_spqr_solve.exe")
declare -a rankDefSolns=("18.3360" "26.5028" "50.8747" "6.12592e-14" "33.2357" "14.0199" "7.31938" "8.99553e-15" 
	"4.23281" "3.46410" "1.82574" "30.9959")
declare -a rankDefMatrices=("lp_ship12l.mtx" "Franz1.mtx" "GL7d26.mtx" "cis-n4c6-b2.mtx" "lp_modszk1.mtx" 
	"rel5.mtx" "ch5-5-b1.mtx" "n3c5-b2.mtx" "ch4-4-b1.mtx" "n3c5-b1.mtx" "n3c4-b1.mtx" "cis-n4c6-b3.mtx")

# Supporting functions
test_dense_solvers(){
	# first argument matrix to test
	for solver in "${denseSolvers[@]}"
	do
		run_and_print_res $solver $1 
	done
}

test_sparse_solvers(){
	# first argument matrix to test
	for solver in "${sparseSolvers[@]}"
	do
		run_and_print_res $solver $1 
	done
}

test_rankDef_solvers(){
	# first argument matrix to test
	for solver in "${rankDefSolvers[@]}"
	do
		run_and_print_res_rankDef $solver $1 
	done
}

run_and_print_res(){
	# first argument solver name; second argument matrix name
	res=$( (ulimit -t 10;$solverDir/$1 $matDir/$2 |grep Residual |awk '{print $NF}'))
	echo $1 gives a solution with residual $res
}

run_and_print_res_rankDef(){
	# first argument solver name; second argument matrix name
	res=$( (ulimit -t 10;$solverDir/$1 $rankDefMatDir/$2 |grep Residual |awk '{print $NF}'))
	echo $1 gives a solution with residual $res
}


# Top layer routine for testing dense solvers
test_dense_solvers_allMat() 
{
	for (( i=0; i<"${#denseMatrices[@]}"; i++));
	do
		echo TESTING ${denseMatrices[$i]}
		# echo We should get ${solns[$i]}
		test_dense_solvers ${denseMatrices[$i]}
		echo
	done
}

# Top layer routine for testing sparse solvers
test_sparse_solvers_allMat() 
{
	for (( i=0; i<"${#sparseMatrices[@]}"; i++));
	do
		echo TESTING ${sparseMatrices[$i]}
		# echo We should get ${solns[$i]}
		test_sparse_solvers ${sparseMatrices[$i]}
		echo
	done
}

# Top layer routine for testing rank deficient inputs
test_rankDef_solvers_allMat() 
{
	for (( i=0; i<"${#rankDefMatrices[@]}"; i++));
	do
		echo TESTING ${rankDefMatrices[$i]}
		echo We should get ${rankDefSolns[$i]}
		test_rankDef_solvers ${rankDefMatrices[$i]}
		echo
	done
}


test_dense_solvers_allMat

test_sparse_solvers_allMat

test_rankDef_solvers_allMat

