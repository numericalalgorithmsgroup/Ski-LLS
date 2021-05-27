# Overview

Ski-LLS (SKetchIng for Linear Least Squares) is a C++ package for finding solutions to large-scale over-determined linear least square problems using modern dimensionality reduction techniques (sketching). Specifically, given a matrix A with n rows and d columns, a vector b with n rows, it finds a vector x that minimises the 2-norm of the vector Ax-b. 

Ski-LLS is faster and more robust than Blendenpik and LAPACK on large over-determined data matrices, e.g. matrices having 40,000 rows and 4,000 columns. Ski-LLS is 10 times faster than Sparse QR and incomplete-Cholesky preconditioned LSQR on large, ill-conditioned sparse data matrices, e.g. matrices with 120,000 rows and 5,000 columns with 1% non-zeros. 

Please see doc/Doc.pdf for more detailed advantages of Ski-LLS, the installation instructions and the C++ library interface. 

# License

Ski-LLS is available under the BSD license.

# Reference

Hashing embeddings of optimal dimension, with applications to linear least squares. C. Cartis, J. Fiala and Z. Shao, Optimization-online preprint, 2021. 
http://www.optimization-online.org/DB_HTML/2021/05/8413.html


Sparse sketching for sparse linear least squares. C. Cartis, J. Fiala and Z. Shao, ICML Workshop “Beyond first order methods for ML systems”, ICML2020.
https://drive.google.com/file/d/1BacyZwtZSKZBBblLC7x4SikuDQaLHuYK/view





 
 


