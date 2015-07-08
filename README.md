# DIRAC_Rcode
This is the R code for DIRAC algorithm (Eddy et al., 2010). Permutation
function has been modified as compared to the original DIRAC algorithm.
permutation is performed at the individual pathway level and then
pvalue is corrected using Benjamini Hochberg correction.  This code
also parallelizes the permutation testing to make the running of the
code fast(cores = 16).

This code has the input as a expression matrix with rownames as human gene names and columns as the samples.
First row of the matrix represents the classes. For example, in the "Main run" there are two conditions "C9" and "T9"


