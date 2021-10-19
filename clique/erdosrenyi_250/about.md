Specifications
--------------

* 250x250 dissimilarity matrix
* each off diagonal entry is an integer between 1 and N = (500^2 - 500)/2
* entries that lie strictly below the diagonal were filled using a random permutation of the integers 1, ..., N
* thus, given any two pairs of distinct integers (i,j) and (k,l), it is equally probably that M[i,j] > M[k,l] or M[i,j] < M[k,l]
* matrix was generated using `erdosrenyi.jl`
