## Big G problem
These codes are to deal with very large dataset, e.g., if you have 1M SNP panel and 50k+ ID.
Typically to construct a G matrix for GBLUP with ordinary GEMM,
i.e., several line coding in many languages,
it might take months to finish.

Here I:
1. divide genotypes into several ID groups.
2. calculate the G blocks one by one with MKL lapack functions.
3. merge the blocks into one big G matrix
4. calculate G inverse

Of 64k ID and 600k SNP, it takes 1+ days for G and around 30m to calculate G^-1.
