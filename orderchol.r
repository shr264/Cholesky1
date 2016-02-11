orderchol <- function(sig,permute){
(p = dim(sig)[1])
(A = Cholesky(sig,perm = permute))
(filledgraph = (abs(expand(A)$L)>0)+t(abs(expand(A)$L)>0)>0)
(order = A@perm)
return(list(ordered = sig[(order+1),(order+1)],filledg = filledgraph))}
