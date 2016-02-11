choleskyalgo1 <- function(Y, threshold) {
    n = dim(Y)[1]
    p = dim(Y)[2]
    data = as.double(Y)
    out <- .C("choleskyalgo1", as.integer(n), as.integer(p), data, as.double(diag(p)) , as.double(threshold), as.double(diag(p)), as.double(diag(p)), as.double(diag(p)))
    return(out[[5]])
}
