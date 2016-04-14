choleskyalgo1 <- function(Y, adjMat) {
    n = dim(Y)[1];
    p = dim(Y)[2];
    data = as.double(Y);
    out <- .C("choleskyalgo1adj", as.integer(n), as.integer(p), data, as.double(diag(p)) , as.integer(adjMat), as.double(diag(p)), as.double(diag(p)), as.double(diag(p)));
    return(out);
}
