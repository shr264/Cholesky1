                                        #Thresholding: Decide on a value such that values on the off-diagonal smaller than it are set to 0
threshS <- function(S,tol){
(S = as.matrix(S))
(T = cov2cor(S))
(p  = dim(S)[1])
for (i in 1:(p-1)){
    (temp = which(abs(T[row(T) == (col(T) - i)])<tol))
    S[row(S) == (col(S) - i)][temp] = 0
    S[row(S) == (col(S) + i)][temp] = 0
}
(adjMat = (abs(S)>0))
return(adjMat)}
