### Function to convert an adjacency matrix to an ordered MISSING
### edge set 
medgeset <- function(A){
##A = test4$adjMat
##A = adjMat
(p = length(A[,1]))
(E = matrix(0,0,2))
for(i in 2:p){
    (zero = which(A[i,1:i]==0))
    if(length(zero)>0){
        (flintrack = cbind(rep(i,length(zero)),zero))
        (E = rbind(E,flintrack))}}
return(E)}
