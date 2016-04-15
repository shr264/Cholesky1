graphfiller <- function(adjMat){
(p = dim(adjMat)[1])
(filledadjacencyMatrix = adjMat)
for(i in p:1){
    for (jt in 1:(i-1)){
        j = i-jt
            for(k in 1:j){
        if((k<j)&&(j<i)){
            if((filledadjacencyMatrix[i,k]==1)&&(filledadjacencyMatrix[j,k]==1)){
                if(filledadjacencyMatrix[i,j]!=1){
                    filledadjacencyMatrix[i,j] <- TRUE
                    filledadjacencyMatrix[j,i] <- TRUE}}
        }
            }
        }
}
return(filledadjacencyMatrix)
}
