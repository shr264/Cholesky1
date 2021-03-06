covargraph3 <- function(Y,adjctMat,order,tol,standard){
time1 = proc.time()
                                        #Step1 of Concgraph. we calculate sample covariance S
#maybe standardize Y first, but this didnt work
if(missing(standard)){Y = Y}
else{if(standard==1){Y = scale(Y, center = TRUE, scale = TRUE)}
    else{Y = Y}}

cat('initial step: calculate S\n')
(n = dim(Y)[1])
(p = dim(Y)[2])
(S = t(Y)%*%Y/n)
                                        #Step4: Threshold to find graph
                                        #aka: Make the adjacency matrix
if(missing(adjctMat)){
    if(missing(tol)){
        cat('Missing Adjacency Matrix. Thresholding using 0.1\n')
        (adjctMat = threshS(S,0.1))
    }
    else{
        cat('Missing Adjacency Matrix. Thresholding using',tol,'\n')
        (adjctMat = threshS(S,tol))
    }}

(E1 = medgeset(adjctMat))
if(dim(E1)[1]==0){
    cat('No missing edges.Calculating Ld\n')
    (Ld1 = cholcalc(S,adjctMat,ldl=FALSE)$Lower)  
    L1 = Ld1}
else{
    cat('Missing edges, Using Algo 1!\n')
    cat('Calculating Ld\n')
    (Ld1 = t(chol(S)))
    (L1obj = choladjest(E1,Ld1))
    (L1 = L1obj$L)}

(omegahat = L1%*%t(L1))
                                        #return to original ordering
time2 = proc.time()
timeelapsed = time2 - time1
return(list(omegahat=omegahat,
            adjMat = adjctMat,
            Lowerd = Ld1,
            Lower = L1,time = timeelapsed, count = L1obj$count))}
