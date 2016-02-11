concentrationgraphChol <- function(Y,adjctMat,order,tol,standard){
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
                                        
                                        #Step2: fill-reducing order or user specified order

#ch <- chol(S)
#dd <- diag(ch)
#L <- t(ch/dd)
#DD <- dd^2
#S = L%*%t(L)

(S =  Matrix(c(S),ncol = p, nrow = p, sparse = TRUE))

if(missing(order)){
    cat('Missing order! Using Cholesky to find filled graph and order \n')
    (orderS = orderchol(S,permute=TRUE))
    (S = orderS$ordered)
    (filledadjacencyMatrix = orderS$filledg)
}
else {
    cat('Order Provided! \n')
    S = S[order,order]
}

S = as.matrix(S)
                                        #Step4: Threshold to find graph
                                        #aka: Make the adjacency matrix
if(missing(adjctMat)){
    if(missing(tol)){
        cat('Missing Adjacency Matrix. Thresholding using 0.1\n')
        (adjctMat = threshS(S,0.1))
    }
    else {
        cat('Missing Adjacency Matrix. Thresholding using',tol,'\n')
        (adjctMat = threshS(S,tol))
    }}


 #Step6: Fill the graph
                                    
if(!(missing(order))){                                       
(filledadjacencyMatrix = graphfiller(adjctMat))
}

if(maxCsize>n){
        cat('Maximal Clique is larger than n. Algorithm will probably fail!\n')}

                                        #step7: calculate lower triangular matrix Ld for filled graph
(Ld1 = cholcalc(S,filledadjacencyMatrix,ldl=FALSE)$Lower)

                                        #Step8: implement algorithm 1.
(E1 = medgeset(adjctMat))
if(dim(E1)[1]==0){
    cat('No missing edges.Calculating Ld\n')
    L = Ld}
else{cat('Missing edges, Using Algo 1!\n')
     (L1obj = choladjest(E1,Ld1))
     (L1 = L1obj$L)}
(omegahat1 = L1%*%t(L1))
time2 = proc.time()
timeelapsed = time2 - time1
return(list(omegahat=omegahat1,
            adjMat = adjctMat,
            filledadjMat = filledadjacencyMatrix,
            Lowerd = Ld1,
            Lower = L1,time = timeelapsed,count = L1obj$count))}
