covargraphBanded <- function(Y,adjctMat,nbands,standard){
    time1 = proc.time()
    if(missing(nbands)){
        cat('Argument nbands missing. using nbands = 1!\n')
        nbands = 1}
                                        #Step1 of Concgraph. we calculate sample covariance S
#maybe standardize Y first, but this didnt work
if(missing(standard)){Y = Y}
else{if(standard==1){Y = scale(Y, center = TRUE, scale = TRUE)}
    else{Y = Y}}

cat('initial step: calculate S\n')
(n = dim(Y)[1])
(p = dim(Y)[2])
(S = t(Y)%*%Y/n)
                                     
cat('Using Algo 1 for Banded Matrices!\n')
cat('Calculating Ld\n')
(Ld1 = t(chol(S)))
(L1 = band(Ld1, -nbands, 0))

(omegahat = as.matrix(L1%*%t(L1)))
if(missing(adjctMat)){adjctMat=omegahat!=0}

time2 = proc.time()
timeelapsed = time2 - time1
return(list(omegahat=omegahat,
            adjMat = adjctMat,
            Lowerd = Ld1,
            Lower = L1,time = timeelapsed))}
