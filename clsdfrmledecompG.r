clsdfrmlefordecompG <- function(S,adjmat,ldl){
(L = S)
(sqD = S)
(D = S)
L[] = 0
sqD[] = 0
D[] = 0
(p = dim(L)[1])

for (i in 1:p){(Snew = S[i:p,i:p])
               if(i<p){(newp = dim(Snew)[1])
                       (Sii = Snew[1,1])
                       (temp = which(adjmat[(i+1):p,i]>0))
                       (dimtemp = length(temp))
                       if(dimtemp>0){
                       (Stempi = S[(i+1):p,i])
                         (Sdoti = Stempi[temp])
                         (Stempgp = S[(i+1):p,(i+1):p])
                         if(length(Stempgp)>1){
                             (Sgp = Stempgp[temp,temp])}
                         else {Sgp = Stempgp}
                       (xi = -solve(Sgp)%*%Sdoti)
                       
                        (dii = sqrt(1/(Sii + t(Sdoti)%*%xi)))
                        
                       (L[i,i] = 1)
                       (L[(temp+i),i] = xi)
                       (sqD[i,i] = dii)
                       (D[i,i] = dii^2)}
                       else{
                       (dii = sqrt(1/(Sii)))
                       (L[i,i] = 1)
                       (L[(temp+i),i] = 0)
                       (sqD[i,i] = dii)
                       (D[i,i] = dii^2)}
                   }
                else {dii = sqrt(1/Snew)
                      L[i,i] = 1
                      sqD[i,i] = dii
                      D[i,i] = dii^2}}
if(missing(adjmat)){
    cat('Error!Missing Adjacency Matrix!')
    return(0)}
return(L%*%D%*%t(L))}
