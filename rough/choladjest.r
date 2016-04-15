### Function to adjust a given Cholesky matrix $L$ so that $LL^t$ 
### reflects zeros based on the ordered missing edge set $E$
choladjest <- function(E,Ld){
(Lnew <- as.matrix(Ld))
(a <- length(E[,1]))
###i = 1
###i = i+1
cnt=0
for(i in 1:a){
    if(E[i,2]>1){
        (Lnew[E[i,1],E[i,2]] = -sum(Lnew[E[i,1],1:(E[i,2]-1)]*Lnew[E[i,2],1:(E[i,2]-1)]))
        if((is.nan((Lnew[E[i,1],E[i,2]])))||(is.infinite((Lnew[E[i,1],E[i,2]])))){
            cnt = cnt + 1
            return(list(L = NULL,count = cnt))}

        (Lnew[E[i,1],E[i,2]] = (Lnew[E[i,1],E[i,2]])/(Lnew[E[i,2],E[i,2]]))
        if((is.nan((Lnew[E[i,1],E[i,2]])))||(is.infinite((Lnew[E[i,1],E[i,2]])))){
            cnt = cnt + 1
            return(list(L = NULL,count = cnt))}
}
    else{       
        Lnew[E[i,1],E[i,2]] = 0}}
return(list(L = Lnew,count = cnt))}
