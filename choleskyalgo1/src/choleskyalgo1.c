/*
 
 cholesky code. to compile use: g++ -o cholesky1 cholesky1.cpp -llapack -lblas 
 
 */
#include <sys/param.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include "choleskyalgo1.h"




void inverse(double* A, int N)
{
    int* IPIV;
    IPIV = (int *)malloc((N+1)*sizeof(int));
    
    int LWORK = N*N;
    
    double* WORK;
    WORK = (double *)malloc(LWORK*sizeof(double));
    
    int INFO;
    
    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

void initvec(double* v, int N){
    for(int i= 0; i<N; ++i){
        v[i]= 0.0;
        }
   }

void matvecprod(double* A, double* v, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, lda= N, incx= 1, incy= 1;
    double* tmp;
    tmp = (double *)malloc(N*sizeof(double));
    initvec(tmp, N);
    dgemv_(&no,&m,&n,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    free(tmp);
}

void transmatvecprod(double* A, double* v, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, lda= N, incx= 1, incy= 1;
    double* tmp;
    tmp = (double *)malloc(N*sizeof(double));
    initvec(tmp, N);
    dgemv_(&tr,&m,&n,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    free(tmp);
}

void vecmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= 1, k= N, lda= N, ldb= N, ldc= N;
    double* tmp;
    tmp = (double *)malloc(N*sizeof(double));
    initvec(tmp, N);
    dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    free(tmp);
}

void matmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, k= N, lda= N, ldb= N, ldc= N;
    double* tmp;
    tmp = (double *)malloc(N*sizeof(double));
    initvec(tmp, N*N);
    dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N*N; ++i){
        u[i]= tmp[i];
    }
    free(tmp);
}

void mattransmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, k= N, lda= N, ldb= N, ldc= N;
    double* tmp;
    tmp = (double *)malloc(N*N*sizeof(double));
    initvec(tmp, N*N);
    dgemm_(&no,&tr,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N*N; ++i){
        u[i]= tmp[i];
    }
    free(tmp);
}


inline int getIndexOfLargestElement(double arr[], int size) {
    int largestIndex = 0;
    for (int index = largestIndex; index < size; index++) {
        if (arr[largestIndex] < arr[index]) {
            largestIndex = index;
        }
    }
    return largestIndex;
}

inline int getIndexOfSmallestElement(double arr[], int size) {
    int largestIndex = 0;
    for (int index = largestIndex; index < size; index++) {
        if (arr[largestIndex] > arr[index]) {
            largestIndex = index;
         }
    }
    return largestIndex;
}

void choleskyalgo1 (int *nIn, int *pIn, double *data, double *SS, double *thresholdIn, double *L, double *D, double *LD) {
  int p = *pIn;
  int n = *nIn;
  double Ssum = 0;
    double threshold = *thresholdIn;
    double S[p][p];
  double Correlation[p][p];
  double newS[p][p];
  double adjMat[p][p];
  double fadjMat[p][p];
  double invS[p*p];
  double inverseS[p][p];
    int medgeset[p][2];
    int medgesetlen = 0;
    int Sswaps[p][2];

  memset(L, 0, sizeof(double) * p * p);
  memset(D, 0, sizeof(double) * p * p);
    memset(LD, 0, sizeof(double) * p * p);
    memset(Sswaps, 0, sizeof(int) * p*2);

   // first we calculate S

  for(int i = 0; i<p; i++){
    for(int j = 0; j<p; j++){
      for (int k = 0; k < n; k++){
	Ssum = Ssum + data[k*p+i]*data[k*p+j];
      }
      S[i][j] = Ssum/n;
        Rprintf("\n S at %d,%d is %f", i , j, S[i][j]);
      Ssum = 0;
    }
  }

    // now we convert the covariance matrix to a correlation matrix

    for(int i = 0; i<p; i++){
      for(int j = 0; j<p; j++){
	Correlation[i][j] = S[i][j]/(pow(S[i][i],0.5)*pow(S[j][j],0.5));
          Rprintf("\n Correlation at %d,%d is %f", i , j, Correlation[i][j]);
    }
  }


    for(int i = 0; i<p; i++){
      for(int j = 0; j<p; j++){
	if(fabs(Correlation[i][j])>threshold){
	  newS[i][j] = S[i][j];
        adjMat[i][j] = 1;
        fadjMat[i][j] = 1;
	}
	else{
        Rprintf("\n Thresholding at %d , %d", i , j);
	  newS[i][j] = 0;
        adjMat[i][j] = 0;
        fadjMat[i][j] = 0;
	}
    }
  }
    
    //this is the sorting portion
 
    
    for(int i = 0; i < p; i++){
        Rprintf("\n Beginning Sorting at %d",i);
        double *zeros;
        zeros = (double *)malloc((p-i)*sizeof(double));
        
        memset(zeros, 0, sizeof(double) * (p-i));
        
        
    for(int c=i; c<p ;c++){
            int zerosum = 0;
            for(int d=i; d<p ;d++){
                if(adjMat[c][d]==0){
                    zerosum += 1;
                
                }
            }
        zeros[c-i] = zerosum;
        }
        
        int max_pos = i;
        
        for (int index = max_pos; index < p; index++) {
            if (zeros[max_pos] < zeros[index]) {
                max_pos= index;
            }
        }
 
        if(max_pos!=i){
            Sswaps[i][0] = i;
            Sswaps[i][1] = max_pos;
            Rprintf("\n Swapping row/column %d with %d",i, max_pos);
            // column-wise swaps
            for(int c = 0; c < p; c++){
                double temp = S[c][i];
                S[c][i] = S[c][max_pos];
                S[c][max_pos] = temp;
                
                double adjtemp = adjMat[c][i];
                adjMat[c][i] = adjMat[c][max_pos];
                adjMat[c][max_pos] = adjtemp;
                
                double fadjtemp = fadjMat[c][i];
                fadjMat[c][i] = fadjMat[c][max_pos];
                fadjMat[c][max_pos] = fadjtemp;
            }

            // row-wise swaps
            for(int c = 0; c < p; c++){
                double temp = S[i][c];
                S[i][c] = S[max_pos][c];
                S[max_pos][c] = temp;
                
                double adjtemp = adjMat[i][c];
                adjMat[i][c] = adjMat[max_pos][c];
                adjMat[max_pos][c]= adjtemp;
                
                double fadjtemp = fadjMat[i][c];
                fadjMat[i][c] = fadjMat[max_pos][c];
                fadjMat[max_pos][c]= fadjtemp;
                
            }
        }
        free(zeros);
    }
    
    
    //filling the graph
    Rprintf("\n Filling graph");
    for(int c = (p-1); c >= 0 ; c--){
        for (int j = (c-1); j >= 0; j--){
            for(int k = 0; k<j; k++){
                if((k<j)&&(j<c)){
                    if((fadjMat[c][k]==1)&&(fadjMat[j][k]==1)){
                        if(fadjMat[c][j]!=1){
                            fadjMat[c][j] = 1;
                            fadjMat[j][c] = 1;
                            medgeset[medgesetlen][0] = c;
                            medgeset[medgesetlen][1] = j;
                            medgesetlen +=1;
                        }}
                }
            }
        }
    }
    
    //choilesky calclulation
    
    
    L[(p-1)*p+(p-1)] = 1;
    D[(p-1)*p+(p-1)] = sqrt(1/S[p-1][p-1]);
    
    for (int i1 = 0; i1 < (p-1); i1++){
        Rprintf("\n Cholesky at i = %d",i1);
        int i = i1;
        double Sii = S[i][i];
        double temp;
        
        int counter = 0;
        
        for (int c = (i+1); c < p; c++){
            if(fadjMat[i][c]>0){
                counter += 1;
            }
        }
        
        Rprintf("\n Number of non-zero entries is %d",counter);
        
        if(counter==0){
            double dii = pow(1/(Sii),0.5);
            L[i*p+i] = 1;
            D[i*p+i] = dii;
        }
        else{
            Rprintf("Counter Non-zero");
            double *Sdoti;
            Sdoti = (double *)malloc(counter*sizeof(double));
            
            double *Sgp;
            Sgp = (double *)malloc(counter*counter*sizeof(double));
            
            double *u;
            u = (double *)malloc(counter*sizeof(double));
            
            memset(u, 0, sizeof(double) * counter);
            
            Rprintf("Memory Allocated");
            
            int count1 = 0;
        
            for (int c = (i+1); c < p; c++){
        
                if(fadjMat[i][c]>0){
                    Sdoti[count1] = S[i][c];
                    
                    int count2 = 0;
                    
                    for(int d = (i+1); d < p; d++){
                        
                        if(fadjMat[d][i]>0){
                            
                            Sgp[count2 + counter*count1]=S[d][c];

                            if(count2<counter){
                                                                count2 += 1;}
                        }
                        
                    }
                   
                    count1 += 1;
                    
                }
                
            }
 
            inverse(Sgp,counter);
            
            matvecprod(Sgp,Sdoti,u,counter);
            
            double tempprod = 0;
            
            for(int c = 0; c < counter; c++){
                tempprod = tempprod - u[c]*Sdoti[c];
            }
            
            double dii;
            dii = sqrt(1/(Sii + tempprod));
            
            L[i*p+i] = 1;
            int count3 = 0;
            for(int c = (i+1); c < p ; c++){
                if(fadjMat[i][c]>0){
                    L[i*p+c] = -u[count3];
                    count3 += 1;
                }
            }
            D[i*p+i] = dii;
            
            free(Sdoti);
            free(Sgp);
            free(u);
            
            Rprintf("Memory Freed");
        }
    }
    //end cholesky
    
    //multiply L by D
    
    Rprintf("\n Multiply L by D");
    matmatprod(D,L,LD,p);
    
    
    //this is algorithm 1 from the paper
    
    Rprintf("\n Beginning Algoritm 1");
    for (int c = 0; c < medgesetlen; c++){
        if(medgeset[c][1]>1){
            double Lsum = 0;
            for(int k = 0; k < (medgeset[c][1]-1); k++){
                Lsum += LD[medgeset[c][0]*p+k]*LD[medgeset[c][1]*p+k];
            }
            LD[medgeset[c][0]*p+medgeset[c][1]] = -Lsum/LD[medgeset[c][1]*p+medgeset[c][1]];
        }
        else{
            LD[medgeset[c][0]*p+medgeset[c][1]] = 0;
        }
    }
    
    // end algorothim 1
        
    //swap back to original order
    
     Rprintf("\n Swapping back to original order");
    for(int i = (p-1); i>=0;i--){
         // row-wise swaps
        for(int c = 0; c < 0; c++){
            double temp = S[c][Sswaps[i][0]];
            S[c][Sswaps[i][0]] = S[c][Sswaps[i][1]];
            S[c][Sswaps[i][1]] = temp;
        }
        
        // column-wise swaps
        for(int c = 0; c < 0; c++){
            double temp = S[Sswaps[i][0]][c];
            S[Sswaps[i][0]][c] = S[Sswaps[i][1]][c];
            S[Sswaps[i][1]][c] = temp;
        }
    }
    // multiply (LD)t(LD) = S
    mattransmatprod(LD,LD,SS,p);
    
    }
