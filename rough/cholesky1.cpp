/*
 
 cholesky code. to compile use: g++ -o cholesky1 cholesky1.cpp -llapack -lblas 
 
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <math.h>
#include <cstdio>

using namespace std;

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    void dsbmv_(const char *uplo,
                       const int *n,
                       const int *k,
                       const double *alpha,
                       const double *a,
                       const int *lda,
                       const double *x,
                       const int *incx,
                       const double *beta,
                       double *y,
                       const int *incy);

      // product C= alphaA.B + betaC                                               
   void dgemm_(char* TRANSA, char* TRANSB, const int* M,
               const int* N, const int* K, double* alpha, double* A,
               const int* LDA, double* B, const int* LDB, double* beta,
               double* C, const int* LDC);
    // product Y= alphaA.X + betaY                                               
   void dgemv_(char* TRANS, const int* M, const int* N,
               double* alpha, double* A, const int* LDA, double* X,
               const int* INCX, double* beta, double* C, const int* INCY);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete IPIV;
    delete WORK;
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
    double* tmp= new double[N];
    initvec(tmp, N);
    dgemv_(&no,&m,&n,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    delete [] tmp;
}

void transmatvecprod(double* A, double* v, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, lda= N, incx= 1, incy= 1;
    double* tmp= new double[N];
    initvec(tmp, N);
    dgemv_(&tr,&m,&n,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    delete [] tmp;
}

void vecmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= 1, k= N, lda= N, ldb= N, ldc= N;
    double* tmp= new double[N];
    initvec(tmp, N);
    dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N; ++i){
        u[i]= tmp[i];
    }
    delete [] tmp;
}

void matmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, k= N, lda= N, ldb= N, ldc= N;
    double* tmp= new double[N*N];
    initvec(tmp, N*N);
    dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N*N; ++i){
        u[i]= tmp[i];
    }
    delete [] tmp;
}

void mattransmatprod(double* v, double* A, double* u, int N){
    double alpha= 1.0, beta= 0.0;
    char no= 'N', tr= 'T';
    int m= N, n= N, k= N, lda= N, ldb= N, ldc= N;
    double* tmp= new double[N*N];
    initvec(tmp, N*N);
    dgemm_(&no,&tr,&m,&n,&k,&alpha,A,&lda,v,&ldb,&beta,tmp,&ldc);
    for(int i= 0; i<N*N; ++i){
        u[i]= tmp[i];
    }
    delete [] tmp;
}


int getIndexOfLargestElement(double arr[], int size) {
    int largestIndex = 0;
    for (int index = largestIndex; index < size; index++) {
        if (arr[largestIndex] < arr[index]) {
            largestIndex = index;
        }
    }
    return largestIndex;
}

int getIndexOfSmallestElement(double arr[], int size) {
    int largestIndex = 0;
    for (int index = largestIndex; index < size; index++) {
        if (arr[largestIndex] > arr[index]) {
            largestIndex = index;
        }
    }
    return largestIndex;
}

int main() {
  int p = 4;
  int n = 5;
  double Ssum = 0;
  double S[p][p];
  double Y[n*p];
  double data[n][p];
  double Correlation[p][p];
  double threshold;
  double newS[p][p];
  double adjMat[p][p];
  double fadjMat[p][p];
  double invS[p*p];
  double inverseS[p][p];
  double L[p][p];
  double D[p][p];
    int medgeset[p][2];
    int medgesetlen = 0;
    int Sswaps[p][2];

  memset(L, 0, sizeof(double) * p * p);
  memset(D, 0, sizeof(double) * p * p);
    memset(Sswaps, 0, sizeof(int) * p*2);

  /*  
   for (int c = 0; c < n; c++)
     for (int d = 0; d < p; d++){
       cout << "Please enter "<< (c+1) << "," << (d+1) << "th element"<< endl;
       cin >> data[c][d];
     }
  */
    data[0][0] = 1.492901;
    data[0][1] = 0.2290852;
    data[0][2] = -2.9856011;
    data[0][3] = -2.258880082;
    data[1][0] = 3.291887;
    data[1][1] = 1.3735888;
    data[1][2] = -0.0656945;
    data[1][3] = -0.006295644;
    data[2][0] = -0.425941;
    data[2][1] = -1.2202329;
    data[2][2] = -3.6724960;
    data[2][3] = -2.223700360;
    data[3][0] = -4.607645;
    data[3][1] = -2.5586426;
    data[3][2] =  0.1473857;
    data[3][3] = -0.081566226;
    data[4][0] = 3.480783;
    data[4][1] = 2.6759905;
    data[4][2] = 1.5341220;
    data[4][3] = 0.186644492;
    

   printf("The data matrix:-\n");
 
    for (int c = 0; c < n; c++)
      for (int d = 0; d < p; d++)
	cout  << (c+1) << "," << (d+1) << ":"<<data[c][d]  << endl;


   

   // first we calculate S

  for(int i = 0; i<n; i++){
    for(int j = 0; j<p; j++){
      for (int k = 0; k < n; k++){
	Ssum = Ssum + data[k][i]*data[k][j];
      }
      S[i][j] = Ssum/n;
      Ssum = 0;
    }
  }
      printf("The sample covariance matrix:-\n");
 
    for (int c = 0; c < p; c++)
      for (int d = 0; d < p; d++)
	cout  << (c+1) << "," << (d+1) << ":"<<S[c][d]  << endl;

    // now we convert the covariance matrix to a correlation matrix

    for(int i = 0; i<p; i++){
      for(int j = 0; j<p; j++){
	Correlation[i][j] = S[i][j]/(pow(S[i][i],0.5)*pow(S[j][j],0.5));
    }
  }

   printf("The sample correlation matrix:-\n");
 
    for (int c = 0; c < p; c++)
      for (int d = 0; d < p; d++)
	cout  << (c+1) << "," << (d+1) << ":"<< Correlation[c][d]  << endl;

    // now we calculate the new S

    cout << "Please enter threshold: "<< endl;
    cin >> threshold;
    

    for(int i = 0; i<p; i++){
      for(int j = 0; j<p; j++){
	if(Correlation[i][j]>threshold){
	  newS[i][j] = S[i][j];
        adjMat[i][j] = 1;
        fadjMat[i][j] = 1;
	}
	else{
	  newS[i][j] = 0;
        adjMat[i][j] = 0;
        fadjMat[i][j] = 0;
	}
    }
  }

    fadjMat[0][0] = 1;
    fadjMat[0][1] = 1;
    fadjMat[0][2] = 0;
    fadjMat[0][3] = 1;
    fadjMat[1][0] = 1;
    fadjMat[1][1] = 1;
    fadjMat[1][2] = 1;
    fadjMat[1][3] = 0;
    fadjMat[2][0] = 0;
    fadjMat[2][1] = 1;
    fadjMat[2][2] = 1;
    fadjMat[2][3] = 1;
    fadjMat[3][0] = 1;
    fadjMat[3][1] = 0;
    fadjMat[3][2] = 1;
    fadjMat[3][3] = 1;
    
    
    
    printf("The new sample covariance matrix:-\n");
 
    for (int c = 0; c < p; c++)
      for (int d = 0; d < p; d++)
	cout  << (c+1) << "," << (d+1) << ":"<< newS[c][d]  << endl;

    for (int c = 0; c < p; c++)
      for (int d = 0; d < p; d++)
	invS[c*p+d] = newS[c][d];
    
    inverse(invS,p);
    
    printf("The inverse sample covariance matrix:-\n");
 
    for (int c = 0; c < p; c++)
      for (int d = 0; d < p; d++){
	cout  << (c+1) << "," << (d+1) << ":"<< invS[c*p+d]  << endl;
	inverseS[c][d] = invS[c*p+d];
      }
    
    /* this is the sorting portion
     
     */
    
    for(int i = 0; i < p; i++){
        double *zeros;
        zeros = (double *)malloc((p-i)*sizeof(double));
        
        memset(zeros, 0, sizeof(double) * (p-i));
        
        printf("The zeros vector:-\n");
        
        for (int c = i; c < p; c++){
            cout  << c << ":" << zeros[c]  << endl;
        }
        
        for(int c=i; c<p ;c++){
            int zerosum = 0;
            for(int d=i; d<p ;d++){
                if(adjMat[c][d]==0){
                    zerosum += 1;
                    cout << "zerosum at i = "<< i << " and " << c  << "," << d << "is: " << zerosum << endl;
                }
            }
        zeros[c-i] = zerosum;
        }
        
        printf("The zeros vector:-\n");
        
        for (int c = i; c < p; c++){
            cout  << c << ":"<< zeros[c]  << endl;
        }
        
        int max_pos = i;
        int min_pos = i;
        
        for (int index = max_pos; index < p; index++) {
            if (zeros[max_pos] < zeros[index]) {
                max_pos= index;
            }
            if (zeros[min_pos] > zeros[index]) {
                min_pos= index;
            }
        }
    
        cout << "Index of max element: "
            << max_pos
            << endl;
    
        cout << "Index of min element: "
        << min_pos
        << endl;

        double break1;
        
        cout << "Please enter break: "<< endl;
        cin >> break1;
        
        if(max_pos!=min_pos){
            Sswaps[i][0] = min_pos;
            Sswaps[i][1] = max_pos;
            // column-wise swaps
            for(int c = 0; c < p; c++){
                double temp = S[c][min_pos];
                S[c][min_pos] = S[c][max_pos];
                S[c][max_pos] = temp;
                
                double adjtemp = adjMat[c][min_pos];
                adjMat[c][min_pos] = adjMat[c][max_pos];
                adjMat[c][max_pos] = adjtemp;
                
                double fadjtemp = fadjMat[c][min_pos];
                fadjMat[c][min_pos] = fadjMat[c][max_pos];
                fadjMat[c][max_pos] = fadjtemp;
            }

            // row-wise swaps
            for(int c = 0; c < p; c++){
                double temp = S[min_pos][c];
                S[min_pos][c] = S[max_pos][c];
                S[max_pos][c] = temp;
                
                double adjtemp = adjMat[min_pos][c];
                adjMat[min_pos][c] = adjMat[max_pos][c];
                adjMat[max_pos][c]= adjtemp;
                
                double fadjtemp = fadjMat[min_pos][c];
                fadjMat[min_pos][c] = fadjMat[max_pos][c];
                fadjMat[max_pos][c]= fadjtemp;
                
            }
        }
        free(zeros);
    }
    
    double break1;
    
    cout << "Please enter break: "<< endl;
    cin >> break1;
    
    printf("The row-column swapped sample covariance matrix:-\n");
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++)
            cout  << (c+1) << "," << (d+1) << ":"<< S[c][d]  << endl;
    
    //filled adjacency matrix, i.e. it creates a decomposable graph
    
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
    
    printf("The adjacency matrix:-\n");
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++)
            cout  << (c+1) << "," << (d+1) << ":"<< adjMat[c][d]  << endl;
    
    cout << "Please enter break: "<< endl;
    cin >> break1;
    
    printf("The filled adjacency matrix:-\n");
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++)
            cout  << (c+1) << "," << (d+1) << ":"<< fadjMat[c][d]  << endl;
    
    cout << "Please enter break: "<< endl;
    cin >> break1;
    
    //this is the cholesky calculation
    
    L[p-1][p-1] = 1;
    D[p-1][p-1] = 1/S[p-1][p-1];
    
    for (int i1 = 0; i1< (p-1); i1++){
        int i = i1;
        cout << "Cholesky Calcluation at i = " << i << endl;
        double Sii = S[i][i];
        double temp;
        
        int counter = 0;
        
        for (int c = (i+1); c < p; c++){
            if(fadjMat[i][c]>0){
                counter += 1;
            }
        }
        
        cout << "Non-zero entries in row " << (i+1) << "i is "<< counter << endl;
        
        if(counter==0){
            cout << "Counter 0"<< endl;
            double dii = pow(1/(Sii),0.5);
            L[i][i] = 1;
            D[i][i] = pow(dii,2);
        }
        else{
            cout << "Counter non-zero" <<endl;
            double *Sdoti;
            Sdoti = (double *)malloc(counter*sizeof(double));
            //Sdoti=new double[counter];
            double *Sgp;
            Sgp = (double *)malloc(counter*counter*sizeof(double));
            //Sgp=new double[counter*counter];
            double *u;
            u = (double *)malloc(counter*sizeof(double));
            //u=new double[counter];
            
            memset(u, 0, sizeof(double) * counter);
            
            for(int c = 0; c < counter ; c++){
                cout << c << ": " << u[c] << endl;
            }
            
            int count1 = 0;
            cout << "count1 set to 0" <<endl;
            for (int c = (i+1); c < p; c++){
                cout << "c = " << c <<endl;
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
            
            cout << "Printing Sgp" << endl;
            
            for (int c = 0; c < counter; c++)
                for (int d = 0; d < counter; d++){
                    cout  << (c+1) << "," << (d+1) << ":"<< Sgp[c*counter+d]  << endl;
                }
            
            cout << "Printing Sdoti" << endl;
            
            for (int c = 0; c < counter; c++){
                cout  << (c+1) << ":"<< Sdoti[c]  << endl;
            }
            
            double break1;
            
            cout << "Please enter break: "<< endl;
            cin >> break1;
            
            printf("Calculating & Printing inverse\n");
            
            inverse(Sgp,counter);
            
            for (int c = 0; c < counter; c++)
                for (int d = 0; d < counter; d++){
                    cout  << (c+1) << "," << (d+1) << ":"<< Sgp[c*counter+d]  << endl;
                }
            
            printf("Calcuating and Printing matrix vec prod:-\n");
            
            matvecprod(Sgp,Sdoti,u,counter);
            
            
            printf("u is:-\n");
            for (int c = 0; c < counter; c++){
                cout  << (c+1) << ":"<< u[c]  << endl;
                
            }
            
            double tempprod = 0;
            
            for(int c = 0; c < counter; c++){
                tempprod = tempprod - u[c]*Sdoti[c];
            }
            
            printf("Calcuating dii:-\n");
            
            double dii;
            dii = sqrt(1/(Sii + tempprod));
            
            cout << "Temp Product is: " << tempprod << endl;
            cout << "Sii: " << Sii << endl;
            cout << "dii: " << dii << endl;
            
            
            L[i][i] = 1;
            int count3 = 0;
            for(int c = (i+1); c < p ; c++){
                if(fadjMat[i][c]>0){
                    L[i][c] = -u[count3];
                    cout << "L" << i <<","<< c  << "," << count3 << ":" << u[count3] << endl;
                    count3 += 1;
                }
            }
            D[i][i] = pow(dii,2);
            cout << "D"<< i <<","<< i  << ":" << 	D[i][i] << endl;
            
            free(Sdoti);
            free(Sgp);
            free(u);
            //delete Sdoti;
            //delete Sgp;
            //delete u;
        }
    }
    //end cholesky
    
    
    //this is algorithm 1 from the paper
    
    for (int c = 0; c < medgesetlen; c++){
        if(medgeset[c][1]>1){
            double Lsum = 0;
            for(int k = 0; k < (medgeset[c][1]-1); k++){
                Lsum += L[medgeset[c][0]][k]*L[medgeset[c][1]][k];
            }
            L[medgeset[c][0]][medgeset[c][1]] = -Lsum/L[medgeset[c][1]][medgeset[c][1]];
        }
        else{
            L[medgeset[c][0]][medgeset[c][1]] = 0;
        }
    }
    
    // end algorothim 1
    
    
    printf("L is:-\n");
      
      for (int c = 0; c < p; c++)
          for (int d = 0; d < p; d++){
	  cout  << (c+1) << "," << (d+1) << ":"<< L[c][d]  << endl;

	}

      printf("D is:-\n");

      for (int c = 0; c < p; c++)
	for (int d = 0; d < p; d++){
	  cout  << (c+1) << "," << (d+1) << ":"<< D[c][d]  << endl;
	}

    double vecL[p*p];
    double vecD[p*p];
    double vecU[p*p];
    double vecprod[p*p];
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            vecD[c*p+d] = D[c][d];
        }

    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            vecL[c*p+d] = L[c][d];
        }
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            vecU[c*p+d] = 0;
        }
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            vecprod[c*p+d] = 0;
        }
    
    matmatprod(vecD,vecL,vecU,p);
    
    mattransmatprod(vecL,vecU,vecprod,p);
    
    printf("LDL^t is:-\n");
    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            cout  << (c+1) << "," << (d+1) << ":"<< vecprod[c*p+d]  << endl;
        }
    
    // column-wise swaps back to original ordering
    
    
    // row-wise swaps back to original ordering

    
    return 0;
}
