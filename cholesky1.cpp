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
    int m= N, n= N, lda= N, incx= N, incy= N;
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
    int m= N, n= N, lda= N, incx= N, incy= N;
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

int main() {
  int p = 3;
  int n = 5;
  double sum = 0;
  double S[p][p];
  double Y[n*p];
  double data[n][p];
  double Correlation[p][p];
  double threshold;
  double newS[p][p];
  double invS[p*p];
  double inverseS[p][p];
  double L[p][p];
  double D[p][p];

  memset(L, 0, sizeof(double) * p * p);
  memset(D, 0, sizeof(double) * p * p); 

  /*  
   for (int c = 0; c < n; c++)
     for (int d = 0; d < p; d++){
       cout << "Please enter "<< (c+1) << "," << (d+1) << "th element"<< endl;
       cin >> data[c][d];
     }
  */
  for(int c=0;c<n*p;c++)
    Y[c] = (c+1)%5;

   
   for (int c = 0; c < n; c++)
     for (int d = 0; d < p; d++){
       data[c][d] = Y[c*p+d];
     }
  

   // first we calculate S

  for(int i = 0; i<n; i++){
    for(int j = 0; j<p; j++){
      for (int k = 0; k < n; k++){
	sum = sum + data[k][i]*data[k][j];
      }
      S[i][j] = sum/n;
      sum = 0;
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
	}
	else{
	  newS[i][j] = 0;
	}
    }
  }

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

    //this is the cholesky calculation
    
    L[p-1][p-1] = 1;
    D[p-1][p-1] = 1/newS[p-1][p-1];
    
    for (int i1 = 0; i1< (p-1); i1++){
      int i = i1;
      cout << "Cholesky Calcluation at i = " << i << endl;
      double Sii = newS[i][i];
      double temp;

      int counter = 0;
    
      for (int c = (i+1); c < p; c++){
	if(fabs(newS[i][c])>0){counter += 1;}
      }
      
      cout << "Non-zero entries in row " << (i+1) << "i s "<< counter << endl;
    
      if(counter==0){
          cout << "Counter 0"<< endl;
	int dii = pow(1/(Sii),0.5);
	L[i][i] = 1;
	D[i][i] = pow(dii,2);
      }
      else{
          cout << "Counter non-zero" <<endl;
	double *Sdoti;
	Sdoti=new double[counter];
	double *Sgp;
	Sgp=new double[counter*counter];
	double *u;
	u=new double[counter];
          
	memset(u, 0, sizeof(double) * counter);

	for (int c = (i+1); c < p; c++){
	  int count1 = 0;
	  if(fabs(newS[i][c])>0){
	    Sdoti[count1] = newS[i][c];
	    cout << count1 << "," << Sdoti[count1] << ","<< newS[i][c] << endl;
	    int count2 = 0;
	    for(int d = (i+1); d < p; d++){
	      if(fabs(newS[c][d])>0){
		Sgp[count1*counter+count2]=newS[c][d];
		cout << count1 << "," << count2 << "," << Sgp[count1*counter+count2] << endl;
		count2 += 1;
	      }
        count1 += 1;
	    }  
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

	printf("Calculating inverse\n");
	
	inverse(Sgp,counter);
	
	printf("Calcuating matrix vec prod:-\n");
	
	matvecprod(Sgp,Sdoti,u,counter);
	
	double tempprod;
    
	for(int c = 0; c < counter; c++){
        tempprod = -u[c]*Sdoti[c];
    }
          
	printf("Calcuating dii:-\n");
	
          double dii;
          dii = sqrt(1/(Sii + tempprod));
          
          cout << "Temp Product is: " << tempprod << endl;
          cout << "Sii: " << Sii << endl;
          cout << "dii: " << dii << endl;

	
	L[i][i] = 1;
	for(int c = (i+1); c < p ; c++){
	  int count1 = 0;
	  if(newS[i][c]>0){
	    L[i][c] = -u[count1];
	    cout << "L" << i <<","<< c  << ":" << u[count1] << endl;
	    count1 += 1;
	  }
	}
	D[i][i] = pow(dii,2);
	cout << "D"<< i <<","<< i  << ":" << 	D[i][i] << endl;
	delete Sdoti;
	delete Sgp;
	delete u;	
      }
    }
    
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
    
    return 0;
}
