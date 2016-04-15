//
//  LUdecomp.cpp
//  
//
//  Created by Syed Rahman on 2/5/16.
//
//

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

int main(){
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
    
    for(int i = 0; i<n; i++){
        for(int j = 0; j<p; j++){
            for (int k = 0; k < n; k++){
                Ssum = Ssum + data[k][i]*data[k][j];
            }
            S[i][j] = Ssum/n;
            Ssum = 0;
        }
    }
    
    fadjMat[0][0] = 1;
    fadjMat[0][1] = 1;
    fadjMat[0][2] = 0;
    fadjMat[0][3] = 1;
    fadjMat[1][0] = 1;
    fadjMat[1][1] = 1;
    fadjMat[1][2] = 1;
    fadjMat[1][3] = 1;
    fadjMat[2][0] = 0;
    fadjMat[2][1] = 1;
    fadjMat[2][2] = 1;
    fadjMat[2][3] = 1;
    fadjMat[3][0] = 1;
    fadjMat[3][1] = 1;
    fadjMat[3][2] = 1;
    fadjMat[3][3] = 1;
    
    
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
            cout << "c = " << c << "counter = "<< counter << endl;
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
        Sdoti=new double[counter];
        double *Sgp;
        Sgp=new double[counter*counter];
        double *u;
        u=new double[counter];
        
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
                cout << count1 << "," << Sdoti[count1] << ","<< S[i][c] << endl;
                int count2 = 0;
                cout << "count2 set to 0" <<endl;
                for(int d = (i+1); d < p; d++){
                    cout << "d = " << d <<endl;
                    if(fadjMat[d][i]>0){
                        
                        Sgp[count2 + counter*count1]=S[d][c];
                        cout << c << "," << d << "," << S[d][c] << endl;
                        cout << count2 << "+" << count1 << "*" << counter << endl;
                        cout << count2 +counter*count1 << "," << Sgp[count2 +counter*count1] << endl;
                        if(count2<counter){
                            cout << "count2 increasing" <<endl;
                            count2 += 1;}
                    }
                    
                    }
cout << "count1 increasing" <<endl;
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
        delete Sdoti;
        delete Sgp;
        delete u;	
    }
    }
    //end cholesky
    
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
    
}
