#ifndef FUNCTION_H_INCLUDED
#define FUNCTION_H_INCLUDED

void inverse(double* A, int N);

void initvec(double* v, int N);

void matvecprod(double* A, double* v, double* u, int N);

void transmatvecprod(double* A, double* v, double* u, int N);

void vecmatprod(double* v, double* A, double* u, int N);

void matmatprod(double* v, double* A, double* u, int N);

void mattransmatprod(double* v, double* A, double* u, int N);

int getIndexOfLargestElement(double arr[], int size);

int getIndexOfSmallestElement(double arr[], int size);

void choleskyalgo1 (int *nIn, int *pIn, double *data, double *SS, double *thresholdIn, double *L, double *D, double *LD);

void choleskyalgo1adj (int *nIn, int *pIn, double *data, double *SS, int *adjMatIn, double *L, double *D, double *LD);

#endif
