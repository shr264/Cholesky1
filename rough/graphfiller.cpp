//
//  graphfiller.cpp
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
#include "graphfiller.hpp"


using namespace std;


//filled adjacency matrix, i.e. it creates a decomposable graph


int main(){
    int p = 4;
    int fadjMat[p][p];
    
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

    for(int c = (p-1); c >= 0 ; c--){
        cout << "c = " << c << endl;
        for (int j = (c-1); j >= 0; j--){
            cout << "j = " << j << endl;
            for(int k = 0; k<j; k++){
                cout << "k = " << k << endl;
                if((k<j)&&(j<c)){
                    cout << "Entries " << c << "," << k << fadjMat[c][k] << endl;
                    cout << "Entries " <<  j << "," << k << fadjMat[j][k] << endl;
                    cout << "Entries " <<  c << "," << j << fadjMat[j][k] << endl;
                    if((fadjMat[c][k]==1)&&(fadjMat[j][k]==1)){
                        if(fadjMat[c][j]!=1){
                            cout << "Filling entry" << c << "," << j << endl;
                            fadjMat[c][j] = 1;
                            fadjMat[j][c] = 1;
                    }}
            }
        }
    }
}

    
    for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++){
            cout << "Please enter "<< (c+1) << "," << (d+1) << "th element:"<< fadjMat[c][d] <<endl;
        }
    
}