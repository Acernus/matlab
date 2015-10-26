#include <iostream>
#include <vector>
#include "mex.h"
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int aRow = static_cast<int>(mxGetM(prhs[0]));
	int aCol = static_cast<int>(mxGetN(prhs[0]));
	int bRow = static_cast<int>(mxGetM(prhs[1]));
	int bCol = static_cast<int>(mxGetN(prhs[1]));
    vector<vector<double> > a(aRow, vector<double>(aCol)), b(aRow, vector<double>(aCol));
    double *d1 = mxGetPr(prhs[0]);
    double *d2 = mxGetPr(prhs[1]);
    vector<double> temp(aRow * aCol);
    for(int i = 0; i < aRow * aCol; ++i) {
        temp[i] = d1[i];
    }
    mexPrintf("temp array :\n");
    for(int i = 0; i < aRow * aCol;++i) {
        mexPrintf("%f ", temp[i]);
    }
    mexPrintf("\n");
    for (int i = 0; i < aRow; ++i) {
		for (int j = 0; j < aCol; ++j) {
			a[i][j] = d1[j * aRow + i];
            b[i][j] = d2[j * aRow + i];
		}
	}
    mexPrintf("array a:\n");
    for (int i = 0; i < aRow; ++i) {
		for (int j = 0; j < aCol; ++j) {
			mexPrintf("%f ", a[i][j]);
		}
        mexPrintf("\n");
	}
    mexPrintf("array b:\n");
    for (int i = 0; i < bRow; ++i) {
		for (int j = 0; j < bCol; ++j) {
			mexPrintf("%f ", b[i][j]);
		}
        mexPrintf("\n");
	}

	double *output;
    plhs[0] = mxCreateDoubleMatrix(aRow, aCol, mxREAL);
    output = mxGetPr(plhs[0]);
    int k = 0;
    for (int i = 0; i < aCol; ++i) {
        for (int j = 0; j < aRow; ++j) {
            output[k++] = a[j][i] + b[j][i]; 
        }
    }
}