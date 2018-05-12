#include "mex.h"

void fillgrid(int Ixyzlength, double *Ixyz1, short A[], float *Ixy){
  for (int i = 0 ; i <  Ixyzlength; i++) {
    *(Ixyz1+i)=Ixy[(int)A[i]+1000];
  }
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
	double *Ixyz1;
  short *A;
  float *Ixy;

	/* Check for proper number of arguments. */
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:invalidNumInputs", "Incorrect number of inputs");
	}

	if (nlhs > 2) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:maxlhs", "Too many output arguments.");
	}

  int Ixyzlength = mxGetNumberOfElements(prhs[0]);
	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix((mwSize) Ixyzlength, 1,mxREAL);

	/* Assign pointers to each input and output. */
	A = (short*) mxGetData(prhs[0]);

	Ixy = (float*) mxGetData(prhs[1]);

  //Output
	Ixyz1 = mxGetPr(plhs[0]);

	/* Call the subroutine. */
	fillgrid(Ixyzlength, Ixyz1, A, Ixy);
}
