#include "mex.h"
#include "math.h"

double *linearto3d(double index,  double xgridsize, double ygridsize);

double* linearto3d(double index,  double xgridsize, double ygridsize){
  static double coordinate[3];
  coordinate[0] = fmod(index-1,xgridsize) + 1;
  coordinate[1] = fmod(((index - coordinate[0])/xgridsize) ,ygridsize) + 1;
  coordinate[2] = ((index - coordinate[0])/xgridsize-coordinate[1] + 1)/ ygridsize + 1;
  return coordinate;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
	double index , xgridsize, ygridsize;

	/* Check for proper number of arguments. */
	if (nrhs != 3) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:invalidNumInputs", "Incorrect number of inputs");
	}

	if (nlhs > 4) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:maxlhs", "Too many output arguments.");
	}
  //Create input arguments

  index = mxGetScalar(prhs[0]);
  xgridsize = mxGetScalar(prhs[1]);
  ygridsize = mxGetScalar(prhs[2]);

  double* coordinate = linearto3d(index, xgridsize, ygridsize);

  plhs[0] = mxCreateDoubleMatrix((mwSize)(int) 1, 1 ,mxREAL);
  plhs[1] = mxCreateDoubleMatrix((mwSize)(int) 1, 1 ,mxREAL);
  plhs[2] = mxCreateDoubleMatrix((mwSize)(int) 1, 1 ,mxREAL);

  *mxGetPr(plhs[0]) = coordinate[0];
  *mxGetPr(plhs[1]) = coordinate[1];
  *mxGetPr(plhs[2]) = coordinate[2];

	/* Call the subroutine. */

}
