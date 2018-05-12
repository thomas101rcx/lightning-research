#include "mex.h"
#include "math.h"

double calcudis(double x1, double y1, double z1, double x2, double y2, double z2, double x, double y, double z);
double *linearto3D(double index,  double xgridsize, double ygridsize);

double* linearto3D(double index,  double xgridsize, double ygridsize){
  static double coordinate[3];
  coordinate[0] = fmod(index-1,xgridsize) + 1;
  coordinate[1] = fmod(((index - coordinate[0])/xgridsize) ,ygridsize) + 1;
  coordinate[2] = ((index - coordinate[0])/xgridsize-coordinate[1] + 1)/ ygridsize + 1;
  return coordinate;
}

double calcudis(double x1, double y1, double z1, double x2, double y2, double z2, double x, double y, double z){
  double c = 299792.458;
  double timedelay;
  timedelay = (sqrt(pow((x - x1),2) + pow((y - y1),2)  + pow((z - z1) ,2)) / c) - (sqrt(pow(x - x2,2)+ pow(y - y2,2) + pow(z - z2,2)) / c);
  return timedelay;
}

void timegrid(double timegridsize, double xgridsize, double ygridsize, double X1Solx, double X1Soly, double Y1Solx, double Y1Soly,
              double Combx, double Comby, double Combz, double Combx2, double Comby2, double Combz2, double *timedifferencematrix){
  for (double i = 1 ; i < timegridsize  ;i++){
      double X1, Y1, Z1;
      double* coordinate = linearto3D(i,xgridsize,ygridsize);
      X1 = (coordinate[0] - X1Solx)*(X1Soly);
      Y1 = (coordinate[1] - Y1Solx)*(Y1Soly);
      Z1 = (coordinate[2] - 1) *0.1;
      *(timedifferencematrix+ (int)i-1) = calcudis(Combx, Comby, Combz, Combx2, Comby2, Combz2, X1, Y1, Z1);
  }
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
	double timegridsize, xgridsize, ygridsize, X1Solx, X1Soly, Y1Solx, Y1Soly, Combx, Comby, Combz, Combx2, Comby2, Combz2;
  double *timedifferencematrix;

	/* Check for proper number of arguments. */
	if (nrhs != 14) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:invalidNumInputs", "Incorrect number of inputs");
	}

	if (nlhs > 2) {
		mexErrMsgIdAndTxt("MATLAB:fillgrid:maxlhs", "Too many output arguments.");
	}

  timegridsize = mxGetScalar(prhs[0]);
	/* Create matrix for the return argument. */
	plhs[0] = mxCreateDoubleMatrix((mwSize)(int) timegridsize, 1 ,mxREAL);
	/* Assign pointers to each input and output. */
  xgridsize = mxGetScalar(prhs[1]);
  ygridsize = mxGetScalar(prhs[2]);
  X1Solx = mxGetScalar(prhs[3]);
  X1Soly = mxGetScalar(prhs[4]);
  Y1Solx = mxGetScalar(prhs[5]);
  Y1Soly = mxGetScalar(prhs[6]);
  Combx = mxGetScalar(prhs[7]);
  Comby = mxGetScalar(prhs[8]);
  Combz = mxGetScalar(prhs[9]);
  Combx2 = mxGetScalar(prhs[10]);
  Comby2 = mxGetScalar(prhs[11]);
  Combz2 = mxGetScalar(prhs[12]);

  //double *coordinate = linearto3D(151*40*101, 151, 151);
  //Output
	timedifferencematrix = mxGetPr(plhs[0]);
	/* Call the subroutine. */
	timegrid(timegridsize, xgridsize, ygridsize, X1Solx, X1Soly, Y1Solx, Y1Soly,Combx, Comby, Combz, Combx2, Comby2, Combz2, timedifferencematrix);
}
