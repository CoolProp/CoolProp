#include "mex.h"   /*--This one is required*/

/* Prototype for the Props function to be called.  Can't use the CoolProp.h header because there are a lot of
   c++ parts in the header that cannot be easily hidden when compiling */
double HAProps(char *Output, char* Name1, double Prop1, char* Name2, double Prop2, char* Name3, double Prop3);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int npol,Nel;
    int m,n,u,v,k,Output_len,Name1_len,Name2_len,Name3_len;
    int *c,status;
    double *d,Prop1,Prop2,Prop3,x,y;
	char *Output,*Name1,*Name2,*Name3,errstr[1000];
    mxArray *cMat[1];
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    
    if (nrhs == 7 && mxIsChar (prhs[0]) && mxIsChar (prhs[1]) && mxIsChar (prhs[3]) && mxIsChar (prhs[5]))
    {
        
        /* Get the output (it is a string) (+1 for the NULL terminator)*/
        Output_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Output = mxCalloc(Output_len, sizeof(char));
        status = mxGetString(prhs[0], Output, Output_len);
        
        /* Get the output (it is a string) (+1 for the NULL terminator)*/
        Name1_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Name1 = mxCalloc(Name1_len, sizeof(char));
        mxGetString(prhs[1], Name1, Name1_len);
        
        Name2_len=(mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
        Name2 = mxCalloc(Name2_len, sizeof(char));
        mxGetString(prhs[3], Name2, Name2_len);
        
        Name3_len=(mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
        Name3 = mxCalloc(Name3_len, sizeof(char));
        mxGetString(prhs[5], Name3, Name3_len);
        
        Prop1 = mxGetScalar(prhs[2]);
        Prop2 = mxGetScalar(prhs[4]);
        Prop3 = mxGetScalar(prhs[6]);

        x = HAProps(Output,Name1,Prop1,Name2,Prop2,Name3,Prop3);
        *mxGetPr(plhs[0])=x;
    }
    else
    {
        mexErrMsgTxt("Bad inputs to HAProps");
    }
}