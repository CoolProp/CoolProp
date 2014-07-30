#include "CoolProp.h"
#include "CoolPropTools.h"
#include "mex.h"   /*--This one is required*/
#include <cstdio>
#include <string>
#include <cstring>
#include "float.h"

extern "C" void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int npol,Nel;
    int m,n,u,v,k,Output_len,Name1_len,Name2_len,Ref_len,Param_len;
    int *c,status;
    double *d,Prop1,Prop2,x,y;
	char *Output,*Name1,*Name2,*Ref,*Param,errstr[1000],errstr2[1000],fluidslist[10000];
    double val;
    mxArray *cMat[1];
    
    if (nrhs == 2 && mxIsChar (prhs[0]) && mxIsChar (prhs[1]))
    {
        /* Get the refrigerant (it is a string) (+1 for the NULL terminator)*/
        Ref_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Ref = (char*)mxCalloc(Ref_len, sizeof(char));
        status = mxGetString(prhs[0], Ref, Ref_len);
        
        /* Get the output (it is a string) (+1 for the NULL terminator)*/
        Output_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Output = (char*)mxCalloc(Output_len, sizeof(char));
        status = mxGetString(prhs[1], Output, Output_len);
        
        /* Try to shortcut to get the strings for the fluid */
        if (!strcmp(Output,"aliases") || !strcmp(Output,"CAS") || !strcmp(Output,"CAS_number") || !strcmp(Output,"ASHRAE34") || !strcmp(Output,"REFPROPName") || !strcmp(Output,"REFPROP_name") || !strcmp(Output,"TTSE_mode"))
        {
            std::string fluidslist = CoolProp::get_fluid_param_string(Ref,Output);
            plhs[0] = mxCreateString(fluidslist.c_str());
            return;
        }
        
        /* Ok, now try to get the values for the fluid */
        
        /* Create matrix for the return argument. */
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

        /* Get the value*/
        val = CoolProp::PropsSI(Ref,"T",0,"P",0,Output);
        /* If it is a good value, return it*/
        if (ValidNumber(val))
        {
            *mxGetPr(plhs[0]) = val; 
        }
        /* Otherwise there was an error, return the CoolProp error*/
        else
        {
            std::string errstr = CoolProp::get_global_param_string("errstring");
            sprintf(errstr2,"CoolProp Error: %s",errstr.c_str());
            mexErrMsgTxt(errstr2);
        }
    }
    else if (nrhs == 6 && mxIsChar (prhs[0]) && mxIsChar (prhs[1])  && mxIsChar (prhs[3])  && mxIsChar (prhs[5])  )
    {
        /* Create matrix for the return argument. */
        plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
        
        /* Get the refrigerant (it is a string) (+1 for the NULL terminator)*/
        Ref_len=(mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
        Ref = (char*)mxCalloc(Ref_len, sizeof(char));
        mxGetString(prhs[5], Ref, Ref_len);
        
        /* Get the output (it is a string) (+1 for the NULL terminator)*/
        Output_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Output = (char*)mxCalloc(Output_len, sizeof(char));
        mxGetString(prhs[0], Output, Output_len);
        
        Name1_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Name1 = (char*)mxCalloc(Name1_len, sizeof(char));
        mxGetString(prhs[1], Name1, Name1_len);
        
        Name2_len=(mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
        Name2 = (char*)mxCalloc(Name2_len, sizeof(char));
        mxGetString(prhs[3], Name2, Name2_len);
        
        Prop1 = mxGetScalar(prhs[2]);
        Prop2 = mxGetScalar(prhs[4]);

        val = CoolProp::PropsSI(Output,Name1,Prop1,Name2,Prop2,Ref);
        
        /* If it is a good value, return it */
        if (ValidNumber(val))
        {
            *mxGetPr(plhs[0]) = val; 
        }
        /* Otherwise there was an error, return the CoolProp error*/
        else
        {
            std::string errstr = CoolProp::get_global_param_string("errstring");
            sprintf(errstr2,"CoolProp Error: %s",errstr.c_str());
            mexErrMsgTxt(errstr2);
        }
    }
    else if (nrhs == 1 && mxIsChar (prhs[0]))
    {
        Ref_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Ref = (char*)mxCalloc(Ref_len, sizeof(char));
        status = mxGetString(prhs[0], Ref, Ref_len);
        std::string s;
        if (!strcmp(Ref,"FluidsList")){
            s = CoolProp::get_global_param_string("FluidsList");
        }
        else if (!strcmp(Ref,"version")){
            s = CoolProp::get_global_param_string("version");
        }
        else if (!strcmp(Ref,"gitrevision")){
            s = CoolProp::get_global_param_string("gitrevision");
        }
        else
        {
            sprintf(errstr2,"single input is invalid: %s",Ref);
            mexErrMsgTxt(errstr2);
        }
        plhs[0] = mxCreateString(s.c_str());
        return;
    }
    else
    {
        mexErrMsgTxt("Props must either receive two strings or the signature Props(Output,Param1,Value1,Param2,Value2,FluidName)");
    }
}
