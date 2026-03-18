/*
 * combiner_adm1_bsm2.c calculates the output concentrations when adding two ADM1 flow  
 * streams together. Output temperature always based on 'heat content' of the
 * influent flows, i.e. parameter TEMPMODEL is not used. If all input flow rates are
 * less or equal to zero then all outputs are zero.
 */

#define S_FUNCTION_NAME combiner_adm1_bsm2

#include "simstruc.h"

/*
 * mdlInitializeSizes - initialize the sizes array
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);   /* number of discrete states             */
    ssSetNumInputs(        S, 66);  /* number of inputs                      */
    ssSetNumOutputs(       S, 33);  /* number of outputs                     */
    ssSetDirectFeedThrough(S, 1);   /* direct feedthrough flag               */
    ssSetNumSampleTimes(   S, 1);   /* number of sample times                */
    ssSetNumSFcnParams(    S, 0);   /* number of input arguments             */
    ssSetNumRWork(         S, 0);   /* number of real work vector elements   */
    ssSetNumIWork(         S, 0);   /* number of integer work vector elements*/
    ssSetNumPWork(         S, 0);   /* number of pointer work vector elements*/
}

/*
 * mdlInitializeSampleTimes - initialize the sample times array
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}


/*
 * mdlInitializeConditions - initialize the states
 */
static void mdlInitializeConditions(double *x0, SimStruct *S)
{
}

/*
 * mdlOutputs - compute the outputs
 */

static void mdlOutputs(double *y, double *x, double *u, SimStruct *S, int tid)
{
    
  if ((u[26] > 0) || (u[59] > 0)) {
    y[0]=(u[0]*u[26] + u[33]*u[59])/(u[26]+u[59]);
    y[1]=(u[1]*u[26] + u[34]*u[59])/(u[26]+u[59]);
    y[2]=(u[2]*u[26] + u[35]*u[59])/(u[26]+u[59]);
    y[3]=(u[3]*u[26] + u[36]*u[59])/(u[26]+u[59]);
    y[4]=(u[4]*u[26] + u[37]*u[59])/(u[26]+u[59]);
    y[5]=(u[5]*u[26] + u[38]*u[59])/(u[26]+u[59]);
    y[6]=(u[6]*u[26] + u[39]*u[59])/(u[26]+u[59]);
    y[7]=(u[7]*u[26] + u[40]*u[59])/(u[26]+u[59]);
    y[8]=(u[8]*u[26] + u[41]*u[59])/(u[26]+u[59]);
    y[9]=(u[9]*u[26] + u[42]*u[59])/(u[26]+u[59]);
    y[10]=(u[10]*u[26] + u[43]*u[59])/(u[26]+u[59]);
    y[11]=(u[11]*u[26] + u[44]*u[59])/(u[26]+u[59]);
    y[12]=(u[12]*u[26] + u[45]*u[59])/(u[26]+u[59]);
    y[13]=(u[13]*u[26] + u[46]*u[59])/(u[26]+u[59]);
    y[14]=(u[14]*u[26] + u[47]*u[59])/(u[26]+u[59]);
    y[15]=(u[15]*u[26] + u[48]*u[59])/(u[26]+u[59]);
    y[16]=(u[16]*u[26] + u[49]*u[59])/(u[26]+u[59]);
    y[17]=(u[17]*u[26] + u[50]*u[59])/(u[26]+u[59]);
    y[18]=(u[18]*u[26] + u[51]*u[59])/(u[26]+u[59]);
    y[19]=(u[19]*u[26] + u[52]*u[59])/(u[26]+u[59]);
    y[20]=(u[20]*u[26] + u[53]*u[59])/(u[26]+u[59]);
    y[21]=(u[21]*u[26] + u[54]*u[59])/(u[26]+u[59]);
    y[22]=(u[22]*u[26] + u[55]*u[59])/(u[26]+u[59]);
    y[23]=(u[23]*u[26] + u[56]*u[59])/(u[26]+u[59]);
    y[24]=(u[24]*u[26] + u[57]*u[59])/(u[26]+u[59]);
    y[25]=(u[25]*u[26] + u[58]*u[59])/(u[26]+u[59]);
    
    y[26]=u[26]+u[59];                               /* Flow rate */
    
    y[27]=(u[27]*u[26] + u[60]*u[59])/(u[26]+u[59]);    /* Temp */
         
    /* Dummy states */
    y[28]=(u[28]*u[26] + u[61]*u[59])/(u[26]+u[59]);
    y[29]=(u[29]*u[26] + u[62]*u[59])/(u[26]+u[59]);
    y[30]=(u[30]*u[26] + u[63]*u[59])/(u[26]+u[59]);
    y[31]=(u[31]*u[26] + u[64]*u[59])/(u[26]+u[59]);
    y[32]=(u[32]*u[26] + u[65]*u[59])/(u[26]+u[59]);

  }
  else {
    y[0]=0;
    y[1]=0;
    y[2]=0;
    y[3]=0;
    y[4]=0;
    y[5]=0;
    y[6]=0;
    y[7]=0;
    y[8]=0;
    y[9]=0;
    y[10]=0;
    y[11]=0;
    y[12]=0;
    y[13]=0;
    y[14]=0;
    y[15]=0;
    y[16]=0;
    y[17]=0;
    y[18]=0;
    y[19]=0;
    y[20]=0;
    y[21]=0;
    y[22]=0;
    y[23]=0;
    y[24]=0;
    y[25]=0;
    y[26]=0;
    y[27]=0;
    y[28]=0;
    y[29]=0;
    y[30]=0;
    y[31]=0;
    y[32]=0;

  }
}

/*
 * mdlUpdate - perform action at major integration time step
 */

static void mdlUpdate(double *x, double *u, SimStruct *S, int tid)
{
}

/*
 * mdlDerivatives - compute the derivatives
 */
static void mdlDerivatives(double *dx, double *x, double *u, SimStruct *S, int tid)
{
}


/*
 * mdlTerminate - called when the simulation is terminated.
 */
static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

