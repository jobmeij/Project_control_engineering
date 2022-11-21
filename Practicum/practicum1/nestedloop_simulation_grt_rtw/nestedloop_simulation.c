/*
 * nestedloop_simulation.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "nestedloop_simulation".
 *
 * Model version              : 1.16
 * Simulink Coder version : 8.12 (R2017a) 16-Feb-2017
 * C source code generated on : Thu Oct 18 14:04:24 2018
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "nestedloop_simulation.h"
#include "nestedloop_simulation_private.h"

/* Block signals (auto storage) */
B_nestedloop_simulation_T nestedloop_simulation_B;

/* Continuous states */
X_nestedloop_simulation_T nestedloop_simulation_X;

/* Block states (auto storage) */
DW_nestedloop_simulation_T nestedloop_simulation_DW;

/* Real-time model */
RT_MODEL_nestedloop_simulation_T nestedloop_simulation_M_;
RT_MODEL_nestedloop_simulation_T *const nestedloop_simulation_M =
  &nestedloop_simulation_M_;

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 24;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  nestedloop_simulation_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  nestedloop_simulation_step();
  nestedloop_simulation_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  nestedloop_simulation_step();
  nestedloop_simulation_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T lo;
  uint32_T hi;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T y;
  real_T sr;
  real_T si;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

/* Model step function */
void nestedloop_simulation_step(void)
{
  if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    /* set solver stop time */
    if (!(nestedloop_simulation_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&nestedloop_simulation_M->solverInfo,
                            ((nestedloop_simulation_M->Timing.clockTickH0 + 1) *
        nestedloop_simulation_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&nestedloop_simulation_M->solverInfo,
                            ((nestedloop_simulation_M->Timing.clockTick0 + 1) *
        nestedloop_simulation_M->Timing.stepSize0 +
        nestedloop_simulation_M->Timing.clockTickH0 *
        nestedloop_simulation_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(nestedloop_simulation_M)) {
    nestedloop_simulation_M->Timing.t[0] = rtsiGetT
      (&nestedloop_simulation_M->solverInfo);
  }

  {
    real_T *lastU;
    int_T ci;
    real_T rtb_Controller1;
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Constant: '<S3>/Constant' */
      nestedloop_simulation_B.Constant = nestedloop_simulation_P.Constant_Value;
    }

    /* S-Function (ref3): '<S3>/S-Function' */

    /* Level2 S-Function Block: '<S3>/S-Function' (ref3) */
    {
      SimStruct *rts = nestedloop_simulation_M->childSfunctions[0];
      sfcnOutputs(rts,0);
    }

    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Gain: '<Root>/Gain2' incorporates:
       *  RandomNumber: '<Root>/Random Number2'
       */
      nestedloop_simulation_B.Gain2 = nestedloop_simulation_P.Gain2_Gain *
        nestedloop_simulation_DW.NextOutput;
    }

    /* ManualSwitch: '<Root>/Manual Switch' */
    if (nestedloop_simulation_P.ManualSwitch_CurrentSetting == 1) {
      nestedloop_simulation_B.r = nestedloop_simulation_B.SFunction[2];
    } else {
      nestedloop_simulation_B.r = nestedloop_simulation_B.Gain2;
    }

    /* End of ManualSwitch: '<Root>/Manual Switch' */

    /* TransferFcn: '<S1>/Plant H21_2' */
    nestedloop_simulation_B.PlantH21_2 = 0.0;
    for (ci = 0; ci < 9; ci++) {
      nestedloop_simulation_B.PlantH21_2 +=
        nestedloop_simulation_P.PlantH21_2_C[ci] *
        nestedloop_simulation_X.PlantH21_2_CSTATE[ci];
    }

    /* End of TransferFcn: '<S1>/Plant H21_2' */
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Gain: '<Root>/Gain' incorporates:
       *  RandomNumber: '<Root>/Random Number'
       */
      nestedloop_simulation_B.Gain = nestedloop_simulation_P.Gain_Gain *
        nestedloop_simulation_DW.NextOutput_h;
    }

    /* Sum: '<Root>/Sum' incorporates:
     *  Sum: '<Root>/Sum1'
     */
    nestedloop_simulation_B.e = nestedloop_simulation_B.r -
      (nestedloop_simulation_B.PlantH21_2 + nestedloop_simulation_B.Gain);
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    }

    /* Derivative: '<Root>/Derivative1' */
    if ((nestedloop_simulation_DW.TimeStampA >=
         nestedloop_simulation_M->Timing.t[0]) &&
        (nestedloop_simulation_DW.TimeStampB >=
         nestedloop_simulation_M->Timing.t[0])) {
      nestedloop_simulation_B.Derivative1 = 0.0;
    } else {
      rtb_Controller1 = nestedloop_simulation_DW.TimeStampA;
      lastU = &nestedloop_simulation_DW.LastUAtTimeA;
      if (nestedloop_simulation_DW.TimeStampA <
          nestedloop_simulation_DW.TimeStampB) {
        if (nestedloop_simulation_DW.TimeStampB <
            nestedloop_simulation_M->Timing.t[0]) {
          rtb_Controller1 = nestedloop_simulation_DW.TimeStampB;
          lastU = &nestedloop_simulation_DW.LastUAtTimeB;
        }
      } else {
        if (nestedloop_simulation_DW.TimeStampA >=
            nestedloop_simulation_M->Timing.t[0]) {
          rtb_Controller1 = nestedloop_simulation_DW.TimeStampB;
          lastU = &nestedloop_simulation_DW.LastUAtTimeB;
        }
      }

      nestedloop_simulation_B.Derivative1 = (nestedloop_simulation_B.PlantH21_2
        - *lastU) / (nestedloop_simulation_M->Timing.t[0] - rtb_Controller1);
    }

    /* End of Derivative: '<Root>/Derivative1' */
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    }

    /* TransferFcn: '<S1>/Plant 1' */
    nestedloop_simulation_B.Plant1 = 0.0;
    for (ci = 0; ci < 7; ci++) {
      nestedloop_simulation_B.Plant1 += nestedloop_simulation_P.Plant1_C[ci] *
        nestedloop_simulation_X.Plant1_CSTATE[ci];
    }

    /* End of TransferFcn: '<S1>/Plant 1' */
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    }

    /* Sum: '<Root>/Sum4' */
    nestedloop_simulation_B.Sum4 = nestedloop_simulation_B.Plant1 -
      nestedloop_simulation_B.PlantH21_2;
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Gain: '<Root>/Gain1' incorporates:
       *  RandomNumber: '<Root>/Random Number1'
       */
      nestedloop_simulation_B.d = nestedloop_simulation_P.Gain1_Gain *
        nestedloop_simulation_DW.NextOutput_he;
    }

    /* TransferFcn: '<Root>/ControllerH21_2' */
    rtb_Controller1 = 0.0;
    for (ci = 0; ci < 5; ci++) {
      rtb_Controller1 += nestedloop_simulation_P.ControllerH21_2_C[ci] *
        nestedloop_simulation_X.ControllerH21_2_CSTATE[ci];
    }

    /* End of TransferFcn: '<Root>/ControllerH21_2' */

    /* Sum: '<Root>/Sum2' */
    nestedloop_simulation_B.u = nestedloop_simulation_B.d + rtb_Controller1;
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Gain: '<S1>/Gain' incorporates:
       *  RandomNumber: '<S1>/Random Number'
       */
      nestedloop_simulation_B.Gain_a = nestedloop_simulation_P.Gain_Gain_m *
        nestedloop_simulation_DW.NextOutput_b;
    }

    /* Sum: '<S1>/Sum' incorporates:
     *  Sum: '<S1>/Sum1'
     */
    nestedloop_simulation_B.Sum = nestedloop_simulation_B.u -
      (nestedloop_simulation_B.Plant1 + nestedloop_simulation_B.Gain_a);
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Gain: '<S1>/Gain1' incorporates:
       *  RandomNumber: '<S1>/Random Number1'
       */
      nestedloop_simulation_B.Gain1 = nestedloop_simulation_P.Gain1_Gain_g *
        nestedloop_simulation_DW.NextOutput_k;
    }

    /* Sum: '<S1>/Sum2' incorporates:
     *  TransferFcn: '<S1>/Controller 1'
     */
    nestedloop_simulation_B.Sum2 = (((nestedloop_simulation_P.Controller1_C[0] *
      nestedloop_simulation_X.Controller1_CSTATE[0] +
      nestedloop_simulation_P.Controller1_C[1] *
      nestedloop_simulation_X.Controller1_CSTATE[1]) +
      nestedloop_simulation_P.Controller1_C[2] *
      nestedloop_simulation_X.Controller1_CSTATE[2]) +
      nestedloop_simulation_P.Controller1_D * nestedloop_simulation_B.Sum) +
      nestedloop_simulation_B.Gain1;
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    }

    /* Saturate: '<S1>/Saturation' */
    if (nestedloop_simulation_B.Sum2 >
        nestedloop_simulation_P.Saturation_UpperSat) {
      nestedloop_simulation_B.Saturation =
        nestedloop_simulation_P.Saturation_UpperSat;
    } else if (nestedloop_simulation_B.Sum2 <
               nestedloop_simulation_P.Saturation_LowerSat) {
      nestedloop_simulation_B.Saturation =
        nestedloop_simulation_P.Saturation_LowerSat;
    } else {
      nestedloop_simulation_B.Saturation = nestedloop_simulation_B.Sum2;
    }

    /* End of Saturate: '<S1>/Saturation' */
  }

  if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(nestedloop_simulation_M->rtwLogInfo,
                        (nestedloop_simulation_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    real_T *lastU;
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Update for RandomNumber: '<Root>/Random Number2' */
      nestedloop_simulation_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
        (&nestedloop_simulation_DW.RandSeed) *
        nestedloop_simulation_P.RandomNumber2_StdDev +
        nestedloop_simulation_P.RandomNumber2_Mean;

      /* Update for RandomNumber: '<Root>/Random Number' */
      nestedloop_simulation_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw_snf
        (&nestedloop_simulation_DW.RandSeed_d) *
        nestedloop_simulation_P.RandomNumber_StdDev +
        nestedloop_simulation_P.RandomNumber_Mean;
    }

    /* Update for Derivative: '<Root>/Derivative1' */
    if (nestedloop_simulation_DW.TimeStampA == (rtInf)) {
      nestedloop_simulation_DW.TimeStampA = nestedloop_simulation_M->Timing.t[0];
      lastU = &nestedloop_simulation_DW.LastUAtTimeA;
    } else if (nestedloop_simulation_DW.TimeStampB == (rtInf)) {
      nestedloop_simulation_DW.TimeStampB = nestedloop_simulation_M->Timing.t[0];
      lastU = &nestedloop_simulation_DW.LastUAtTimeB;
    } else if (nestedloop_simulation_DW.TimeStampA <
               nestedloop_simulation_DW.TimeStampB) {
      nestedloop_simulation_DW.TimeStampA = nestedloop_simulation_M->Timing.t[0];
      lastU = &nestedloop_simulation_DW.LastUAtTimeA;
    } else {
      nestedloop_simulation_DW.TimeStampB = nestedloop_simulation_M->Timing.t[0];
      lastU = &nestedloop_simulation_DW.LastUAtTimeB;
    }

    *lastU = nestedloop_simulation_B.PlantH21_2;

    /* End of Update for Derivative: '<Root>/Derivative1' */
    if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
      /* Update for RandomNumber: '<Root>/Random Number1' */
      nestedloop_simulation_DW.NextOutput_he = rt_nrand_Upu32_Yd_f_pw_snf
        (&nestedloop_simulation_DW.RandSeed_b) *
        nestedloop_simulation_P.RandomNumber1_StdDev +
        nestedloop_simulation_P.RandomNumber1_Mean;

      /* Update for RandomNumber: '<S1>/Random Number' */
      nestedloop_simulation_DW.NextOutput_b = rt_nrand_Upu32_Yd_f_pw_snf
        (&nestedloop_simulation_DW.RandSeed_l) *
        nestedloop_simulation_P.RandomNumber_StdDev_k +
        nestedloop_simulation_P.RandomNumber_Mean_k;

      /* Update for RandomNumber: '<S1>/Random Number1' */
      nestedloop_simulation_DW.NextOutput_k = rt_nrand_Upu32_Yd_f_pw_snf
        (&nestedloop_simulation_DW.RandSeed_o) *
        nestedloop_simulation_P.RandomNumber1_StdDev_k +
        nestedloop_simulation_P.RandomNumber1_Mean_o;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(nestedloop_simulation_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(nestedloop_simulation_M)!=-1) &&
          !((rtmGetTFinal(nestedloop_simulation_M)-
             (((nestedloop_simulation_M->Timing.clockTick1+
                nestedloop_simulation_M->Timing.clockTickH1* 4294967296.0)) *
              0.00048828125)) > (((nestedloop_simulation_M->Timing.clockTick1+
              nestedloop_simulation_M->Timing.clockTickH1* 4294967296.0)) *
            0.00048828125) * (DBL_EPSILON))) {
        rtmSetErrorStatus(nestedloop_simulation_M, "Simulation finished");
      }

      if (rtmGetStopRequested(nestedloop_simulation_M)) {
        rtmSetErrorStatus(nestedloop_simulation_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&nestedloop_simulation_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++nestedloop_simulation_M->Timing.clockTick0)) {
      ++nestedloop_simulation_M->Timing.clockTickH0;
    }

    nestedloop_simulation_M->Timing.t[0] = rtsiGetSolverStopTime
      (&nestedloop_simulation_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.00048828125s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.00048828125, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      nestedloop_simulation_M->Timing.clockTick1++;
      if (!nestedloop_simulation_M->Timing.clockTick1) {
        nestedloop_simulation_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void nestedloop_simulation_derivatives(void)
{
  int_T is;
  XDot_nestedloop_simulation_T *_rtXdot;
  _rtXdot = ((XDot_nestedloop_simulation_T *) nestedloop_simulation_M->derivs);

  /* Derivatives for TransferFcn: '<S1>/Plant H21_2' */
  for (is = 0; is < 9; is++) {
    _rtXdot->PlantH21_2_CSTATE[is] = 0.0;
    _rtXdot->PlantH21_2_CSTATE[0] += nestedloop_simulation_P.PlantH21_2_A[is] *
      nestedloop_simulation_X.PlantH21_2_CSTATE[is];
  }

  for (is = 0; is < 8; is++) {
    _rtXdot->PlantH21_2_CSTATE[is + 1] +=
      nestedloop_simulation_X.PlantH21_2_CSTATE[is];
  }

  _rtXdot->PlantH21_2_CSTATE[0] += nestedloop_simulation_B.Plant1;

  /* End of Derivatives for TransferFcn: '<S1>/Plant H21_2' */

  /* Derivatives for TransferFcn: '<S1>/Plant 1' */
  for (is = 0; is < 7; is++) {
    _rtXdot->Plant1_CSTATE[is] = 0.0;
    _rtXdot->Plant1_CSTATE[0] += nestedloop_simulation_P.Plant1_A[is] *
      nestedloop_simulation_X.Plant1_CSTATE[is];
  }

  for (is = 0; is < 6; is++) {
    _rtXdot->Plant1_CSTATE[is + 1] += nestedloop_simulation_X.Plant1_CSTATE[is];
  }

  _rtXdot->Plant1_CSTATE[0] += nestedloop_simulation_B.Saturation;

  /* End of Derivatives for TransferFcn: '<S1>/Plant 1' */

  /* Derivatives for TransferFcn: '<Root>/ControllerH21_2' */
  for (is = 0; is < 5; is++) {
    _rtXdot->ControllerH21_2_CSTATE[is] = 0.0;
    _rtXdot->ControllerH21_2_CSTATE[0] +=
      nestedloop_simulation_P.ControllerH21_2_A[is] *
      nestedloop_simulation_X.ControllerH21_2_CSTATE[is];
  }

  _rtXdot->ControllerH21_2_CSTATE[1] +=
    nestedloop_simulation_X.ControllerH21_2_CSTATE[0];
  _rtXdot->ControllerH21_2_CSTATE[2] +=
    nestedloop_simulation_X.ControllerH21_2_CSTATE[1];
  _rtXdot->ControllerH21_2_CSTATE[3] +=
    nestedloop_simulation_X.ControllerH21_2_CSTATE[2];
  _rtXdot->ControllerH21_2_CSTATE[4] +=
    nestedloop_simulation_X.ControllerH21_2_CSTATE[3];
  _rtXdot->ControllerH21_2_CSTATE[0] += nestedloop_simulation_B.e;

  /* End of Derivatives for TransferFcn: '<Root>/ControllerH21_2' */

  /* Derivatives for TransferFcn: '<S1>/Controller 1' */
  _rtXdot->Controller1_CSTATE[0] = 0.0;
  _rtXdot->Controller1_CSTATE[0] += nestedloop_simulation_P.Controller1_A[0] *
    nestedloop_simulation_X.Controller1_CSTATE[0];
  _rtXdot->Controller1_CSTATE[1] = 0.0;
  _rtXdot->Controller1_CSTATE[0] += nestedloop_simulation_P.Controller1_A[1] *
    nestedloop_simulation_X.Controller1_CSTATE[1];
  _rtXdot->Controller1_CSTATE[2] = 0.0;
  _rtXdot->Controller1_CSTATE[0] += nestedloop_simulation_P.Controller1_A[2] *
    nestedloop_simulation_X.Controller1_CSTATE[2];
  _rtXdot->Controller1_CSTATE[1] += nestedloop_simulation_X.Controller1_CSTATE[0];
  _rtXdot->Controller1_CSTATE[2] += nestedloop_simulation_X.Controller1_CSTATE[1];
  _rtXdot->Controller1_CSTATE[0] += nestedloop_simulation_B.Sum;
}

/* Model initialize function */
void nestedloop_simulation_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)nestedloop_simulation_M, 0,
                sizeof(RT_MODEL_nestedloop_simulation_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&nestedloop_simulation_M->solverInfo,
                          &nestedloop_simulation_M->Timing.simTimeStep);
    rtsiSetTPtr(&nestedloop_simulation_M->solverInfo, &rtmGetTPtr
                (nestedloop_simulation_M));
    rtsiSetStepSizePtr(&nestedloop_simulation_M->solverInfo,
                       &nestedloop_simulation_M->Timing.stepSize0);
    rtsiSetdXPtr(&nestedloop_simulation_M->solverInfo,
                 &nestedloop_simulation_M->derivs);
    rtsiSetContStatesPtr(&nestedloop_simulation_M->solverInfo, (real_T **)
                         &nestedloop_simulation_M->contStates);
    rtsiSetNumContStatesPtr(&nestedloop_simulation_M->solverInfo,
      &nestedloop_simulation_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&nestedloop_simulation_M->solverInfo,
      &nestedloop_simulation_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&nestedloop_simulation_M->solverInfo,
      &nestedloop_simulation_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&nestedloop_simulation_M->solverInfo,
      &nestedloop_simulation_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&nestedloop_simulation_M->solverInfo,
                          (&rtmGetErrorStatus(nestedloop_simulation_M)));
    rtsiSetRTModelPtr(&nestedloop_simulation_M->solverInfo,
                      nestedloop_simulation_M);
  }

  rtsiSetSimTimeStep(&nestedloop_simulation_M->solverInfo, MAJOR_TIME_STEP);
  nestedloop_simulation_M->intgData.y = nestedloop_simulation_M->odeY;
  nestedloop_simulation_M->intgData.f[0] = nestedloop_simulation_M->odeF[0];
  nestedloop_simulation_M->intgData.f[1] = nestedloop_simulation_M->odeF[1];
  nestedloop_simulation_M->intgData.f[2] = nestedloop_simulation_M->odeF[2];
  nestedloop_simulation_M->contStates = ((X_nestedloop_simulation_T *)
    &nestedloop_simulation_X);
  rtsiSetSolverData(&nestedloop_simulation_M->solverInfo, (void *)
                    &nestedloop_simulation_M->intgData);
  rtsiSetSolverName(&nestedloop_simulation_M->solverInfo,"ode3");
  nestedloop_simulation_M->solverInfoPtr = (&nestedloop_simulation_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = nestedloop_simulation_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    nestedloop_simulation_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    nestedloop_simulation_M->Timing.sampleTimes =
      (&nestedloop_simulation_M->Timing.sampleTimesArray[0]);
    nestedloop_simulation_M->Timing.offsetTimes =
      (&nestedloop_simulation_M->Timing.offsetTimesArray[0]);

    /* task periods */
    nestedloop_simulation_M->Timing.sampleTimes[0] = (0.0);
    nestedloop_simulation_M->Timing.sampleTimes[1] = (0.00048828125);

    /* task offsets */
    nestedloop_simulation_M->Timing.offsetTimes[0] = (0.0);
    nestedloop_simulation_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(nestedloop_simulation_M, &nestedloop_simulation_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = nestedloop_simulation_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    nestedloop_simulation_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(nestedloop_simulation_M, 10.0);
  nestedloop_simulation_M->Timing.stepSize0 = 0.00048828125;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    nestedloop_simulation_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(nestedloop_simulation_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(nestedloop_simulation_M->rtwLogInfo, (NULL));
    rtliSetLogT(nestedloop_simulation_M->rtwLogInfo, "tout");
    rtliSetLogX(nestedloop_simulation_M->rtwLogInfo, "");
    rtliSetLogXFinal(nestedloop_simulation_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(nestedloop_simulation_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(nestedloop_simulation_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(nestedloop_simulation_M->rtwLogInfo, 0);
    rtliSetLogDecimation(nestedloop_simulation_M->rtwLogInfo, 1);
    rtliSetLogY(nestedloop_simulation_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(nestedloop_simulation_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(nestedloop_simulation_M->rtwLogInfo, (NULL));
  }

  nestedloop_simulation_M->solverInfoPtr = (&nestedloop_simulation_M->solverInfo);
  nestedloop_simulation_M->Timing.stepSize = (0.00048828125);
  rtsiSetFixedStepSize(&nestedloop_simulation_M->solverInfo, 0.00048828125);
  rtsiSetSolverMode(&nestedloop_simulation_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) memset(((void *) &nestedloop_simulation_B), 0,
                sizeof(B_nestedloop_simulation_T));

  /* states (continuous) */
  {
    (void) memset((void *)&nestedloop_simulation_X, 0,
                  sizeof(X_nestedloop_simulation_T));
  }

  /* states (dwork) */
  (void) memset((void *)&nestedloop_simulation_DW, 0,
                sizeof(DW_nestedloop_simulation_T));

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &nestedloop_simulation_M->NonInlinedSFcns.sfcnInfo;
    nestedloop_simulation_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, (&rtmGetErrorStatus(nestedloop_simulation_M)));
    rtssSetNumRootSampTimesPtr(sfcnInfo,
      &nestedloop_simulation_M->Sizes.numSampTimes);
    nestedloop_simulation_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (nestedloop_simulation_M)[0]);
    nestedloop_simulation_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr
      (nestedloop_simulation_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,
                   nestedloop_simulation_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(nestedloop_simulation_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(nestedloop_simulation_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (nestedloop_simulation_M));
    rtssSetStepSizePtr(sfcnInfo, &nestedloop_simulation_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested
      (nestedloop_simulation_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &nestedloop_simulation_M->derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &nestedloop_simulation_M->zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo, &nestedloop_simulation_M->blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &nestedloop_simulation_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &nestedloop_simulation_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &nestedloop_simulation_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &nestedloop_simulation_M->solverInfoPtr);
  }

  nestedloop_simulation_M->Sizes.numSFcns = (1);

  /* register each child */
  {
    (void) memset((void *)
                  &nestedloop_simulation_M->NonInlinedSFcns.childSFunctions[0],
                  0,
                  1*sizeof(SimStruct));
    nestedloop_simulation_M->childSfunctions =
      (&nestedloop_simulation_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    nestedloop_simulation_M->childSfunctions[0] =
      (&nestedloop_simulation_M->NonInlinedSFcns.childSFunctions[0]);

    /* Level2 S-Function Block: nestedloop_simulation/<S3>/S-Function (ref3) */
    {
      SimStruct *rts = nestedloop_simulation_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod =
        nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset =
        nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap =
        nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &nestedloop_simulation_M->
                         NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, nestedloop_simulation_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts,
                           &nestedloop_simulation_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts,
                           &nestedloop_simulation_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory of model methods 4 */
      {
        ssSetModelMethods4(rts,
                           &nestedloop_simulation_M->NonInlinedSFcns.methods4[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts,
                         &nestedloop_simulation_M->NonInlinedSFcns.statesInfo2[0]);
        ssSetPeriodicStatesInfo(rts,
          &nestedloop_simulation_M->NonInlinedSFcns.periodicStatesInfo[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 1);
        ssSetPortInfoForInputs(rts,
          &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &nestedloop_simulation_B.Constant;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 3);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            nestedloop_simulation_B.SFunction));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts, "nestedloop_simulation/Subsystem1/S-Function");
      ssSetRTModel(rts,nestedloop_simulation_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 1);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)
                       nestedloop_simulation_P.SFunction_P1_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &nestedloop_simulation_DW.SFunction_RWORK[0]);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &nestedloop_simulation_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 1);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 50);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &nestedloop_simulation_DW.SFunction_RWORK[0]);
      }

      /* registration */
      ref3(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 0.0);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 0;

      /* set compiled values of dynamic vector attributes */
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
    }
  }

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(nestedloop_simulation_M->rtwLogInfo, 0.0,
    rtmGetTFinal(nestedloop_simulation_M),
    nestedloop_simulation_M->Timing.stepSize0, (&rtmGetErrorStatus
    (nestedloop_simulation_M)));

  /* Start for Constant: '<S3>/Constant' */
  nestedloop_simulation_B.Constant = nestedloop_simulation_P.Constant_Value;

  {
    uint32_T tseed;
    int32_T r;
    int32_T t;
    real_T tmp;

    /* InitializeConditions for S-Function (ref3): '<S3>/S-Function' */
    /* Level2 S-Function Block: '<S3>/S-Function' (ref3) */
    {
      SimStruct *rts = nestedloop_simulation_M->childSfunctions[0];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for RandomNumber: '<Root>/Random Number2' */
    tmp = floor(nestedloop_simulation_P.RandomNumber2_Seed);
    if (rtIsNaN(tmp) || rtIsInf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-tmp : (uint32_T)tmp;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    nestedloop_simulation_DW.RandSeed = tseed;
    nestedloop_simulation_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf
      (&nestedloop_simulation_DW.RandSeed) *
      nestedloop_simulation_P.RandomNumber2_StdDev +
      nestedloop_simulation_P.RandomNumber2_Mean;

    /* End of InitializeConditions for RandomNumber: '<Root>/Random Number2' */

    /* InitializeConditions for TransferFcn: '<S1>/Plant H21_2' */
    memset(&nestedloop_simulation_X.PlantH21_2_CSTATE[0], 0, 9U * sizeof(real_T));

    /* InitializeConditions for RandomNumber: '<Root>/Random Number' */
    tmp = floor(nestedloop_simulation_P.RandomNumber_Seed);
    if (rtIsNaN(tmp) || rtIsInf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-tmp : (uint32_T)tmp;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    nestedloop_simulation_DW.RandSeed_d = tseed;
    nestedloop_simulation_DW.NextOutput_h = rt_nrand_Upu32_Yd_f_pw_snf
      (&nestedloop_simulation_DW.RandSeed_d) *
      nestedloop_simulation_P.RandomNumber_StdDev +
      nestedloop_simulation_P.RandomNumber_Mean;

    /* End of InitializeConditions for RandomNumber: '<Root>/Random Number' */

    /* InitializeConditions for Derivative: '<Root>/Derivative1' */
    nestedloop_simulation_DW.TimeStampA = (rtInf);
    nestedloop_simulation_DW.TimeStampB = (rtInf);

    /* InitializeConditions for TransferFcn: '<S1>/Plant 1' */
    for (r = 0; r < 7; r++) {
      nestedloop_simulation_X.Plant1_CSTATE[r] = 0.0;
    }

    /* End of InitializeConditions for TransferFcn: '<S1>/Plant 1' */

    /* InitializeConditions for RandomNumber: '<Root>/Random Number1' */
    tmp = floor(nestedloop_simulation_P.RandomNumber1_Seed);
    if (rtIsNaN(tmp) || rtIsInf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-tmp : (uint32_T)tmp;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    nestedloop_simulation_DW.RandSeed_b = tseed;
    nestedloop_simulation_DW.NextOutput_he = rt_nrand_Upu32_Yd_f_pw_snf
      (&nestedloop_simulation_DW.RandSeed_b) *
      nestedloop_simulation_P.RandomNumber1_StdDev +
      nestedloop_simulation_P.RandomNumber1_Mean;

    /* End of InitializeConditions for RandomNumber: '<Root>/Random Number1' */

    /* InitializeConditions for TransferFcn: '<Root>/ControllerH21_2' */
    for (r = 0; r < 5; r++) {
      nestedloop_simulation_X.ControllerH21_2_CSTATE[r] = 0.0;
    }

    /* End of InitializeConditions for TransferFcn: '<Root>/ControllerH21_2' */

    /* InitializeConditions for RandomNumber: '<S1>/Random Number' */
    tmp = floor(nestedloop_simulation_P.RandomNumber_Seed_a);
    if (rtIsNaN(tmp) || rtIsInf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-tmp : (uint32_T)tmp;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    nestedloop_simulation_DW.RandSeed_l = tseed;
    nestedloop_simulation_DW.NextOutput_b = rt_nrand_Upu32_Yd_f_pw_snf
      (&nestedloop_simulation_DW.RandSeed_l) *
      nestedloop_simulation_P.RandomNumber_StdDev_k +
      nestedloop_simulation_P.RandomNumber_Mean_k;

    /* End of InitializeConditions for RandomNumber: '<S1>/Random Number' */

    /* InitializeConditions for RandomNumber: '<S1>/Random Number1' */
    tmp = floor(nestedloop_simulation_P.RandomNumber1_Seed_p);
    if (rtIsNaN(tmp) || rtIsInf(tmp)) {
      tmp = 0.0;
    } else {
      tmp = fmod(tmp, 4.294967296E+9);
    }

    tseed = tmp < 0.0 ? (uint32_T)-(int32_T)(uint32_T)-tmp : (uint32_T)tmp;
    r = (int32_T)(tseed >> 16U);
    t = (int32_T)(tseed & 32768U);
    tseed = ((((tseed - ((uint32_T)r << 16U)) + t) << 16U) + t) + r;
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else {
      if (tseed > 2147483646U) {
        tseed = 2147483646U;
      }
    }

    nestedloop_simulation_DW.RandSeed_o = tseed;
    nestedloop_simulation_DW.NextOutput_k = rt_nrand_Upu32_Yd_f_pw_snf
      (&nestedloop_simulation_DW.RandSeed_o) *
      nestedloop_simulation_P.RandomNumber1_StdDev_k +
      nestedloop_simulation_P.RandomNumber1_Mean_o;

    /* End of InitializeConditions for RandomNumber: '<S1>/Random Number1' */

    /* InitializeConditions for TransferFcn: '<S1>/Controller 1' */
    nestedloop_simulation_X.Controller1_CSTATE[0] = 0.0;
    nestedloop_simulation_X.Controller1_CSTATE[1] = 0.0;
    nestedloop_simulation_X.Controller1_CSTATE[2] = 0.0;
  }
}

/* Model terminate function */
void nestedloop_simulation_terminate(void)
{
  /* Terminate for S-Function (ref3): '<S3>/S-Function' */
  /* Level2 S-Function Block: '<S3>/S-Function' (ref3) */
  {
    SimStruct *rts = nestedloop_simulation_M->childSfunctions[0];
    sfcnTerminate(rts);
  }
}
