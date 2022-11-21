/*
 * nestedloop_simulation.h
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

#ifndef RTW_HEADER_nestedloop_simulation_h_
#define RTW_HEADER_nestedloop_simulation_h_
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef nestedloop_simulation_COMMON_INCLUDES_
# define nestedloop_simulation_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "rt_logging.h"
#endif                                 /* nestedloop_simulation_COMMON_INCLUDES_ */

#include "nestedloop_simulation_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->blkStateChange = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals (auto storage) */
typedef struct {
  real_T Constant;                     /* '<S3>/Constant' */
  real_T SFunction[3];                 /* '<S3>/S-Function' */
  real_T Gain2;                        /* '<Root>/Gain2' */
  real_T r;                            /* '<Root>/Manual Switch' */
  real_T PlantH21_2;                   /* '<S1>/Plant H21_2' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T e;                            /* '<Root>/Sum' */
  real_T Derivative1;                  /* '<Root>/Derivative1' */
  real_T Plant1;                       /* '<S1>/Plant 1' */
  real_T Sum4;                         /* '<Root>/Sum4' */
  real_T d;                            /* '<Root>/Gain1' */
  real_T u;                            /* '<Root>/Sum2' */
  real_T Gain_a;                       /* '<S1>/Gain' */
  real_T Sum;                          /* '<S1>/Sum' */
  real_T Gain1;                        /* '<S1>/Gain1' */
  real_T Sum2;                         /* '<S1>/Sum2' */
  real_T Saturation;                   /* '<S1>/Saturation' */
} B_nestedloop_simulation_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T NextOutput;                   /* '<Root>/Random Number2' */
  real_T NextOutput_h;                 /* '<Root>/Random Number' */
  real_T TimeStampA;                   /* '<Root>/Derivative1' */
  real_T LastUAtTimeA;                 /* '<Root>/Derivative1' */
  real_T TimeStampB;                   /* '<Root>/Derivative1' */
  real_T LastUAtTimeB;                 /* '<Root>/Derivative1' */
  real_T NextOutput_he;                /* '<Root>/Random Number1' */
  real_T NextOutput_b;                 /* '<S1>/Random Number' */
  real_T NextOutput_k;                 /* '<S1>/Random Number1' */
  real_T SFunction_RWORK[50];          /* '<S3>/S-Function' */
  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope10_PWORK;                     /* '<Root>/Scope10' */

  struct {
    void *LoggedData;
  } Scope11_PWORK;                     /* '<Root>/Scope11' */

  struct {
    void *LoggedData;
  } Scope13_PWORK;                     /* '<Root>/Scope13' */

  struct {
    void *LoggedData;
  } Scope14_PWORK;                     /* '<Root>/Scope14' */

  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<Root>/Scope2' */

  struct {
    void *LoggedData;
  } Scope9_PWORK;                      /* '<Root>/Scope9' */

  struct {
    void *LoggedData;
  } Scope1_PWORK_a;                    /* '<S1>/Scope1' */

  struct {
    void *LoggedData;
  } Scope2_PWORK_n;                    /* '<S1>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<S1>/Scope3' */

  uint32_T RandSeed;                   /* '<Root>/Random Number2' */
  uint32_T RandSeed_d;                 /* '<Root>/Random Number' */
  uint32_T RandSeed_b;                 /* '<Root>/Random Number1' */
  uint32_T RandSeed_l;                 /* '<S1>/Random Number' */
  uint32_T RandSeed_o;                 /* '<S1>/Random Number1' */
} DW_nestedloop_simulation_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T PlantH21_2_CSTATE[9];         /* '<S1>/Plant H21_2' */
  real_T Plant1_CSTATE[7];             /* '<S1>/Plant 1' */
  real_T ControllerH21_2_CSTATE[5];    /* '<Root>/ControllerH21_2' */
  real_T Controller1_CSTATE[3];        /* '<S1>/Controller 1' */
} X_nestedloop_simulation_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T PlantH21_2_CSTATE[9];         /* '<S1>/Plant H21_2' */
  real_T Plant1_CSTATE[7];             /* '<S1>/Plant 1' */
  real_T ControllerH21_2_CSTATE[5];    /* '<Root>/ControllerH21_2' */
  real_T Controller1_CSTATE[3];        /* '<S1>/Controller 1' */
} XDot_nestedloop_simulation_T;

/* State disabled  */
typedef struct {
  boolean_T PlantH21_2_CSTATE[9];      /* '<S1>/Plant H21_2' */
  boolean_T Plant1_CSTATE[7];          /* '<S1>/Plant 1' */
  boolean_T ControllerH21_2_CSTATE[5]; /* '<Root>/ControllerH21_2' */
  boolean_T Controller1_CSTATE[3];     /* '<S1>/Controller 1' */
} XDis_nestedloop_simulation_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Parameters (auto storage) */
struct P_nestedloop_simulation_T_ {
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S3>/S-Function'
                                        */
  real_T ref_part[18];                 /* Variable: ref_part
                                        * Referenced by: '<S3>/S-Function'
                                        */
  real_T Constant_Value;               /* Expression: 1
                                        * Referenced by: '<S3>/Constant'
                                        */
  real_T RandomNumber2_Mean;           /* Expression: 0
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  real_T RandomNumber2_StdDev;         /* Computed Parameter: RandomNumber2_StdDev
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  real_T RandomNumber2_Seed;           /* Expression: 0
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  real_T Gain2_Gain;                   /* Expression: 0.1
                                        * Referenced by: '<Root>/Gain2'
                                        */
  real_T PlantH21_2_A[9];              /* Computed Parameter: PlantH21_2_A
                                        * Referenced by: '<S1>/Plant H21_2'
                                        */
  real_T PlantH21_2_C[9];              /* Computed Parameter: PlantH21_2_C
                                        * Referenced by: '<S1>/Plant H21_2'
                                        */
  real_T RandomNumber_Mean;            /* Expression: 0
                                        * Referenced by: '<Root>/Random Number'
                                        */
  real_T RandomNumber_StdDev;          /* Computed Parameter: RandomNumber_StdDev
                                        * Referenced by: '<Root>/Random Number'
                                        */
  real_T RandomNumber_Seed;            /* Expression: 0
                                        * Referenced by: '<Root>/Random Number'
                                        */
  real_T Gain_Gain;                    /* Expression: 1e-3
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Plant1_A[7];                  /* Computed Parameter: Plant1_A
                                        * Referenced by: '<S1>/Plant 1'
                                        */
  real_T Plant1_C[7];                  /* Computed Parameter: Plant1_C
                                        * Referenced by: '<S1>/Plant 1'
                                        */
  real_T RandomNumber1_Mean;           /* Expression: 0
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  real_T RandomNumber1_StdDev;         /* Computed Parameter: RandomNumber1_StdDev
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  real_T RandomNumber1_Seed;           /* Expression: 0
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  real_T Gain1_Gain;                   /* Expression: 0
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T ControllerH21_2_A[5];         /* Computed Parameter: ControllerH21_2_A
                                        * Referenced by: '<Root>/ControllerH21_2'
                                        */
  real_T ControllerH21_2_C[5];         /* Computed Parameter: ControllerH21_2_C
                                        * Referenced by: '<Root>/ControllerH21_2'
                                        */
  real_T RandomNumber_Mean_k;          /* Expression: 0
                                        * Referenced by: '<S1>/Random Number'
                                        */
  real_T RandomNumber_StdDev_k;        /* Computed Parameter: RandomNumber_StdDev_k
                                        * Referenced by: '<S1>/Random Number'
                                        */
  real_T RandomNumber_Seed_a;          /* Expression: 0
                                        * Referenced by: '<S1>/Random Number'
                                        */
  real_T Gain_Gain_m;                  /* Expression: 1e-4
                                        * Referenced by: '<S1>/Gain'
                                        */
  real_T RandomNumber1_Mean_o;         /* Expression: 0
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  real_T RandomNumber1_StdDev_k;       /* Computed Parameter: RandomNumber1_StdDev_k
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  real_T RandomNumber1_Seed_p;         /* Expression: 0
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: 0
                                        * Referenced by: '<S1>/Gain1'
                                        */
  real_T Controller1_A[3];             /* Computed Parameter: Controller1_A
                                        * Referenced by: '<S1>/Controller 1'
                                        */
  real_T Controller1_C[3];             /* Computed Parameter: Controller1_C
                                        * Referenced by: '<S1>/Controller 1'
                                        */
  real_T Controller1_D;                /* Computed Parameter: Controller1_D
                                        * Referenced by: '<S1>/Controller 1'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 2.5
                                        * Referenced by: '<S1>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -2.5
                                        * Referenced by: '<S1>/Saturation'
                                        */
  uint8_T ManualSwitch_CurrentSetting; /* Computed Parameter: ManualSwitch_CurrentSetting
                                        * Referenced by: '<Root>/Manual Switch'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_nestedloop_simulation_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[1];
    SimStruct *childSFunctionPtrs[1];
    struct _ssBlkInfo2 blkInfo2[1];
    struct _ssSFcnModelMethods2 methods2[1];
    struct _ssSFcnModelMethods3 methods3[1];
    struct _ssSFcnModelMethods4 methods4[1];
    struct _ssStatesInfo2 statesInfo2[1];
    ssPeriodicStatesInfo periodicStatesInfo[1];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[1];
      real_T const *UPtrs0[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[1];
      mxArray *params[1];
      struct _ssDWorkRecord dWork[1];
      struct _ssDWorkAuxRecord dWorkAux[1];
    } Sfcn0;
  } NonInlinedSFcns;

  X_nestedloop_simulation_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T blkStateChange;
  real_T odeY[24];
  real_T odeF[3][24];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T options;
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_nestedloop_simulation_T nestedloop_simulation_P;

/* Block signals (auto storage) */
extern B_nestedloop_simulation_T nestedloop_simulation_B;

/* Continuous states (auto storage) */
extern X_nestedloop_simulation_T nestedloop_simulation_X;

/* Block states (auto storage) */
extern DW_nestedloop_simulation_T nestedloop_simulation_DW;

/* Model entry point functions */
extern void nestedloop_simulation_initialize(void);
extern void nestedloop_simulation_step(void);
extern void nestedloop_simulation_terminate(void);

/* Real-time Model object */
extern RT_MODEL_nestedloop_simulation_T *const nestedloop_simulation_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'nestedloop_simulation'
 * '<S1>'   : 'nestedloop_simulation/Motor to mass1'
 * '<S2>'   : 'nestedloop_simulation/Subsystem'
 * '<S3>'   : 'nestedloop_simulation/Subsystem1'
 */
#endif                                 /* RTW_HEADER_nestedloop_simulation_h_ */
