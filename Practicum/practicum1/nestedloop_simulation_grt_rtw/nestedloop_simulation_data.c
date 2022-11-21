/*
 * nestedloop_simulation_data.c
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

/* Block parameters (auto storage) */
P_nestedloop_simulation_T nestedloop_simulation_P = {
  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S3>/S-Function'
   */
  { 3.0, 6.0 },

  /*  Variable: ref_part
   * Referenced by: '<S3>/S-Function'
   */
  { 0.0, 26.0, 0.0, 0.0, 0.89767, -1.0, 26.0, 0.0, 0.0, 30.0, 30.0, 0.0, 0.031,
    0.031, 0.0, 1.0E+6, 1.0E+6, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S3>/Constant'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  1.0,                                 /* Computed Parameter: RandomNumber2_StdDev
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number2'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<Root>/Gain2'
                                        */

  /*  Computed Parameter: PlantH21_2_A
   * Referenced by: '<S1>/Plant H21_2'
   */
  { -6688.09523809524, -1.886147770219199E+7, -2.2308875451415142E+10,
    -4.1680977576215684E+12, -4.2118585705887295E+15, -2.0380168528316675E+17,
    -1.8289894833104712E+20, -0.0, -0.0 },

  /*  Computed Parameter: PlantH21_2_C
   * Referenced by: '<S1>/Plant H21_2'
   */
  { -7.1428571428571432, -19149.659863945566, 3.1094482237339383E+8,
    -1.0856756529772404E+12, 1.397827608409619E+15, -1.2541642171271803E+17,
    1.8289894833104712E+20, 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number'
                                        */
  1.0,                                 /* Computed Parameter: RandomNumber_StdDev
                                        * Referenced by: '<Root>/Random Number'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number'
                                        */
  0.001,                               /* Expression: 1e-3
                                        * Referenced by: '<Root>/Gain'
                                        */

  /*  Computed Parameter: Plant1_A
   * Referenced by: '<S1>/Plant 1'
   */
  { -6680.9523809523826, -1.8747089947089948E+7, -2.1729570840681953E+10,
    -2.7630805408583184E+12, -2.7434842249657065E+15, -0.0, -0.0 },

  /*  Computed Parameter: Plant1_C
   * Referenced by: '<S1>/Plant 1'
   */
  { 0.0, -4761.9047619047615, 3.1712018140589572E+7, -8.8274124464600662E+10,
    9.9468099997200531E+13, -5.179026343047507E+15, 6.5321052975373967E+18 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  1.0,                                 /* Computed Parameter: RandomNumber1_StdDev
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Random Number1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Gain1'
                                        */

  /*  Computed Parameter: ControllerH21_2_A
   * Referenced by: '<Root>/ControllerH21_2'
   */
  { -254.46900494077329, -900897.48973143671, -5.8639070457783014E+7,
    -1.7533636386120436E+8, -0.0 },

  /*  Computed Parameter: ControllerH21_2_C
   * Referenced by: '<Root>/ControllerH21_2'
   */
  { 0.0, 40118.833050459209, 457514.42319036252, 2.6313362777747126E+9,
    1.6525062978405193E+10 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/Random Number'
                                        */
  1.4142135623730951,                  /* Computed Parameter: RandomNumber_StdDev_k
                                        * Referenced by: '<S1>/Random Number'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/Random Number'
                                        */
  0.0001,                              /* Expression: 1e-4
                                        * Referenced by: '<S1>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  1.0,                                 /* Computed Parameter: RandomNumber1_StdDev_k
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/Random Number1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S1>/Gain1'
                                        */

  /*  Computed Parameter: Controller1_A
   * Referenced by: '<S1>/Controller 1'
   */
  { -3148.3068654090453, -87780.656327464792, -2.0950445052451703E+8 },

  /*  Computed Parameter: Controller1_C
   * Referenced by: '<S1>/Controller 1'
   */
  { -630227.56238842639, 9.4447520687587969E+6, -4.14553810068265E+10 },
  204.3735100998498,                   /* Computed Parameter: Controller1_D
                                        * Referenced by: '<S1>/Controller 1'
                                        */
  2.5,                                 /* Expression: 2.5
                                        * Referenced by: '<S1>/Saturation'
                                        */
  -2.5,                                /* Expression: -2.5
                                        * Referenced by: '<S1>/Saturation'
                                        */
  1U                                   /* Computed Parameter: ManualSwitch_CurrentSetting
                                        * Referenced by: '<Root>/Manual Switch'
                                        */
};
