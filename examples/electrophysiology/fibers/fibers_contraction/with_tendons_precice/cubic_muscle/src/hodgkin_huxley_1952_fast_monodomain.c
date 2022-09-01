#include <math.h>
#include <vc_or_std_simd.h>  // this includes <Vc/Vc> or a Vc-emulating wrapper of <experimental/simd> if available
#include <iostream> 
/*
   There are a total of 9 entries in the algebraic variable array.
   There are a total of 4 entries in each of the rate and state variable arrays.
   There are a total of 9 entries in the constant variable array.
 */
/*
 * VOI is time in component environment (millisecond).
 * STATES[0] is V in component membrane (millivolt).
 * CONSTANTS[0] is E_R in component membrane (millivolt).
 * CONSTANTS[1] is Cm in component membrane (microF_per_cm2).
 * ALGEBRAIC[0] is i_Na in component sodium_channel (microA_per_cm2).
 * ALGEBRAIC[4] is i_K in component potassium_channel (microA_per_cm2).
 * ALGEBRAIC[8] is i_L in component leakage_current (microA_per_cm2).
 * CONSTANTS[2] is i_Stim in component membrane (microA_per_cm2).
 * CONSTANTS[3] is g_Na in component sodium_channel (milliS_per_cm2).
 * CONSTANTS[6] is E_Na in component sodium_channel (millivolt).
 * STATES[1] is m in component sodium_channel_m_gate (dimensionless).
 * STATES[2] is h in component sodium_channel_h_gate (dimensionless).
 * ALGEBRAIC[1] is alpha_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[5] is beta_m in component sodium_channel_m_gate (per_millisecond).
 * ALGEBRAIC[2] is alpha_h in component sodium_channel_h_gate (per_millisecond).
 * ALGEBRAIC[6] is beta_h in component sodium_channel_h_gate (per_millisecond).
 * CONSTANTS[4] is g_K in component potassium_channel (milliS_per_cm2).
 * CONSTANTS[7] is E_K in component potassium_channel (millivolt).
 * STATES[3] is n in component potassium_channel_n_gate (dimensionless).
 * ALGEBRAIC[3] is alpha_n in component potassium_channel_n_gate (per_millisecond).
 * ALGEBRAIC[7] is beta_n in component potassium_channel_n_gate (per_millisecond).
 * CONSTANTS[5] is g_L in component leakage_current (milliS_per_cm2).
 * CONSTANTS[8] is E_L in component leakage_current (millivolt).
 * RATES[0] is d/dt V in component membrane (millivolt).
 * RATES[1] is d/dt m in component sodium_channel_m_gate (dimensionless).
 * RATES[2] is d/dt h in component sodium_channel_h_gate (dimensionless).
 * RATES[3] is d/dt n in component potassium_channel_n_gate (dimensionless).
 */

using Vc::double_v; 

/* This file was created by opendihu at 2022/8/26 11:40:05.
 * It is designed for the FastMonodomainSolver.
  */

// helper functions
Vc::double_v exponential(Vc::double_v x);
Vc::double_v pow2(Vc::double_v x);
Vc::double_v pow3(Vc::double_v x);
Vc::double_v pow4(Vc::double_v x);

Vc::double_v exponential(Vc::double_v x)
{
  //return Vc::exp(x);
  // it was determined the x is always in the range [-12,+12] for the Hodgkin-Huxley model

  // exp(x) = lim n→∞ (1 + x/n)^n, we set n=1024
  x = 1.0 + x / 1024.;
  for (int i = 0; i < 10; i++)
  {
    x *= x;
  }
  return x;

  // relative error of this implementation:
  // x    rel error
  // 0    0
  // 1    0.00048784455634225593
  // 3    0.0043763626896140342
  // 5    0.012093715791500804
  // 9    0.038557535762274039
  // 12   0.067389808619653505
}

Vc::double_v pow2(Vc::double_v x)
{
  return x*x;
}
Vc::double_v pow3(Vc::double_v x)
{
  return x*(pow2(x));
}

Vc::double_v pow4(Vc::double_v x)
{
  return pow2(pow2(x));
}

// set initial values for all states
#ifdef __cplusplus
extern "C"
#endif

void initializeStates(Vc::double_v states[]) 
{
  states[0] = -75;
  states[1] = 0.05;
  states[2] = 0.6;
  states[3] = 0.325;
}

// compute one Heun step
#ifdef __cplusplus
extern "C"
#endif

void compute0DInstance(Vc::double_v states[], std::vector<Vc::double_v> &parameters, double currentTime, double timeStepWidth, bool stimulate,
                       bool storeAlgebraicsForTransfer, std::vector<Vc::double_v> &algebraicsForTransfer, const std::vector<int> &algebraicsForTransferIndices, double valueForStimulatedPoint) 
{
  // assert that Vc::double_v::size() is the same as in opendihu, otherwise there will be problems
  if (Vc::double_v::size() != 4)
  {
    std::cout << "Fatal error in compiled library of source file \"src/hodgkin_huxley_1952_fast_monodomain.c\", size of SIMD register in compiled code (" << Vc::double_v::size() << ") does not match opendihu code (4)." << std::endl;
    std::cout << "Delete library such that it will be regenerated with the correct compile options!" << std::endl;
    exit(1);
  }

  // define constants
  const double constant0 = -75;
  const double constant1 = 1;
  const double constant2 = 0;
  const double constant3 = 120;
  const double constant4 = 36;
  const double constant5 = 0.3;
  const double constant6 = constant0+115.000;
  const double constant7 = constant0 - 12.0000;
  const double constant8 = constant0+10.6130;

  // compute new rates, rhs(y_n)
  const double_v algebraic1 = ( - 0.100000*(states[0]+50.0000))/(exponential(- (states[0]+50.0000)/10.0000) - 1.00000);
  const double_v algebraic5 =  4.00000*exponential(- (states[0]+75.0000)/18.0000);
  const double_v rate1 =  algebraic1*(1.00000 - states[1]) -  algebraic5*states[1];
  const double_v algebraic2 =  0.0700000*exponential(- (states[0]+75.0000)/20.0000);
  const double_v algebraic6 = 1.00000/(exponential(- (states[0]+45.0000)/10.0000)+1.00000);
  const double_v rate2 =  algebraic2*(1.00000 - states[2]) -  algebraic6*states[2];
  const double_v algebraic3 = ( - 0.0100000*(states[0]+65.0000))/(exponential(- (states[0]+65.0000)/10.0000) - 1.00000);
  const double_v algebraic7 =  0.125000*exponential((states[0]+75.0000)/80.0000);
  const double_v rate3 =  algebraic3*(1.00000 - states[3]) -  algebraic7*states[3];
  const double_v algebraic0 =  constant3*pow3(states[1])*states[2]*(states[0] - constant6);
  const double_v algebraic4 =  constant4*pow4(states[3])*(states[0] - constant7);
  const double_v algebraic8 =  constant5*(states[0] - constant8);
  const double_v rate0 = - (- parameters[0]+algebraic0+algebraic4+algebraic8)/constant1;

  // algebraic step
  // compute y* = y_n + dt*rhs(y_n), y_n = state, rhs(y_n) = rate, y* = algebraicState
  double_v algebraicState0 = states[0] + timeStepWidth*rate0;
  const double_v algebraicState1 = states[1] + timeStepWidth*rate1;
  const double_v algebraicState2 = states[2] + timeStepWidth*rate2;
  const double_v algebraicState3 = states[3] + timeStepWidth*rate3;



  // if stimulation, set value of Vm (state0)
  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::size()); i++)
    {
      algebraicState0[i] = valueForStimulatedPoint;
    }
  }
  // compute new rates, rhs(y*)
  const double_v algebraicAlgebraic1 = ( - 0.100000*(algebraicState0+50.0000))/(exponential(- (algebraicState0+50.0000)/10.0000) - 1.00000);
  const double_v algebraicAlgebraic5 =  4.00000*exponential(- (algebraicState0+75.0000)/18.0000);
  const double_v algebraicRate1 =  algebraicAlgebraic1*(1.00000 - algebraicState1) -  algebraicAlgebraic5*algebraicState1;
  const double_v algebraicAlgebraic2 =  0.0700000*exponential(- (algebraicState0+75.0000)/20.0000);
  const double_v algebraicAlgebraic6 = 1.00000/(exponential(- (algebraicState0+45.0000)/10.0000)+1.00000);
  const double_v algebraicRate2 =  algebraicAlgebraic2*(1.00000 - algebraicState2) -  algebraicAlgebraic6*algebraicState2;
  const double_v algebraicAlgebraic3 = ( - 0.0100000*(algebraicState0+65.0000))/(exponential(- (algebraicState0+65.0000)/10.0000) - 1.00000);
  const double_v algebraicAlgebraic7 =  0.125000*exponential((algebraicState0+75.0000)/80.0000);
  const double_v algebraicRate3 =  algebraicAlgebraic3*(1.00000 - algebraicState3) -  algebraicAlgebraic7*algebraicState3;
  const double_v algebraicAlgebraic0 =  constant3*pow3(algebraicState1)*algebraicState2*(algebraicState0 - constant6);
  const double_v algebraicAlgebraic4 =  constant4*pow4(algebraicState3)*(algebraicState0 - constant7);
  const double_v algebraicAlgebraic8 =  constant5*(algebraicState0 - constant8);
  const double_v algebraicRate0 = - (- parameters[0]+algebraicAlgebraic0+algebraicAlgebraic4+algebraicAlgebraic8)/constant1;


  // final step
  // y_n+1 = y_n + 0.5*[rhs(y_n) + rhs(y*)]
  states[0] += 0.5*timeStepWidth*(rate0 + algebraicRate0);
  states[1] += 0.5*timeStepWidth*(rate1 + algebraicRate1);
  states[2] += 0.5*timeStepWidth*(rate2 + algebraicRate2);
  states[3] += 0.5*timeStepWidth*(rate3 + algebraicRate3);

  if (stimulate)
  {
    for (int i = 0; i < std::min(3,(int)Vc::double_v::size()); i++)
    {
      states[0][i] = valueForStimulatedPoint;
    }
  }
  // store algebraics for transfer
  if (storeAlgebraicsForTransfer)
  {
    for (int i = 0; i < algebraicsForTransferIndices.size(); i++)
    {
      const int algebraic = algebraicsForTransferIndices[i];
      switch (algebraic)
      {
        case 0:
          algebraicsForTransfer[i] = algebraicAlgebraic0;
          break;
        case 1:
          algebraicsForTransfer[i] = algebraicAlgebraic1;
          break;
        case 2:
          algebraicsForTransfer[i] = algebraicAlgebraic2;
          break;
        case 3:
          algebraicsForTransfer[i] = algebraicAlgebraic3;
          break;
        case 4:
          algebraicsForTransfer[i] = algebraicAlgebraic4;
          break;
        case 5:
          algebraicsForTransfer[i] = algebraicAlgebraic5;
          break;
        case 6:
          algebraicsForTransfer[i] = algebraicAlgebraic6;
          break;
        case 7:
          algebraicsForTransfer[i] = algebraicAlgebraic7;
          break;
        case 8:
          algebraicsForTransfer[i] = algebraicAlgebraic8;
          break;

      }
    }
  }
}

// compute the rhs for a single instance, this can be used for computation of the equilibrium values of the states
#ifdef __cplusplus
extern "C"
#endif
void computeCellMLRightHandSideSingleInstance(void *context, double t, double *states, double *rates, double *algebraics, double *parameters)
{
  double VOI = t;   /* current simulation time */

  /* define constants */
  double CONSTANTS[9];
  CONSTANTS[0] = -75;
  CONSTANTS[1] = 1;
  CONSTANTS[2] = 0;
  CONSTANTS[3] = 120;
  CONSTANTS[4] = 36;
  CONSTANTS[5] = 0.3;
  CONSTANTS[6] = CONSTANTS[0]+115.000;
  CONSTANTS[7] = CONSTANTS[0] - 12.0000;
  CONSTANTS[8] = CONSTANTS[0]+10.6130;

  algebraics[1] = ( - 0.100000*(states[0]+50.0000))/(exp(- (states[0]+50.0000)/10.0000) - 1.00000);
  algebraics[5] =  4.00000*exp(- (states[0]+75.0000)/18.0000);
  rates[1] =  algebraics[1]*(1.00000 - states[1]) -  algebraics[5]*states[1];
  algebraics[2] =  0.0700000*exp(- (states[0]+75.0000)/20.0000);
  algebraics[6] = 1.00000/(exp(- (states[0]+45.0000)/10.0000)+1.00000);
  rates[2] =  algebraics[2]*(1.00000 - states[2]) -  algebraics[6]*states[2];
  algebraics[3] = ( - 0.0100000*(states[0]+65.0000))/(exp(- (states[0]+65.0000)/10.0000) - 1.00000);
  algebraics[7] =  0.125000*exp((states[0]+75.0000)/80.0000);
  rates[3] =  algebraics[3]*(1.00000 - states[3]) -  algebraics[7]*states[3];
  algebraics[0] =  CONSTANTS[3]*pow(states[1], 3.00000)*states[2]*(states[0] - CONSTANTS[6]);
  algebraics[4] =  CONSTANTS[4]*pow(states[3], 4.00000)*(states[0] - CONSTANTS[7]);
  algebraics[8] =  CONSTANTS[5]*(states[0] - CONSTANTS[8]);
  rates[0] = - (- parameters[0]+algebraics[0]+algebraics[4]+algebraics[8])/CONSTANTS[1];
}
