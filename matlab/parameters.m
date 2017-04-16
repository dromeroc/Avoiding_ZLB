function [SIGG,PHI,EPSILON,THETA,BETA_LR,ZETA,PHI_PI,PHI_Y,RHO_B,SIGG_B,...
    PI_LR,R_LR,Y_LR] = parameters
% This function: set the parameterization of the NKM assuming a zero
% inflation rate in steady state

% Elasticity of consumption
SIGG    = 2;

% Elasticity of labor
PHI     = 1;

% Elasticity of varieties
EPSILON = 6;

% Frequency of price adjustment on Calvo
THETA   = 0.75;

% Discount factor in the long-run
BETA_LR = 0.994;

% Adjustment price parameter on Rotemberg
ZETA    = (EPSILON-1)*THETA/((1-THETA)*(1-BETA_LR*THETA));

% Taylor rule--reaction to inflation
PHI_PI  = 1.5;

% Taylor rule--reaction to output
PHI_Y   = 0.25;

% Taylor rule--reaction to past output
%PHI_FG  = 2;

% Persistence of preference shock
RHO_B   = 0.8;

% Volatility of preference shock
SIGG_B  = 0.0025;

% Steady state value of inflation
PI_LR   = 1.005;

% Steady state value of interest rate
R_LR    = PI_LR/BETA_LR;

% Steady state value of output
PSI     = 1/(1-0.5*ZETA*((PI_LR-1)^2));
WP      = ((1-BETA_LR)*ZETA*(PI_LR-1)*PI_LR-(1-EPSILON))/EPSILON;
C       = (WP/(PSI^PHI))^(1/(SIGG+PHI));
Y_LR    = PSI*C;
