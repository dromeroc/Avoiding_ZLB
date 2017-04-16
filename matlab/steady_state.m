function [WP,C,PI,RMC,Y,N,R,BETA] = steady_state
% This function: computes the steady state of the model

% Call parameters
[SIGG,PHI,EPSILON,~,BETA_LR,ZETA,~,~,~,~,PI_LR,R_LR,Y_LR] = parameters;

% Preferences
BETA = log(BETA_LR);

% Inflation
PI  = PI_LR;

% Interest rate
R   = PI/BETA_LR;

% Real marginal cost
RMC = ((1-BETA_LR)*ZETA*(PI-1)*PI-(1-EPSILON))/EPSILON;

% Marginal cost
WP  = RMC;

% Auxiliar variable (connecting output with consumption in the aggregate)
PSI = 1/(1-0.5*ZETA*((PI-1)^2));

% Consumption
C   = (WP/(PSI^PHI))^(1/(SIGG+PHI));

% Output
Y   = PSI*C;

% Labor
N   = Y;

% Sanity check
if any([R-R_LR,Y-Y_LR]~=0);
    warning('Caution with parameters and steady state')
end