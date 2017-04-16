function res = nk_model_p(x,pars)
% This function compute the residuals of the system of equations that
% characterize the model.
%
% Input
%
%   -x:         (2*nvar+nshocks+nEerror)x1 vector of variables
%   -params:    11x1 vector of parameters
%
% Output
%
%   -res:       nvarx1 vector of residual equations

% Read parameters
SIGG    = pars(1);
PHI     = pars(2);
EPSILON = pars(3);
THETA   = pars(4);
BETA_LR = pars(5);
ZETA    = pars(6);
PHI_PI  = pars(7);
PHI_Y   = pars(8);
RHO_B   = pars(9);
SIGG_B  = pars(10);
PI_LR   = pars(11);
R_LR    = pars(12);
Y_LR    = pars(13);
PHI_FG  = pars(14);

% Decomposing contemporaneous variables
b       = x(1);
c       = x(2);
n       = x(3);
y       = x(4);
Pi      = x(5);
r       = x(6);
wp      = x(7);
rmc     = x(8);

% Decomposing lagged variables
b_lag   = x(9);
c_lag   = x(10);
n_lag   = x(11);
y_lag   = x(12);
Pi_lag  = x(13);
r_lag   = x(14);
wp_lag  = x(15);
rmc_lag = x(16);

% Decomposing exogenous shocks
eps_b   = x(17);

% Decomposing expectational errors
ERR_C   = x(18); % Error in the Euler equation
ERR_P   = x(19); % Error in the Euler of pricing

% Model
res = zeros(8,1);

% Law of movement of preference shock
res(1) = -b+log(BETA_LR)*(1-RHO_B)+RHO_B*b_lag+eps_b;

% Intratemporal equation
res(2) = -(c^SIGG)*(n^PHI)+wp;

% Euler equation for consumption
res(3) = -1/r_lag+exp(b)*((c_lag/c)^SIGG)*(1/Pi)+ERR_C;

% Euler equation for pricing
res(4) = -((1-EPSILON)+EPSILON*rmc_lag-ZETA*(Pi_lag-1)*Pi_lag)-...
    exp(b)*((c_lag/c)^SIGG)*(ZETA*(Pi-1)*Pi*y/y_lag)+ERR_P;

% Real marginal cost
res(5) = -rmc+wp;

% Taylor rule
res(6) = -(r/R_LR)+((Pi/PI_LR)^PHI_PI)*((y/Y_LR)^PHI_Y)*((Pi_lag/PI_LR)^PHI_FG);

% Production function
res(7) = -y+n;

% Aggregation
res(8) = -y+(1/(1-0.5*ZETA*(Pi-1)^2))*c;
