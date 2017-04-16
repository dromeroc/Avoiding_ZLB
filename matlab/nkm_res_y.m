function [res,C,N,Y,PI,R,WP,RMC,cres] = nkm_res_y(guess,fspace,Grid,B1,QW,PHI_FG,zlb)
% This function: compute the residuals of the NKM using collocation method
%
% Inputs:
%     - guess:    2*nA*nBx1 vector of guesses
%     - fspace:   defines interpolation scheme
%     - Grid:     nA*nBx1 vector of TFP and discount factor (in logs)
%     - Yn:       nA*nBx1 vector of natural output (pre-allocated)
%     - A1:       nA*nB*nQx1 vector of TFP in period t+1 (in logs)
%     - B1:       nA*nB*nQx1 vector of discount factor in period t+1 (in logs)
%     - QW:       nQx1 vector of weights for expectations
%     - zlb:      1: solve the model with the ZLB constraint
%
% Outputs:
%     - res:      2*nA*nBx1 vector of residuals
%     - C:        nA*nBx1 vector of consumption
%     - N:        nA*nBx1 vector of labor
%     - Y:        nA*nBx1 vector of output
%     - PI:       nA*nBx1 vector of inflation
%     - R:        nA*nBx1 vector of nominal interest rate
%     - WP:       nA*nBx1 vector of real wage
%     - RMC:      nA*nBx1 vector of real marginal cost
%     - Yg:       nA*nBx1 vector of output gap
%     - cres:     2*nA*nBx1 vector of consumption residuals (from Euler)

% Read parameters
[SIGG,PHI,EPSILON,~,~,ZETA,PHI_PI,PHI_Y,~,~,PI_LR,R_LR,Y_LR] = parameters;

% State dimensions
bN  = fspace.n(1);
fgN = fspace.n(2);

% States in levels
FG = exp(Grid(:,2));

% Compute guess
guess   = reshape(guess,bN*fgN,2);
C       = funeval(guess(:,1),fspace,Grid); % Consumption
PI      = funeval(guess(:,2),fspace,Grid); % Inflation

% Compute other variables in period t
PSI     = 1./(1-(ZETA/2)*(PI-1).^2);
Y       = PSI.*C;
N       = Y;
WP      = (C.^SIGG).*(N.^PHI);
RMC     = WP;
R       = R_LR*((PI/PI_LR).^PHI_PI).*((Y/Y_LR).^PHI_Y).*((FG/Y_LR).^PHI_FG);

if zlb==1
    R   = max(R,1);
end

% Compute next period guess
C1      = funeval(guess(:,1),fspace,[B1 kron(log(Y),ones(length(QW),1))]);
PI1     = funeval(guess(:,2),fspace,[B1 kron(log(Y),ones(length(QW),1))]);

% Compute other variables in period t+1
PSI1    = 1./(1-(ZETA/2)*(PI1-1).^2);
Y1      = PSI1.*C1;

% Compute expected values
% When length(Grid)==order of approximation (standard case), aN_aux=aN.
% This version is more general because allows to compute the solution for a
% larger grid in order to evaluate Euler residuals
N_aux   = length(Grid);
P_aux   = kron(ones(N_aux,1),QW);
EC      = sum(reshape(exp(B1).*(C1.^-SIGG).*(1./PI1).*P_aux,length(QW),N_aux))';
EP      = sum(reshape(exp(B1).*(C1.^-SIGG).*(PI1-1).*PI1.*Y1.*P_aux,length(QW),N_aux))';

% Compute residuals
res1    = 1-R.*(C.^SIGG).*EC;
res2    = Y.*(1-EPSILON+EPSILON*RMC-ZETA*(PI-1).*PI)+ZETA*(C.^SIGG).*EP;
res     = [res1;res2];

% Compute consumption residual (to verify accuracy)
cres    = 1-((R.*EC).^(-1/SIGG))./C;