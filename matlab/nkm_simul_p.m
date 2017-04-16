function [C,N,Y,PI,R,WP,RMC,Y_exp,R_exp,PI_exp,C_exp,RMC_exp,U] = ...
    nkm_simul_p(fspace,B_proc,par,QP,QW,PHI_FG,zlb,stst,hor)
% This function: simulates the model given a proces for TFP and discount 
% factor (externaly given) and the solution of the model. 
% Valid for simulation and IRF
%
% Inputs:
%     - fspace:   defines interpolation scheme
%     - A_proc:   simulated stochastic process for log TFP
%     - B_proc:   simulated stochastic process for log discount factor
%     - par:      parameters that solve the NKM approximation
%     - QP:       vector of possible shocks in t+1
%     - QW:       vector of weights for expectations
%     - zlb:      1 if simulates subject to the ZLB
%     - hor:      horizon for expectations
%         
% Outputs:
%     - C:        vector of simulated consumption
%     - N:        vector of simulated labor
%     - Y:        vector of simulated output
%     - PI:       vector of simulated inflation
%     - R:        vector of simulated interest rate
%     - WP:       vector of simulated real wage
%     - RMC:      vector of simulated real marginal cost
%     - Yn:       vector of simulated natural output
%     - Yg:       vector of simulated output gap
%     - Y_exp:    vector of simulated expected output
%     - R_exp:    vector of simulated expected interest rate
%     - PI_exp:   vector of simulated expected inflation
%     - C_exp:    vector of simulated expected consumption
%     - RMC_exp:  vector of simulated expected consumption
%     - Yg_exp:   vector of simulated expected consumption

% Length of simulation
T  = length(B_proc);

% State dimensions
bN  = fspace.n(1);
fgN = fspace.n(2);

% Call parameters that solve the model
par = reshape(par,bN*fgN,2);

% Call parameters
[SIGG,PHI,~,~,BETA_LR,ZETA,PHI_PI,PHI_Y,RHO_B,~,PI_LR,R_LR,Y_LR] = parameters;

% Pre-allocate matrices to save results (because we want to compute also
% expectations, we need to do it recursively)
C       = zeros(T,1);
N       = zeros(T,1);
Y       = zeros(T,1);
PI      = zeros(T,1);
R       = zeros(T,1);
WP      = zeros(T,1);
RMC     = zeros(T,1);
U       = zeros(T,1);
FG      = zeros(T,1);   FG(1) = stst(5); % Inflation
Y_exp   = zeros(T,hor);
R_exp   = zeros(T,hor);
PI_exp  = zeros(T,hor);
C_exp   = zeros(T,hor);
RMC_exp = zeros(T,hor);
for t=2:T
    
    % Compute variables in t
    C(t)    = funeval(par(:,1),fspace,[B_proc(t) log(FG(t-1))]);
    PI(t)   = funeval(par(:,2),fspace,[B_proc(t) log(FG(t-1))]);
    PSI     = 1./(1-(ZETA/2)*(PI(t)-1).^2);
    Y(t)    = PSI.*C(t);
    N(t)    = Y(t);
    WP(t)   = (C(t).^SIGG).*(N(t).^PHI);
    RMC(t)  = WP(t);
    R(t)    = R_LR*((PI(t)/PI_LR).^PHI_PI).*((Y(t)/Y_LR).^PHI_Y).*((FG(t-1)/PI_LR).^PHI_FG);
    U(t)    = (C(t)^(1-SIGG))/(1-SIGG)-(N(t)^(1+PHI))/(1+PHI);
    if zlb;R(t)=max(R(t),1);end
    FG(t)   = PI(t);
    
    % Compute expectations h periods ahead
    B_aux = B_proc(t);  % Auxiliar variable for discount factor tomorrow
    F_aux = FG(t);      % Auxiliar variable for forward guidance tomorrow
    for j=1:hor
        % Compute possible values for preferences tomorrow
        B_next = kron(log(BETA_LR)*(1-RHO_B)+RHO_B*B_aux,ones(length(QP),1))+QP;
        B_next = min(max(B_next,fspace.a(1)),fspace.b(1));

        % Compute possible values for variables tomorrow
        C1    = funeval(par(:,1),fspace,[B_next log(F_aux)*ones(length(QP),1)]);
        PI1   = funeval(par(:,2),fspace,[B_next log(F_aux)*ones(length(QP),1)]);
        PSI1  = 1./(1-(ZETA/2)*(PI1-1).^2);
        Y1    = PSI1.*C1;
        N1    = Y1;
        WP1   = (C1.^SIGG).*(N1.^PHI);
        RMC1  = WP1;
        R1    = R_LR*((PI1/PI_LR).^PHI_PI).*((Y1/Y_LR).^PHI_Y).*...
            ((F_aux*ones(length(QP),1)/PI_LR).^PHI_FG);
        if zlb;R1=max(R1,1);end

        % Compute expected values
        Y_exp(t,j)    = QW'*Y1;
        R_exp(t,j)    = QW'*R1; if zlb;R_exp(t,j) = max(R_exp(t,j),1);end
        PI_exp(t,j)   = QW'*PI1;
        C_exp(t,j)    = QW'*C1;
        RMC_exp(t,j)  = QW'*RMC1;
        
        % Update B_aux with their expected values
        B_aux = (1-RHO_B)*log(BETA_LR)+RHO_B*B_aux;
        F_aux = PI_exp(t,j);
    end
end

% Adjust output
C       = C(2:end);
N       = N(2:end);
Y       = Y(2:end);
PI      = PI(2:end);
R       = R(2:end);
WP      = WP(2:end);
RMC     = RMC(2:end);
U       = U(2:end);
Y_exp	= Y_exp(2:end,:);
R_exp	= R_exp(2:end,:);
PI_exp	= PI_exp(2:end,:);
C_exp	= C_exp(2:end,:);
RMC_exp = RMC_exp(2:end,:);