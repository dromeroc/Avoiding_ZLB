%% Expectations at the Zero Lower Bound
%   Damian Romero
%   December 2016
%
% This version: Solves the New Keynesian Model with preferences and forward
% guidance (past inflation)

% Housekeeping
clear
close all
clc

% Forward guidance parameters
FGpar = [0 0.1 0.5 0.7 1 1.2 1.499 1.7 2];
for j=1:length(FGpar)
    % Parameters
    PHI_FG = FGpar(j);
    [SIGG,PHI,EPSILON,THETA,BETA_LR,ZETA,PHI_PI,PHI_Y,RHO_B,...
        SIGG_B,PI_LR,R_LR,Y_LR,PHI_FG] = parameters2(PHI_FG);
    pars = [SIGG,PHI,EPSILON,THETA,BETA_LR,ZETA,PHI_PI,PHI_Y,RHO_B,SIGG_B,...
        PI_LR,R_LR,Y_LR,PHI_FG];
    
    % Compute steady state
    [WP,C,PI,RMC,Y,N,R,BETA] = steady_state;
    stst = [BETA,C,N,Y,PI,R,WP,RMC];
    
    % Amplied vector of steady state values (to include lagged and current
    % state values, fundamental shock and expectational error shocks)
    stst2 = [stst stst 0 0 0];
    
    % Solving the linearized model
    nvar    = 8;    % # of variables
    nshocks = 1;    % # of shocks
    nEerr   = 2;    % # of expectational errors
    [Hy,Hz,eu] = mysolvemodel(@nk_model_p,stst2,pars,nvar,nshocks,nEerr);
    
    % Space for state variable
    bN      = 20;
    Bs      = sqrt(SIGG_B^2/(1-RHO_B^2));
    Bprec   = 5;
    Bmin    = -Bprec*Bs+log(BETA_LR);
    Bmax    = Bprec*Bs+log(BETA_LR);
    
    fgN     = 20;
    pct     = 0.1;
    fgmin   = log(stst(5)*(1-pct)); % inflation stst-pct
    fgmax   = log(stst(5)*(1+pct)); % inflation stst+pct
    
    % Preferences for expectations
    % Our task is to characterize carefully expectations. Therefore, we
    % requiere a lot of points to compute them.
    qBN     = 30;
    [QP,QW] = qnwnorm(qBN,0,SIGG_B^2);
    
    % Define settings for the nonlinear solution
    LB      = [Bmin fgmin];
    UB      = [Bmax fgmax];
    Order   = [bN fgN];
    aptype  = 'spli';
    spliord = 1;
    fspace  = fundefn(aptype,Order,LB,UB,spliord);
    nodes   = funnode(fspace);
    nodes{1}(lookup(nodes{1},stst(1),3))= stst(1);          % Put ss of B in nodes
    nodes{2}(lookup(nodes{2},log(stst(5)),3))= log(stst(5));% Put ss of Pi in nodes
    Grid    = gridmake(nodes);
    Basis   = funbas(fspace);
    
    % Define guess (we want to approximate consumption and inflation, which are
    % variables 2 and 5 in the linear system, respectively. Beta is the first
    % variable in the system and inflation is the fifth variable.)
    c_guess     = stst(2)+Hy(2,1)*(Grid(:,1)-stst(1))+Hy(2,5)*(Grid(:,2)-log(stst(5)));
    pi_guess    = stst(5)+Hy(5,1)*(Grid(:,1)-stst(1))+Hy(5,5)*(Grid(:,2)-log(stst(5)));
    guess       = [Basis\c_guess;Basis\pi_guess];
    
    % Optimization options
    optset('broyden','maxit',100)
    optset('broyden','tol',1e-8)
    optset('broyden','maxsteps',50)
    optset('broyden','showiters',1)
    
    % Pre-allocate variables
    [B1,bound2] = preallocated(Grid,fspace,bN,fgN,QP);
    
    % Define a larger grid to evaluate residuals
    bN_large    = 100;
    fgN_large   = 100;
    Grid_large  = gridmake(linspace(Bmin,Bmax,bN_large)',linspace(fgmin,fgmax,fgN_large)');
    B1_large    = preallocated(Grid_large,fspace,bN_large,fgN_large,QP);
    
    % Solve the model and get new parameters
    
    % Solve the model without the ZLB
    zlb = 0;
    par1 = broyden('nkm_res_p',guess,fspace,Grid,B1,QW,PHI_FG,zlb);
    
    % Solve the model with the ZLB (use as initial guess previous solution)
    zlb = 1;
    par2 = broyden('nkm_res_p',par1,fspace,Grid,B1,QW,PHI_FG,zlb);
    
    % Recover policy functions without the ZLB
    zlb = 0;
    [~,C1,N1,Y1,PI1,R1,WP1,RMC1] = nkm_res_p(par1,fspace,Grid,B1,QW,PHI_FG,zlb);
    
    % Check Euler residuals without the ZLB
    [EulerRes1,~,~,~,~,~,~,~,CRes1] = nkm_res_p(par1,fspace,Grid_large,...
        B1_large,QW,PHI_FG,zlb);
    
    % Recover policy functions with the ZLB
    zlb = 1;
    [~,C2,N2,Y2,PI2,R2,WP2,RMC2] = nkm_res_p(par2,fspace,Grid,B1,QW,PHI_FG,zlb);
    
    % Check Euler residuals with the ZLB
    [EulerRes2,~,~,~,~,R2_large,~,~,CRes2] = nkm_res_p(par2,fspace,Grid_large,...
        B1_large,QW,PHI_FG,zlb);
    
    % Table with analysis of residuals
    Resids      = zeros(2,4);
    Resids(1,:) = [mean(EulerRes1) max(abs(EulerRes1)) mean(CRes1) max(abs(CRes1))];
    Resids(2,:) = [mean(EulerRes2) max(abs(EulerRes2)) mean(CRes2) max(abs(CRes2))];
    %disp(Resids)
    
    % Impulse-response
    
    T = 20;
    B_proc_irf = zeros(T+1,1);
    Bshock = zeros(T+1,1); Bshock(2) = RHO_B*log(BETA_LR)/SIGG_B+Bmax/SIGG_B;
    for i=2:T+1
        B_proc_irf(i) = log(BETA_LR)*(1-RHO_B)+RHO_B*B_proc_irf(i-1)+Bshock(i)*SIGG_B;
        B_proc_irf(i) = max(min(B_proc_irf(i),fspace.b(1)),fspace.a(1));
    end
    
    % Set some properties for plots
%     set(groot, 'defaultAxesTickLabelInterpreter','latex')
%     set(groot, 'defaultLegendInterpreter','latex')
    
    % Compute IRFs without ZLB
    zlb         = 0;            % ZLB or not
    hor         = 8;            % Horizon for expectations
    irf_h_sel   = [4 hor];      % Selection matrix to show expectations
    [C,~,Y,PI,R,~,RMC,Y_exp,R_exp,PI_exp,C_exp,RMC_exp,U] = ...
        nkm_simul_p(fspace,B_proc_irf,par1,QP,QW,PHI_FG,zlb,stst,hor);
    irf.noZ = struct('C',C,'Y',Y,'PI',PI,'R',R,'RMC',RMC,'Y_exp',Y_exp,...
        'R_exp',R_exp,'PI_exp',PI_exp,'C_exp',C_exp,'RMC_exp',RMC_exp,'U',U);
    
    % Compute IRFs with ZLB
    zlb = 1;
    [Cz,~,Yz,PIz,Rz,~,RMCz,Yz_exp,Rz_exp,PIz_exp,Cz_exp,RMCz_exp,Uz] = ...
        nkm_simul_p(fspace,B_proc_irf,par2,QP,QW,PHI_FG,zlb,stst,hor);
    irf.Z = struct('C',Cz,'Y',Yz,'PI',PIz,'R',Rz,'RMC',RMCz,'Y_exp',Yz_exp,...
        'R_exp',Rz_exp,'PI_exp',PIz_exp,'C_exp',Cz_exp,'RMC_exp',RMCz_exp,'U',Uz);
    
    % Simulation and distribution
    
    rng('default')
    T           = 50000;
    B_proc_sim  = zeros(T+1,1); Bshock = randn(T+1,1); Btrunc = zeros(T+1,1);
    for i=2:T+1
        B_proc_sim(i)   = log(BETA_LR)*(1-RHO_B)+RHO_B*B_proc_sim(i-1)+Bshock(i)*SIGG_B;
        Btrunc(i)       = (B_proc_sim(i)>fspace.b(1)) + (B_proc_sim(i)<fspace.a(1));
        B_proc_sim(i)   = max(min(B_proc_sim(i),fspace.b(1)),fspace.a(1));
    end
    %B_proc_sim = B_proc_sim(2:end);
    if mean(Btrunc)>0.05;warning('Too much truncation');end
    
    % Compute simulations without ZLB
    zlb = 0;
    [Cs,~,Ys,PIs,Rs,~,RMCs,Ys_exp,Rs_exp,PIs_exp,Cs_exp,RMCs_exp,Us] = ...
        nkm_simul_p(fspace,B_proc_sim,par1,QP,QW,PHI_FG,zlb,stst,hor);
    sim.noZ = struct('C',Cs,'Y',Ys,'PI',PIs,'R',Rs,'RMC',RMCs,'Y_exp',Ys_exp,...
        'R_exp',Rs_exp,'PI_exp',PIs_exp,'C_exp',Cs_exp,'RMC_exp',RMCs_exp,'U',Us);
    
    % Compute simulations with ZLB
    zlb = 1;
    [Csz,~,Ysz,PIsz,Rsz,~,RMCsz,Ysz_exp,Rsz_exp,PIsz_exp,Csz_exp,RMCsz_exp,Usz] = ...
        nkm_simul_p(fspace,B_proc_sim,par2,QP,QW,PHI_FG,zlb,stst,hor);
    sim.Z = struct('C',Csz,'Y',Ysz,'PI',PIsz,'R',Rsz,'RMC',RMCsz,'Y_exp',Ysz_exp,...
        'R_exp',Rsz_exp,'PI_exp',PIsz_exp,'C_exp',Csz_exp,'RMC_exp',RMCsz_exp,'U',Usz);
    
    spentZLB = ...
        [mean((PIs_exp(:,irf_h_sel)<stst(5))) mean((PIsz_exp(:,irf_h_sel)<stst(5)));
        mean((Ys_exp(:,irf_h_sel)<stst(4))) mean((Ysz_exp(:,irf_h_sel)<stst(4)));
        mean((RMCs_exp(:,irf_h_sel)<stst(8))) mean((RMCsz_exp(:,irf_h_sel)<stst(8)));
        mean((Rs_exp(:,irf_h_sel)<stst(6))) mean((Rsz_exp(:,irf_h_sel)<stst(6)))];
    spentZLB = spentZLB(:,[1 3 2 4]);
    
%     figure(1)
%     subplot(2,2,1)
%     plot(1:T,PI,1:T,PIz,'--','LineWidth',2),hold on,plot(1:T,stst(5)*ones(T,1),'k')
%     title('Inflation','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,2)
%     plot(1:T,Y,1:T,Yz,'--','LineWidth',2),hold on,plot(1:T,stst(4)*ones(T,1),'k')
%     title('Output','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,3)
%     plot(1:T,RMC,1:T,RMCz,'--','LineWidth',2),hold on,plot(1:T,stst(8)*ones(T,1),'k')
%     title('Real marginal cost','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,4)
%     plot(1:T,R,1:T,Rz,'--','LineWidth',2),hold on,plot(1:T,stst(6)*ones(T,1),'k')
%     title('Interest rate','Interpreter','latex')
%     legend('No ZLB','ZLB','Location','SouthEast')
%     set(gca,'FontSize',13)
%     
%     figure(2)
%     subplot(2,2,1)
%     plot(1:T,PI_exp(:,irf_h_sel),1:T,PIz_exp(:,irf_h_sel),'--','LineWidth',2),
%     hold on,plot(1:T,stst(5)*ones(T,1),'--k')
%     title('Expected inflation','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,2)
%     plot(1:T,Y_exp(:,irf_h_sel),1:T,Yz_exp(:,irf_h_sel),'--','LineWidth',2),
%     hold on,plot(1:T,stst(4)*ones(T,1),'k')
%     title('Expected output','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,3)
%     plot(1:T,RMC_exp(:,irf_h_sel),1:T,RMCz_exp(:,irf_h_sel),'--','LineWidth',2),
%     hold on,plot(1:T,stst(8)*ones(T,1),'k')
%     title('Expected real marginal cost','Interpreter','latex')
%     set(gca,'FontSize',13)
%     
%     subplot(2,2,4)
%     plot(1:T,R_exp(:,irf_h_sel),1:T,Rz_exp(:,irf_h_sel),'--','LineWidth',2),
%     hold on,plot(1:T,stst(6)*ones(T,1),'k')
%     title('Expected interest rate','Interpreter','latex')
%     hL=legend('No ZLB, $h=4$','No ZLB, $h=8$','ZLB, $h=4$','ZLB, $h=8$');
%     newPosition = [0 -0.02 1 0.1];
%     newUnits = 'normalized';
%     set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
%         ,'FontSize',9);
%     set(gca,'FontSize',13)

    % Name of outputs
    out_name                = strcat('PI_',num2str(j));
    results(1).(out_name)   = par1;
    results(2).(out_name)   = par2;
    results(3).(out_name)   = Resids;
    results(4).(out_name)   = irf;
    results(5).(out_name)   = PHI_FG;
    results(6).(out_name)   = pars;
    results(7).(out_name)   = B_proc_irf;
    results(8).(out_name)   = sim;
    results(9).(out_name)   = spentZLB;
    disp(j)
    
end

% %% Simulation and distribution
% 
% rng('default')
% T           = 1000;
% B_proc_sim  = zeros(T+1,1); Bshock = randn(T+1,1); Btrunc = zeros(T+1,1);
% for i=2:T+1
%     B_proc_sim(i)   = log(BETA_LR)*(1-RHO_B)+RHO_B*B_proc_sim(i-1)+Bshock(i)*SIGG_B;
%     Btrunc(i)       = (B_proc_sim(i)>fspace.b) + (B_proc_sim(i)<fspace.a);
% 	B_proc_sim(i)   = max(min(B_proc_sim(i),fspace.b),fspace.a);
% end
% B_proc_sim = B_proc_sim(2:end);
% if mean(Btrunc)>0.05;warning('Too much truncation');end
% 
% % Compute simulations without ZLB
% zlb = 0;
% [Cs,~,Ys,PIs,Rs,~,RMCs,Ys_exp,Rs_exp,PIs_exp,Cs_exp,RMCs_exp] = ...
%     nkm_simul(fspace,B_proc_sim,par1,QP,QW,zlb,hor);
% 
% % Compute simulations with ZLB
% zlb = 1;
% [Csz,~,Ysz,PIsz,Rsz,~,RMCsz,Ysz_exp,Rsz_exp,PIsz_exp,Csz_exp,RMCsz_exp] = ...
%     nkm_simul(fspace,B_proc_sim,par2,QP,QW,zlb,hor);
% 
% % Produce figures
% simul_figures;
% 
% % Produce table with time spent below the steady-state
% %[stst(5) stst(4) stst(8) stst(6)]; % inflation, output, RMC, interest rate
% spentZLB = ...
% [mean((PIs_exp(:,irf_h_sel)<stst(5))) mean((PIsz_exp(:,irf_h_sel)<stst(5)));
%     mean((Ys_exp(:,irf_h_sel)<stst(4))) mean((Ysz_exp(:,irf_h_sel)<stst(4)));
%     mean((RMCs_exp(:,irf_h_sel)<stst(8))) mean((RMCsz_exp(:,irf_h_sel)<stst(8)));
%     mean((Rs_exp(:,irf_h_sel)<stst(6))) mean((Rsz_exp(:,irf_h_sel)<stst(6)))];
% spentZLB = spentZLB(:,[1 3 2 4]);
% disp(spentZLB)