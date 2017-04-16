%% Expectations at the Zero Lower Bound
%   Damian Romero
%   December 2016
%
% This version: Analyze the solutions for forward guidance cases

% Housekeeping
clear
close all
clc

% Load solutions
res_y = load('results_spli1_y'); res_y = res_y.results;
res_p = load('results_spli1_p'); res_p = res_p.results;

% Compute steady state
[WP,C,PI,RMC,Y,N,R,BETA] = steady_state;
stst = [BETA,C,N,Y,PI,R,WP,RMC];

% Set some properties for plots
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
irf_h_sel   = [4 8];
T           = 20;

%% 1) Comparison between the model including the ZLB and not including it 
%       and without forward guidance
% Here we present the comparison in the responses of the model with and
% without the ZLB when the economy is subject to a large negative demand
% shock that induces the ZLB to be binding. Compare what happens in the
% economy by considering or not the restriction.

figure(1) % Level
subplot(2,2,1)
plot(1:T,100*(res_y(4).Yg_1.noZ.PI/stst(5)-1),...
    1:T,100*(res_y(4).Yg_1.Z.PI/stst(5)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(res_y(4).Yg_1.noZ.Y/stst(4)-1),...
    1:T,100*(res_y(4).Yg_1.Z.Y/stst(4)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(res_y(4).Yg_1.noZ.RMC/stst(8)-1),...
    1:T,100*(res_y(4).Yg_1.Z.RMC/stst(8)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,res_y(4).Yg_1.noZ.R,...
    1:T,res_y(4).Yg_1.Z.R,...
    '--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Interest rate','Interpreter','latex')
% legend('No ZLB','ZLB','Location','SouthEast')
% set(gca,'FontSize',13)
hL=legend('No ZLB','ZLB');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',12);
set(gca,'FontSize',13)

figure(2) % Expectations
subplot(2,2,1)
plot(1:T,100*(res_y(4).Yg_1.noZ.PI_exp(:,irf_h_sel)/stst(5)-1),...
    1:T,100*(res_y(4).Yg_1.Z.PI_exp(:,irf_h_sel)/stst(5)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(res_y(4).Yg_1.noZ.Y_exp(:,irf_h_sel)/stst(4)-1),...
    1:T,100*(res_y(4).Yg_1.Z.Y_exp(:,irf_h_sel)/stst(4)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(res_y(4).Yg_1.noZ.RMC_exp(:,irf_h_sel)/stst(8)-1),...
    1:T,100*(res_y(4).Yg_1.Z.RMC_exp(:,irf_h_sel)/stst(8)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,res_y(4).Yg_1.noZ.R_exp(:,irf_h_sel),...
    1:T,res_y(4).Yg_1.Z.R_exp(:,irf_h_sel),...
    '--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Expected interest rate','Interpreter','latex')
hL=legend('No ZLB, $h=4$','No ZLB, $h=8$','ZLB, $h=4$','ZLB, $h=8$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

%% 2) Comparison between the model including the ZLB and with forward guidance
% Here we present the comparison between the model subject to the ZLB but
% with no forward guidance policy and the case where forward guidance
% policy is active. To make comparison fair and illustrative enough, the
% former case present the responses of relevant variables when the forward
% guidance channel exists but is not relevant (i.e. the parameter that
% governs the forward guidance is zero), while the latter has a value for
% the parameter equal to 0.5 to make it strong enough to make the point.

figure(3) % Level
subplot(2,2,1)
plot(1:T,100*(res_y(4).Yg_1.Z.PI/stst(5)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.PI/stst(5)-1),...
    1:T,100*(res_p(4).PI_5.Z.PI/stst(5)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(res_y(4).Yg_1.Z.Y/stst(4)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.Y/stst(4)-1),...
    1:T,100*(res_p(4).PI_5.Z.Y/stst(4)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(res_y(4).Yg_1.Z.RMC/stst(8)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.RMC/stst(8)-1),...
    1:T,100*(res_p(4).PI_5.Z.RMC/stst(8)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,res_y(4).Yg_1.Z.R,':o',...
    1:T,res_y(4).Yg_5.Z.R,...
    1:T,res_p(4).PI_5.Z.R,...
    '--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Interest rate','Interpreter','latex')
%legend('No ZLB','ZLB','Location','SouthEast')
set(gca,'FontSize',13)
hL=legend('ZLB no FG','ZLB FG $\pi$','ZLB FG $y$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);

figure(4) % Expectations, h=4
subplot(2,2,1)
plot(1:T,100*(res_y(4).Yg_1.Z.PI_exp(:,4)/stst(5)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.PI_exp(:,4)/stst(5)-1),...
    1:T,100*(res_p(4).PI_5.Z.PI_exp(:,4)/stst(5)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(res_y(4).Yg_1.Z.Y_exp(:,4)/stst(4)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.Y_exp(:,4)/stst(4)-1),...
    1:T,100*(res_p(4).PI_5.Z.Y_exp(:,4)/stst(4)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(res_y(4).Yg_1.Z.RMC_exp(:,4)/stst(8)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.RMC_exp(:,4)/stst(8)-1),...
    1:T,100*(res_p(4).PI_5.Z.RMC_exp(:,4)/stst(8)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,res_y(4).Yg_1.Z.R_exp(:,4),':o',...
    1:T,res_y(4).Yg_5.Z.R_exp(:,4),...
    1:T,res_p(4).PI_5.Z.R_exp(:,4),...
    '--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Expected interest rate','Interpreter','latex')
hL=legend('ZLB no FG','ZLB FG $\pi$','ZLB FG $y$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

figure(5) % Expectations, h=8
subplot(2,2,1)
plot(1:T,100*(res_y(4).Yg_1.Z.PI_exp(:,8)/stst(5)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.PI_exp(:,8)/stst(5)-1),...
    1:T,100*(res_p(4).PI_5.Z.PI_exp(:,8)/stst(5)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(res_y(4).Yg_1.Z.Y_exp(:,8)/stst(4)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.Y_exp(:,8)/stst(4)-1),...
    1:T,100*(res_p(4).PI_5.Z.Y_exp(:,8)/stst(4)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(res_y(4).Yg_1.Z.RMC_exp(:,8)/stst(8)-1),':o',...
    1:T,100*(res_y(4).Yg_5.Z.RMC_exp(:,8)/stst(8)-1),...
    1:T,100*(res_p(4).PI_5.Z.RMC_exp(:,8)/stst(8)-1),...
    '--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,res_y(4).Yg_1.Z.R_exp(:,8),':o',...
    1:T,res_y(4).Yg_5.Z.R_exp(:,8),...
    1:T,res_p(4).PI_5.Z.R_exp(:,8),...
    '--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Expected interest rate','Interpreter','latex')
hL=legend('ZLB no FG','ZLB FG $\pi$','ZLB FG $y$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

%% 3) Comparison between the model including the ZLB and with forward
%       guidance checking different levels and types
% Here we rescue responses of the economy to the shocks under different
% configurations of the forward guidance policy. In particular, we extract
% responses for the different levels and generate plots that shows the
% percentual differences with respect to steady state of each variable

% Summarize information
R_y     = zeros(20,9); Re_y     = zeros(20,8,9);
Y_y     = zeros(20,9); Ye_y     = zeros(20,8,9);
PI_y    = zeros(20,9); PIe_y    = zeros(20,8,9);
RMC_y   = zeros(20,9); RMCe_y   = zeros(20,8,9);

for i=1:9
    % Name
    out_name = strcat('Yg_',num2str(i));
    
    % FG with output in levels
    R_y(:,i)    = res_y(4).(out_name).Z.R;
    PI_y(:,i)   = res_y(4).(out_name).Z.PI;
    Y_y(:,i)    = res_y(4).(out_name).Z.Y;
    RMC_y(:,i)  = res_y(4).(out_name).Z.RMC;
    
    % FG with output in expectations
    Re_y(:,:,i)   = res_y(4).(out_name).Z.R_exp;
    PIe_y(:,:,i)  = res_y(4).(out_name).Z.PI_exp;
    Ye_y(:,:,i)   = res_y(4).(out_name).Z.Y_exp;
    RMCe_y(:,:,i) = res_y(4).(out_name).Z.RMC_exp;
    
end

R_p     = zeros(20,9); Re_p     = zeros(20,8,9);
Y_p     = zeros(20,9); Ye_p     = zeros(20,8,9);
PI_p    = zeros(20,9); PIe_p    = zeros(20,8,9);
RMC_p   = zeros(20,9); RMCe_p   = zeros(20,8,9);
for i=1:9
    % Name
    out_name = strcat('PI_',num2str(i));
    
    % FG with inflation in levels
    R_p(:,i)    = res_p(4).(out_name).Z.R;
    PI_p(:,i)   = res_p(4).(out_name).Z.PI;
    Y_p(:,i)    = res_p(4).(out_name).Z.Y;
    RMC_p(:,i)  = res_p(4).(out_name).Z.RMC;
    
    % FG with inflation in expectations
    Re_p(:,:,i)   = res_p(4).(out_name).Z.R_exp;
    PIe_p(:,:,i)  = res_p(4).(out_name).Z.PI_exp;
    Ye_p(:,:,i)   = res_p(4).(out_name).Z.Y_exp;
    RMCe_p(:,:,i) = res_p(4).(out_name).Z.RMC_exp;
    
end

% Plot in levels
% stst = [BETA,C,N,Y,PI,R,WP,RMC];
FGsel = [3 9];

figure(6)
subplot(2,2,1)
plot(1:T,100*(PI_y(:,FGsel)/stst(5)-1),...
    1:T,100*(PI_p(:,FGsel)/stst(5)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
plot(1:T,100*(Y_y(:,FGsel)/stst(4)-1),...
    1:T,100*(Y_p(:,FGsel)/stst(4)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
plot(1:T,100*(RMC_y(:,FGsel)/stst(8)-1),...
    1:T,100*(RMC_p(:,FGsel)/stst(8)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
plot(1:T,R_y(:,FGsel),...
    1:T,R_p(:,FGsel),'--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Interest rate','Interpreter','latex')
hL=legend('ZLB FG $y=0.5$','ZLB FG $y=2$','ZLB FG $\pi=0.5$','ZLB FG $\pi=2$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

% Plot in expectations for h=4
figure(7)
subplot(2,2,1)
aux1(:,:) = PIe_y(:,4,FGsel);
aux2(:,:) = PIe_p(:,4,FGsel);
plot(1:T,100*(aux1/stst(5)-1),...
    1:T,100*(aux2/stst(5)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
aux1(:,:) = Ye_y(:,4,FGsel);
aux2(:,:) = Ye_p(:,4,FGsel);
plot(1:T,100*(aux1/stst(4)-1),...
    1:T,100*(aux2/stst(4)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
aux1(:,:) = RMCe_y(:,4,FGsel);
aux2(:,:) = RMCe_p(:,4,FGsel);
plot(1:T,100*(aux1/stst(8)-1),...
    1:T,100*(aux2/stst(8)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
aux1(:,:) = Re_y(:,4,FGsel);
aux2(:,:) = Re_p(:,4,FGsel);
plot(1:T,aux1,...
    1:T,aux2,'--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Expected interest rate','Interpreter','latex')
hL=legend('ZLB FG $y=0.5$','ZLB FG $y=2$','ZLB FG $\pi=0.5$','ZLB FG $\pi=2$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

% Plot in expectations for h=8
figure(8)
subplot(2,2,1)
aux1(:,:) = PIe_y(:,8,FGsel);
aux2(:,:) = PIe_p(:,8,FGsel);
plot(1:T,100*(aux1/stst(5)-1),...
    1:T,100*(aux2/stst(5)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,2)
aux1(:,:) = Ye_y(:,8,FGsel);
aux2(:,:) = Ye_p(:,8,FGsel);
plot(1:T,100*(aux1/stst(4)-1),...
    1:T,100*(aux2/stst(4)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,3)
aux1(:,:) = RMCe_y(:,8,FGsel);
aux2(:,:) = RMCe_p(:,8,FGsel);
plot(1:T,100*(aux1/stst(8)-1),...
    1:T,100*(aux2/stst(8)-1),'--','LineWidth',2),
hold on,plot(1:T,zeros(T,1),'k')
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(2,2,4)
aux1(:,:) = Re_y(:,8,FGsel);
aux2(:,:) = Re_p(:,8,FGsel);
plot(1:T,aux1,...
    1:T,aux2,'--','LineWidth',2),
hold on,plot(1:T,stst(6)*ones(T,1),'k')
title('Expected interest rate','Interpreter','latex')
hL=legend('ZLB FG $y=0.5$','ZLB FG $y=2$','ZLB FG $\pi=0.5$','ZLB FG $\pi=2$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

%% 4) Comparison between the model including the ZLB and with forward
%       guidance checking different levels and types
% Same exercise as before but showing long-run distributions.

% Summarize information
Rs_y     = zeros(50000,9); Res_y     = zeros(50000,8,9);
Ys_y     = zeros(50000,9); Yes_y     = zeros(50000,8,9);
PIs_y    = zeros(50000,9); PIes_y    = zeros(50000,8,9);
RMCs_y   = zeros(50000,9); RMCes_y   = zeros(50000,8,9);
for i=1:9
    % Name
    out_name = strcat('Yg_',num2str(i));
    
    % FG with output in levels
    Rs_y(:,i)    = res_y(8).(out_name).Z.R;
    PIs_y(:,i)   = res_y(8).(out_name).Z.PI;
    Ys_y(:,i)    = res_y(8).(out_name).Z.Y;
    RMCs_y(:,i)  = res_y(8).(out_name).Z.RMC;
    
    % FG with output in expectations
    Res_y(:,:,i)   = res_y(8).(out_name).Z.R_exp;
    PIes_y(:,:,i)  = res_y(8).(out_name).Z.PI_exp;
    Yes_y(:,:,i)   = res_y(8).(out_name).Z.Y_exp;
    RMCes_y(:,:,i) = res_y(8).(out_name).Z.RMC_exp;
    
end

nbins   = 100;

% Plot in levels
% stst = [BETA,C,N,Y,PI,R,WP,RMC];
FGsel = [1 3 9];
grey1 = [0 0 0]+.5;
grey2 = [0 0 0]+.75;

% Distribution of variables in levels
figure(9)
subplot(2,2,1)
[c1,b1] = hist(PIs_y(:,FGsel),nbins);
plot(b1,c1(:,1),'k','LineWidth',2),hold on,
plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
yL = get(gca,'YLim');line([stst(5) stst(5)],yL,'Color','k');
title('Inflation','Interpreter','latex'),axis tight
set(gca,'FontSize',13)
 
subplot(2,2,2)
[c1,b1] = hist(Ys_y(:,FGsel),nbins);
plot(b1,c1(:,1),'k','LineWidth',2),hold on,
plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
yL = get(gca,'YLim');line([stst(4) stst(4)],yL,'Color','k');
title('Output','Interpreter','latex'),axis tight
set(gca,'FontSize',13)
 
subplot(2,2,3)
[c1,b1] = hist(RMCs_y(:,FGsel),nbins);
plot(b1,c1(:,1),'k','LineWidth',2),hold on,
plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
yL = get(gca,'YLim');line([stst(8) stst(8)],yL,'Color','k');
title('Real marginal cost','Interpreter','latex'),axis tight
set(gca,'FontSize',13)
 
subplot(2,2,4)
[c1,b1] = hist(Rs_y(:,FGsel),nbins);
plot(b1,c1(:,1),'k','LineWidth',2),hold on,
plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
yL = get(gca,'YLim');line([stst(6) stst(6)],yL,'Color','k');
title('Interest rate','Interpreter','latex'),axis tight
hL=legend('ZLB no FG','ZLB FG $y=0.5$','ZLB FG $y=2$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)
 
% % Distribution of variables in expectations (h=4)
% figure(10)
% subplot(2,2,1)
% [c1,b1] = hist(PIes_y(:,4,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(5) stst(5)],yL,'Color','k');
% title('Expected inflation','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,2)
% [c1,b1] = hist(Yes_y(:,4,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(4) stst(4)],yL,'Color','k');
% title('Expected output','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,3)
% [c1,b1] = hist(RMCes_y(:,4,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(8) stst(8)],yL,'Color','k');
% title('Expected real marginal cost','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,4)
% [c1,b1] = hist(Res_y(:,4,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(6) stst(6)],yL,'Color','k');
% title('Expected interest rate','Interpreter','latex')
% hL=legend('ZLB no FG','ZLB FG $y=0.5$','ZLB FG $y=2$');
% newPosition = [0 -0.02 1 0.1];
% newUnits = 'normalized';
% set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
%     ,'FontSize',9);
% set(gca,'FontSize',13)
% 
% % Distribution of variables in expectations (h=8)
% figure(11)
% subplot(2,2,1)
% [c1,b1] = hist(PIes_y(:,8,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(5) stst(5)],yL,'Color','k');
% title('Expected inflation','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,2)
% [c1,b1] = hist(Yes_y(:,8,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(4) stst(4)],yL,'Color','k');
% title('Expected output','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,3)
% [c1,b1] = hist(RMCes_y(:,8,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(8) stst(8)],yL,'Color','k');
% title('Expected real marginal cost','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,4)
% [c1,b1] = hist(Res_y(:,8,FGsel),nbins);
% plot(b1,c1(:,1),'k','LineWidth',2),hold on,
% plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
% plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
% yL = get(gca,'YLim');line([stst(6) stst(6)],yL,'Color','k');
% title('Expected interest rate','Interpreter','latex')
% hL=legend('ZLB no FG','ZLB FG $y=0.5$','ZLB FG $y=2$');
% newPosition = [0 -0.02 1 0.1];
% newUnits = 'normalized';
% set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
%     ,'FontSize',9);
% set(gca,'FontSize',13)

%% 5) Welfare

% Reproduce the shock simulation used in previous exercise
[~,~,~,~,BETA_LR,~,~,~,RHO_B,SIGG_B] = parameters;
bN      = 20;
Bs      = sqrt(SIGG_B^2/(1-RHO_B^2));
Bprec   = 5;
Bmin    = -Bprec*Bs+log(BETA_LR);
Bmax    = Bprec*Bs+log(BETA_LR);
    
rng('default')
T           = 50000;
B_proc_sim  = zeros(T+1,1); Bshock = randn(T+1,1); Btrunc = zeros(T+1,1);
for i=2:T+1
    B_proc_sim(i)   = log(BETA_LR)*(1-RHO_B)+RHO_B*B_proc_sim(i-1)+Bshock(i)*SIGG_B;
    B_proc_sim(i)   = max(min(B_proc_sim(i),Bmax),Bmin);
end
B_proc_sim = B_proc_sim(2:end);

W = zeros(50000,9);
FGpar = [0 0.1 0.5 0.7 1 1.2 1.5 1.7 2];
for i=1:9
    % Name
    out_name = strcat('Yg_',num2str(i));
    
    % FG with output in levels
    W(:,i)     = res_y(8).(out_name).Z.U;
end

figure(11)
% subplot(1,2,1)
% plot(FGpar,cumprod(exp(B_proc_sim))'*W,'--o','LineWidth',2)
% %xlabel('$\phi_{FG}$','Interpreter','latex'),
% set(gca,'FontSize',13)

%subplot(1,2,2)
[c1,b1] = hist(W(:,FGsel),nbins);
plot(b1,c1(:,1),'k','LineWidth',2),hold on,
plot(b1,c1(:,2),'--','Color',grey1,'LineWidth',2),
plot(b1,c1(:,3),'-.','Color',grey2,'LineWidth',2),hold off
axis([-1.51 -1.5058 0 21000])
legend('ZLB no FG','ZLB FG $y=0.5$','ZLB FG $y=2$','Location','NorthWest');