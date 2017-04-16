% Figure 3: ergodic density for level variables
nbins   = 100;
figure(3)
subplot(3,2,1)
[c1,b1] = hist(PIs,nbins);
[c2,b2] = hist(PIsz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(6) stst(6)],yL,'Color','k');
title('Inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,2)
[c1,b1] = hist(Ys,nbins);
[c2,b2] = hist(Ysz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(5) stst(5)],yL,'Color','k');
title('Output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,3)
[c1,b1] = hist(Cs,nbins);
[c2,b2] = hist(Csz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(3) stst(3)],yL,'Color','k');
title('Consumption','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,4)
[c1,b1] = hist(RMCs,nbins);
[c2,b2] = hist(RMCsz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(9) stst(9)],yL,'Color','k');
title('Real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,5)
[c1,b1] = hist(Ygs,nbins);
[c2,b2] = hist(Ygsz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(11) stst(11)],yL,'Color','k');
title('Output gap','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,6)
[c1,b1] = hist(Rs,nbins);
[c2,b2] = hist(Rsz,nbins);
a1 = min([min(find(c1(:)>0)) min(find(c2(:)>0))]);
a2 = max([max(find(c1(:)>0)) max(find(c2(:)>0))]);
plot(b1(a1:a2),100*c1(a1:a2)/T,b2(a1:a2),100*c2(a1:a2)/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(7) stst(7)],yL,'Color','k');axis tight
title('Interest rate','Interpreter','latex')
legend('No ZLB','ZLB','Location','NorthEast')
set(gca,'FontSize',13)

% % Figure 4: ergodic cumulative distribution for level variables
% figure(4)
% subplot(2,2,1)
% [c1,b1] = hist(PIs,nbins);
% [c2,b2] = hist(PIsz,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(6),3)+1),aux2(lookup(b2,stst(6),3)+1));
% plot(b1,aux1,b2,aux2,'--','LineWidth',2)
% line([stst(6) stst(6)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(6)],[aux1(lookup(b1,stst(6),3)+1) aux1(lookup(b1,stst(6),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(6)],[aux2(lookup(b2,stst(6),3)+1) aux2(lookup(b2,stst(6),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Inflation','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,2)
% [c1,b1] = hist(Ys,nbins);
% [c2,b2] = hist(Ysz,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(5),3)+1),aux2(lookup(b2,stst(5),3)+1));
% plot(b1,aux1,b2,aux2,'--','LineWidth',2)
% line([stst(5) stst(5)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux1(lookup(b1,stst(5),3)+1) aux1(lookup(b1,stst(5),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux2(lookup(b2,stst(5),3)+1) aux2(lookup(b2,stst(5),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Output','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,3)
% [c1,b1] = hist(Cs,nbins);
% [c2,b2] = hist(Csz,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(3),3)+1),aux2(lookup(b2,stst(3),3)+1));
% plot(b1,aux1,b2,aux2,'--','LineWidth',2)
% line([stst(3) stst(3)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux1(lookup(b1,stst(3),3)+1) aux1(lookup(b1,stst(3),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux2(lookup(b2,stst(3),3)+1) aux2(lookup(b2,stst(3),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Consumption','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,4)
% [c1,b1] = hist(Rs,nbins);
% [c2,b2] = hist(Rsz,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(7),3)+1),aux2(lookup(b2,stst(7),3)+1));
% plot(b1,aux1,b2,aux2,'--','LineWidth',2)
% line([stst(7) stst(7)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(7)],[aux1(lookup(b1,stst(7),3)+1) aux1(lookup(b1,stst(7),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(7)],[aux2(lookup(b2,stst(7),3)+1) aux2(lookup(b2,stst(7),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Interest rate','Interpreter','latex')
% legend('No ZLB','ZLB','Location','NorthWest')
% set(gca,'FontSize',13)

% Figure 5: ergodic density for expected variables
figure(5)
subplot(3,2,1)
[c3,b3] = hist(PIs_exp,nbins);
[c4,b4] = hist(PIsz_exp,nbins);
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(6) stst(6)],yL,'Color','k');axis tight
title('Expected inflation','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,2)
[c3,b3] = hist(Ys_exp,nbins);
[c4,b4] = hist(Ysz_exp,nbins);
%plot(b3,100*c3(:,irf_h_sel)/T,b4,100*c4(:,irf_h_sel)/T,'--','LineWidth',2)
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(5) stst(5)],yL,'Color','k');axis tight
title('Expected output','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,3)
[c3,b3] = hist(Cs_exp,nbins);
[c4,b4] = hist(Csz_exp,nbins);
%plot(b3,100*c3(:,irf_h_sel)/T,b4,100*c4(:,irf_h_sel)/T,'--','LineWidth',2)
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(3) stst(3)],yL,'Color','k');axis tight
title('Expected consumption','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,4)
[c3,b3] = hist(RMCs_exp,nbins);
[c4,b4] = hist(RMCsz_exp,nbins);
%plot(b3,100*c3(:,irf_h_sel)/T,b4,100*c4(:,irf_h_sel)/T,'--','LineWidth',2)
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(9) stst(9)],yL,'Color','k');axis tight
title('Expected real marginal cost','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,5)
[c3,b3] = hist(Ygs_exp,nbins);
[c4,b4] = hist(Ygsz_exp,nbins);
%plot(b3,100*c3(:,irf_h_sel)/T,b4,100*c4(:,irf_h_sel)/T,'--','LineWidth',2)
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(11) stst(11)],yL,'Color','k');axis tight
title('Expected output gap','Interpreter','latex')
set(gca,'FontSize',13)

subplot(3,2,6)
[c3,b3] = hist(Rs_exp,nbins);
[c4,b4] = hist(Rsz_exp,nbins);
%plot(b3,100*c3(:,irf_h_sel)/T,b4,100*c4(:,irf_h_sel)/T,'--','LineWidth',2)
a1 = zeros(length(irf_h_sel),1); a2 = a1;
for i=1:length(irf_h_sel)
    a1(i) = min([min(find(c3(:,irf_h_sel(i))>0)) min(find(c4(:,irf_h_sel(i))>0))]);
    a2(i) = max([max(find(c3(:,irf_h_sel(i))>0)) max(find(c4(:,irf_h_sel(i))>0))]);
end
plot(b3(a1(1):a2(1)),100*c3(a1(1):a2(1),irf_h_sel(1))/T,...
    b3(a1(2):a2(2)),100*c3(a1(2):a2(2),irf_h_sel(2))/T,...
    b4(a1(1):a2(1)),100*c4(a1(1):a2(1),irf_h_sel(1))/T,'--',...
    b4(a1(2):a2(2)),100*c4(a1(2):a2(2),irf_h_sel(2))/T,'--','LineWidth',2)
yL = get(gca,'YLim');line([stst(7) stst(7)],yL,'Color','k');axis tight
title('Expected interest rate','Interpreter','latex')
%legend('No ZLB','ZLB','Location','NorthEast')
hL=legend('No ZLB, $h=4$','No ZLB, $h=8$','ZLB, $h=4$','ZLB, $h=8$');
newPosition = [0 -0.02 1 0.1];
newUnits = 'normalized';
set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
    ,'FontSize',9);
set(gca,'FontSize',13)

% % Figure 6: ergodic cumulative distribution for expected variables
% figure(6)
% subplot(2,2,1)
% [c1,b1] = hist(PIs_exp,nbins);
% [c2,b2] = hist(PIsz_exp,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(6),3)+1),aux2(lookup(b2,stst(6),3)+1));
% plot(b1,aux1(:,irf_h_sel),b2,aux2(:,irf_h_sel),'--','LineWidth',2)
% line([stst(6) stst(6)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(6)],[aux1(lookup(b1,stst(6),3)+1) aux1(lookup(b1,stst(6),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(6)],[aux2(lookup(b2,stst(6),3)+1) aux2(lookup(b2,stst(6),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Expected inflation','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,2)
% [c1,b1] = hist(Ys_exp,nbins);
% [c2,b2] = hist(Ysz_exp,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(5),3)+1),aux2(lookup(b2,stst(5),3)+1));
% plot(b1,aux1(:,irf_h_sel),b2,aux2(:,irf_h_sel),'--','LineWidth',2)
% line([stst(5) stst(5)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux1(lookup(b1,stst(5),3)+1) aux1(lookup(b1,stst(5),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux2(lookup(b2,stst(5),3)+1) aux2(lookup(b2,stst(5),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Expected output','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,3)
% [c1,b1] = hist(Cs_exp,nbins);
% [c2,b2] = hist(Csz_exp,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(3),3)+1),aux2(lookup(b2,stst(3),3)+1));
% plot(b1,aux1(:,irf_h_sel),b2,aux2(:,irf_h_sel),'--','LineWidth',2)
% line([stst(3) stst(3)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux1(lookup(b1,stst(3),3)+1) aux1(lookup(b1,stst(3),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(5)],[aux2(lookup(b2,stst(3),3)+1) aux2(lookup(b2,stst(3),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Expected consumption','Interpreter','latex')
% set(gca,'FontSize',13)
% 
% subplot(2,2,4)
% [c1,b1] = hist(Rs_exp,nbins);
% [c2,b2] = hist(Rsz_exp,nbins);
% aux1 = cumsum(c1)/T;
% aux2 = cumsum(c2)/T;
% aux3 = max(aux1(lookup(b1,stst(7),3)+1),aux2(lookup(b2,stst(7),3)+1));
% plot(b1,aux1(:,irf_h_sel),b2,aux2(:,irf_h_sel),'--','LineWidth',2)
% line([stst(7) stst(7)],[aux3 0],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(7)],[aux1(lookup(b1,stst(7),3)+1) aux1(lookup(b1,stst(7),3)+1)],'Color',[0 0 0],'LineStyle','--');
% line([min([b1(:);b2(:)]) stst(7)],[aux2(lookup(b2,stst(7),3)+1) aux2(lookup(b2,stst(7),3)+1)],'Color',[0 0 0],'LineStyle','--');
% axis tight
% title('Expected interest rate','Interpreter','latex')
% %legend('No ZLB','ZLB','Location','NorthWest')
% hL=legend('No ZLB, $h=4$','No ZLB, $h=8$','ZLB, $h=4$','ZLB, $h=8$');
% newPosition = [0 -0.02 1 0.1];
% newUnits = 'normalized';
% set(hL,'Position',newPosition,'Units',newUnits,'Orientation','horizontal'...
%     ,'FontSize',9);
% set(gca,'FontSize',13)