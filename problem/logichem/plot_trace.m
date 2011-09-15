% % ===========
% % script for 1 feedback
% % ===========
% load paperresult_1f_tr_k1
% tr = paperresult_1f_tr_k1;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% %xlim([0, 600])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,8), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'total')
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_1f_k1.eps')
% 
% % ===========
% % script for 1 feedback bar-plot
% % ===========
% mycmap = [1     0     0
%      0     0     1
%      0     1     0];
% colormap(mycmap)
% load paperresult_1f_tr_k07
% load paperresult_1f_tr_k85
% load paperresult_1f_tr_k1
% tr3 = paperresult_1f_tr_k07;
% tr2 = paperresult_1f_tr_k85;
% tr1 = paperresult_1f_tr_k1;
% 
% xpmean(1,:) = mean(tr1(600:1000,2:4));
% xpmean(2,:) = mean(tr2(600:1000,2:4));
% xpmean(3,:) = mean(tr3(600:1000,2:4));
% bar(xpmean,1)
% xlabel('k_1','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC')
% %set(get(h,'axis'),'FontSize',24)
% set(gca,'XTickLabel',{'1', '0.85', '0.7'},'fontsize', 14);
% box on;
% print('-depsc','fig_logi_1f_bar.eps')


% % ===========
% % script for 2 feedback, woundhealing
% % ===========
% load paperresult_tr_2f_wondhealing 
% tr = paperresult_tr_2f_wondhealing;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 2)
% xlim([0, 200])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 2)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC')
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_2f_woundhealing.eps')
% figure
% plot(tr(:,1), tr(:,7), '-r','linewidth', 2)
% hold on
% xlim([0, 200])
% plot(tr(:,1), tr(:,10), '-g','linewidth', 1.5)
% %plot(tr(:,1), tr(:,11), '-b','linewidth', 1.1)
% legend('p_0','v_0')
% xlabel('t','fontsize', 18);
% ylabel('Stem cell property','fontsize', 18)
% print('-depsc','fig_logi_2f_woundhealing_supp.eps')

% % ===========
% % script for 2 feedback, self-recovery
% % ===========
% load paperresult_tr_2f_selfrecover
% tr = paperresult_tr_2f_selfrecover;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,8), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_2f_selfrecover.eps')


% % ===========
% % script for 3 feedback, woundhealing
% % ===========
% load paperresult_tr_3f_wondhealing
% tr = paperresult_tr_3f_wondhealing;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 2)
% xlim([0, 200])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 2)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC')
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_3f_woundhealing.eps')
% figure
% plot(tr(:,1), tr(:,7), '-r','linewidth', 2)
% hold on
% xlim([0, 200])
% plot(tr(:,1), tr(:,10), '-g','linewidth', 2)
% plot(tr(:,1), tr(:,11), '-b','linewidth', 2)
% legend('p_0','v_0', 'p_{sym}')
% xlabel('t','fontsize', 18);
% ylabel('Stem cell property','fontsize', 18)
% print('-depsc','fig_logi_3f_woundhealing_supp.eps')


% % ===========
% % script for 3 feedback, self-recovery
% % ===========
% load paperresult_tr_3f_selfrecover
% tr = paperresult_tr_3f_selfrecover;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,8), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_3f_selfrecover.eps')

% 
% % ===========
% % script for 1 feedback, mutation
% % ===========
% load paperresult_geneticstability_tr1f
% tr = paperresult_geneticstability_tr1f;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.0)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.0)
% plot(tr(:,1), tr(:,5), ':k','linewidth', 1.0)
% plot(tr(:,1), tr(:,6), '-k','linewidth', 2.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC','TDC', 'MC1', 'MC2',2)
% %plot(tr(length(tr),1), tr(length(tr),3), 'k*', 'markersize', 8)
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_1f_mutation.eps')


% % ===========
% % script for 2 feedback, mutation
% % ===========
% load paperresult_geneticstability_tr2f_k09 
% tr = paperresult_geneticstability_tr2f_k09;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.0)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.0)
% plot(tr(:,1), tr(:,5), ':k','linewidth', 1.0)
% plot(tr(:,1), tr(:,6), '-k','linewidth', 1.0)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC')
% plot(tr(length(tr),1), tr(length(tr),2), 'k*', 'markersize', 8)
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_2f_mutation.eps')


% 
% % ===========
% % script for 3 feedback, mutation
% % ===========
% load paperresult_geneticstability_tr3f
% tr = paperresult_geneticstability_tr3f;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.0)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.0)
% plot(tr(:,1), tr(:,5), ':k','linewidth', 1.0)
% plot(tr(:,1), tr(:,6), '-k','linewidth', 1.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC')
% plot(tr(length(tr),1), tr(length(tr),2), 'k*', 'markersize', 8)
% ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_3f_mutation.eps')


% % ===========
% % script for srI_short
% % ===========
% load tr
% tr = tr;
% plot(tr(:,1), tr(:,2), 'r')
% xlim([350, 550])
% hold on
% plot(tr(:,1), tr(:,3), 'm')
% plot(tr(:,1), tr(:,4), 'g')
% plot(tr(:,1), tr(:,5), 'k')
% plot(tr(:,1), tr(:,8), 'b')
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'Total')
% print('-depsc','srI_short.eps')

% % ===========
% % script for srI_pv
% % ===========
% load tr
% tr = tr;
% plot(tr(:,1), tr(:,7), 'b')
% xlim([00, 600])
% hold on
% plot(tr(:,1), tr(:,10), 'r')
% xlabel('t','fontsize', 18);
% ylabel('value','fontsize',18)
% legend('p_0', 'v_0')
% print('-depsc','srI_pv.eps')

% % ===========
% % plot the average MC vs pm
% % ===========
% hold on
% load xbar_pm
% pm = xbar_pm
% l = length(pm)
% plot(xbar_pm(1:5:l,1), (pm(1:5:l,5)+pm(1:5:l,6))./pm(1:5:l,7), '-*')
% plot(p1m(1:5:l), r(1:5:l), '-r^')
% xlabel('p_m','fontsize', 18)
% ylim([0,1])
% box on;
% ylabel('MCs ratio','fontsize', 18)
% legend('sto', 'det', 'Location', 'northwest')
% print('-depsc','fig_pm.eps')

% ===========
% script for sym and asym division
% ===========
% load tr
% tr = tr;
% plot(tr(:,1), tr(:,2), 'r')
% %xlim([0, 600])
% hold on
% plot(tr(:,1), tr(:,3), 'm')
% plot(tr(:,1), tr(:,4), 'g')
% plot(tr(:,1), tr(:,5), 'k')
% plot(tr(:,1), tr(:,8), 'b')
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'Total')
% print('-depsc','asym.eps')

% ===========
% script pm
% ===========
load paperresult_2f_pm
load paperresult_3f_pm
datapm1 = paperresult_2f_pm;
datapm2 = paperresult_3f_pm;
h = plot(datapm1(:,1), (datapm1(:,5)+datapm1(:,6))./datapm1(:,7), '-r*','markersize', 10, 'linewidth', 1.1)
hold on
plot(datapm2(:,1), (datapm2(:,5)+datapm2(:,6))./datapm2(:,7), '-o','markersize', 10, 'linewidth', 1.1)
xlabel('p_m','fontsize', 18);
ylabel('ratio of mutant type','fontsize',18)
ylim([0,1]);
legend('feedback scheme II', 'feedback scheme III', 2)
%set(get(h,'axis'),'FontSize',24)
print('-depsc','fig_logi_pm.eps')