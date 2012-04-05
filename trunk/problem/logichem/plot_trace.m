% ===========
% test
% ===========
load tr
clf;
tr = tr;
plot(tr(:,1), tr(:,2), '-r')
%xlim([150, 350])
hold on
plot(tr(:,1), tr(:,3), '-b')
plot(tr(:,1), tr(:,4), '-g')
plot(tr(:,1), tr(:,5), '-k')
% figure
plot(tr(:,1), tr(:,7), '-m')
xlabel('t','fontsize', 18);
ylabel('cell population','fontsize',18)
legend('SC', 'TAC', 'TDC', 'MC', 'Total')

% figure

% % ===========
% % test
% % ===========
% load test2
% data = test2;
% vm = 0.1:.1:3.009;
% pm = 0.0:0.02:1.0001;
% [X,Y] = meshgrid(pm, vm);
% wg_barmap = reshape(data(:,3), length(vm), length(pm));
% pcolor(X, Y, wg_barmap);
% xlabel('p_m','fontsize', 18);
% ylabel('v_m','fontsize',18);
% zlabel('ratio of mutant cell','fontsize',18);
% print('-depsc','pmvm_response_nf_ayms.eps')
% 
% figure
% wg_barmap = reshape(data(:,6), length(vm), length(pm));
% pcolor(X, Y, wg_barmap);
% xlabel('p_m','fontsize', 18);
% ylabel('v_m','fontsize',18);
% zlabel('ratio of mutant cell','fontsize',18);
% print('-depsc','pmvm_response_nf_ayms.eps')


% load test2
% data = test2;
% hold on
% for i=1: length(data)
%     plot(data(i, 2)*(data(i,1)-0.5), data(i, 6),'*')
% end

% % ===========
% % script for 2 feedback, self-recovery, asym/asym
% % ===========
% load pr2_slefrecover_logi_aa
% tr = pr2_slefrecover_logi_aa;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,7), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_selfrecover_aa.eps')
% 
% figure
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,9), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,10), '-r','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('p_0, v_0','fontsize',18)
% legend('p_0', 'v_0')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% box on
% print('-depsc','fig_logi_selfrecover_pvplot_aa.eps')

% % ===========
% % script for 2 feedback, self-recovery, sym/sym
% % ===========
% load pr2_slefrecover_logi_ss
% tr = pr2_slefrecover_logi_ss;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,7), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_selfrecover_ss.eps')
% 
% figure
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,9), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,10), '-r','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('p_0, v_0','fontsize',18)
% legend('p_0', 'v_0')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% box on
% print('-depsc','fig_logi_selfrecover_pvplot_ss.eps')

% % ===========
% % script for 2 feedback, self-recovery, asym/sym
% % ===========
% load pr2_slefrecover_logi_as
% tr = pr2_slefrecover_logi_as;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,7), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_selfrecover_as.eps')
% 
% figure
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,9), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,10), '-r','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('p_0, v_0','fontsize',18)
% legend('p_0', 'v_0')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% box on
% print('-depsc','fig_logi_selfrecover_pvplot_as.eps')

% % ===========
% % script for 2 feedback, self-recovery, sym/asym
% % ===========
% load pr2_slefrecover_logi_sa
% tr = pr2_slefrecover_logi_sa;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,7), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_selfrecover_sa.eps')
% 
% figure
% xlim([0, 450])
% hold on
% plot(tr(:,1), tr(:,9), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,10), '-r','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('p_0, v_0','fontsize',18)
% legend('p_0', 'v_0')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% box on
% print('-depsc','fig_logi_selfrecover_pvplot_sa.eps')
% 
% % ===========
% % script for mutation system1
% % ===========
% load pr2_system1_trace 
% tr = pr2_system1_trace;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.0)f
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.0)
% plot(tr(:,1), tr(:,5), '-k','linewidth', 1.0)
% plot(tr(:,1), tr(:,6), ':k','linewidth', 1.0)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MTAC','MSC','location','northwest')
% %plot(tr(length(tr),1), tr(length(tr),2), 'k*', 'markersize', 8)
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_mutation_system1.eps')



% % ===========
% % script for mutation system2
% % ===========
% load pr2_system2_trace 
% tr = pr2_system2_trace;
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.0)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.0)
% plot(tr(:,1), tr(:,5), '-k','linewidth', 1.0)
% plot(tr(:,1), tr(:,6), ':k','linewidth', 1.0)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MTAC','MSC','location','northwest')
% %plot(tr(length(tr),1), tr(length(tr),2), 'k*', 'markersize', 8)
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_mutation_system2.eps')


% % ===========
% % script for pie
% % ===========
% load pr2_sample_system1
% sample = pr2_sample_system1;
% counter = zeros(1,4);
% for i=1:length(sample)
%     if ( sum(sample(i,2:6)) == 0 )
%         counter(4) = counter(4) + 1;
%     end
%     if ( sample(i,5) == 0 && sample(i,6) ~= 0 )
%         counter(1) = counter(1) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) == 0 )
%         counter(3) = counter(3) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) ~= 0)
%         counter(2) = counter(2) + 1;
%     end
% end
% subplot(1, 2, 1), pie(counter(1:3))
% %colormap cool
% %print('-depsc','fig_logi_mpie_system1.eps')
% 
% % ===========
% % script for pie
% % ===========
% load pr2_sample_system2
% sample = pr2_sample_system2;
% counter = zeros(1,4);
% for i=1:length(sample)
%     if ( sum(sample(i,2:6)) == 0 )
%         counter(4) = counter(4) + 1;
%     end
%     if ( sample(i,5) == 0 && sample(i,6) ~= 0 )
%         counter(1) = counter(1) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) == 0 )
%         counter(3) = counter(3) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) ~= 0)
%         counter(2) = counter(2) + 1;
%     end
% end
% subplot(1, 2, 2), pie(counter(1:3))
% %colormap cool
% print('-depsc','fig_logi_mpie.eps')

% % ===========
% % script for stable solution
% % ===========
% load pr2_stable_ss06
% clf;
% tr = pr2_stable_ss06;
% plot(tr(:,1), tr(:,2), 'r','linewidth', 2)
% xlim([0, 200])
% ylim([0, 250])
% hold on
% plot(tr(:,1), tr(:,3), 'b','linewidth', 2)
% plot(tr(:,1), tr(:,4), 'g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 2)
% % figure
% plot(tr(:,1), tr(:,7), 'm','linewidth', 2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'Total')
% print('-depsc','fig_logi_stable_ss06.eps')

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

% % ===========
% % script pm
% % ===========
% load paperresult_2f_pm
% load paperresult_3f_pm
% datapm1 = paperresult_2f_pm;
% datapm2 = paperresult_3f_pm;
% h = plot(datapm1(:,1), (datapm1(:,5)+datapm1(:,6))./datapm1(:,7), '-r*','markersize', 10, 'linewidth', 1.1)
% hold on
% plot(datapm2(:,1), (datapm2(:,5)+datapm2(:,6))./datapm2(:,7), '-o','markersize', 10, 'linewidth', 1.1)
% xlabel('p_m','fontsize', 18);
% ylabel('ratio of mutant type','fontsize',18)
% ylim([0,1]);
% legend('feedback scheme II', 'feedback scheme III', 2)
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_pm.eps')
% 
% % ===========
% % script pm ultra
% % ===========
% load pr2_p0v0_data2
% data = pr2_p0v0_data2
% ymap = 0.1:.1:3.009;
% xmap = 0.0:0.02:1.0001;
% [X,Y] = meshgrid(xmap, ymap);
% wg_barmap = reshape((data(:,6)+data(:,7))./(data(:,3)+data(:,4)+data(:,5)+data(:,6)+data(:,7)), length(ymap), length(xmap))
% surf(X, Y, wg_barmap)
% xlabel('p_m','fontsize', 18);
% ylabel('v_m','fontsize',18);
% zlabel('ratio of mutant cell','fontsize',18);
% xlim([0.5,1]);
% %print('-depsc','fig_logi_pmvm.eps')
% figure
% plot(xmap, wg_barmap(10, :),'linewidth', 2)
% hold on
% plot(xmap, wg_barmap(7, :),'g', 'linewidth', 2)
% plot(xmap, wg_barmap(15, :),'r', 'linewidth', 2)
% %plot(xmap, wg_barmap(25, :),'r', 'linewidth', 2)
% xlim([0.5,1]);
% legend('v_m = 0.7', 'v_m=1', 'v_m=1.5','location', 'northwest')
% xlabel('p_m','fontsize', 18);
% ylabel('ratio of mutant cell','fontsize',18);
% print('-depsc','fig_logi_pmvm_slicev.eps')
