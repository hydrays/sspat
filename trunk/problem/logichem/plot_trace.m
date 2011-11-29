% % ===========
% % test
% % ===========
% load tr
% clf;
% tr = tr;
% plot(tr(:,1), tr(:,2), 'r')
% %xlim([150, 350])
% hold on
% plot(tr(:,1), tr(:,3), 'b')
% plot(tr(:,1), tr(:,4), 'g')
% plot(tr(:,1), tr(:,5), 'k')
% % figure
% plot(tr(:,1), tr(:,7), 'm')
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'Total')
% print('-depsc','asym.eps')

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
% % script for mutation
% % ===========
% load pr_mutation 
% tr = pr_mutation;
% h = plot(tr(1:5000,1), tr(1:5000,2), '-r','linewidth', 1.0)
% %xlim([0, 250])
% hold on
% plot(tr(1:5000,1), tr(1:5000,3), '-b','linewidth', 1.0)
% plot(tr(1:5000,1), tr(1:5000,4), '-g','linewidth', 1.0)
% plot(tr(1:5000,1), tr(1:5000,5), '-k','linewidth', 1.0)
% %plot(tr(1:5000,1), tr(1:5000,6), '-k','linewidth', 1.0)
% %plot(tr(:,1), tr(:,8), '-m','linewidth', 1.0)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC')
% %plot(tr(length(tr),1), tr(length(tr),2), 'k*', 'markersize', 8)
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% print('-depsc','fig_logi_mutation.eps')
% 
% counter = 0;
% for i=1:length(sample_aa)
%     if ( sample_aa(i, 2) == 0 )
%         counter = counter + 1;
%     end
% end

% ===========
% script for pie
% ===========
load sample_survive
counter = zeros(1,4);
for i=1:length(sample_survive)
    if ( sample_survive(i,12) ~= 0 )
        counter(4) = counter(4) + 1;
    end
    if ( sample_survive(i,11) ~= 0 && sample_survive(i,12) == 0 && sample_survive(i, 2)~=0)
        counter(2) = counter(2) + 1;
    end
    if ( sample_survive(i,11) == 0 && sample_survive(i,12) == 0 )
        counter(1) = counter(1) + 1;
    end
    if ( sample_survive(i,11) ~= 0 && sample_survive(i,12) == 0 && sample_survive(i, 2)==0)
        counter(3) = counter(3) + 1;
    end
end
pie(counter)
colormap cool
print('-depsc','fig_logi_mpie.eps')

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