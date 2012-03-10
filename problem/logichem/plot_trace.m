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

% % ===========
% % test for mutation
% % ===========
% load pr2_random_mutation_sample1
% load pr2_random_mutation_sample2
% load pr2_random_mutation_sample3;
% load sam4
% sam1 = pr2_random_mutation_sample1;
% sam2 = pr2_random_mutation_sample2;
% sam3 = pr2_random_mutation_sample3;
% load sam1
% load sam2
% load sam3
% load sam4
% clf;
% 
% plot(sam2(1:1000,19), sam2(1:1000,20), '*') 
% xlim([0.5,1])
% title('v_0 change, MTDCs feedback'), xlabel('p_m'), ylabel('v_m')
% subplot(2,2,1),...
% plot(sam3(1:1000,19), sam3(1:1000,20), '*');
% xlim([0.5,1]);
% text(0.51, 0.6, '$\mathbf{\bar{T}=665}$', 'color', 'r', 'fontsize', 12, 'Interpreter', 'latex');
% title('v_0 constant, MTDCs do not feedback'), xlabel('p_m'), ylabel('v_0')
% 
% subplot(2,2,2), plot(sam4(1:1000,19), sam4(1:1000,20), '*'), xlim([0.5,1]), text(0.51, 0.6, '$\mathbf{\bar{T}=22700}$', 'color', 'r', 'fontsize', 12, 'Interpreter', 'latex'), title('v_0 change, MTDCs do not feedback'), xlabel('p_m'), ylabel('v_0')
% 
% subplot(2,2,3), plot(sam1(1:1000,19), sam1(1:1000,20), '*'), xlim([0.5,1]), text(0.51, 0.6, '$\mathbf{\bar{T}=483}$', 'color', 'r', 'fontsize', 12, 'Interpreter', 'latex'), title('v_0 constant, MTDCs feedback'), xlabel('p_m'), ylabel('v_0')
% 
% subplot(2,2,4), plot(sam2(1:1000,19), sam2(1:1000,20), '*'), xlim([0.5,1]), text(0.51, 0.6, '$\mathbf{\bar{T}=678}$', 'color', 'r', 'fontsize', 12, 'Interpreter', 'latex'), title('v_0 change, MTDCs feedback'), xlabel('p_m'), ylabel('v_0')
% 
% print('-depsc','fig_logi_random_mutation.eps')

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
% h = plot(tr(:,1), tr(:,2), '-r','linewidth', 1.1)
% xlim([100, 250])
% hold on
% plot(tr(:,1), tr(:,3), '-b','linewidth', 1.1)
% plot(tr(:,1), tr(:,4), '-g','linewidth', 1.1)
% plot(tr(:,1), tr(:,5), 'k','linewidth', 1.1)
% plot(tr(:,1), tr(:,7), '-m','linewidth', 1.1)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'total')
% print('-depsc','fig_logi_selfrecover_aa_zoomin.eps')
% 
% figure
% hl1 = line(tr(:,1),tr(:,9),'Color','b','linewidth', 1.1);
% xlabel('t','fontsize', 18);
% ylabel('p_0','fontsize',18)
% legend('p_0')
% %ylim([0, 600])
% %set(get(h,'axis'),'FontSize',24)
% box on
% hold on
% ax1 = gca;
% xlim([0, 450])
% ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% line(1000,0, 'Color','b','linewidth', 1.1);
% hl2 = line(tr(:,1),tr(:,10),'Color','r','Parent',ax2,'linewidth', 1.1);
% xlim([0, 450])
% ylabel('v_0','fontsize',18)
% legend('p_0', 'v_0')
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


% ===========
% script for pie
% ===========
% load pr2_sample_system1
% sample = pr2_sample_system1;
% counter = zeros(1,4);
% for i=1:length(sample)
%     if ( sum(sample(i,2:6)) == 0 )
%         counter(4) = counter(4) + 1;
%     end
%     if ( sample(i,5) == 0 && sample(i,6) ~= 0 )
%         counter(3) = counter(3) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) == 0 )
%         counter(1) = counter(1) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) ~= 0)
%         counter(2) = counter(2) + 1;
%     end
% end
% subplot(1, 2, 2), pie(counter(1:3))
% title('System 2')
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
%         counter(3) = counter(3) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) == 0 )
%         counter(1) = counter(1) + 1;
%     end
%     if ( sample(i,5) ~= 0 && sample(i,6) ~= 0)
%         counter(2) = counter(2) + 1;
%     end
% end
% subplot(1, 2, 1), pie(counter(1:3))
% %colormap cool
% title('System 1')
% print('-depsc','fig_logi_mpie.eps')

% % ===========
% % script for stable solution
% % ===========
% load pr2_stable_aa06
% clf;
% tr = pr2_stable_aa06;
% plot(tr(:,1), tr(:,2), 'r','linewidth', 2)
% xlim([0, 200])
% ylim([0, 250])
% partition(1,1) = mean(tr(50:200,2))
% partition(1,2) = mean(tr(50:200,3))
% partition(1,3) = mean(tr(50:200,4))
% hold on
% plot(tr(:,1), tr(:,3), 'b','linewidth', 2)
% plot(tr(:,1), tr(:,4), 'g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 2)
% % figure
% plot(tr(:,1), tr(:,7), 'm','linewidth', 2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'Total')
% print('-depsc','fig_logi_stable_aa06.eps')
% 
% load pr2_stable_aa08
% clf;
% tr = pr2_stable_aa08;
% plot(tr(:,1), tr(:,2), 'r','linewidth', 2)
% xlim([0, 200])
% ylim([0, 250])
% partition(2,1) = mean(tr(50:200,2))
% partition(2,2) = mean(tr(50:200,3))
% partition(2,3) = mean(tr(50:200,4))
% hold on
% plot(tr(:,1), tr(:,3), 'b','linewidth', 2)
% plot(tr(:,1), tr(:,4), 'g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 2)
% % figure
% plot(tr(:,1), tr(:,7), 'm','linewidth', 2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'Total')
% print('-depsc','fig_logi_stable_aa08.eps')
% 
% load pr2_stable_aa10
% clf;
% tr = pr2_stable_aa10;
% plot(tr(:,1), tr(:,2), 'r','linewidth', 2)
% xlim([0, 200])
% ylim([0, 250])
% partition(3,1) = mean(tr(50:200,2))
% partition(3,2) = mean(tr(50:200,3))
% partition(3,3) = mean(tr(50:200,4))
% hold on
% plot(tr(:,1), tr(:,3), 'b','linewidth', 2)
% plot(tr(:,1), tr(:,4), 'g','linewidth', 2)
% %plot(tr(:,1), tr(:,5), 'k','linewidth', 2)
% % figure
% plot(tr(:,1), tr(:,7), 'm','linewidth', 2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'Total')
% print('-depsc','fig_logi_stable_aa10.eps')
% 
% figure
% map(1, :) = [1, 0, 0];
% map(2, :) = [0, 0, 1];
% map(3, :) = [0, 1, 0];
% colormap(map)
% bar(partition)
% legend('SC', 'TAC', 'TDC')
% set(gca,'xticklabel',{'k=0.6','k=0.8','k=1'}, 'fontsize', 12)
% print('-depsc','fig_logi_stable_aa_partition.eps')

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

% ===========
% script pm ultra
% ===========
load pr2_p0v0_data2
data = pr2_p0v0_data2
ymap = 0.1:.1:3.009;
xmap = 0.0:0.02:1.0001;
[X,Y] = meshgrid(xmap, ymap);
wg_barmap = reshape((data(:,6)+data(:,7))./(data(:,3)+data(:,4)+data(:,5)+data(:,6)+data(:,7)), length(ymap), length(xmap))
wg_barmap2 = reshape((data(:,6))./(data(:,3)+data(:,4)+data(:,5)+data(:,6)+data(:,7)), length(ymap), length(xmap))
xmap = 0.5:0.02:1.0001;
[X,Y] = meshgrid(xmap, ymap);
wg_barmap = wg_barmap(:, 26:51)
wg_barmap2 = wg_barmap2(:, 26:51)
pcolor(Y, X, wg_barmap)
%surf(X, Y, wg_barmap)

colorbar;
xlabel('v_m','fontsize', 18);
ylabel('p_m','fontsize',18);
zlabel('ratio of mutant cell','fontsize',18);
%xlim([0.5,1]);
print('-depsc','fig_logi_pmvm.eps')

% figure
% plot(xmap, wg_barmap2(7, :),'g', 'linewidth', 2)
% hold on
% plot(xmap, wg_barmap2(10, :),'linewidth', 2)
% plot(xmap, wg_barmap2(15, :),'r', 'linewidth', 2)
% %plot(xmap, wg_barmap(25, :),'r', 'linewidth', 2)
% xlim([0.5,1]);
% legend('v_m = 0.7', 'v_m=1', 'v_m=1.5','location', 'northwest')
% 
% xlabel('p_m','fontsize', 18);
% ylabel('ratio of cancer stem cell','fontsize',18);
% print('-depsc','fig_logi_pmvm_slicev.eps')
figure

for index = 1:length(ymap)
    pmax(index) = 0;
    prob(index) = 0;
    for temp_j = 1:length(xmap)
        if ( wg_barmap(index, temp_j) > pmax(index) )
            pmax(index) = wg_barmap(index, temp_j)
            %pmax(index) = xmap(temp_j)
            prob(index) = xmap(temp_j);
            csc(index) = wg_barmap2(index, temp_j)
            ctt(index) = wg_barmap(index, temp_j)
        end
    end
end
hold on
%plot(ymap, prob, 'o', 'markersize', 10, 'LineWidth',2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor',[.49 1 .63])
box on
plot(ymap, csc,'ks', 'LineWidth',2, 'markersize', 12, 'markerfacecolor', 'r')
% plot(ymap, ctt,'ks', 'LineWidth',2, 'markersize', 12, 'markerfacecolor', 'g')

xlabel('v_m','fontsize', 18);
ylabel('ratio of cancer stem cell','fontsize',18);
print('-depsc','fig_logi_pmvm_csc.eps')
    
