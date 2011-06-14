% ===========
% script for srI_long
% ===========
% load tr
% tr = tr;
% plot(tr(:,1), tr(:,2), 'r')
% xlim([0, 600])
% hold on
% plot(tr(:,1), tr(:,3), 'm')
% plot(tr(:,1), tr(:,4), 'g')
% plot(tr(:,1), tr(:,5), 'k')
% plot(tr(:,1), tr(:,8), 'b')
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC', 'Total')
% print('-depsc','srI_long.eps')

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

% ===========
% plot the average MC vs pm
% ===========
hold on
load xbar_pm
pm = xbar_pm
l = length(pm)
plot(xbar_pm(1:5:l,1), (pm(1:5:l,5)+pm(1:5:l,6))./pm(1:5:l,7), '-*')
plot(p1m(1:5:l), r(1:5:l), '-r^')
xlabel('p_m','fontsize', 18)
ylim([0,1])
box on;
ylabel('MCs ratio','fontsize', 18)
legend('sto', 'det', 'Location', 'northwest')
print('-depsc','fig_pm.eps')