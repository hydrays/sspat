% ============
% Response Time
% ============

% load switch1
% [n1, x] = hist(switch1(:,1), 100);
% plot(x, n1, 'g');
% hold on;
% 
% load switch2
% [n2, x] = hist(switch2(:,1), 100);
% plot(x, n2, 'b');
% 
% load switch3
% [n3, x] = hist(switch3(:,1), 100);
% plot(x, n3, 'r');
% 
% xlim([1000, 1100]);
% 
% legend('s1', 's2', 's3')
% print('-depsc','fig_switch_dist.eps')


% ============
% PolII Distribution
% ============
clear all
% load re_hist01
% re_hist = re_hist01;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:10) = sum(re_hist(:, 4))/7;
% plot(x);
hold on;
% 
% load re_hist02
% re_hist = re_hist02;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:10) = sum(re_hist(:, 4))/7;
% plot(x,'g');

load re_hist
%re_hist = re_hist001;
x(1) = 0;
x(2) = sum(re_hist(:, 3));
x(3:10) = sum(re_hist(:, 4))/5;
plot(x,'b');
%print('-depsc','fig_polII_dist.eps')