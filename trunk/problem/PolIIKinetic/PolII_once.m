% ============
% Response Time
% ============

% clear all
% 
% hold on;
% 
% load switch2
% [n2, x] = hist(switch2(:,1), 100);
% plot(x, n2, 'b', 'linewidth', 2);
% 
% 
% load switch3
% [n3, x] = hist(switch3(:,1), 100);
% plot(x, n3, 'r', 'linewidth', 2);
% 
% 
% legend('s1', 's2')
% box on;
% print('-depsc','fig_switch_dista.eps')
% 
% 
% figure
% hold on;
% 
% load switch1
% [n1, x] = hist(switch1(:,1), 100);
% %n1 = n1*5;
% plot(x, n1, 'g', 'linewidth', 2);
% 
% load switch3
% [n3, x] = hist(switch3(:,1), 100);
% plot(x, n3, 'r', 'linewidth', 2);
% 
% load switch2
% [n2, x] = hist(switch2(:,1), 100);
% plot(x, n2, 'b', 'linewidth', 2);
% 
% %xlim([1000, 1200]);
% 
% legend('s1', 's2', 's3')
% box on
% print('-depsc','fig_switch_distb.eps')


% % % ============
% % % PolII Distribution
% % % ============
% clear all
% 
% t = 0:15;
% 
load re_hist0
re_hist = re_hist0;
x(1) = 0;
x(2) = sum(re_hist(:, 2));
x(3:16) = sum(re_hist(:, 3))/7;
plot(t, x,'r','linewidth', 1);
hold on;
% 
load re_hist01
re_hist = re_hist01;
x(1) = 0;
x(2) = sum(re_hist(:, 2));
x(3:16) = sum(re_hist(:, 3))/7;
plot(t, x,'b','linewidth', 2);
% hold on;
% 
% load re_hist01
% re_hist = re_hist01;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'g','linewidth', 1);
% 
% load re_hist015
% re_hist = re_hist015;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'b','linewidth', 1);
% 
% load re_hist02
% re_hist = re_hist02;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'m','linewidth', 1);
% 
% load re_hist025
% re_hist = re_hist025;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'r','linewidth', 1);
% 
% load re_hist03
% re_hist = re_hist03;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'c','linewidth', 1);
% 
% 
% xlabel('DNA position')
% ylabel('number counted')
% legend('s_3=0', 's_3=0.005','s_3=0.01','s_3=0.015');
% print('-depsc','fig_polII_dist.eps')

% load re_hist0005
% re_hist = re_hist0005;
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:10) = sum(re_hist(:, 4))/5;
% plot(x,'g','linewidth', 2);
% 

% t = 0:15;
% hold on
% load re_hist
% x(1) = 0;
% x(2) = sum(re_hist(:, 3));
% x(3:16) = sum(re_hist(:, 4))/7;
% plot(t, x,'g');
% 

% % ============
% % Response Time
% % ============
% clear all
% load re_expression_p01
% n = size(re_expression_p01);
% x = zeros(n(2),1);
% y = zeros(n(2),1);
% for i=1: n(2) 
%     x(i) = sum(re_expression_p01(:,i));
%     %y(i) = 50*abs(sin(0.01*10*i));
% end
% plot(10*(1:n(2)), x,'b');
% hold on;
% load re_expression_p02
% n = size(re_expression_p02);
% x = zeros(n(2),1);
% y = zeros(n(2),1);
% for i=1: n(2) 
%     x(i) = sum(re_expression_p02(:,i));
%     %y(i) = 50*abs(sin(0.01*10*i));
% end
% plot(10*(1:n(2)), x,'r');
% load re_expression_p005
% n = size(re_expression_p005);
% x = zeros(n(2),1);
% y = zeros(n(2),1);
% for i=1: n(2) 
%     x(i) = sum(re_expression_p005(:,i));
%     %y(i) = 50*abs(sin(0.01*10*i));
% end
% plot(10*(1:n(2)), x,'g');
% xlabel('time')
% ylabel('number of E')
% legend('p=0.01', 'p=0.02', 'p=0.005')
% % %plot(y, 'k')
% %text(3000, 500, 'Brg off','fontsize', 18)
% annotation('textarrow',[0.4,0.331],[0.85, 0.8],'String','Brg off','FontSize',14)
% print('-depsc','fig_brg_off.eps')

%print('-depsc','fig_express.eps')
% 
% % ============
% % Signal S1
% % ============
% clear all
% load re_signal_s1a
% re_signal=re_signal_s1a;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 100*re_signal(1,2*i);
% end
% plot(10*(1:m), x,'b');
% hold on;
% plot(10*(1:m), y, 'k', 'linewidth', 2)
% text(1000, 120, 'T=1600','fontsize', 18)
% print('-depsc','fig_signal_s1_large.eps')
% 
% clear all
% figure
% load re_signal_s1b
% re_signal=re_signal_s1b;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 100*re_signal(1,2*i);
% end
% plot(10*(1:m), x,'b');
% hold on;
% plot(10*(1:m), y, 'k', 'linewidth', 2)
% text(1000, 120, 'T=200','fontsize', 18)
% print('-depsc','fig_signal_s1_small.eps')



% % ============
% % Signal S3
% % ============
% clear all
% load re_signal_s3a
% re_signal=re_signal_s3a;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 500*re_signal(1,2*i);
% end
% plot(10*(1:m), x,'b');
% hold on;
% plot(10*(1:m), y, 'k', 'linewidth', 2)
% text(200, 1520, 'T=200, p=0.01','fontsize', 18)
% print('-depsc','fig_signal_s3a.eps')
% 
% clear all
% figure
% load re_signal_s3b
% re_signal=re_signal_s3b;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 500*re_signal(1,2*i);
% end
% plot(10*(1:m), x,'b');
% hold on;
% ylim([0, 8000]);
% plot(10*(1:m), y, 'k', 'linewidth', 2)
% text(200, 1520, 'T=200, p=0.1','fontsize', 18)
% print('-depsc','fig_signal_s3b.eps')


% % ============
% % Signal S3
% % ============
% clear all
% load re_signal
% %re_signal=re_signal_s3_1000;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 100*re_signal(1,2*i);
% end
% plot((1:m), x,'b');
% hold on;
% plot((1:m), y, 'k', 'linewidth', 2)
%xlim([0, 30000])
%text(13000, 120, 'T=6000','fontsize', 18)
%print('-depsc','fig_signal_s1_large.eps')

% clear all
% figure
% load re_signal1000
% re_signal=re_signal1000;
% n = size(re_signal);
% m = n(2)/2;
% x = zeros(m,1);
% y = zeros(m,1);
% for i=1:m 
%     i
%     x(i) = sum(re_signal(:,2*i-1));
%     y(i) = 100*re_signal(1,2*i);
% end
% plot(50*(1:m), x,'b');
% hold on;
% plot(50*(1:m), y, 'k', 'linewidth', 2)
% xlim([0, 30000])
% text(13000, 120, 'T=1000','fontsize', 18)
% print('-depsc','fig_signal_s1_small.eps')

% % % ============
% % % Plot r value for the simplifed model.
% % % ============
% x=0:0.01:1;
% y=0.02;
% z = (x.*x+x-y*y)./(x+y).^2;
% plot(x, z, 'linewidth', 2)
% hold on;
% xlabel('\lambda');
% ylabel('output-input-ratio');
% print('-depsc','rvalue.eps')