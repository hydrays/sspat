% ===========
% plot ode trajectory, three species
% ===========
% y0 = zeros(4,1);
% y0 = [0 10 0 1];
% sum(y0);
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5]);
% [T,Y] = ode45(@ode_onefeed,[0 20],y0,options);
% plot(T,Y(:,1),'r-',T,Y(:,2),'-.m',T,Y(:,3),'-g',T,Y(:,1)+Y(:,2)+Y(:,3),'-b');

% ===========
% plot ode trajectory, four species
% ===========
% y0 = zeros(4,1);
% y0 = [100 300 1600 1];
% sum(y0);
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5]);
% [T,Y] = ode45(@ode_onefeed,[0 400],y0,options);
% plot(T,Y(:,1),'r-',T,Y(:,2),'-.m',T,Y(:,3),'-g',T,Y(:,4),'-k');

% ===========
% plot ode trajectory, 2 feedback self-recover
% ===========
y0 = zeros(3,1);
y0 = [50 100 350];
options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5]);
[T1,Y1] = ode45(@funcs2, [0 50],y0,options);

y0 = zeros(5,1);
y0(1:3) = Y1(length(Y1),:);
y0(4) = 1;
y0(5) = 0;
y0(2) = y0(2) - 1;
options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5 1e-5]);
[T2,Y2] = ode45(@funcs, [0 200],y0,options);
T2 = T2 + 50;
S2 = Y2(:,1) + Y2(:,2) + Y2(:,3) + Y2(:,4) + Y2(:,5);
S1 = Y1(:,1) + Y1(:,2) + Y1(:,3);
plot(T2,Y2(:,1),'r-',T2,Y2(:,2),'-b',T2,Y2(:,3),'-g',T2,Y2(:,4),'-k',T2, S2,'-m','linewidth',2)
hold on
plot(T1,Y1(:,1),'r-',T1,Y1(:,2),'-b',T1,Y1(:,3),'-g',T1,S1,'-m','linewidth',2)
xlabel('t','fontsize', 18);
ylabel('cell population','fontsize',18)
legend('SC', 'TAC', 'TDC', 'MC','total');
print('-depsc','fig_logi_ode_2f_selfrecover.eps')

% % ===========
% % plot ode trajectory, 3 feedback self-recover
% % ===========
% y0 = zeros(3,1);
% y0 = [50 100 350];
% options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5]);
% [T1,Y1] = ode45(@funcs2, [0 50],y0,options);
% 
% y0 = zeros(5,1);
% y0(1:3) = Y1(length(Y1),:);
% y0(4) = 1;
% y0(5) = 0;
% y0(2) = y0(2) - 1;
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5 1e-5]);
% [T2,Y2] = ode45(@funcs, [0 200],y0,options);
% T2 = T2 + 50;
% S2 = Y2(:,1) + Y2(:,2) + Y2(:,3) + Y2(:,4) + Y2(:,5);
% S1 = Y1(:,1) + Y1(:,2) + Y1(:,3);
% plot(T2,Y2(:,1),'r-',T2,Y2(:,2),'-b',T2,Y2(:,3),'-g',T2,Y2(:,4),'-k',T2, S2,'-m','linewidth',2)
% hold on
% plot(T1,Y1(:,1),'r-',T1,Y1(:,2),'-b',T1,Y1(:,3),'-g',T1,S1,'-m','linewidth',2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC','total');
% print('-depsc','fig_logi_ode_3f_selfrecover.eps')

% % ===========
% % plot ratio vs pm
% % ===========
% p1m = 0:0.001:1;
% r = zeros(length(p1m),1);
% for i = 1:length(p1m)
%     [i,length(p1m)]
%     r(i) = funcs(p1m(i));
% end
% plot(p1m, r, '-r*')