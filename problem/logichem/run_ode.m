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
% plot ode trajectory, five species
% ===========
% y0 = zeros(5,1);
% y0 = [200 300 600 0.00000000000001 0];
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5 1e-5]);
% [T,Y] = ode45(@funcs, [0 600],y0,options);
% plot(T,Y(:,1),'r-',T,Y(:,2),'-.m',T,Y(:,3),':g',T,Y(:,4),'--k','linewidth',2)
% xlabel('t','fontsize', 18);
% ylabel('cell population','fontsize',18)
% legend('SC', 'TAC', 'TDC', 'MC');
% print('-depsc','fig_det_trace.eps')

% ===========
% plot ratio vs pm
% ===========
p1m = 0:0.001:1;
r = zeros(length(p1m),1);
for i = 1:length(p1m)
    [i,length(p1m)]
    r(i) = funcs(p1m(i));
end
plot(p1m, r, '-r*')