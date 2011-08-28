
% ===========
% ode, three species
% ===========
% function dy = ode_onefeed(t,y)
% dy = zeros(3,1);    % a column vector
% N = 2000;
% gain = 2.2/N; %5.65 for hill function; 2 for exp1; 
% %gain2 = 3/N;
% %lambda_max = 4;
% p0 = 1/(1.01+gain*y(3));
% %p0 = exp(-gain*y(3));
% %lambda = lambda_max /(1+gain2*xp(3));
% p1 = 0.4;
% %lambda = lambda_max*exp(-gain2*y(3));
% lambda = 0.01*(sum(y)-N);
% dy(1) = p0*y(1) - (1-p0)*y(1) - lambda*y(1);
% dy(2) = 2*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
% dy(3) = -0.1*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
% end

% ===========
% ode, four species
% ===========
% function dy = ode_onefeed(t,y)
% dy = zeros(3,1);    % a column vector
% N = 2000;
% gain = 1.8/N;
% gain2 = 4/N;
% vm = 2.5;
% p0 = 1/(1.01+gain*y(3));
% v = vm/(1+gain2*y(3));
% p1 = 0.4;
% lambda = 0.01*(sum(y)-N);
% dy(1) = v*p0*y(1) - v*(1-p0)*y(1) - lambda*y(1);
% dy(2) = 2*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
% dy(3) = -0.1*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
% dy(4) = y(4) - lambda*y(4);
% end


% ===========
% ode, five species 2 feedback
% ===========
function dy = funcs(t,y)
dy = zeros(5,1);    % a column vector
N = 500;
k1 = 0.9/N;
k2 = 4/N;
vm = 2.5;
p0 = 1/(1.01+k1*(y(3)+y(5)));
v = vm/(1+k2*(y(3)+y(5)));
p1 = 0.45;
pm = 1.0;
lambda = max(0, 0.04*(sum(y)-N));
dy(1) = v*p0*y(1) - v*(1-p0)*y(1) - lambda*y(1);
dy(2) = 2*v*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
dy(3) = -0.2*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
dy(4) = pm*y(4) -(1-pm)*y(4) - lambda*y(4);
dy(5) = -0.2*y(5) + 2*(1-pm)*y(4) - lambda*y(5);
end

% % ===========
% % ode, five species 3 feedback
% % ===========
% function dy = funcs(t,y)
% dy = zeros(5,1);    % a column vector
% N = 500;
% k1 = 0.6/N;
% k2 = 4/N;
% k3 = 3/N;
% vm = 2.5;
% p0 = 1/(1.01+k1*(y(3)+y(5)));
% v = vm/(1+k2*(y(3)+y(5)));
% p1 = 0.4;
% psym = 1/(1+(k3*(y(3)+y(5)))^2);
% pm = 1.0;
% lambda = max(0, 0.04*(sum(y)-N));
% dy(1) = v*psym*p0*y(1) - v*psym*(1-p0)*y(1) - lambda*y(1);
% dy(2) = v*(1-psym)*y(1) + 2*v*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
% dy(3) = -0.2*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
% dy(4) = pm*y(4) -(1-pm)*y(4) - lambda*y(4);
% dy(5) = -0.2*y(5) + 2*(1-pm)*y(4) - lambda*y(5);
% end


% % ===========
% % function to get ratio for a ceitain pm
% % ===========
% 
% function y = funcs(pm)
% 
% function dy = ode(t,y)
% dy = zeros(5,1);    % a column vector
% N = 2000;
% gain = 0.8/N;
% gain2 = 4/N;
% vm = 2.5;
% p0 = 1/(1.01+gain*(y(3)+y(5)));
% v = vm/(1+gain2*(y(3)+y(5)));
% p1 = 0.45;
% lambda = 0.01*(sum(y)-N);
% dy(1) = v*p0*y(1) - v*(1-p0)*y(1) - lambda*y(1);
% dy(2) = 2*v*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
% dy(3) = -0.1*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
% dy(4) = pm*y(4) -(1-pm)*y(4) - lambda*y(4);
% dy(5) = -0.1*y(5) + 2*(1-pm)*y(4) - lambda*y(5);
% end
% 
% y0 = zeros(5,1);
% y0 = [200 300 600 0.00000001 0];
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5 1e-5 1e-5]);
% [T,Y] = ode45(@ode, [0 1000],y0,options);
% % plot(T,Y(:,1),'r-',T,Y(:,2),'-.b',T,Y(:,4),'--k',T,Y(:,3),':g',T,Y(:,5),'linewidth',2)
% % %xlim([600,700])
% % pause(1)
% 
% yratio = 0;
% t_temp = 0;
% for i = 1:length(Y)-1
%     if (T(i) > 200)
%         yratio = yratio + (T(i+1)-T(i))*(Y(i,4)+Y(i,5))/(sum(Y(i,:)));
%         t_temp = t_temp + (T(i+1)-T(i));
%     end
% end
% yratio = yratio / t_temp
% y = yratio;
% end

