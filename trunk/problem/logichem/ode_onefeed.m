function dy = ode_onefeed(t,y)
dy = zeros(3,1);    % a column vector
N = 2000;
gain = 1.2/N; %5.65 for hill function; 2 for exp1; 
%gain2 = 3/N;
%lambda_max = 4;
p0 = 1/(1.01+gain*y(3));
%p0 = exp(-gain*y(3));
%lambda = lambda_max /(1+gain2*xp(3));
p1 = 0.45;
%lambda = lambda_max*exp(-gain2*y(3));
lambda = 0.01*(sum(y)-N);
dy(1) = p0*y(1) - (1-p0)*y(1) - lambda*y(1);
dy(2) = 2*(1-p0)*y(1) + p1*y(2) - (1-p1)*y(2) - lambda*y(2);
dy(3) = -0.1*y(3) + 2*(1-p1)*y(2) - lambda*y(3);
end