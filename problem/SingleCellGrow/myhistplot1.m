x = csvread('x1.csv');
[n, z] = hist(x, 10);
n = n./sum(n);
plot(z, n)
hold on

x = csvread('x2.csv');
[n, z] = hist(x, 10);
n = n./sum(n);
plot(z, n, 'c')
hold on

x = csvread('x3.csv');
[n, z] = hist(x, 10);
n = n./sum(n);
plot(z, n, 'r')
hold on

% 
% x = csvread('cellsum1.csv');
% [n, z] = hist(x, 10);
% n = n./sum(n);
% plot(z, n,'r')
% hold on
% 
% x = csvread('cellsum2.csv');
% [n, z] = hist(x, 10);
% n = n./sum(n);
% plot(z, n,'c')
% % hold on
% 
% x = csvread('cellsum3.csv');
% [n, z] = hist(x, 10);
% n = n./sum(n);
% plot(z, n,'k')
% hold on

% x = csvread('cellsuma2.csv');
% [n, z] = hist(x, 10);
% n = n./sum(n);
% plot(z, n, 'c')
% hold on
% 
% 
% x = csvread('cellsuma.csv');
% [n, z] = hist(x, 10);
% n = n./sum(n);
% plot(z, n, 'r')
% hold on
