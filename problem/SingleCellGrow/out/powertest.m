% x=csvread('xE8.csv');
% p=csvread('paraE8.csv');
% plplot(x, p(1), p(2));
% title('E8', 'fontsize', 30)
% text(10, 0.1, ['xmin=', num2str(p(1))], 'fontsize', 20)
% text(10, 0.06, ['alpha=', num2str(p(2))], 'fontsize', 20)
% print('-depsc', 'powertestE8.eps')

x=csvread('xFor02.csv');
% p=csvread('paraFor02.csv');
% plplot(x, p(1), p(2));
% title('E8', 'fontsize', 30)
% text(10, 0.1, ['xmin=', num2str(p(1))], 'fontsize', 20)
% text(10, 0.06, ['alpha=', num2str(p(2))], 'fontsize', 20)
% print('-depsc', 'powertestFor02.eps')
n = length(x);
c = [sort(x) (n:-1:1)'./n];
loglog(c(:,1),c(:,2), 'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]);
print('-depsc', 'loglogplot.eps')

figure
hist(x, 1000)
print('-depsc', 'histplot.eps')