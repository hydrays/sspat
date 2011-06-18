ms = 6;
lw = 1;
mycmap = [1     1    1
     1     0     0
     0     0     1
     0     1     0];
colormap(mycmap)

for index = 700:715
clf
s1 = ['0000000' num2str(index)];
s = s1((length(s1)-4):length(s1));
s2 = s1((length(s1)-4):length(s1));
IN  = ['m' s2 '.dat'];
gIN  = ['g' s2 '.dat']; 
popline = load(IN);
popgene = load(gIN);
[L, D] = size(popline);
D = D - 3;
i = index + 1
pcolor(popline(1:L,1:D)')
xlim([1,L])
%ylim([1,D])
axis equal
shading flat
hold on
plot(popgene(:,1), 0.5*D*(popgene(:,2)-1),'-+')
% plot(1:L, D*popline(1:L, D+1),'m', 'linewidth', lw)
% plot(1:L, 0.1*popline(1:L, D+2), 'g', 'linewidth', lw)
% plot(1:L, .5*D*popline(1:L, D+3), 'y', 'linewidth', lw)
s1 = ['0000000' num2str(index)];
s = s1((length(s1)-4):length(s1));
s
%eval(['print(''-dpng'',''slice',s,''');'])
%eval(['print(''-depsc'',''slice',s,''');'])
pause(0.5)
end