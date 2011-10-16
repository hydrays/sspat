ms = 6;
lw = 1;
mycmap = [1     1    1
     1     0     0
     0     0     1
     0     1     0];
colormap(mycmap)

j = 1;
for index =0:1:100
clf
s1 = ['0000000' num2str(index)];
s = s1((length(s1)-4):length(s1));
s2 = s1((length(s1)-4):length(s1));
IN  = ['m' s2 '.dat'];
gIN  = ['g' s2 '.dat']; 
popline = load(IN);
%popgene = load(gIN);
[L, D] = size(popline);
D = D - 1;
i = index + 1
hold on
%box on;
% xlim([1,L])
%ylim([0,100])

pcolor(popline(1:L,1:200)')
axis equal
shading flat
%plot(popgene(:,1), 0.5*D*(popgene(:,2)-1),'-+')

plot(1:L, 0.1*popline(1:L, D+1),'m', 'linewidth', lw)
% plot(1:L, 0.1*popline(1:L, D+2), 'g', 'linewidth', lw)
% plot(1:L, popline(1:L, D+1), 'm', 'linewidth', lw)
% s1 = ['0000000' num2str(index)];
% s = s1((length(s1)-4):length(s1));
% s = s;
% eval(['print(''-dpng'',''slice',s,''');'])
% eval(['print(''-depsc'',''select',s,''');'])
% genemean(j) = mean(popgene(:,2));
pause(0.5)
j = j + 1;
end

% plot(5*(1:length(genemean)-1), genemean(2:length(genemean)))
% xlabel('t','fontsize', 18);
% ylabel('mean','fontsize',18)
