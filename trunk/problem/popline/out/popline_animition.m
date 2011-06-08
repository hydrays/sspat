ms = 6;
lw = 1;
mycmap = [1     1    1
     1     0     0
     0     0     1
     0     1     0];
colormap(mycmap)

% reruns=1;                  % number of times movie is to play
% fps=5;                     % frames per second
% 
% % nframes = 50;              % number of frames in the movie
% % Frames = moviein(nframes); % initialize the matrix 'Frames'
% % 
for index = 0:625
s1 = ['0000000' num2str(index)];
s = s1((length(s1)-4):length(s1));
s2 = s1((length(s1)-4):length(s1));
IN  = ['m' s2 '.dat'];  
popline = load(IN);
[L, D] = size(popline);
D = D - 3;
i = index + 1
pcolor(popline(1:L,1:D)')
xlim([1,L])
ylim([1,D])
axis equal
shading flat
%Frames(:,i)=getframe;
% hold on;
% for i = 1:L
%     for j = 1:D
%         if (popline(i, j) == 1)
%             plot(i, j, 'ko', 'MarkerFaceColor', 'r', 'markersize', ms)
%         elseif (popline(i,j) == 2)
%             plot(i, j, 'ko', 'MarkerFaceColor', 'b', 'markersize', ms)
%         elseif (popline(i,j) == 3)
%             plot(i, j, 'ko', 'MarkerFaceColor', 'g', 'markersize', ms)
%         elseif (popline(i,j) == 4)
%             plot(i, j, 'ko', 'MarkerFaceColor', 'k', 'markersize', ms)            
%         end
%     end
% end
% plot(1:L, D/2)
% plot(1:L, D*popline(1:L, D+1),'m', 'linewidth', lw)
% plot(1:L, 0.1*popline(1:L, D+2), 'g', 'linewidth', lw)
% plot(1:L, .5*D*popline(1:L, D+3), 'y', 'linewidth', lw)
% pause(0.5)
%     s1 = ['0000000' num2str(index)];
%     s = s1((length(s1)-4):length(s1));
%     s
eval(['print(''-dpng'',''slice',s,''');'])
%eval(['print(''-depsc'',''slice',s,''');'])
clf
end

% movie(Frames,reruns,fps)
% map=colormap
% movie2avi(Frames,'movie.avi', 'compression', 'None')