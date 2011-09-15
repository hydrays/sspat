load tr200_2
tr = tr200_2;
figure
plot(tr(:,1), tr(:, 2),'r')
hold on
plot(tr(:,1), tr(:, 3),'b')
plot(tr(:,1), tr(:, 4),'g')