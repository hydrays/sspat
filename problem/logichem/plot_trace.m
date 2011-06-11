load tr
tr = tr;
plot(tr(:,1), tr(:,2), 'r')
xlim([350, 550])
hold on
plot(tr(:,1), tr(:,3), 'm')
plot(tr(:,1), tr(:,4), 'g')
plot(tr(:,1), tr(:,5), 'k')
plot(tr(:,1), tr(:,8), 'b')
legend('SC', 'TAC', 'TDC', 'MC', 'Total')
