load switch1
[n1, x] = hist(switch1(:,1), 100);
plot(x, n1, 'g');
hold on;

load switch2
[n2, x] = hist(switch2(:,1), 100);
plot(x, n2, 'b');

load switch3
[n3, x] = hist(switch3(:,1), 100);
plot(x, n3, 'r');

xlim([1000, 1100]);

legend('s1', 's2', 's3')
print('-depsc','fig_switch_dist.eps')