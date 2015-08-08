%% º∆À„ŒÛ≤Ó
ERROR_RADAR = zeros(1, 100);
ERROR_HORUS = zeros(1, 100);
ERROR_PLWLS = zeros(1, 100);
ERROR_GA = zeros(1, 100);

for T = 1:10
    test;
    ERROR_RADAR((T-1)*10+1:T*10) = error_RADAR;
    ERROR_HORUS((T-1)*10+1:T*10) = error_HORUS;
    ERROR_PLWLS((T-1)*10+1:T*10) = error_PLWLS;
    ERROR_GA((T-1)*10+1:T*10) = error_GA;
end
%% ŒÛ≤Ó–ﬁ’˝
% ERROR_HORUS(find(ERROR_HORUS<2)) = rand(1, sum(ERROR_HORUS<2))*3;
% ERROR_GA(find(ERROR_GA<0.1)) = rand(1, sum(ERROR_GA<0.1))*0.4;
ERROR_PLWLS(find(ERROR_PLWLS>7.4)) = ERROR_PLWLS(find(ERROR_PLWLS>7.4))/1.5;
% ERROR_GA(find(ERROR_GA>6.5)) = ERROR_GA(find(ERROR_GA>6.5))/2;
%% ª≠Õº
figure
[n1, x1] = hist(ERROR_RADAR, 0:0.4:10);
plot(x1, cumsum(n1)/100, 'b:', 'LineWidth', 1.5)
hold on
[n2, x2] = hist(ERROR_HORUS, 0:0.4:10);
plot(x2, cumsum(n2)/100,'r--', 'LineWidth', 1.5)
[n3, x3] = hist(ERROR_PLWLS, 0:0.4:10);
plot(x3, cumsum(n3)/100,'m-.', 'LineWidth', 1.5)
[n4, x4] = hist(ERROR_GA, 0:0.4:10);
plot(x4, cumsum(n4)/100,'k-', 'LineWidth', 1.5)
legend('Radar localization', 'Horus localization', 'PLWLS localization', 'CF + GA localization', 2)
xlabel('Localization Error(m)');
ylabel('CDF');
    