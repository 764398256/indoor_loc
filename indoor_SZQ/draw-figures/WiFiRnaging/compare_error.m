%% 各种方法误差对比
indoor_test;
num = sum(n(1, :));
errorIndoor= {error_Radar, error_Indoor_WifiRPCAVA, error_PLWLS,...
    error_Horus, error_Indoor_WifiLR_ALL, error_Indoor_WifiRPCA, Error_TSLRL_RF};

% % 50%误差
error50 = zeros(7, 1);
error80 = zeros(7, 1);
Mean_error = zeros(7, 1);
for i = 1:7
    error50(i) = stepVec(find_error(n(i, :), num, 0.5));
    error80(i) = stepVec(find_error(n(i, :), num, 0.8));
    Mean_error(i) = stepVec(mean(errorIndoor{i}(:)));
end

xx = 1:3;
yy = [error50  Mean_error error80];

figure 
h = bar(xx', yy');
hp = findobj(h,'type','patch'); 
hatch(hp(1),45,'k','-',4,2);
hatch(hp(2),45,'k','--',4,2);
hatch(hp(3),45,'k','-.',4,2);
hatch(hp(4),45,'k',':',4,2);
hatch(hp(5),135,'k','-',4,2);
hatch(hp(6),135,'k','--',4,2);
hatch(hp(7),135,'k','-.',4,2);

grid on;
legend('Radar', 'indoor WiFiRPCAVA', 'PLWLS',...
    'Horus', 'indoor WifiLR', 'indoor WifiRPCA', 'ZY TSLRLRF', 2);
set(gca,'xticklabel',{'50% Error',  'Mean Error', '80% Error'});
ylabel('Localization Error(m)');

set(gcf,'Position',[100,100,800,600]);

% set(gcf,'PaperPositionMode','auto');
% saveas(gcf, 'E:\金山快盘\sharebox\717zhj@163.com\WSNs_localization\our-scheme\indoor-wifi-ranging\indoor_SZQ\draw-figures\WiFiRnaging\compareError.eps')
% print(gcf,'-dpsc2','indoortest');