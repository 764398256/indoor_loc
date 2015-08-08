%% 运行‘indoor’算法
close all
clc

% Square_Speed_Control_Song;
%% load data
load('Square_Speed_Error_Radar_ALL.mat')
load('Square_Speed_Error_PLWLS_ALL.mat')
load('Square_Speed_Error_Horus_ALL.mat')
load('Square_Speed_Error_Indoor_WifiRPCAVA_ALL.mat')
load('Square_Speed_Error_Centaur_ALL.mat')
error_Indoor_WifiRPCAVA = error_Indoor;
load('Square_Speed_Error_Indoor_WifiLR_ALL.mat')
error_Indoor_WifiLR_ALL = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCA_ALL.mat')
error_Indoor_WifiRPCA = error_Indoor;
load('Square_Speed_Error_TSLRL_RF_ALL.mat')

%% set figure
stepVec = 0:0.5:8;
x = repmat(stepVec,8,1);
n = zeros(8, length(stepVec));

[n(1, :), ~] = hist(error_Radar(:), stepVec);
[n(2, :), ~] = hist(error_Indoor_WifiRPCAVA(:), stepVec);
[n(3, :), ~] = hist(error_PLWLS(:), stepVec);
[n(4, :), ~] = hist(error_Horus(:), stepVec);
[n(5, :), ~] = hist(error_Indoor_WifiLR_ALL(:), stepVec);
[n(6, :), ~] = hist(error_Indoor_WifiRPCA(:), stepVec);
[n(7, :), ~] = hist(Error_TSLRL_RF(:), stepVec);
[n(8, :), ~] = hist(error_Centaur(:), stepVec);


%% 画图
figure
plot(x(1, :), cumsum(n(1, :))/Network.T/Network.N,'bx-', 'MarkerSize', 5, 'LineWidth', 1)
hold on
plot(x(2, :), cumsum(n(2, :))/Network.T/Network.N,'r.-', 'MarkerSize', 3, 'LineWidth', 1)
hold on
plot(x(3, :), cumsum(n(3, :))/Network.T/Network.N,'ko-', 'MarkerSize', 3, 'LineWidth', 1)
hold on
plot(x(4, :), cumsum(n(4, :))/Network.T/Network.N,'y*-', 'MarkerSize', 3, 'LineWidth', 1)
hold on
plot(x(5, :), cumsum(n(5, :))/Network.T/Network.N,'m>-', 'MarkerSize', 3, 'LineWidth', 1)
% [n(6, :), x(6, :)] = hist(error_Indoor_WifiRPCA(:), 0:0.3:8);
plot(x(6, :), cumsum(n(6, :))/Network.T/Network.N,'g<-', 'MarkerSize', 3, 'LineWidth', 1)
hold on
plot(x(7, :), cumsum(n(7, :))/Network.T/Network.N,'cs-', 'MarkerSize', 3, 'LineWidth', 1)
hold on
plot(x(8, :), cumsum(n(8, :))/Network.T/Network.N,'kx-', 'MarkerSize', 5, 'LineWidth', 1)
hold on

set(gcf,'Position',[100,100,600,450]);

legend('Radar', 'indoor WiFiRPCAVA', 'PLWLS',...
    'Horus', 'indoor WifiLR', 'indoor WifiRPCA', 'ZY TSLRLRF', 'Centaur', 4);

% box off
% ax(2, :) = axes('Position',get(gca,'Position'),...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none',...
%     'XColor','k','YColor','k');
%
% set(ax2,'YTick', []);
% set(ax2,'XTick', []);
% box on
% set(gca,'xtick',0:0.3:8);
ylabel('CDF');
xlabel('localization error (m)');
% set(gcf,'PaperPositionMode','auto');
% saveas(gcf, 'E:\金山快盘\sharebox\717zhj@163.com\WSNs_localization\our-scheme\indoor-wifi-ranging\indoor_SZQ\draw-figures\WiFiRnaging\localizationError.eps')
% print(gcf,'-dpsc2','indoortest');
