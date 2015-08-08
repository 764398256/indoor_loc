load('Square_Speed_Error_Indoor_WifiRPCAVA_5_ALL.mat');
error_5 = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCAVA_6_ALL.mat');
error_6 = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCAVA_7_ALL.mat');
error_7 = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCAVA_8_ALL.mat');
error_8 = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCAVA_9_ALL.mat');
error_9 = error_Indoor;
load('Square_Speed_Error_Indoor_WifiRPCAVA_10_ALL.mat');
error_10 = error_Indoor;

%% set figure
stepVec = 0:0.3:8;
n = zeros(6, length(stepVec));

[n(1, :), ~] = hist(error_5(:), stepVec);
[n(2, :), ~] = hist(error_6(:), stepVec);
[n(3, :), ~] = hist(error_7(:), stepVec);
[n(4, :), ~] = hist(error_8(:), stepVec);
[n(5, :), ~] = hist(error_9(:), stepVec);
[n(6, :), ~] = hist(error_10(:), stepVec);

% info --setting information about the figure, containing the following fields:
%       LineStyle       --line style of the curves to be plot, a string cell with size of 1*n;
%       MarkerSize      --size of the markers of the data points;
%       LineWidth       --width of the lines to be plot;
%       MarkerFaceColor --face colors of the makers, a string cell with size of 1*n;
%       LegendText      --legend text of the legends, a string cell with size of 1*n;
%       row1num         --the number of legends on the first row;
%       xlabel          --label of the x axis;
%       ylabel          --label of the y axis;
%       title           --title of the figure;
%% call discoDrawFigure draw figure
info.plotType = 1;
info.ytick = 0:.2:1;
info.xtick = 0:1:9;
info.axis = [0 8 0 1];
info.xlabel = 'localization error (m)';
info.ylabel = 'CDF';
info.title = 'localization error for different communication distance';
info.LengendText = {'longestCommDis=5', 'longestCommDis=6', 'longestCommDis=7',...
    'longestCommDis=8', 'longestCommDis=9', 'longestCommDis=10'};
stepY = cumsum(n, 2)/Network.T/Network.N;
[handle]=discoDrawFigure(stepVec, stepY, info);

% saveas(gcf, 'E:\½ðÉ½¿ìÅÌ\sharebox\717zhj@163.com\WSNs_localization\our-scheme\indoor-wifi-ranging\indoor_SZQ\draw-figures\WiFiRnaging\testComDis.eps')
% % print(gcf,'-dpsc2','indoortest');