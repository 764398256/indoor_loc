%% 定位主函数

%% 加载数据
global numtrain numap numdev numT
load('new_trace.mat')

numdev = 10;%设置测试设备个数

get_test_data;

%% radar定位
k = 4;
[ G_RADAR, error_RADAR ] = radar_loc( xtrain, ytrain, rss_train_mean, rss_test, G_real, k );

%% horus定位
l = 4;
[ para_m, para_v ] = get_para(rss_train);
[ G_HORUS, error_HORUS ] = horus_loc( G_real, para_m, para_v, rss_test, xtrain, ytrain, l );

%% PLWLS定位
[ G_PLWLS, error_PLWLS ] = Peer_assist( dist, rss_test, rss_train_mean, G_real, G_RADAR, xtrain, ytrain );

%% GA 定位
k = 4;
[ dist_m ] = dist_map( xtrain, ytrain);
[ G_GA, error_GA ] = ga_loc( G_real, rss_test, dist, xtrain, ytrain, dist_m, para_m, para_v, k, 6);

% %% 画图
% figure
% [n1, x1] = hist(error_RADAR, 1:0.1:16);
% plot(x1, cumsum(n1)/numdev,'b*-')
% hold on
% [n2, x2] = hist(error_HORUS, 1:0.1:16);
% plot(x2, cumsum(n2)/numdev,'r+-')
% [n3, x3] = hist(error_PLWLS, 1:0.1:16);
% plot(x3, cumsum(n3)/numdev,'k.-')
% [n4, x4] = hist(error_GA, 1:0.1:16);
% plot(x4, cumsum(n4)/numdev,'yo-')
% legend('RADAR localization', 'HORUS localization', 'PLWLS localization', 'GA localization', 2)

