function [ error_Centaur, MError_Centaur, MMError_Centaur ] = centaur_loc( Network,X,Y,ErrorName)
% CENTAUR ��λ�㷨

k = 6;%����HORUS��λ�����������������
% k0 = 4;%CENTAUR��λK����K��ֵ
% dist_var = 1;%����ֲ��Ĳ�������
numtrain = Network.numtrain;
numdev = Network.rx;
numap = Network.numap;
numT = Network.T;
dist = Network.dist;
rss_test = Network.rssTest;
rss_train = Network.rssTrain;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
xcoor = X.data;
ycoor = Y.data;
G_real(:, 1, :) = xcoor;
G_real(:, 2, :) = ycoor;
[ para_m, para_v ] = get_para(rss_train, numap, numtrain);
[ dist_m ] = dist_map( xtrain, ytrain, numtrain);
error_Centaur = zeros(numdev, numT);
stepLength = 9;
for i = 1:numdev/stepLength
    fprintf('CENTAUR start. step: %d\n', i);

    [~, error_Centaur((i-1)*stepLength+1:i*stepLength, :) ] = ...
        centaur_loc_core( G_real((i-1)*stepLength+1:i*stepLength, :, :), ...
        rss_test(((i-1)*stepLength+1:i*stepLength),:, :), dist(((i-1)*stepLength+1:i*stepLength),...
        ((i-1)*stepLength+1:i*stepLength), :), xtrain, ytrain, dist_m, para_m, para_v, k);
end
% G_CENTAUR = LOC_CENTAUR;
MError_Centaur = mean(error_Centaur);
MMError_Centaur = mean(MError_Centaur);

