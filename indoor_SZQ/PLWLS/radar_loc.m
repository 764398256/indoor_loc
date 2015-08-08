function [ error_Radar, MError_Radar, MMError_Radar, G_Radar ] = radar_loc( Network,X,Y,ErrorName )
%WiFi localization
%test the rssi singal at test position and repeat 'numrep' times

numT = Network.T;
numdev = Network.rx;
k = 4;
rss_test = Network.rssTest;
rss_train_mean = Network.rssTrainMean;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
xcoor = X.data;
ycoor = Y.data;
G_real(:, 1, :) = xcoor;
G_real(:, 2, :) = ycoor;

error_Radar = zeros(numdev, numT);
G_Radar=zeros(numdev, 2, numT);

for i = 1:numT
    rss_testT = rss_test(:, :,i);
    [idx, ~] = knnsearch(rss_train_mean, rss_testT,'k',k);
    
    x = zeros(numdev,1);
    y = zeros(numdev,1);
    
    for j=1:numdev
        x(j,1) = mean(xtrain(idx(j,:)));
        y(j,1) = mean(ytrain(idx(j,:)));
    end
    
    G_Radar(:, :, i) = [x y];
    error_Radar(:, i) = sqrt(sum((G_Radar(:,:,i) - G_real(:,:,i)).^2, 2));
    
end
MError_Radar = mean(error_Radar);
MMError_Radar = mean(MError_Radar);
end

