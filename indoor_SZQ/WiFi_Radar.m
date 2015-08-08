function [PIni] =WiFi_Radar(Network)

numtrain = Network.numtrain;
numac = Network.numap;
numdev = Network.N;
SnapShots = Network.T;
xtrain = Network.xtrain;
ytrain = Network.ytrain;

%wifi
rssTrainMean = Network.rssTrainMean; %ѵ���õ���������rss����
xalg = zeros(numdev, SnapShots);
yalg = zeros(numdev, SnapShots);
for t = 1:SnapShots
    rssTest = Network.rssTest(:,:,t);%iʱ�̵��������ж�λ��rss����
    
    k = 4;
    [idx, dist] = knnsearch(rssTrainMean, rssTest,'k',k);
    
    
    for j=1:numdev
        xalg(j,t) = mean(xtrain(idx(j,:)));
        yalg(j,t) = mean(ytrain(idx(j,:)));
    end
end
PIni(:, 1, :) = xalg;
PIni(:, 2, :) = yalg;