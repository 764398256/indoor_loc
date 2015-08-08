%插值构造训练点接收信号强度
numtrain = Network.numtrain;
numdev = Network.rx;
numT = Network.T;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
rssTrainMean = Network.rssTrainMean;
numap = Network.numap;
k = 4;

distAll = zeros(numdev, numT, numtrain);
for i = 1:numtrain
    distAll(:, :, i) = sqrt((X.data - xtrain(i)).^2 + (Y.data - ytrain(i)).^2);
end
[distSorted, idx] = sort(distAll, 3);

rssTest = zeros(numdev, numap, numT);

for i = 1:numdev
    for j = 1:numT
%         a = rssTrainMean(reshape(idx(i, j, 1:k), 1, k), :)
%         repmat(reshape(distSorted(i, j, 1:k), 1, k), numap, 1)'
        rssTest(i, :, j) = sum(rssTrainMean(reshape(idx(i, j, 1:k), 1, k),...
            :)./repmat(reshape(distSorted(i, j, 1:k), 1, k), numap, 1)')/sum(1/distSorted(i, j, 1:k)); 
    end
end