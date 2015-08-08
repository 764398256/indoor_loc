%% indoor our scheme
Square_Speed_Control;
load('Square_Speed_Error_Indoor_WifiRpcaCo_ALL.mat')

%% wifi localization
rsstrain = Network.rsstrain; %训练得到的样本点rss矩阵
numdev = Network.numdev;
numT = Network.numT;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
G_WiFi = zeros(numdev, 2, numT);
Error_WiFi = zeros(numdev, numT);
xcoor = Network.X.xcoor;
ycoor = Network.Y.ycoor;


k = 4;

for i = 1:numT

    rssloc = Network.rssloc(:,:,i);%i时刻的用来进行定位的rss矩阵
    [idx, dist] = knnsearch(rsstrain, rssloc,'k',k);
    for j=1:numdev
        G_WiFi(j, 1, i) = mean(xtrain(idx(j,:)));
        G_WiFi(j, 2, i) = mean(ytrain(idx(j,:)));
    end
    Error_WiFi(:, i) = sqrt((G_WiFi(:, 1, i) - xcoor(:, 1, i)).^2 + ...
            (G_WiFi(:, 2, i) - ycoor(:, 1, i)).^2);
    
end

%% 画图
figure
[n1, x1] = hist(mean(Error_WiFi), 1:0.1:6);
plot(x1, cumsum(n1)/numT,'b')
hold on
[n2, x2] = hist(MError_Indoor, 1:0.1:6);
plot(x2, cumsum(n2)/numT,'r')
legend('localization error WiFi', 'localization error ourscheme', 2) 

