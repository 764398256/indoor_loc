function [x,y] = runLocPLWLS(Network, i)

numtrain = Network.numtrain;
numac = Network.numac;
numdev = Network.numdev;
xtrain = Network.xtrain;
ytrain = Network.ytrain;

rsstrain = Network.rsstrain; %训练得到的样本点rss矩阵
rss = Network.rssloc(:,:,i);%i时刻的用来进行定位的rss矩阵

k = 4;
[idx, dist] = knnsearch(rsstrain, rss,'k',k);

x = zeros(numdev,1);
y = zeros(numdev,1);

for j=1:numdev
    x(j,1) = mean(xtrain(idx(j,:)));
    y(j,1) = mean(ytrain(idx(j,:)));
end

%x,y为根据wifi利用knn得到的定位坐标结果，可在这里进行进一步处理，添加函数，提高定位精度。
%....
%.....
dist=Network.dis(:,:,i);
G_ini=[x y];
G_real = [Network.X.xcoor(:, i) Network.Y.ycoor(:, i)];
%PLWLS算法
[ G_est ] = Peer_assist( dist, rss, rsstrain, G_ini, G_real, 3, 3, 4 );
x = G_est(:,1);
y = G_est(:,2);
