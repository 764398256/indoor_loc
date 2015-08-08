function [x,y] = runPLWLS(Network,i,x,y)


rsstrain = Network.rsstrain; %训练得到的样本点rss矩阵
rss = Network.rssloc(:,:,i);%i时刻的用来进行定位的rss矩阵

dist=Network.dis(:,:,i);
G_ini=[x y];
G_real = [Network.X.xcoor(:, i) Network.Y.ycoor(:, i)];
%PLWLS算法
[ G_est ] = Peer_assist( dist, rss, rsstrain, G_ini, G_real, 3, 3, 4 );
x = G_est(:,1);
y = G_est(:,2);