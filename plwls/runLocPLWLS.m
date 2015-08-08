function [x,y] = runLocPLWLS(Network, i)

numtrain = Network.numtrain;
numac = Network.numac;
numdev = Network.numdev;
xtrain = Network.xtrain;
ytrain = Network.ytrain;

rsstrain = Network.rsstrain; %ѵ���õ���������rss����
rss = Network.rssloc(:,:,i);%iʱ�̵��������ж�λ��rss����

k = 4;
[idx, dist] = knnsearch(rsstrain, rss,'k',k);

x = zeros(numdev,1);
y = zeros(numdev,1);

for j=1:numdev
    x(j,1) = mean(xtrain(idx(j,:)));
    y(j,1) = mean(ytrain(idx(j,:)));
end

%x,yΪ����wifi����knn�õ��Ķ�λ������������������н�һ��������Ӻ�������߶�λ���ȡ�
%....
%.....
dist=Network.dis(:,:,i);
G_ini=[x y];
G_real = [Network.X.xcoor(:, i) Network.Y.ycoor(:, i)];
%PLWLS�㷨
[ G_est ] = Peer_assist( dist, rss, rsstrain, G_ini, G_real, 3, 3, 4 );
x = G_est(:,1);
y = G_est(:,2);
