function [x,y] = runPLWLS(Network,i,x,y)


rsstrain = Network.rsstrain; %ѵ���õ���������rss����
rss = Network.rssloc(:,:,i);%iʱ�̵��������ж�λ��rss����

dist=Network.dis(:,:,i);
G_ini=[x y];
G_real = [Network.X.xcoor(:, i) Network.Y.ycoor(:, i)];
%PLWLS�㷨
[ G_est ] = Peer_assist( dist, rss, rsstrain, G_ini, G_real, 3, 3, 4 );
x = G_est(:,1);
y = G_est(:,2);