function [xalg,yalg] = runLocIndoor(Network,i,xalg,yalg,hopDistance)
%%%%%

numtrain = Network.numtrain;
numac = Network.numac;
numdev = Network.numdev;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
pro = char(Network.pro);
X = Network.X;
Y = Network.Y;

%wifi
rsstrain = Network.rsstrain; %训练得到的样本点rss矩阵
rssloc = Network.rssloc(:,:,i);%i时刻的用来进行定位的rss矩阵

k = 5;
[idx, dist] = knnsearch(rsstrain, rssloc,'k',k);


for j=1:numdev
    xalg(j,i) = mean(xtrain(idx(j,:)));
    yalg(j,i) = mean(ytrain(idx(j,:)));
end

%wifi end

switch pro
    case 'Wifi'
        
    case 'PLWLS'
        [xalg(:,i),yalg(:,i)] = runPLWLS(Network,i,xalg(:,i),yalg(:,i));
    case 'WifiRpca'
        Network.tw = 8;
        Network.tr = 3;
        [xalg,yalg] = runWifiRpca(Network,i,xalg,yalg);
    case 'WifiRpcaCo'
        %%%WifiRpca + or not
        Network.tw = 8;
        Network.tr = 3;
        %%%[xalg,yalg] = runWifiRpca(Network,i,xalg,yalg);
        [xalg,yalg] = runWifiCo(Network,i,xalg,yalg,hopDistance);
    otherwise
        disp('There is no this algorithm, you must have made a mistake !!');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end