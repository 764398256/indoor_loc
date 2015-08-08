%% ��'Square_Speed_Trace_All.mat'�л�ȡ���ݣ����С�indoor������
%���� X,Y����8��10��struct��8����1 2 5 10 20 30 40 50��8���ٶȣ�10����10�������trace
%ÿ��trace��55��50�ľ��󣬷ֱ����55��node��device����50��ʱ��Ƭ

load('Square_Speed_Trace_All.mat')
numdev = 10;
numT = 50;
xcoor = X(4, 1).data;
xcoor = xcoor(1:10, :);
ycoor = Y(4, 1).data;
ycoor = ycoor(1:10, :);
clear X Y
% load('E:\workspace\MATLAB\wifi_rpca_dev\traces\trace1.0.mat')
Network.numdev = numdev;
Network.numT = numT;
Network.numac = 5;
Network.numtrain = 99;
X.xcoor = xcoor;
Y.ycoor = ycoor;
%����������
dis = zeros(numdev, numdev, numT);
for i = 1:numT
    dis(:, :, i) = sqrt((repmat(X.xcoor(:, i), 1, numdev) - ...
        repmat(X.xcoor(:, i), 1, numdev)').^2 + (repmat(Y.ycoor(:, i), 1, numdev)...
        - repmat(Y.ycoor(:, i), 1, numdev)').^2);
end
Network.dis = dis;
Network.G_RADAR = G_RADAR;
% %%
% load('indoortrace.mat')

%% ���С�indoor���㷨

AlsName={'MSL','MCL','TSLRL','TSLRL_RF','CODEC','CODEC_RF','DISCO1','DISCO2','DISCO_RF1','DISCO_RF2','IMCL','AENL','Orbit','Indoor'};
A = 14;
ErrorName='Square_Speed_Error'; 
Network.pro = 'WifiRpcaCo'; % Wifi,PLWLS,WifiRpca,WifiRpcaCo
diary(cattime([ErrorName,'_',AlsName{A},'_running_log'],'month2hour','.txt'));


Square_Speed_ErrorCalculate(Network,X,Y,AlsName{A},ErrorName); %run error calculating function
diary off;

%% ��ͼ
figure
[n1, x1] = hist(mean(error_RADAR), 1:0.1:6);
plot(x1, cumsum(n1)/numT,'b')
hold on
[n2, x2] = hist(MError_Indoor, 1:0.1:6);
plot(x2, cumsum(n2)/numT,'r')
legend('localization error_ini', 'localization indoor', 2) 