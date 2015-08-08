load('E:\½ðÉ½¿ìÅÌ\sharebox\717zhj@163.com\WSNs_localization\our-scheme\wifi_rpca_dev\run-error-vs-speed\Square_Speed_Error_Indoor_WifiRpcaCo_ALL.mat')
MError_Indoor1 = MError_Indoor;
load('E:\workspace\MATLAB\wifi_rpca_dev\run-error-vs-speed\Square_Speed_Error_Indoor_WifiRpcaCo_ALL.mat');
figure
[n1, x1] = hist(MError_Indoor1, 1:0.3:10);
plot(x1, cumsum(n1)/30,'b*-')
hold on
[n2, x2] = hist(MError_Indoor, 1:0.3:10);
plot(x2, cumsum(n2)/30,'r+-')