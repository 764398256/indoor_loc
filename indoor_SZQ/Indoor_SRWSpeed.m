function [errDT,errT,error]=Indoor_SRWSpeed(Network,X,Y,ErrorName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Author: Liang Zhang
%               Modified by Junhong Ye on 27 Feb. 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to calculate the avrage error of speed varying SRW scenarios
% with algorithm Obrit.
disp([' ']);
disp(['     Enter Indoor_', char(Network.pro),'_SRWSpeed --> NumAP:',num2str(Network.numac),' || Numtrain ',num2str(Network.numtrain),' || NumDevice ',num2str(Network.numdev),' || NumT ',num2str(Network.numT),'.']);

if nargin==0
    clear all
    close all
    clc
    changepath(mfilename('fullpath'));
    %´ý²¹³ä
end

%%
alname='Indoor';
errDT=zeros(Network.numdev,Network.numT);
errT = zeros(1,Network.numT);
xcoor = X.xcoor;
ycoor = Y.ycoor;
xalg = zeros(Network.numdev,Network.numT);
yalg = zeros(Network.numdev,Network.numT);
Network.R = 1;%connect radius

connectMtx = Network.dis;
connectMtx(Network.dis<=Network.R)=1;
% connectMtx(connectMtx~=1) = 0;
% connectMtx
hopDistance=zeros(Network.numdev,Network.numdev,Network.numT);   
for s=1:Network.numT
    hopDistance(:,:,s)=graphallshortestpaths(sparse(connectMtx(:,:,s)));
%     D(:,:,s) = graphallshortestpaths(sparse(Network.dis(:,:,s)));
end
% index = Network.dis==inf;
% Network.dis(index) = 0;
% D(Network.dis~=0) = 0;
% Network.shdis = D;

%
% a = Network.dis(:, :, 4)
% b = hopDistance(:, :, 4)


for i=1:Network.numT
%for i=1
%     DISCO_Initialization;
     [xalg,yalg] = runLocIndoor(Network,i,xalg,yalg,hopDistance);
     
     [mean, median, ninety] = Terror_computer(xcoor(:,i),ycoor(:,i),xalg(:,i),yalg(:,i),Network.numdev);
     showtime(['         Current Time: ',num2str(i),'||Tcdfmean: ',num2str(mean),'||Tcdfmedian: ',num2str(median),'||Tcdfninety: ',num2str(ninety),' || Algorithm: ',alname,'_',char(Network.pro),' || Finished at:']);
end

LocCoor = [xalg;yalg];
save (['LocCoor_',char(Network.pro),'.mat'], 'LocCoor');  
%save LocCoor.mat LocCoor;
%% mark
[errDT,errT,error] =error_computer(xcoor,ycoor,xalg,yalg,Network.numdev,Network.numT); %error after rpca
% [amean,amedian,aninety] = aTerror_computer(errDT,Network.tw);
[amean,amedian,aninety] = aTerror_computer(errDT,8);
showtime(['         NumT: ',num2str(Network.numT),'||cdfmean: ',num2str(amean),'||cdfmedian: ',num2str(amedian),'||cdfninety: ',num2str(aninety),' || Algorithm: ',alname,'_',char(Network.pro),' || Finished at:']);
end




