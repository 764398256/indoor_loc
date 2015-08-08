function [errDT,errT,error]=PLWLS_SRWSpeed(Network,X,Y,ErrorName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Author: Liang Zhang
%               Modified by Junhong Ye on 27 Feb. 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to calculate the avrage error of speed varying SRW scenarios
% with algorithm Obrit.
disp([' ']);
disp(['     Enter PLWLS_SRWSpeed --> NumAP:',num2str(Network.numac),' || Numtrain ',num2str(Network.numtrain),' || NumDevice ',num2str(Network.numdev),' || NumT ',num2str(Network.numT),'.']);

if nargin==0
    clear all
    close all
    clc
    changepath(mfilename('fullpath'));
    %´ý²¹³ä
end

%%
alname='PLWLS';
errDT=zeros(Network.numdev,Network.numT);
errT = zeros(1,Network.numT);
xcoor = X.xcoor;
ycoor = Y.ycoor;
xalg = zeros(Network.numdev,Network.numT);
yalg = zeros(Network.numdev,Network.numT);
for i=1:Network.numT
%for i=1
%     DISCO_Initialization;
     [xalg(:,i), yalg(:,i)] = runLocPLWLS(Network,i);
     %err{k}(i,:)=Network.currentError;      %time*T matrix
     showtime(['         Current Time: ',num2str(i),' || Algorithm: ',alname,' || Finished at:']);
end

LocCoor = struct('X',{xalg},'Y',{yalg});
save LocCoor.mat LocCoor;
[errDT,errT,error] =error_computer(xcoor,ycoor,xalg,yalg,Network.numdev,Network.numT);
% err1(k,:)=mean(err{k}(:,end-9:end),2);   %step*T matrix
% err2=mean(err1,2);       %step*1 vector
% varargout{1}=err;
end




