function [error_Indoor,err2,varargout]=ZY_SRWSpeed(Network,X,Y,ErrorName)
%This function is to caculate the error of swati core in SRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Modified by Junhong Ye on 6 Jan. 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' ']);
disp(['     Enter ZY_SRWSpeed --> Hole:',num2str(Network.time),' || Speed Step ',num2str(Network.speedStep),' || Repeat Times ',num2str(Network.time),'.']);
if nargin==0
    addpath ..\SRW ..\minFunc ..\Matlab end
    clc;
    clear all;
    close all;
    changepath;
    if 0
        path1 = mfilename('fullpath');  %get the path of the current script
        i=strfind(path1,'\');
        path1=path1(1:i(end));
        cd(path1);                  %set the working path to be the path of this script
    end
    
    Network.rx=50;
    Network.ra=5;
    asscript=1;

    Network.N=Network.rx+Network.ra;
    
    Network.NewData=1;          %switch for using new data
    Network.RefreshIniData=1;   %switch for refresh data
    Network.TracePlot=0;        %switch for node trace sample plotting
    
    Network.speed=[1,2,5,10:10:50];
    Network.minSpeed=0;         %the minimux Network.speed
    Network.speedStep=numel(Network.speed);%the number points in x-axis.
    Network.time=10;            %the total round times.
    Network.T=30;               %the snapshort.
    
    Network.length=200;
    Network.width=200;

    Network.R=50;
    err1=zeros(Network.speedStep,Network.time);
    
    Network.noise=0;                    %the percentage of the noise in the distance.
    Network.Type.Trans='TOA';
    Network.Perct.Loss=1;
    Network.Perct.TOA=1;
    
    fileName = strcat('srw_iniPos_','topology_','s');       %iniPos?srw_iniPos_topology_s.mat??????SRW_generate_iniPos.m??
    load(fileName);
    
    if Network.NewData
        for k=1:Network.speedStep
            Network.maxSpeed=Network.speed(k);
            for i=1:Network.time
                [X(k,i).data,Y(k,i).data]=SRW(Network);
            end
        end
        save('TSLRL_SRWSpeed.mat');
        disp('New data generated.');
    else
        load TSLRL_SRWSpeed.mat;
    end
else
    if  nargin==3
        disp(['"ErrorName" not assigned !']);
        ErrorName='';
    end
    asscript=0;
end
Network.N=Network.rx+Network.ra;

tic1=tic;
%% 
if Network.Perct.TOA==0
    alname='TSLRL_RF';
else
    alname='TSLRL';
end

error_Indoor = [];
for k=1:Network.speedStep
    Network.currentSpeed = Network.speed(k);
    err{k}=zeros(Network.time,Network.T);
    for i=1:Network.time
        Network.round = i;
        DISCO_Initialization;
        
        P=zeros(Network.N,2,Network.T);
        P(:,1,:)=X(k,i).data;
        P(:,2,:)=Y(k,i).data;
        
        %paramenter changed by song
%         PIni = WiFi_Horus(Network);
        PIni = WiFi_Horus(Network) + 5*(rand(Network.N, 2, Network.T) - ...
            rand(Network.N, 2, Network.T));
        PIni(PIni < 0) = 0;
        last = 30;%????
        usedSnapshots = 3;%????
        
        numOfSeg = last / usedSnapshots;
%         tempErr = [];% changed by Song
        for j = 1:numOfSeg
            PP = P(:,:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
            Network.PP = PP;
            PPIni = PIni(:,:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
            Network.PIni = PPIni;
            [errmat temp currentToc]=swatiCore(Network);
%             tempErr = [tempErr ; temp];%changed by Song
            err{k}(i,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots) = temp;%changed by Song
            error_Indoor = [error_Indoor errmat];
        end
        %   save(num2str(currentToc),'currentToc');
        showtime(['         Repeat Time: ',num2str(i),' || Speed: ',num2str(Network.speed(k)),' || Algorithm: ',alname,' || Finished at:']);       
    end
    err1(k,:)=mean(err{k}(:,end-9:end),2);
    showtime(['------------> Speed: ',num2str(Network.speed(k)),' || Algorithm: ',alname,' ||  Finished at:']);
end
err2=mean(err1,2);       %step*1 vector
varargout{1}=err;

if asscript
    plot(Network.speed,err2,'r-V');
end
speed=Network.speed;

% % try
% %     save(cattime([ErrorName,'_',alname,'_',num2str(Network.HoleNum),'_holes']),'speed','err1','err2','err');
% %     save([ErrorName,'_',alname,'_',num2str(Network.HoleNum),'_holes.mat'],'speed','err1','err2','err');
% %     disp(' ');
% %     showtime(['========== Finished ---> Hole: ',num2str(Network.HoleNum),' || Algorithm: ',alname,' || Scenario: ',ErrorName,' || Used Time: ',num2str(toc(tic1)),'s || Exitting ',alname,'_Speed at']);
% % catch
% %     disp('Fuck! Error occurs in saving files in ',alname,'_Speed.');
% % end

end
