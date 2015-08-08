function [err1,err2,varargout]=ZY_VR(Network,X,Y,ErrorName)
%This function is to caculate the error of swati core in SRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Modified by Junhong Ye on 6 Jan. 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([' ']);
disp(['     Enter ZY_SRWSpeed --> Hole:',num2str(Network.HoleNum),' || Speed Step ',num2str(Network.speedStep),' || Repeat Times ',num2str(Network.time),'.']);
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
        ErrorName='Temp';
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

%save disk space
allRatiosInfo = Network.allRatiosInfo;
Network.allRatiosInfo = [];
dirName = num2str(Network.currentSpeed);
mkdir(dirName);
for k=1:length(Network.usedRatios)
    Network.currentRatio = Network.usedRatios(k);
    err{k}=zeros(Network.time,Network.T);
    for i=1:Network.time
        resultFileName = [dirName,'/k-',num2str(k),'-times-',num2str(i),'.mat'];
        if exist(resultFileName,'file')
            load(resultFileName);
            Network.currentError = currentErr;
        else
            
            Network.round = i;
            DISCO_Initialization;
            
            P=zeros(Network.N,2,Network.T);
            P(:,1,:)=X(Network.currentSpeed==Network.speed,i).data;
            P(:,2,:)=Y(Network.currentSpeed==Network.speed,i).data;
            curMeaDisEnt = allRatiosInfo(Network.currentSpeed==Network.speed,i).meaDisEnt;
            curNoiseCon = allRatiosInfo(Network.currentSpeed==Network.speed,i).noiseCon;
            curNoiseDis = allRatiosInfo(Network.currentSpeed==Network.speed,i).noiseDis;
            curMeaDis = allRatiosInfo(Network.currentSpeed==Network.speed,i).meaDis;
            
            last = 15;
            usedSnapshots = 15;
            numOfSeg = last / usedSnapshots;
            tempErr = [];
            for j = 1:numOfSeg
                Network.PP = P(:,:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
                Network.meaDisEnt = curMeaDisEnt(:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
                Network.noiseCon = curNoiseCon(:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
                Network.noiseDis = curNoiseDis(:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
                Network.meaDis = curMeaDis(:,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots);
                [temp currentToc]= swatiCore(Network);
                %             [temp currentToc]=swatiCore(Network.rx,Network.ra,usedSnapshots,Network.noise,Network.R,PP,Network.Type,Network.Perct);
                tempErr = [tempErr ; temp];
                err{k}(i,end-j*usedSnapshots+1:end-(j-1)*usedSnapshots) = tempErr;
            end
            
        end
    end
    err1(k,:)=mean(err{k}(:,end-9:end),2);
    showtime(['         Repeat Time: ',num2str(i),' ||Ratio: ',...
        num2str(Network.currentRatio),' || Algorithm: ',alname,' || Finished at:']);
end
err2=mean(err1,2);       %step*1 vector
varargout{1}=err;

% % try
% %     save(cattime([ErrorName,'_',alname,'_',num2str(Network.HoleNum),'_holes']),'speed','err1','err2','err');
% %     save([ErrorName,'_',alname,'_',num2str(Network.HoleNum),'_holes.mat'],'speed','err1','err2','err');
% %     disp(' ');
% %     showtime(['========== Finished ---> Hole: ',num2str(Network.HoleNum),' || Algorithm: ',alname,' || Scenario: ',ErrorName,' || Used Time: ',num2str(toc(tic1)),'s || Exitting ',alname,'_Speed at']);
% % catch
% %     disp('Fuck! Error occurs in saving files in ',alname,'_Speed.');
% % end
