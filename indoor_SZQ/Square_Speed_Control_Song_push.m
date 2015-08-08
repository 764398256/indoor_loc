%============================For Experiment of Varying Normal Node Density==========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Author: Junhong Ye
%                   Date: 26 Mar. 2012
%                   @NEST, HUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script is used for calculating the error for all 6 algorithms in all 7 extensive topologies 
%whose holes number ranges from 1 to 16 ,with all speeds ranging from 1 to 50.
%To get your work done, please assign a number from 1 to 6 to variable "A" and run this script by pressing "F5".
%
%You can run this script on different computers which has synchronized files with MANFEN-PC, namely,NEST-2, nethu-home-2 and nethu-home-3.
%It is OK if you run with A=1 on one computer and then A=2 or any other number not exceeding 6 on another computer referred above.
%All files will be synchronized, just as if you have run all values of A on the same computer.
%
%Recommended operation: run with A=0 first and then run with other values on another computer, 
%don't repeat with the same value of A for times, it will be a waste of time.
%
%Note that: A of 1 to 6 is corresponding to MSL, MCL, TSLRL, TSLRL-RF, CODEC and CODEC-RF.
%
%Caution: Once the data is generated ("HoleVaryingScript" is run), the "HoleVaryingScript" should be annotated 
%         to avoid overwritting the trace data, unless you want to use new trace instead of the stored one.
%

close all;
% clear;
warning off all;
str=pwd;
index_dir=strfind(str,'\');
str_temp=str(1:index_dir(end)-1);
addpath(genpath(str_temp));
%changepath(mfilename('fullpath'));


ParFlag=0;              %open parallel computing

TraceName='Square_Speed_Contril_G703testbed_all.mat';
% TraceName='PLWLS.mat';
% TraceName='trace1.0.mat';    %_test
% TraceName = 'Square_Speed_Trace_All.mat';
% TraceName = 'SRW_Connected_1_Topology_Square_N_50_static_1_NumOfSpeed_1_T_50_Time_1';
ErrorName='Square_Speed_Error';    %_test
AlsName={'MSL','MCL','TSLRL','TSLRL_RF','CODEC','CODEC_RF','DISCO1','DISCO2','DISCO_RF1',...
    'DISCO_RF2','IMCL','AENL','Orbit','Indoor','PLWLS','Radar','Horus','Centaur'};

%A=1;    %corresponding to the following cell array, e.g., 1 for MSL, 2 for MCL, 3 for TSLRL etc.
%        %Choose a number and run this script by pressing F5
if ParFlag
    openparallel;
end

for A = [15]%AlsIndex
    if A==0                 %A=0 means generating new data
        Square_Speed_Script;  %generate the needed traces with size of 8*10, that is totally 80 traces/topologies,
                            %once new trace data is generated, this line should be annotated.
    else
        load(TraceName);
%         TraceName
        for i = 10
%         Network.dist = Network.dist + 0.5*(rand(Network.rx, Network.rx, Network.T)-rand(Network.rx, Network.rx, Network.T));
        Network.longestCommDis = i;
%         Network.G_RADAR = G_RADAR;
%         Network.pro = 'WifiLR'; % Wifi,PLWLS,WifiRpca,WifiRpcaCo
%         Network.pro = 'WifiRPCA'; % Wifi,PLWLS,WifiRpca,WifiRpcaCo
        Network.pro = 'WifiRPCAVA'; % Wifi,PLWLS,WifiRpca,WifiRpcaCo
%         Network.pro = 'WifiRpcaCo'; % Wifi,PLWLS,WifiRpca,WifiRpcaCo
        diary(cattime([ErrorName,'_',AlsName{A},'_running_log'],'month2hour', '.txt'));
%         Network.locType = 'DISCO';
%         Network.locType = 'DEMO';

        %Square_Speed_ErrorCalculate(Network,X,Y,AlsName{A},ErrorName); %run error calculating function
        Square_Speed_ErrorCalculate(Network,X,Y,AlsName{A},ErrorName); %run error calculating function
        diary off;
        end
    end
end

if ParFlag
    matlabpool close;
end
