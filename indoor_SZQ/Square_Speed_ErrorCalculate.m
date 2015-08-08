function Square_Speed_ErrorCalculate(Network,X,Y,Alname,ErrorName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Author: Junhong Ye
%                   Date: 26 Mar. 2012
%                   @NEST, HUST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run this function, you don't need to input actual trace data since it will load the data itself.
% The inputs of this function are: algorithm name as a string, the index of topologies and the repeated times for coresponding topologies.
% The algorithms are chosen from Alnames={'MSL','MCL','TSLRL','TSLRL-RF','CODEC','CODEC-RF'};
% The index of topologies is a vector of numbers from 1 to 7, with holes 1,2,4,6,9,12,16, e.g., HolesIndex=[1,3,7] or 1:7;
% The times for repeating is a vector of numbers from 1 to 10, which has the same length as the topology index.
% There are two optional parameters: P and ps, ps will only make sense when P=0.
% P: control the parallel mode of this function, P=1 will enable the parallel computing among holes appointed for an algorithm
% ps: control the parallel mode of inner function parallel computing, it will only be useable when P=0,
%     ps=1 will enable parallel computing between different speeds,
%     ps=0 will enable parallel computing between different times of repeated running.
%
% Example: HoleVaryingErrorCal('MSL',1:7,10*ones(1,7),0,1)

if nargin<4
    ErrorName='Square_Speed_Error';
end
ErrorV=['Error_',Alname];
MErrorV=['MError_',Alname];
OErrorV=['OError_',Alname];
disp(' ');
disp(['>>>>>>>>>>  Begin --> Algorithm: ',Alname,' || Scenario: ',ErrorName,' || Time:',time2str(' '),' <<<<<<<<<<']);
tic0=tic;
switch Alname
    case 'MSL'
        [Error_MSL,MError_MSL,OError_MSL]                          = MyMSL_SRWSpeed(Network,X,Y,ErrorName); %MyMSL
    case 'MCL'
        [Error_MCL ,MError_MCL ,OError_MCL]                        = MCL_SRWSpeed(Network,X,Y,ErrorName); %MCL
    case 'TSLRL'
        [Error_TSLRL ,MError_TSLRL ,OError_TSLRL]                  = ZY_SRWSpeed(Network,X,Y,ErrorName); %TSLRL
    case 'TSLRL_RF'
        Network.Perct.TOA=0;             %RF version of TSLRL
        [Error_TSLRL_RF ,MError_TSLRL_RF ,OError_TSLRL_RF]         = ZY_SRWSpeed(Network,X,Y,ErrorName); %TSLRL-RF
    case 'CODEC'
        [Error_CODEC ,MError_CODEC ,OError_CODEC]                  = CODEC_SRWSpeed(Network,X,Y,ErrorName); %CODEC
    case 'CODEC_RF'
        Network.Type.Loc='RF';
        [Error_CODEC_RF ,MError_CODEC_RF ,OError_CODEC_RF]         = CODEC_SRWSpeed(Network,X,Y,ErrorName); %CODEC-RF
    case 'DISCO1'
        Network.Rlm=1;
        Network.rangeFree = 0;
        %disp(['T: ',num2str(Network.T)]);
        [Error_DISCO1,MError_DISCO1,OError_DISCO1]                 = DISCO_Square_Speed(Network,X,Y,ErrorName); %DISCO1
    case 'DISCO2'
        Network.rangeFree = 0;
        [Error_DISCO2 ,MError_DISCO2,OError_DISCO2 ]               = DISCO_Square_Speed(Network,X,Y,ErrorName); %DISCO2
    case 'DISCO_RF1'
        Network.Rlm=1;
        Network.rangeFree = 1;
        [Error_DISCO_RF1 ,MError_DISCO_RF1,OError_DISCO_RF1 ]      = DISCO_Square_Speed(Network,X,Y,ErrorName); %DISCO_RF1
    case 'DISCO_RF2'
        Network.rangeFree = 1;
        [Error_DISCO_RF2 ,MError_DISCO_RF2,OError_DISCO_RF2 ]      = DISCO_Square_Speed(Network,X,Y,ErrorName); %DISCO_RF2
    case 'IMCL'
        [Error_IMCL,MError_IMCL,OError_IMCL]                       = IMCL_SRWSpeed(Network,X,Y,ErrorName);%IMCL
    case 'Orbit'
        [Error_Orbit,MError_Orbit,OError_Orbit]                    = Orbit_SRWSpeed(Network,X,Y,ErrorName); %Orbit
    case 'Indoor'
        Alname = [Alname,'_', char(Network.pro)];
        [error_Indoor,MError_Indoor,OError_Indoor]                 = Song_SRWSpeed(Network,X,Y,ErrorName);%indoor
    case 'PLWLS'
        [ error_PLWLS, MError_PLWLS, MMError_PLWLS ]               = Peer_assist( Network,X,Y,ErrorName);%PLWLS
    case 'Radar'
        [ error_Radar, MError_Radar, MMError_Radar,~ ]             = radar_loc( Network,X,Y,ErrorName );%Radar
    case 'Horus'
        [ error_Horus, MError_Horus, MMError_Horus ]               = horus_loc( Network,X,Y,ErrorName );%Horus
    case 'Centaur'
        [ error_Centaur, MError_Centaur, MMError_Centaur ]         = centaur_loc( Network,X,Y,ErrorName );%Centaur
    otherwise
        disp('There is no other algorithms, you must have made a mistake !!');
end
%preparing for including holes number in the file name

% save(cattime([ErrorName,'_',Alname,'_ALL']),ErrorV,MErrorV,OErrorV);        %store all error data of this running with hole numbers and time stamp
save([ErrorName,'_',Alname, '_ALL.mat']);
% save([ErrorName,'_',Alname, '_', num2str(Network.longestCommDis), '_ALL.mat']);
% [ErrorName,'_',Alname, '_ALL.mat']
disp(' ');
disp(['>>>>>>>>>>  Finished --->  Algorithm: ',Alname,' || Scenario: ',ErrorName,' || Used Time: ',num2str(toc(tic0)),'s || Finished at',time2str(' '),' <<<<<<<<<<']);
end
