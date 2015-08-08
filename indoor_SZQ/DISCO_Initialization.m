%DISCO initialization
%%%%%%%%%%%%%%%%%%%%%%%%%Network settings%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(Network,'scenario')
    Network.scenario='Square_Speed_All';
end

if ~isfield(Network,'isStatic')
    Network.isStatic = 0;
end
%%%%%%%%%%%%%%%if using static topology,just run once at getTop and iniLoc
if Network.isStatic
    Network.runTime = 1;
else
    Network.runTime = Network.T;
end

%%%%%%%%%%%%%%%%%%%%%%%%parameters for MCL MSL
if ~isfield(Network,'NS')
    Network.NS = 50;
end

%%%%%%%%%%%%%%%%%%%%%%%%%parameters for DISCO CODEC DEMO MICRO CELL
if ~isfield(Network,'longestCommDis')
    Network.longestCommDis = 8;
end

if ~isfield(Network,'measureType')
    Network.measureType = 'RSS';
end
if ~isfield(Network,'Rlm')
    Network.Rlm = 2;
end
if ~isfield(Network,'rangeFree')
    Network.rangeFree = 0;
end
if ~isfield(Network,'debug')
   Network.debug=1; 
end

if ~isfield(Network,'cons_weight')
   Network.cons_weight = [1 0.1 0.01]; 
end
if ~isfield(Network,'stru')
   Network.stru=1; 
end
if ~isfield(Network,'d')
   Network.d=2; 
end
if ~isfield(Network,'disT')
   Network.disT=15; 
end
if ~isfield(Network,'winLen')
   Network.winLen=15; 
end
if ~isfield(Network,'wd')
   Network.wd=1; 
end
if ~isfield(Network,'cons_rank')
   Network.cons_rank = 3;
end
if ~isfield(Network,'iter')
   Network.iter=20; 
end
if ~isfield(Network,'iter_alt_fit')
   Network.iter_alt_fit=20;
end

%%%%%%%%%%% support multi algorithm combinition
if ~isfield(Network,'iniType')
    Network.iniType = 2;
end
if ~isfield(Network,'locType')
    Network.locType = 'DISCO';
end
if ~isfield(Network,'errorType')
    Network.errorType = 1;
end
if ~isfield(Network,'coorOfFla')
    Network.coorOfFla = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%percent attributes
if ~isfield(Network,'ratioType')
   Network.ratioType = 'NONE';  %% errors in distance
end
if ~isfield(Network,'currentRatio')
   Network.currentRatio = 1;  %% errors in distance
end
if ~isfield(Network,'mesDisDir')
    Network.mesDisDir = [];
end
if ~isfield(Network,'disMeaDir')
    Network.disMeaDir = [];
end


%here MRatio may be: the percent of nodes with disntance measurement
%capacity; the ratio of noise; the probability of loss of distance
%measurements
Network.dirName = strcat(Network.scenario,'_static_',num2str(Network.isStatic),...
    '_N_',num2str(Network.N),'_T_',num2str(Network.T));

str1 = [Network.dirName,'/'];

str=strcat('speed_',num2str(Network.currentSpeed),'_random_',num2str(Network.round),...
    '_MRatio_',num2str(Network.currentRatio),'_Rlm_',num2str(Network.Rlm),...
    '_coorOfFla_',num2str(Network.coorOfFla),'_RF_',num2str(Network.rangeFree),...
    '_MType_',Network.measureType);

Network.locFileName = strcat(str1,str,'_loc.mat');

Network.baseFileName = strcat(str1,str,'_base.mat');

Network.topFileName = strcat(str1,str,'_top.mat');

Network.iniFileName = strcat(str1,str,'_ini.mat');

Network.resFileName = strcat(str1,str,'_res.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
