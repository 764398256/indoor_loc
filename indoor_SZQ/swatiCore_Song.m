function [errmat err currentToc]=swatiCore_Song(Network)
%This function is to caculate the error of swati core in SRW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Modified by Zeqi Song on 28 April. 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx = Network.rx;
ra = Network.ra;
noise1 = Network.noise;
longestCommDis = Network.longestCommDis;
PP = Network.PP;
Type = Network.Type;
Perct = Network.Perct;
numSnapshots = size(PP,3);
Y1 = Network.PIni;
pro = Network.pro;

TC_weight = 1;        % Weight for the Temporal Stability constraint, parameter 'alpha' in the paper
LRank_weight = 1;   % Weight for the Low Rank constraint, parameter 'beta' in the paper
cons.LRank = 3;       % Desired Low Rank, parameter 'r' in the paper
d=2;
X=PP(1:rx,:,:);
A=PP((rx+1):end,:,:);


% error = Y1 - X;
% Y1 = X + error * 0.75;


[EqMatrx,UEqMatrx,MM]=dataDeal_Song(Network);

cons.EQx = EqMatrx.Dx.^2;
cons.EWx = zeros(size(cons.EQx));
cons.EWx(cons.EQx ~= 0) = 1;
cons.EQa = EqMatrx.Da.^2;
cons.EWa = zeros(size(cons.EQa));
cons.EWa(cons.EQa ~= 0) = 1;
cons.UBx = UEqMatrx.UBx.^2;
cons.LBx = UEqMatrx.LBx.^2;
cons.UBa = UEqMatrx.UBa.^2;
cons.LBa = UEqMatrx.LBa.^2;


% Toeplitz matrix as discussed in the paper
% TC = sparse(numSnapshots-2,numSnapshots);
% for i = 1:(numSnapshots-2)
%   TC(i,i) = 1;
%   TC(i,i+1) = -2;
%   TC(i,i+2) = 1;
% end
TC=diag(ones(1,numSnapshots),0)+diag(-2*ones(1,numSnapshots-1),1)+diag(ones(1,numSnapshots-2),2);
TC=TC(1:numSnapshots-2,:);
TC=sparse(TC);

cons.TC = TC*TC_weight;       % Temporal Stability constraint
cons.gamma = LRank_weight;    % Low rank constraint

tic
% numSnapshots = T;
% A = coordinates of anchors; structure is A_dd
% cons = parameters
% Y1 is initial values

W0 = ones(rx,d,numSnapshots);
switch pro
    case 'WifiLR'
        Y2 = LocalizeLowRankSDP_3D_MICRO_enhanced(rx,ra,d,numSnapshots,A,cons,Y1,W0,'off',0);
    case {'WifiRPCA','WifiRPCAVA'}
        Y2 = LocalizeLowRankSDP_3D_RPCA(rx,ra,d,numSnapshots,A,cons,Y1,W0,'off',1,1);
end

switch pro
    case 'WifiRPCAVA'
        temp1 = abs(reshape(X - Y2,size(X,1)*size(X,2),size(X,3)));
        temp2 = abs(reshape(X - Y1,size(X,1)*size(X,2),size(X,3)));
        temp3 = abs(reshape(Y1 - Y2,size(X,1)*size(X,2),size(X,3)));
        
        %%% decide the virtual anchors
        temperr = Y1 - Y2;
%         tempdiserrten = sqrt(sum(temperr.^2 , 2));
        tempdiserrten = sqrt(sum(temperr.^2 , 2)) * 0.1 + abs(temperr(:,1,:) - temperr(:,2,:));
        tempdiserrmat = reshape(tempdiserrten,size(X,1),size(X,3));
        tempdiserrvec = sort(tempdiserrmat(:));
        tempthresh = tempdiserrvec( ceil(length(tempdiserrvec)*0.3));
        
        %%% find the entries of the virtual anchors
        temppos = zeros(size(tempdiserrten));
        temppos( tempdiserrten < tempthresh ) = 1;
        tempposd = repmat(temppos,[1,2,1]);
        
        %%% prepare the constraints
        conssec = cons;
        A = zeros(size(X));
        A(tempposd==1) = Y2(tempposd==1);
        oriA = A;
        
        %%% tricks here: let the virtual anchors are stick to themselves
        conssec.EQa = cons.EQx;
        conssec.UBa = cons.UBx;
        conssec.LBa = cons.LBx;
        for c = 1:size(X,3) % for time slots
            conssec.EQa(:,:,c) = conssec.EQa(:,:,c) + diag(temppos(:,:,c)) * 1e-7;
            conssec.UBa(:,:,c) = conssec.UBa(:,:,c) + diag(temppos(:,:,c)) * 1e-7;
            conssec.LBa(:,:,c) = conssec.LBa(:,:,c) + diag(temppos(:,:,c)) * 1e-7;
            
            conssec.EQa(:,temppos(:,:,c) ~= 1,c) = 0;
            conssec.UBa(:,temppos(:,:,c) ~= 1,c) = 0;
            conssec.LBa(:,temppos(:,:,c) ~= 1,c) = 0;
        end
        conssec.EWa = zeros(size(conssec.EQa));
        conssec.EWa(conssec.EQa ~= 0) = 1;
        
        %%% cut off the unused virtual anchors
        tempdiserrmatrow = ~any(tempdiserrmat < tempthresh,2);
        A(tempdiserrmatrow,:,:) = [];
        ra = rx - sum(tempdiserrmatrow);
        
        conssec.EQa(:,tempdiserrmatrow,:) = [];
        conssec.UBa(:,tempdiserrmatrow,:) = [];
        conssec.LBa(:,tempdiserrmatrow,:) = [];
        conssec.EWa(:,tempdiserrmatrow,:) = [];
        
        %%% copy the other constraints
        conssec.EQx = cons.EQx;
        conssec.UBx = cons.UBx;
        conssec.LBx = cons.LBx;
        conssec.EWx = cons.EWx;
        
        % Y2 = LocalizeLowRankSDP_3D_MICRO_enhanced(rx,ra,d,numSnapshots,A,conssec,Y2,W0,'off',0);
        Y2 = LocalizeLowRankSDP_3D_RPCA(rx,ra,d,numSnapshots,A,conssec,Y2,W0,'off',1,1);
        Y2(oriA~=0) = oriA(oriA~=0);
        
        % DD3 = zeros(rx,rx,numSnapshots);
        % for i=1:numSnapshots
        %   DD3(:,:,i) = PairDist(Y3(:,:,i)).^2;
        % end
        % DD3(cons.EQx~=0) = cons.EQx(cons.EQx~=0);
        % % % DD3(MM==1) = DD(MM==1); %% Liang want to give back the measured distance, but the
        % MM has been deleted here, so use this instead
        
        % Y4 = zeros(rx,d,numSnapshots);
        % % Use MDS to map back the distance matrix obtained into desired dimension, in this case 2
        % for i=1:numSnapshots
        %   Y4(:,:,i) = MDS(DD3(:,:,i),A(:,:,i));
        % end
        
        %tic
        % Use above solution as initial solution for localization in desired dimension, in this case 2
        % numSnapshots = T;
        % A = coordinates of anchors; structure is A_dd
        % cons = parameters
        % Y3 is initial values
        % Y2 = LocalizeLowRankSDP_3D(rx,ra,d,numSnapshots,A,conssec,Y3,'off',0);
end

currentToc = toc;
toc
% err = mean(mean(sqrt( sum((X(:,:,end-9:end) - Y2(:,:,end-9:end)).^2 , 2) )))/R;

temp = sqrt(sum((X - Y2).^2 , 2));
errmat = reshape(temp,size(X,1),size(X,3));
temp = mean(errmat);
err = temp(:);
%err
