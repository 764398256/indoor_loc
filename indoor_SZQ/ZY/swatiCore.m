function [errmat err currentToc]=swatiCore(Network)

rx = Network.rx;
ra = Network.ra;
noise1 = Network.noise;
R = Network.longestCommDis;
PP = Network.PP;
Type = Network.Type;
Perct = Network.Perct;
numSnapshots = size(PP,3);
PIni = Network.PIni;

% addpath ..\SRW ..\minFunc ..\Matlab end
% clear
% clc
% length=200;
% width=200;
% speed=10;
% minSpeed=0;
% N=55;
% T=30;
% R=50;
% Type.Trans='TOA';
% Perct.Loss=1;
% Perct.TOA=0;
% numSnapshots=30;
% rx=50;
% ra=5;
% noise=0;
% [XX,YY]=SRW(N,length,width,speed,minSpeed,T);
% save('temp11.mat','XX','YY');
% % load temp11.mat
% PP=zeros(N,2,T);
% PP(:,1,:)=XX;
% PP(:,2,:)=YY;

TC_weight = 1;        % Weight for the Temporal Stability constraint, parameter 'alpha' in the paper
LRank_weight = 0.1;   % Weight for the Low Rank constraint, parameter 'beta' in the paper 
cons.LRank = 3;       % Desired Low Rank, parameter 'r' in the paper
d=2;
X=PP(1:rx,:,:);
A=PP((rx+1):end,:,:);

[EqMatrx,UEqMatrx,MM]=dataDeal_ZY(Network);

cons.EQx = EqMatrx.Dx.^2;
cons.EQa = EqMatrx.Da.^2;
cons.UBx = UEqMatrx.UBx.^2;
cons.LBx = UEqMatrx.LBx.^2;
cons.UBa = UEqMatrx.UBa.^2;
cons.LBa = UEqMatrx.LBa.^2;

Dx= EqMatrx.Dx;
Da= EqMatrx.Da;

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

tic%%improve
dd = 2;
A_dd = zeros(ra,d+dd,numSnapshots);
% A_dd(1:ra, 1:d, :) = PIni;%changed by song
A_dd(:,1:d,:) = A;
%Localize in higher dimension in this case 4, as discussed in the paper
X3 = LocalizeLowRankSDP_3D(rx,ra,d+dd,numSnapshots,A_dd,cons,[],'off',0);

DD3 = zeros(rx+ra,rx+ra,numSnapshots);
DD=zeros(rx+ra,rx+ra,numSnapshots);
for i=1:numSnapshots
  DD3(:,:,i) = PairDist([X3(:,:,i); A(:,:,i) zeros(ra,dd)]).^2;
  DD(:,:,i)= [Dx(:,:,i),Da(:,:,i);Da(:,:,i)',PairDist(A(:,:,i))].^2;
end
DD3(DD~=0) = DD(DD~=0);
% % % DD3(MM==1) = DD(MM==1); %% Liang want to give back the measured distance, but the
% MM has been deleted here, so use this instead

%changed by song
Y3 = zeros(rx,d,numSnapshots);
Y3 = PIni;

% Use MDS to map back the distance matrix obtained into desired dimension, in this case 2
% for i=1:numSnapshots
%   Y3(:,:,i) = MDS(DD3(:,:,i),A(:,:,i));
% end

%tic
% Use above solution as initial solution for localization in desired dimension, in this case 2
% numSnapshots = T;
% A = coordinates of anchors; structure is A_dd
% cons = parameters
% Y3 is initial values
Y2 = LocalizeLowRankSDP_3D(rx,ra,d,numSnapshots,A,cons,Y3,'off',0);
currentToc = toc;
toc

% err = mean(mean(sqrt( sum((X(:,:,end-9:end) - Y2(:,:,end-9:end)).^2 , 2) )))/R;
% temp = mean(sqrt(     sum((X - Y2).^2 , 2)    ))/R;%??
% err = temp(:);
temp = sqrt(sum((X - Y2).^2 , 2));
errmat = reshape(temp,size(X,1),size(X,3));
temp = mean(errmat);
err = temp(:);
%err
