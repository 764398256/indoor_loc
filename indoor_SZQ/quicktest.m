% Here is a simple script to illustrate how to make the localization calls

% Below packages will need to be installed for the optimization code to work
addpath /u/lili/MATLAB/SDPLR-1.03-beta /u/lili/MATLAB/SeDuMi_1_21 /u/lili/MATLAB/matlab /u/lili/MATLAB/minFunc -end

rand('seed',sum(100*clock));

clear all;

rx = 50;              % # sensors
ra = ceil(sqrt(rx));  % # anchors 
d  = 2;               % # dimensions
R = 0.2^2;            % communication range^2
numSnapshots = 16;    % no of snapshots to be considered together for localization 
TC_weight = 1;        % Weight for the Temporal Stability constraint, parameter 'alpha' in the paper
LRank_weight = 0.1;   % Weight for the Low Rank constraint, parameter 'beta' in the paper 
cons.LRank = 3;       % Desired Low Rank, parameter 'r' in the paper
noise = 0.10;         % noise can varied, right now it is set to 10%

% Below the regular and anchor node locations across snapshots 
% are initialized, just for a demo of how to use our localization 
% code. Note that such random locations are not realistic in any 
% mobility model, not even random waypoint. For real evaluation
% we need to feed in X and A by generating locations over time
% using some mobility model we choose to evaluate.

% Randomly initializing all regular node locations 
X = rand(rx,d,numSnapshots)-0.5;

% Randomly initializing all anchor node locations
A = rand(ra,d,numSnapshots)-0.5;

% normalization: scale everything such that R = 1
scale = sqrt(R);
R = R/scale^2;
X = X/scale;
A = A/scale;

PP = zeros(rx+ra,d,numSnapshots);
PP(1:rx,:,:) = X;
PP((rx+1):end,:,:) = A;

DD = zeros(rx+ra,rx+ra,numSnapshots);
for i=1:numSnapshots
  DD(:,:,i) = PairDist(PP(:,:,i)).^2;
end
Dx = DD(1:rx,1:rx,:);
Da = DD(1:rx,rx+1:rx+ra,:);

% Adding noise to distance matrix from normal distribution, for the evaluation, 
% as discussed in paper, we add noise to communication range and perturb distances
% accordingly. To set upper bounds and lower bounds we use noisy communication range

Dx_noisy = Dx.*max(eps,1+noise*randn(rx,rx,numSnapshots)).^2;
Da_noisy = Da.*max(eps,1+noise*randn(rx,ra,numSnapshots)).^2;

for i=1:numSnapshots
  Dx_noisy(:,:,i) = (Dx_noisy(:,:,i)+Dx_noisy(:,:,i)')/2;
  DD_noisy(:,:,i) = [Dx_noisy(:,:,i) Da_noisy(:,:,i); Da_noisy(:,:,i)' PairDist(A(:,:,i)).^2];
end

% Keep only distance values >= communication range R
MM = 1.0*(DD<R);
MM((rx+1):end,(rx+1):end,:) = 1;
Mx = 1.0*(Dx<R);
Ma = 1.0*(Da<R);

non_missing_frac = nnz(MM)/prod(size(MM))

% Below we set equality constraints. We do not show how to set 
% upper bounds and lower bounds (cons.UBx, cons.UBa, cons.LBx, cons.LBa)
% in this simple wrapper, but they can be set as discussed in the paper.
cons.EQx = zeros(rx,rx,numSnapshots);
cons.EQa = zeros(rx,ra,numSnapshots);
cons.EQx(Mx==1) = Dx_noisy(Mx==1);
cons.EQa(Ma==1) = Da_noisy(Ma==1);

% Toeplitz matrix as discussed in the paper
TC = sparse(numSnapshots-2,numSnapshots);
for i = 1:(numSnapshots-2)
  TC(i,i) = 1;
  TC(i,i+1) = -2;
  TC(i,i+2) = 1;
end

cons.TC = TC*TC_weight;       % Temporal Stability constraint
cons.gamma = LRank_weight;    % Low rank constraint

tic
dd = 2;
A_dd = zeros(ra,d+dd,numSnapshots);
A_dd(:,1:d,:) = A;
% Localize in higher dimension in this case 4, as discussed in the paper
X3 = LocalizeLowRankSDP_3D(rx,ra,d+dd,numSnapshots,A_dd,cons,[],'off',0);
toc
DD3 = zeros(rx+ra,rx+ra,numSnapshots);
for i=1:numSnapshots
  DD3(:,:,i) = PairDist([X3(:,:,i); A(:,:,i) zeros(ra,dd)]).^2;
end
DD3(MM==1) = DD_noisy(MM==1);
Y3 = zeros(rx,d,numSnapshots);

% Use MDS to map back the distance matrix obtained into desired dimension, in this case 2
for i=1:numSnapshots
  Y3(:,:,i) = MDS(DD3(:,:,i),A(:,:,i));
end
DD3 = zeros(rx+ra,rx+ra,numSnapshots);
for i=1:numSnapshots
  DD3(:,:,i) = PairDist([Y3(:,:,i); A(:,:,i)]).^2;
end

tic
% Use above solution as initial solution for localization in desired dimension, in this case 2
Y2 = LocalizeLowRankSDP_3D(rx,ra,d,numSnapshots,A,cons,Y3,'off',0);
toc
DD2 = zeros(rx+ra,rx+ra,numSnapshots);
for i=1:numSnapshots
  DD2(:,:,i) = PairDist([Y2(:,:,i); A(:,:,i)]).^2;
end

% Calculation of mean absolute error, X: Original known coordinates, Y2: Estimated coordinates
err = mean(mean(sqrt(sum((X - Y2).^2,2))))/sqrt(R)

quit
