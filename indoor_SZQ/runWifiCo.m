function [xalg,yalg]= runWifiCo(Network,i,xalg,yalg,hopDistance)
% profile on%
%
% warning off
% disp(['================ In runWifiCo ===============']);

tw = Network.tw; % 8
tr = Network.tr; % 3

rx = Network.rx;
ra = Network.ra;
noise1 = Network.noise;
R = Network.longestCommDis;
PP = Network.PP;
Type = Network.Type;
Perct = Network.Perct;
numSnapshots = size(PP,3);

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

numdev = Network.numdev;
d = 2;
R = Network.R;


if i>tr
     if i>=tw
         nt = tw;
     else 
         nt = i;
     end

     tmpx = zeros(numdev,1);
     tmpy = zeros(numdev,1);
     
     for j = 1:1  
        HopN = j;
        tmp1 = 1;
        
        for k = i-nt+1:i
            hopD = hopDistance(:,:,k);
%             hop = 1;
            hop = max(hopD(:));
            tNei = findNei(j,hopD(:,:),hop);
            nei{tmp1}.t = k;
            nei{tmp1}.hop = [j,tNei];
            HopN = union(tNei,HopN);
            tmp1 = tmp1+1;
        end
%         HopN
        index = strfind(HopN,j);
        
        rx = size(HopN,2);
        ra = 0;
        A = zeros(ra,d,nt);
        
        
        cons.EQx = zeros(rx,rx,nt);
        cons.EQa = zeros(rx,ra,nt);
        cons.UBx = zeros(rx,rx,nt);
        cons.UBa = zeros(rx,ra,nt);
        cons.LBx = zeros(rx,rx,nt);
        cons.LBa = zeros(rx,ra,nt);
        
        for k = 1:nt
            cons.EQx(:,:,k) = Network.dis(HopN,HopN,nei{k}.t);
            cons.UBx(:,:,k) = R*hopDistance(HopN,HopN,nei{k}.t);
            cons.LBx(:,:,k) = R*(hopDistance(HopN,HopN,nei{k}.t)-1);
        end
        
%         cons.UBx(cons.UBx==inf) = 0; cons.LBx(cons.LBx==inf) = 0;
        cons.LBx(cons.LBx<0) = 0;
        cons.EQx = (cons.EQx).^2;
        cons.UBx = (cons.UBx).^2;
        cons.LBx = (cons.LBx).^2;



        X0 = zeros(rx,d,nt);
        X0(:,1,:) = xalg(HopN,i-nt+1:i);
        X0(:,2,:) = yalg(HopN,i-nt+1:i);

        cons.LRank = 3;
        
%%% DEMO
%         Y = [xalg(HopN,i-ctw+1:i-1); yalg(HopN,i-ctw+1:i-1)]; X0 =
%         zeros(rx,2); X0(:,1) = xalg(HopN,i); X0(:,2) = yalg(HopN,i);
%          X3 = LocalizeLowRankSDP_3D_DEMO(rx,ra,d,nt,A,cons,Y,X0,'off',0);
%%% DEMO

        W0 = ones(rx,d,nt);
%%mark
        X3 = LocalizeLowRankSDP_3D_MICRO_enhanced(rx,ra,d,nt,A,cons,X0,W0,'off',0);
        
        tmpx = X3(:,1,end);
        tmpy = X3(:,2,end);
     end
    
    xalg(:,i) = tmpx;
    yalg(:,i) = tmpy;
end
