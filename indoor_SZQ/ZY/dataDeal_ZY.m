function [EqMatrx,UEqMatrx,MM]=dataDeal_ZY(Network)

rx = Network.rx;
ra = Network.ra;
noise1 = Network.noise;
R = Network.longestCommDis;
PP = Network.PP;
Type = Network.Type;
Perct = Network.Perct;
numSnapshots = size(PP,3);
longestCommDis = R;

% % % % the function is to transform the input data to the form of the input of
% % % % ZY's core.
% % % TransType=Type.Trans;
TOAPct=Perct.TOA;
MM = [];
% % % LossPct=Perct.Loss;
% % % 
% % % %to construct the miss matrix according to the loss rate.
% % % %1:the distance is not missed
% % % %0:the distance is missed.
% % % tempRand=rand(ra+rx,ra+rx,numSnapshots);
% % % tempMiss=zeros(ra+rx,ra+rx,numSnapshots)+1;
% % % tempMiss(tempRand>LossPct)=0;
% % % for i=1:numSnapshots
% % %     tempMiss1=tempMiss(:,:,i)+tempMiss(:,:,i)';
% % %     tempMiss1(tempMiss1==1)=2;
% % %     tempMiss1=tempMiss1/2;
% % %     tempMiss(:,:,i)=tempMiss1;
% % %     clear tempMiss1;
% % % end
% % % 
% % % %add the miss data from the TOA or RSS.
% % % nodeSeq=randperm(rx+ra);
% % % SelctNum=floor((rx+ra)*TOAPct);
% % % SelctBoor=zeros(rx+ra,1);
% % % SelctBoor(nodeSeq(1:SelctNum))=1;
% % % TOASelctNode=find(SelctBoor>0);
% % % if strcmp(TransType,'TOA')
% % %     tempTOA=zeros(rx+ra,rx+ra,numSnapshots);
% % %     tempTOA(TOASelctNode,TOASelctNode,:)=1;
% % %     tempTOA=1-tempTOA;
% % %     tempMiss1=tempMiss.*tempTOA;
% % %     tempMiss(tempMiss1>0)=0;
% % %     tempMiss((rx+1):end,(rx+1):end,:)=1;
% % % else if strcmp(TransType,'RSS')
% % %         tempTOA=zeros(rx+ra,rx+ra,numSnapshots);
% % %         tempTOA(TOASelctNode,:,:)=1;
% % %         tempTOA(:,TOASelctNode,:)=1;
% % %         tempTOA=1-tempTOA;
% % %         tempMiss1=tempMiss.*tempTOA;
% % %         tempMiss(tempMiss1>0)=0;
% % %         tempMiss((rx+1):end,(rx+1):end,:)=1;
% % %     end
% % % end

%add the noise to the distance data.
DD=zeros(rx+ra,rx+ra,numSnapshots);
DD1 = zeros(rx+ra,rx+ra,numSnapshots);
DNoise = zeros(rx+ra,rx+ra,numSnapshots);
% % % noise=max(eps,1+noise1*randn(rx+ra,rx+ra,numSnapshots));
% % % noise((rx+1):end,(rx+1):end,:)=1;
for i=1:numSnapshots
% % %   noise(:,:,i)=(noise(:,:,i)+noise(:,:,i)')/2;
  DD(:,:,i) = PairDist(PP(:,:,i));
end

DD1(DD<=R) = DD(DD<=R);
DNoise(DD<=R) = DD(DD<=R);
% % % RNoise=R*noise;
% % % DNoise=DD./noise;
% % % DNoise(DNoise>RNoise)=0;
% % % MM=DNoise<=RNoise&DNoise~=0;
% % % MM1=MM.*tempMiss;
% % % DD1=DNoise.*MM1;
% % % 
% % % %get the shortest path matrix according to the distance matrix.
% % % DNoise1=DNoise-DD1;
% % % DNoise(DNoise1>0)=RNoise(DNoise1>0);
% % % DNoise=DNoise/10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shortestDisMtx=zeros(rx+ra,rx+ra,numSnapshots);
% EqMatrx.Dx=DD1(1:rx,1:rx,:);
% EqMatrx.Da=DD1(1:rx,(rx+1):end,:);
Dx = zeros(rx,rx,numSnapshots);
Da = zeros(rx,ra,numSnapshots);

%find the upper bound and lower bound.
UBx=zeros(rx,rx,numSnapshots);
UBa=zeros(rx,ra,numSnapshots);
LBx=zeros(rx,rx,numSnapshots);
LBa=zeros(rx,ra,numSnapshots);

%for DNoise, I divide it by 10 and then time 10. That is becasue the
%function graphallshortestpaths() can only deal with small data.
for i=1:numSnapshots
    
    switch Network.ratioType
        case 'NONE'
        case 'measureRatio'
            DNoise(:,:,i) = DNoise(:,:,i) .* Network.meaDisEnt...
                {Network.currentRatio==Network.ratiosMeaDis,i};
        case 'measuredDistance'
            disp(['at i=',num2str(i)]);
            temp = full(Network.meaDis{Network.currentRatio==Network.ratiosMea,i});
            globalConnHop = double(temp <= longestCommDis);
            temp(temp == Inf) = 0;
            DNoise(:,:,i) = temp .* globalConnHop;
            DD1(:,:,i) = DNoise(:,:,i);
        case 'noisedDistance'
            disp(['at i=',num2str(i)]);
            temp = full(Network.meaDis{Network.currentRatio==Network.ratiosMea,i});
            temp(temp == Inf) = 0;
            DNoise(:,:,i) = temp;
            DD1(:,:,i) = DNoise(:,:,i);
    end
    
    if ~TOAPct
        DD1(:,:,i) = zeros(rx+ra,rx+ra);
        temp = DNoise(:,:,i);
        temp(temp>0) = longestCommDis; 
        DNoise(:,:,i) = temp;
    end
    
    shortestDisMtx1=graphallshortestpaths(sparse(DNoise(:,:,i)));   
    if any(any(shortestDisMtx1==inf))
        msg = strcat('Topology not connected at time slot',num2str(i));
        disp(msg);
    end
    shortestDisMtx(:,:,i)=shortestDisMtx1;
    Dx(:,:,i) = DD1(1:rx,1:rx,i);
    Da(:,:,i) = DD1(1:rx,(rx+1):end,i);
    
    DN=DNoise(:,:,i);
    HOP=DN>0;
    NHOP=(1-HOP)&(1-eye(rx+ra));
    LBx1=LBx(:,:,i);
    LBx1(NHOP(1:rx,1:rx))=R;
    LBx(:,:,i)=LBx1;
    
    LBa1=LBa(:,:,i);
    LBa1(NHOP(1:rx,(rx+1):end))=R;
    LBa(:,:,i)=LBa1;
    UBx(:,:,i)=shortestDisMtx(1:rx,1:rx,i);
    UBa(:,:,i)=shortestDisMtx(1:rx,(rx+1):end,i);
    
    UBx1=UBx(:,:,i);
    UBx1(HOP(1:rx,1:rx))=R;
    UBx(:,:,i)=UBx1;
    
    UBa1=UBa(:,:,i);
    UBa1(HOP(1:rx,(rx+1):end))=R;
    UBa(:,:,i)=UBa1;
end

EqMatrx.Dx = Dx;
EqMatrx.Da = Da;
UEqMatrx.UBx=UBx;
UEqMatrx.UBa=UBa;
UEqMatrx.LBx=LBx;
UEqMatrx.LBa=LBa;
