function datapro()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%ѵ�����ݴ���
%numacΪAC����Ŀ��numtrainΪѵ�������ص�ĸ���
%trainΪnumtrain*1��cell
%ÿ��cell�����������󣬵�һ������coorΪѵ���������������󣬵ڶ�������rss���ź�ǿ�Ⱦ���ά����numac
[number,traindata] = xlsread('Main file.xls','h2:i100');
numtrain = size(traindata,1);
numac = 3;
xtrain = zeros(numtrain,1);
ytrain = zeros(numtrain,1);
rsstrain = zeros(numtrain,numac);
%train = cell(numtrain,2);
for i=1:numtrain
    coor = str2num(char(traindata(i,1)));
    xtrain(i,1) = coor(1,1);
    ytrain(i,1) = coor(1,2);
    rss = str2num(char(traindata(i,2)));
    rsstrain(i,:)=rss(:);   
end
%save network.mat train numtrain numac;

%��λ���ݴ���
%numdevΪ�豸��Ŀ��numTΪʱ��Ƭ��Ŀ
%sdno = xlsread('Main file.xls','a2:b148');
numdev = 7;
numT = 21;
[number,locdata] = xlsread('Main file.xls','c2:d148');
rssloc = zeros(numdev,numac,numT);
for i=1:numdev
    for j=1:numT
        coor = str2num(char(locdata((i-1)*numT+j,1)));
        xcoor(i,j) = coor(1,1);
        ycoor(i,j) = coor(1,2);
        rss = str2num(char(locdata((i-1)*numT+j,2)));
        rssloc(i,:,j) = rss(:);
    end
end

%���������

dis = inf*ones(numdev,numdev,numT);
[devdis] = xlsread('Main file.xls','p2:q401');
sizeDis = size(devdis,1);
meadis = zeros(sizeDis,1);
caldis = zeros(sizeDis,1);
meadis(:,1) = devdis(:,1);
caldis(:,1) = devdis(:,2);
caldis = caldis/100;
% for i=1:sizeDis
%     devnum = str2num(char(devdis(i)));
%     dis(devnum(1),devnum(2),numdis(i,1)+1) = numdis(i,3);
% end

for i=1:numT
    for j=1:numdev
        for k=j:numdev
            odis = sqrt((xcoor(j,i)-xcoor(k,i))*(xcoor(j,i)-xcoor(k,i))+(ycoor(j,i)-ycoor(k,i))*(ycoor(j,i)-ycoor(k,i)));
            if(odis>0.3 && odis<=4)
                [IDX,D] = knnsearch(caldis,odis);
                dis(j,k,i) = meadis(IDX);
            elseif(odis<=0.3)
                dis(j,k,i) = odis*(1-rand(1,1)/5+rand(1,1)/5); 
            else
                dis(j,k,i)=inf;
            end
            dis(k,j,i) = dis(j,k,i);
        end
    end
end

tw = 8;
tr = 3;
%save����Ϊ.mat
Network = struct('xtrain',{xtrain},'ytrain',{ytrain},'rsstrain',{rsstrain},'numtrain',{numtrain},...
                    'numac',{numac},'numdev',{numdev},'numT',{numT},'rssloc', {rssloc},'dis',{dis},'tw',{tw},'tr',{tr},'pro',{'Wifi'});
X = struct('xcoor',{xcoor});
Y = struct('ycoor',{ycoor});
%save network.mat train numtrain numac numdev numT rssloc xcoor ycoor;

save trace.mat Network X Y;
end

