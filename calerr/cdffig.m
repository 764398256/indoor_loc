clear;
close all;
warning off all;
clc


FileName = 'Square_Speed_Error_Indoor_';
lenName = length(FileName);

str=pwd;
index_dir=strfind(str,'\');
str_temp=str(1:index_dir(end)-1);
addpath(genpath(str_temp));
if strcmp(str(index_dir(end)+1:end),'calerr')
    str=str(1:index_dir(end)-1);
end
strtmp = [str,'\run-error-vs-speed'];
addpath (strtmp);
x = dir(strtmp);
len = length(dir(strtmp));
j = 0;
str = ['IndoorError'];

for i=1:len
    lenFile = length(x(i).name);
    if lenFile>lenName && strcmp(x(i).name(1:lenName),FileName) && strcmp(x(i).name(lenFile-3:lenFile),'.mat')
        j = j+1;
        load (char(x(i).name))
        label = strfind(x(i).name(lenName+1:end),'_');
        strName = x(i).name(lenName+1:lenName+label-1);
        c{1} = strName;
        str = [str,'_',strName];
        c{2} = Error_Indoor;
        cel{j} = c;
    end
end

lin = ['-',':','+.','*','.'];
col = ['b','r','g','c''y'];
gcf = figure;
set(gcf,'NumberTitle','off','Name',str);

for i=1:size(cel,2)
    c = cel{i};
    strName = c{1};
    tmpError = c{2};
    tmpError = reshape(tmpError(:,8:end),1,[]);
    [h,stats] = cdfplot(tmpError);
    ah{i} = h;
    ast{i} = stats;
    stats.Name = strName;
    astats{i} = stats;
    set(h,'color',col(i),'linestyle',lin(i),'linewidth',1.5);
    mylgd{i} = [strName];
    hold on;
end
xlabel('Error(meter)');
ylabel('Probability');
legend(mylgd,2);
saveas(gcf,str,'fig');

strtmp = [str,'_T'];
mkdir(strtmp);
cd (char(strtmp));
for i=1:size(cel{1}{2},2)
    gcf = figure;
    strT = [str,'_T',num2str(i)];
    set(gcf,'NumberTitle','off','Name',strT);
    for j=1:size(cel,2)
        c = cel{j};
        size(c);
        strName = c{1};
        tmpError = c{2}(:,i);
        [h,stats] = cdfplot(tmpError);
        th{j,i} = h;
        tst{j,i} = stats;
        stats.Name = [strName,'_T',num2str(i)];
        tstats{j}{i} = stats;
        set(h,'color',col(j),'linestyle',lin(j),'linewidth',1.5);
        mylgd{j} = [strName];
        hold on;
    end
    xlabel('Error(meter)');
    ylabel('Probability');
    legend(mylgd,2);
    saveas(gcf,strT,'fig');
    close;
end
cd ..;
% save stats.mat astats tstats;
clear;




