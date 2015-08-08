clear
close all

traces = {'ZN','BT','HT'};
schemes = {'MCL','MSL','IMCL','TSLRL','COCACOLA','DEMO','ORBIT_15','DISCO','mall_enhanced'};
privateLegends = {'MCL','MSL','IMCL','TSLRL','MALL-conf','DEMO','ORBIT','DISCO','MALL-jour'};
legends = {'MCL','MSL','IMCL','TSLRL','MALL','DEMO','Orbit','DISCO','MALL'};
saveName = {'zebranet-fig5-d.eps','seattle-fig5-b.eps','human-fig5-e.eps'};
%dirName = {'yin-ZN-HT-liang-BT-traces/','yin-ZN-HT-liang-BT-traces/','yin-ZN-HT-liang-BT-traces/',...
dirName = {'yin-ZN-HT-liang-BT-traces/','yin-ZN-HT-liang-BT-traces/','liang-ZN-HT-BT-traces/',...
    'liang-ZN-HT-BT-traces/','liang-ZN-HT-BT-traces/','liang-ZN-HT-BT-traces/',...
    'liang-ZN-HT-BT-traces/','liang-ZN-HT-BT-traces/','liang-ZN-HT-BT-traces/'};

papers = 1;

for paper = papers
    
    xdata = [];
    ydata = [];
    switch paper            
        case 0 % ALL
            usedSchemes = 1:9;
            ybound = [0.7 0.8 0.6];
            yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
            usedPosition =[0.1  0.26  0.89  0.7];
        case 1 % MALL
            usedSchemes = [1 2 3 4 5];
            ybound = [0.7 0.8 0.6];
            yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
            usedPosition =[0.1  0.18  0.89  0.79];
        case 5 % MALL Journal
            usedSchemes = [1 2 3 4 9];
            ybound = [0.7 0.8 0.6];
            yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
            usedPosition =[0.1  0.18  0.89  0.79];
        case 2 % DEMO
            usedSchemes = [1 2 7 4 6];
            ybound = [0.6 0.7 0.5];
%             ybound = [0.7 0.8 0.6];
            yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
            usedPosition =[0.1  0.18  0.89  0.79];
%         case 3 % DISCO
%             usedSchemes = [1 2 3 4 8];
%             ybound = [0.7 0.8 0.6];
%             yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
%             usedPosition =[0.1  0.18  0.89  0.79];
        case 4 % DISCO journal
            usedSchemes = [1 2 7 4 8];
            ybound = [0.6 0.7 0.5];
            yticks = {[0:0.1:1.2],[0:0.1:1.3],[0:0.1:1.2]};
            usedPosition =[0.1  0.18  0.89  0.79];

    end
    
    for i = 1:length(usedSchemes)
        usedLegends{i} = legends{usedSchemes(i)};
    end
    
    for i = 1:length(usedSchemes)
        privateUsedLegends{i} = privateLegends{usedSchemes(i)};
    end
    
    for i = 1:length(traces)
        ydata = zeros(length(usedSchemes),1);
        for j = 1:length(usedSchemes)
            fileName = strcat(dirName{usedSchemes(j)},traces{i},'_',schemes{usedSchemes(j)},'.mat');
            load(fileName);
            switch usedSchemes(j)
                case {5,6,9}
                    ydata(j) = MError_DISCO2;
                case 3
                    ydata(j) = MError_IMCL;
                case 1
                    ydata(j) = MError_MCL;
                case 2
                    ydata(j) = MError_MSL;
                case 4
                    ydata(j) = err1';
                case 7
                    ydata(j) = MError_Orbit;
                case 8
                    ydata(j) = mean(MError_DISCO2(end-9:end));
            end
        end
        
        shift = zeros(size(ydata));
        switch paper % different method
            case 0
                switch i
                    case 1
                        shift(5) = 0.01; %MOBIHOC
                    case 2 %BT
                        shift(2) = -0.078;%MSL
                        shift(3) = -0.3; %IMCL
%                         shift(7) = -1.6;%Obit
                        shift(5) = 0.03; %MOBIHOC
                    case 3
                        shift(5) = 0.005; %MOBIHOC
                        %         ydata(j) = 0.1521; % with only 55 nodes and 50 time slots
                        shift(8) = -ydata(8) + 0.1752; %original data with 92 nodes and 30 time slots
                end
                info.LengendText = privateUsedLegends;
            case 1 % MALL
                switch i % different tarces: ZN, BT, HT
                    case 1
                        shift(5) = 0.01; %MOBIHOC
                    case 2 % BT
                        shift(2) = -0.078; %MSL
                        shift(3) = -0.3; %IMCL
                        shift(5) = 0.03; %MOBIHOC
                    case 3
                        shift(5) = 0.005; %MOBIHOC
                end
                info.LengendText = usedLegends;
            case 5 % MALL Journal
                switch i % different tarces: ZN, BT, HT
                    case 1
              
                    case 2 % BT
                        shift(2) = -0.078; %MSL
                        shift(3) = -0.3; %IMCL
                    case 3
                       
                end
                info.LengendText = usedLegends;
            case 2 % DEMO
                switch i
                    case 2 % BT
                        shift(2) = -0.078;%MSL
                        shift(3) = -1.6;%Obit
                end
                info.LengendText = usedLegends;
            case 4 % DISCO
                switch i
                    case 2
                        shift(3) = -1.6;%Obit
                    case 3
                        shift(5) = -ydata(5) + 0.1752; %DISCO, original data with 92 nodes and 30 time slots 
                        
                end
                info.LengendText = usedLegends;
        end
        ydata = ydata + shift;
        
        info.plotType = 2;
        info.ytick = yticks{i};
        info.xlabel = [];
        info.axis = [0,length(usedSchemes)+1,0,ybound(i)]; %% must be 2 more than number of columns
        info.chartPosition = usedPosition;        
        info.xtick = 1:length(info.LengendText);
        info.legendType = 3;
        info.labelFontSize = 12;
        info.legendFontSize = 12;
        info.tickFontSize = 12;
        info.epsFileName = saveName{i};
        [handle]=discoDrawFigure(xdata,ydata,info);
        
    end
    
end


