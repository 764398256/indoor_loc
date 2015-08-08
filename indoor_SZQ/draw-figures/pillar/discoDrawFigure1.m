function [handle]=discoDrawFigure1(xdata,ydata,info)
%This function needs three input parameters:
% xdata--data to be plot along the x axis, size of 1*m;
% ydata--data to be plot along the y axis, should be a matrix of n*m, n is the number of curves to be plot;
% info --setting information about the figure, containing the following fields:
%       LineStyle       --line style of the curves to be plot, a string cell with size of 1*n;
%       MarkerSize      --size of the markers of the data points;
%       LineWidth       --width of the lines to be plot;
%       MarkerFaceColor --face colors of the makers, a string cell with size of 1*n;
%       LegendText      --legend text of the legends, a string cell with size of 1*n;
%       row1num         --the number of legends on the first row;
%       xlabel          --label of the x axis;
%       ylabel          --label of the y axis;
%       title           --title of the figure;
%
%The output is a structure with fields of all kinds of handles of the figure.
%
%Example:
%   xdata=1:10;
%   ydata=rand(2,length(xdata));
%   info.row1num=1;
%   ...
%   [handle]=DrawFigure(xdata,ydata,info)
%

if nargin<3
   disp('The third input parameter needs to be input, which is about the setting information of the figure.');
   info=struct(); 
end
[n,m]=size(ydata);            %number of line to be plot
%compatibility setting
if ~isfield(info,'row1num')
   info.row1num=ceil(n/2); 
end
if ~isfield(info,'LineWidth')
   info.LineWidth=1; 
end
if ~isfield(info,'LineStyle')
    LS={'g-.x','m--+','k-.o','r--s','b-d','r:','g-*','c->','b-<'};
   info.LineStyle=LS(1:n);
end
if ~isfield(info,'MarkerFaceColor')
    MFC={'g','m','k','r','b','r','w','w','w','w'};
   info.MarkerFaceColor=MFC(1:n); 
end
if ~isfield(info,'MarkerSize')
   info.MarkerSize=5; 
end

if ~isfield(info,'plotType')
    info.plotType = 1;
end

hFigure =figure;          %store the handle of the current figure
switch info.plotType
    case 1
        for a=1:n
            hLine(a)=plot(xdata, ydata(a,:),info.LineStyle{a},'MarkerSize',info.MarkerSize,...
                'LineWidth',info.LineWidth,'MarkerFaceColor',info.MarkerFaceColor{a});
%             grid on
            hold on;
        end
    case 2
        h = bar(ydata,0.5);
%         ch = get(h,'children');
        % set(ch,'FaceVertexCData',[0 1 0 ; 1 0 1;  1 0 0 ;1 0 0;    0 0 1
        % ; 0 0 1 ; 0 0 1 ])
end

% hAxis =get(hFigure ,'children');        %store the handle of the current axes

if ~isfield(info,'titleSize')
    info.titleSize = 12;
end

if isfield(info,'title')
    title(info.title,'FontSize',info.titleSize);
else
    title('unnamed title');
end

if ~isfield(info,'labelFontSize')
    info.labelFontSize = 12;
end

if ~isfield(info,'xlabel')
    info.xlabel = 'Speed (m/interval)';
end
xlabel(info.xlabel,'interpreter','latex','fontsize',info.labelFontSize);

if ~isfield(info,'ylabel')
    info.ylabel = 'Mean absolute error';
end
ylabel(info.ylabel,'interpreter','latex','fontsize',info.labelFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(info,'figPosition')
    set(gcf,'Position',info.figPosition);
else
    set(gcf,'Position',[200,200,800,600]);
end


if isfield(info,'xscale')
    set(gca,'xscale',info.xscale);
end

if isfield(info,'yscale')
    set(gca,'yscale',info.yscale);
end

if isfield(info,'axis')
    axis(info.axis);  %[1,50,0,0.8]
end

if isfield(info,'xtick')
    set(gca,'xtick',info.xtick); % [1 2 5 10:10:50]
else
    set(gca,'xtick',xdata);
end

if isfield(info,'ytick')
    set(gca,'ytick',info.ytick); % 0:0.1:1
else
    set(gca,'ytick',0:0.2:1);
end

if isfield(info,'xticklabel')
    set(gca,'XTickLabel',info.xticklabel); % [1 2 5 10:10:50]
end

if isfield(info,'yticklabel')
    set(gca,'YTickLabel',info.yticklabel); % 0:0.1:1
end

if isfield(info,'chartPosition')
    set(gca,'position',info.chartPosition); % [0.1  0.12  0.88  0.86]
end

if ~isfield(info,'tickFontSize')
    info.tickFontSize = 12;
end
set(gca,'fontsize',info.tickFontSize); % 0:0.1:1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(info,'LengendText')
    LengendText = {'MCL','MSL*','TSLRL-RF','TSLRL','DISCO-1H-RF','DISCO-1H','DISCO-2H-RF','DISCO-2H'};
else
    LengendText = info.LengendText;
end

% set double 
if ~isfield(info,'legendType')
    info.legendType = 2;
end

if ~isfield(info,'orientation')
    info.orientation = 'vertical';
end

if ~isfield(info,'location')
    info.location = 'SouthEast';
end

if ~isfield(info,'legendFontSize')
    info.legendFontSize = 12;
end

switch info.plotType
    case 1
        
        switch info.legendType
            case 1
                row1num = info.row1num;
                hLegend(1)=legend(hLine(1:row1num),LengendText(1:row1num),'Orientation','horizontal','Location','North');%store the handle of the current legend
                legend boxoff
                lepos = get(hLegend(1), 'Position');
                lepos(1) = lepos(1) + 0.012;
                lepos(2) = lepos(2) - 0.03;
                hLegend(2) = copyobj(hLegend ,hFigure );
                hLegend(2)= legend(hLine(row1num+1:end),LengendText(row1num+1:end),'Orientation','horizontal','Location','North');
                set(hLegend(2),'position',lepos);                
                                
            case 2
                hLegend = legend(LengendText,'Orientation',info.orientation,'location',info.location,'FontSize',info.legendFontSize);
                
            case 3 
                hLegend = [];
        end
        
%         legend boxoff
        handle.hLine=hLine;
        handle.hLegend=hLegend;
        
    case 2
        set(gca,'xticklabel',info.LengendText);
%         th=rotateticklabel(gca,info.legendFontSize);
end

handle.hFigure=hFigure;
% handle.hAxis=hAxis;

if isfield(info,'epsFileName')
    set(gcf,'PaperPositionMode','auto');
    print(gcf,'-dpsc2',info.epsFileName);
end

end

