X.data = round(X.data*100)/100;
Y.data = round(Y.data*100)/100;

%% draw trace figure
for i = 1:18
    figure
    plot(X.data(i, :), Y.data(i, :), '.');
    for j = 1:30
        text(X.data(i, j) + .1, Y.data(i, j) + .1, [num2str(j) '(' ...
            num2str(X.data(i, j)) ',' num2str(Y.data(i, j)) ')']...
            , 'FontSize', 8, 'FontName', 'Times New Roman');
    end
    title([num2str(i) 'th trace'])
    set(gcf,'outerposition',get(0,'screensize'));
    grid on
    saveas(gcf, ['E:\金山快盘\sharebox\717zhj@163.com\WSNs_localization\our-scheme\indoor-wifi-ranging\indoor_SZQ\traceFigure\' num2str(i) 'th trace.pdf'])
end

%% write trace file
format = repmat('%4.2f ', 1, 30);
fid = fopen('E:\金山快盘\sharebox\717zhj@163.com\WSNs_localization\our-scheme\indoor-wifi-ranging\indoor_SZQ\traceFigure\traceFile.txt', 'w');
for i = 1:18
    fprintf(fid, ['trace no: ' num2str(i) '\n']);
    for j = 1:30
        fprintf(fid, [num2str(j) '(' num2str(X.data(i, j)) ', ' num2str(Y.data(i, j)) ')  ']);
    end
end
fclose(fid);   