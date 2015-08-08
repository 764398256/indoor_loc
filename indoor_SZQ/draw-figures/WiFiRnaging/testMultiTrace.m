figure
hold on
col = {'b', 'b', 'k', 'g', 'y'};
error = cell(1, 4);
for ii = 2:5
    TraceName=['Square_Speed_Contril_G703_', num2str(ii), '.mat'];
    Square_Speed_Control_Song;
    load('Square_Speed_Error_Indoor_WifiRPCAVA_ALL.mat');
    error{ii} = error_Indoor;
    [n1, x1] = hist(error_Indoor(:), 0:0.1:6);
    plot(x1, cumsum(n1)/Network.T/Network.N, col{ii})
    hold on
    
end
%%
figure

for ii = 2:5
    error_Indoor = error{ii};
    [n1, x1] = hist(error_Indoor(:), 0:0.1:6);
    plot(x1, cumsum(n1)/Network.T/Network.N, col{ii})
    hold on
    
end