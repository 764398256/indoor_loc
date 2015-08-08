function [error_Horus, MError_Horus, MMError_Horus ] = horus_loc( Network,X,Y,ErrorName )
% Horus: to estimate the user location based on a received signal strength vector, using the constructed radio map.
% rss_train: train daa; rss_test: rss tested of one time second.

l = 4;
numdev = Network.rx;
numap = Network.numap;
numtrain = Network.numtrain;
rss_train = Network.rssTrain;
numT = Network.T;
rss_test = Network.rssTest;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
xcoor = X.data;
ycoor = Y.data;
G_real(:, 1, :) = xcoor;
G_real(:, 2, :) = ycoor;
[ para_m, para_v ] = get_para(rss_train, numap, numtrain);

G_Horus = zeros(numdev, 2, numT);
error_Horus = zeros(numdev, numT);
for t = 1:numT
    rss_testT = rss_test(:, :, t);
    
    prob = ones(numdev, numtrain);
    loc = zeros(numdev, 2);
    
    for i = 1:numdev
        for j = 1:numtrain
            for k = 1:numap
                prob(i, j) = prob(i, j) * (normcdf(rss_testT(i, k) + 0.5, ...
                    para_m(j, k), para_v(j, k))-normcdf(rss_testT(i, k) - ...
                    0.5, para_m(j, k), para_v(j, k)));
            end
        end
    end
    
    %knn search
    [prob_sorted, idx] = sort(prob, 2, 'descend');
    % save probbb prob
    %arithmetic mean
    % for i = 1:numdev
    %     loc(i, 1) = mean(xtrain(idx(i,1:k)));
    %     loc(i, 2) = mean(ytrain(idx(i,1:k)));
    % end
    
    %weighted arithmetic mean
    for i = 1:numdev
        temp = sum(prob_sorted(i, 1:l));
        loc(i, 1) = sum(prob_sorted(i, 1:l).*xtrain(idx(i, 1:l))')/temp;
        loc(i, 2) = sum(prob_sorted(i, 1:l).*ytrain(idx(i, 1:l))')/temp;
    end
     
    G_Horus(:, :, t) = loc;
    error_Horus(:, t) = sqrt(sum((G_real(:, :, t)-loc).^2, 2));
end

MError_Horus = mean(error_Horus);
MMError_Horus = mean(MError_Horus);
end



