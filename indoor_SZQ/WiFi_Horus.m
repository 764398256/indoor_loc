function [G_Horus] =WiFi_Horus(Network)

numtrain = Network.numtrain;
numac = Network.numap;
numdev = Network.N;
SnapShots = Network.T;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
rss_train = Network.rssTrain;
[ para_m, para_v ] = get_para(rss_train, numac, numtrain);
G_Horus = zeros(numdev, 2, SnapShots);
l = 4;
for t = 1:SnapShots

    rss_testT = Network.rssTest(:,:,t);%i时刻的用来进行定位的rss矩阵
    
    prob = ones(numdev, numtrain);
    loc = zeros(numdev, 2);
    
    for i = 1:numdev
        for j = 1:numtrain
            for k = 1:numac
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
    
end

end