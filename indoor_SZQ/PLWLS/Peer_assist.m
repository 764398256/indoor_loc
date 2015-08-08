function [ error_PLWLS, MError_PLWLS, MMError_PLWLS ] = Peer_assist( Network,X,Y,ErrorName )
% main fun of PLWLS
%  L select the L longest edges
%  P the number of peers assisted
%  K knnsearch the number of neibors

L = 3;%select the L longest edges
P = 3;%the number of peers assisted
K = 4;
numdev = Network.rx;
numT = Network.T;
longestCommDis = Network.longestCommDis;
G_PLWLS = zeros(numdev, 2, numT);
error_PLWLS = zeros(numdev, numT);
dist = Network.dist;
rss_test = Network.rssTest;
rss_train_mean = Network.rssTrainMean;
xtrain = Network.xtrain;
ytrain = Network.ytrain;
xcoor = X.data;
ycoor = Y.data;
G_real(:, 1, :) = xcoor;
G_real(:, 2, :) = ycoor;
[ ~, ~, ~, G_ini] = radar_loc( Network,X,Y,ErrorName);

for l = 1:numT
    fprintf('time slot %d\n', l);
    real_loc = G_real(:, :, l);
    est_loc = G_ini(:, :, l);
    distT = dist(:, :, l);
    error_loc = sqrt(sum((real_loc - est_loc).^2, 2));
    for i = 1:numdev
%         fprintf('time slot %d\n', i);
        error_loc(i) = inf;
        distI = distT(i, :);
        error_loc(distI > longestCommDis) = 10;
        [~, idx] = sort(error_loc);
%         idx = randperm(numdev)';
        idx_assit =[i; idx(1:P)];
        est_loc_assit = est_loc(idx_assit, :);
        real_loc_assit = real_loc(idx_assit, :);
        dist_assit = distT(idx_assit, idx_assit);
        rss_test_assit = rss_test(idx_assit, :, l);
        
        G_constr = MDS_MAP( dist_assit );
        [ loc_rotated, dist_constr ] = rotate_longest( G_constr );
%         disp('MDS error');
%         mean(sqrt((dist_constr(:) - dist_assit(:)).^2))
        [ dire_vec, dire_vec0, dire_vec1, loc_rotated0 ] = edge_direct( loc_rotated, est_loc_assit,L );
        loc_orient = orientation( dire_vec,dire_vec0,dire_vec1,loc_rotated,loc_rotated0 );
        loc_orient=loc_orient + repmat(mean(est_loc_assit(:,:) -loc_orient(:,:)), P+1, 1);
        
%         figure
%         max_text = {'1', '2', '3', '4'};
%         plot(est_loc_assit(:,1),est_loc_assit(:,2), 'k+')
%         hold on
%         plot(real_loc_assit(:,1),real_loc_assit(:,2), 'r*')
%         plot(loc_orient(:,1),loc_orient(:,2), 'ks')
%         plot(G_constr(:,1),G_constr(:,2), 'bo')
%         text(est_loc_assit(:,1)+0.1, est_loc_assit(:,2)+0.1, max_text);
%         text(real_loc_assit(:,1)+0.1, real_loc_assit(:,2)+0.1, max_text)
%         text(loc_orient(:,1)+0.1, loc_orient(:,2)+0.1, max_text)
%         text(G_constr(:,1)+0.1, G_constr(:,2)+0.1, max_text)
        
        est_loc_assit = location_est( loc_orient, est_loc_assit, rss_train_mean, rss_test_assit, xtrain, ytrain, K);
        
%         plot(est_loc_assit(:,1),est_loc_assit(:,2), 'bo')
%         text(est_loc_assit(:,1)+0.1, est_loc_assit(:,2)+0.1, max_text);
%         legend( 'WiFi est location', 'Real localization', 'PLWLS est localization 1', 'PLWLS est localization 2', 2)
        
        G_PLWLS(i, :, l) = est_loc_assit(1, :);
        error_PLWLS(i, l) = sqrt(sum((G_PLWLS(i,:,l) - G_real(i,:,l)).^2, 2));
    end
    
end
MError_PLWLS = mean(error_PLWLS);
MMError_PLWLS = mean(MError_PLWLS);

