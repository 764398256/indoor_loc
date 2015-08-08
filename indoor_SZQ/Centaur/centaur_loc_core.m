function [ G_CENTAUR, error_CENTAUR ] = centaur_loc_core( G_real, rss_test, dist, xtrain, ytrain, dist_m, para_m, para_v, k)
% CENTAUR 定位算法
[numdev, numap, numT] = size(rss_test);
numtrain = length(xtrain);

% k = 6;%设置HORUS定位结果保留的网格点个数
k0 = 4;%CENTAUR定位K近邻K的值
dist_var = 1;%距离分布的参数方差

LOC_CENTAUR = zeros(numdev, 2, numT);
% LOC_HORUS = G_real;
error_Centaur = zeros(numdev, numT);
% error_HORUS = error_CENTAUR;
[ bnet ] = create_bnt( k, numap, numdev );

fprintf('CENTAUR start. numT: %d, numdev %d\n', numT, numdev);

for l = 1:numT;
    tic
    fprintf('time slot %d\n', l);
    [ prob_HORUS, idx_HORUS ] = get_loc( para_m, para_v, rss_test(:, :, l),numdev, numap, numtrain );
    [ bnet ] = set_bnt( bnet, prob_HORUS, idx_HORUS, para_m, para_v, dist_m, dist_var, k, numap, numdev );
    
    %推断引擎设置
    engine = jtree_inf_engine(bnet);
%     engine = cond_gauss_inf_engine(bnet);
%     engine = pearl_inf_engine(bnet, 'max_iter', 30);
    
    
    evidence = cell(1,(numap+1)*numdev+numdev*(numdev-1)/2);
    
    %设置接收信号强度证据
    for i = 1:numdev
        for j = 1:numap
            evidence{i+j*numdev} = rss_test(i, j, l);
        end
    end
    
    %设置设备之间距离证据
    D = tril(ones(numdev),-1);
    for i = 1:numdev
        for j = i+1:numdev
            evidence{(numap+1)*numdev+sum(D(1:(i-1)*numdev + j))} = dist(i, j, l);
        end
    end
    
    %推理
    [engine, ll] = enter_evidence(engine, evidence);
    
    prob  = zeros(numdev, k);
    for i = 1:numdev
        marg = marginal_nodes(engine, i);
        prob(i, :) = marg.T;
    end
    
    [prob_CENTAUR, idx_CENTAUR] = sort(prob, 2, 'descend');
    
    loc_CENTAUR = zeros(numdev, 2);
%     loc_HORUS = zeros(numdev, 2);
    for i = 1:numdev
        temp = sum(prob_CENTAUR(i, 1:k0));
        loc_CENTAUR(i, 1) = sum(prob_CENTAUR(i, 1:k0).*xtrain(idx_HORUS(i, idx_CENTAUR(i, 1:k0)))')/temp;
        loc_CENTAUR(i, 2) = sum(prob_CENTAUR(i, 1:k0).*ytrain(idx_HORUS(i, idx_CENTAUR(i, 1:k0)))')/temp;
%         temp0 = sum(prob_HORUS(i, 1:k0));
%         loc_HORUS(i, 1) = sum(prob_HORUS(i, 1:k0).*xtrain(idx_HORUS(i, 1:k0))')/temp0;
%         loc_HORUS(i, 2) = sum(prob_HORUS(i, 1:k0).*ytrain(idx_HORUS(i, 1:k0))')/temp0;
    end
    
    LOC_CENTAUR(:, :, l) = loc_CENTAUR;
%     LOC_HORUS(:, :, l) = loc_HORUS;
    error_CENTAUR(:, l) = sqrt(sum((G_real(:, :, l)-loc_CENTAUR).^2, 2));
%     error_HORUS(:, l) = sqrt(sum((G_real(:, :, l)-loc_HORUS).^2, 2));
end
G_CENTAUR = LOC_CENTAUR;
toc
