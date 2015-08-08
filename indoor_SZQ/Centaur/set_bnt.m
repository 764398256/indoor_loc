function [ bnet ] = set_bnt( bnet, prob_sorted, idx, para_m, para_v, dist_m, v, k, numap, numdev )
%设置贝叶斯网络参数
%输入参数   bnet: 已经建立的贝叶斯网络
%          prob_sorted: WiFi 定位确定的位置概率
%          idx: 位置概率所对应的网格点序号
%          para_m， para_v: 训练点接收到信号强度的均值，方差
%          dist_m: 网格点之间的距离信息
%          v: 贝叶斯网络证据的参数， 相关系数


%离散节点
for i = 1:numdev
    %离散节点概率用条件概率表表示，并进行归一化处理
    bnet.CPD{i} = tabular_CPD(bnet, i, prob_sorted(1:k)/sum(prob_sorted(1:k)));
%     bnet.CPD{i} = tabular_CPD(bnet, i, 1./k*ones(1,k));
end

%连续节点（高斯节点）
for i = 1:numdev
    for j = 1:numap
        bnet.CPD{i+numdev*j} = gaussian_CPD(bnet, i+numdev*j, 'mean', ...
            para_m(idx(i, 1:k), j)', 'cov', para_v(idx(i, 1:k), j).^2');
    end
end
%连续节点
D = tril(ones(numdev),-1);
for i = 1:numdev
    for j = i+1:numdev
        Mean = dist_m(idx(i, 1:k), idx(j, 1:k));
        Var = v*ones(k);
        bnet.CPD{(numap+1)*numdev+sum(D(1:(i-1)*numdev + j))} = gaussian_CPD(bnet,...
        (numap+1)*numdev+sum(D(1:(i-1)*numdev + j)), 'mean', Mean(:), 'cov', Var(:).^2);
    end
end

end

