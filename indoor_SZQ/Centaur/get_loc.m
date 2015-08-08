function [ prob_sorted, idx ] = get_loc( para_m, para_v, rss_testT, numdev, numap, numtrain)
%计算设备可能的位置以及概率
%输入 para_m: 接收信号强度的均值
%     para_v: 接收信号强度分布的标准差
%     rss_testT: 某个时间片接收到的信号强度值
%     numdev, numap, numtrain: 设备个数， ap个数， 训练结点个数
%输出 prob_sorted: 定位在各个网格点的概率值
%     idx: 对应网格点序号


prob = ones(numdev, numtrain);

 
for i = 1:numdev
    for j = 1:numtrain
        for k = 1:numap
            prob(i, j) = prob(i, j) * (normcdf(rss_testT(i, k)+0.5, ...
                para_m(j, k), para_v(j, k))-normcdf(rss_testT(i, k)-0.5, ...
                para_m(j, k), para_v(j, k)));
        end
    end
end

%knn search
[prob_sorted, idx] = sort(prob, 2, 'descend');

end

