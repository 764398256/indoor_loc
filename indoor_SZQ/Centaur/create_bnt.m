function [ bnet ] = create_bnt( k, numap, numdev )
%构造贝叶斯网络框架
%输入参数 numdev：设备数
%        numap: 基站数目
%        k: 贝叶斯网络离散节点可取值得数目
%输出 bnet: 构造的贝叶斯网络


N = (numap+1)*numdev + numdev*(numdev-1)/2;
dag = zeros(N, N);
D = tril(ones(numdev),-1);
dag0 = zeros(numdev, numdev*(numdev-1)/2);

%构建贝叶斯网络的图
for i = 1:numdev
    for j = i+1:numdev
        dag0(i, sum(D(1:(i-1)*numdev + j))) = 1;
        dag0(j, sum(D(1:(i-1)*numdev + j))) = 1;
    end
end

for i = 1:numdev
    dag(i, numdev+i:numdev:numdev*(numap+1)) = 1;
end

dag(1:numdev, (numap+1)*numdev+1:end) = dag0;
% view(biograph(dag));

%构建贝叶斯网络框架
discrete_nodes = 1:numdev;
node_sizes = ones(1,N);
node_sizes(1:numdev) = k;
onodes = numdev+1:N;%可以被观察的节点
bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes, 'observed', onodes);

end

