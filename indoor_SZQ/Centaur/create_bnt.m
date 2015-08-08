function [ bnet ] = create_bnt( k, numap, numdev )
%���챴Ҷ˹������
%������� numdev���豸��
%        numap: ��վ��Ŀ
%        k: ��Ҷ˹������ɢ�ڵ��ȡֵ����Ŀ
%��� bnet: ����ı�Ҷ˹����


N = (numap+1)*numdev + numdev*(numdev-1)/2;
dag = zeros(N, N);
D = tril(ones(numdev),-1);
dag0 = zeros(numdev, numdev*(numdev-1)/2);

%������Ҷ˹�����ͼ
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

%������Ҷ˹������
discrete_nodes = 1:numdev;
node_sizes = ones(1,N);
node_sizes(1:numdev) = k;
onodes = numdev+1:N;%���Ա��۲�Ľڵ�
bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes, 'observed', onodes);

end

