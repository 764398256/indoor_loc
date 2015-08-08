function [ bnet ] = set_bnt( bnet, prob_sorted, idx, para_m, para_v, dist_m, v, k, numap, numdev )
%���ñ�Ҷ˹�������
%�������   bnet: �Ѿ������ı�Ҷ˹����
%          prob_sorted: WiFi ��λȷ����λ�ø���
%          idx: λ�ø�������Ӧ����������
%          para_m�� para_v: ѵ������յ��ź�ǿ�ȵľ�ֵ������
%          dist_m: �����֮��ľ�����Ϣ
%          v: ��Ҷ˹����֤�ݵĲ����� ���ϵ��


%��ɢ�ڵ�
for i = 1:numdev
    %��ɢ�ڵ�������������ʱ��ʾ�������й�һ������
    bnet.CPD{i} = tabular_CPD(bnet, i, prob_sorted(1:k)/sum(prob_sorted(1:k)));
%     bnet.CPD{i} = tabular_CPD(bnet, i, 1./k*ones(1,k));
end

%�����ڵ㣨��˹�ڵ㣩
for i = 1:numdev
    for j = 1:numap
        bnet.CPD{i+numdev*j} = gaussian_CPD(bnet, i+numdev*j, 'mean', ...
            para_m(idx(i, 1:k), j)', 'cov', para_v(idx(i, 1:k), j).^2');
    end
end
%�����ڵ�
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

