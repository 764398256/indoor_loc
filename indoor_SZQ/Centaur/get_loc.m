function [ prob_sorted, idx ] = get_loc( para_m, para_v, rss_testT, numdev, numap, numtrain)
%�����豸���ܵ�λ���Լ�����
%���� para_m: �����ź�ǿ�ȵľ�ֵ
%     para_v: �����ź�ǿ�ȷֲ��ı�׼��
%     rss_testT: ĳ��ʱ��Ƭ���յ����ź�ǿ��ֵ
%     numdev, numap, numtrain: �豸������ ap������ ѵ��������
%��� prob_sorted: ��λ�ڸ��������ĸ���ֵ
%     idx: ��Ӧ��������


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

