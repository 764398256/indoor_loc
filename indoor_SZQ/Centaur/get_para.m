function [ para_m, para_v ] = get_para( rss_train, numap, numtrain)
%offline train phrase to build the radio map signal strength distributions
%at any location are represented by parametric distributions
%input: rss_train: train data; numtrain; numap
%output: para_m: estimate  mean value; para: data var


para_m = zeros(numtrain, numap);
para_v = zeros(numtrain, numap);


for i = 1:numtrain
    %�ж�ѵ�����ݸ�ʽ
    if isa(rss_train, 'cell')
        temp = rss_train{i};
    else
        temp = rss_train(:, :, i);
    end
    for j = 1:numap
        para_m(i, j) = mean(temp(:, j));
        para_v(i, j) = std(temp(:, j));
    end
end
%��׼��Ϊ0��Ԫ������Ϊ0.1����������
para_v(para_v == 0) = 0.1;

