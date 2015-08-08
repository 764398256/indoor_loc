function [ dist_m ] = dist_map( xtrain, ytrain, numtrain)
%�����������Ŷ�Ӧ�ľ������


dist_m = sqrt((repmat(xtrain,1, numtrain)-repmat(xtrain,1, numtrain)').^2+...
    (repmat(ytrain,1, numtrain)-repmat(ytrain,1, numtrain)').^2);

end

