function [ dist_m ] = dist_map( xtrain, ytrain)
%�����������Ŷ�Ӧ�ľ������

global numtrain

dist_m = sqrt((repmat(xtrain,1, numtrain)-repmat(xtrain,1, numtrain)').^2+...
    (repmat(ytrain,1, numtrain)-repmat(ytrain,1, numtrain)').^2);

end

