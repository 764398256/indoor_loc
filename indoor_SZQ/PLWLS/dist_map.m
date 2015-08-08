function [ dist_m ] = dist_map( xtrain, ytrain)
%构造网格点序号对应的距离矩阵

global numtrain

dist_m = sqrt((repmat(xtrain,1, numtrain)-repmat(xtrain,1, numtrain)').^2+...
    (repmat(ytrain,1, numtrain)-repmat(ytrain,1, numtrain)').^2);

end

