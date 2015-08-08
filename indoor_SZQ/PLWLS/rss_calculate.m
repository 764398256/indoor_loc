function [ rss_loc ] = rss_calculate( loc, rss_train, xtrain, ytrain, L )
%Given a loc matrix, return the rss vector of the loc matrix
%

[m, ~] = size(loc);
[~, n] = size(rss_train);
rss_loc = zeros(m, n);

idx = knnsearch([xtrain ytrain], loc, 'k', L);
for i = 1:m
    rss_loc(i, :) = mean(rss_train(idx(i, :), :));
end