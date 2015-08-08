function D = PairDist(X, dist)

%
% D = PAIRDIST(X,dist) returns a square matrix that gives the pairwise
% distance between points in X using distance function 'dist'
%
% file:        PairDist.m

  if nargin < 2, dist = 'euclidean'; end
  
  D = squareform(pdist(X,dist));
