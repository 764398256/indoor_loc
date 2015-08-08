function X = MDS(D,A)
%
% [X] = MDS(D,A) applies classic multi-dimensional scaling to estimate
% the sensor location matrix X
%
% Input:
%
%       D:    PairDist([X;A])^2 gives the pair-wise squared Euclidean 
%             between all the sensors and anchors (size: (rx+ra)-by-(rx+ra))
%
%       A:    known location of the anchors (size: ra-by-d)
%
% Output:
%
%       X:    estimated location of the sensrs (size: rx-by-d)
%
% file:        MDS.m

  % get size information
  n      = size(D,1);
  [ra,d] = size(A);
  rx     = n - ra;
  
  %
  % Double center D in a memory efficient manner.  Functionally the code
  % below is equivalent to the following
  %
  % B = -0.5 * (D - repmat(D_avg,n,1) - repmat(D_avg',1,n) + mean(D_avg));
  %
  B = 0.5*(D+D');
  D_avg = mean(B,1);
  D_avg = D_avg - 0.5*mean(D_avg);
  for i = 1:n
    B(i,:) = B(i,:) - D_avg;
  end
  D_avg = D_avg';
  for i = 1:n
    B(:,i) = B(:,i) - D_avg;
  end
  B = -0.5*B;
  
  % ensure symmetry
  B = 0.5*(B + B');
  
  % eigendecomposition
  [U,S] = eigs(B,d,'LA');
  clear B;
  
  % form Z
  Z = U*sqrt(max(S,0));
  
  % no alignment is needed
  if (ra == 0)
    X = Z(1:rx,:);
    return
  end
  
  %
  % find rotation matrix R and center vector w (of size 1xd) s.t.
  %
  %    Z*R + ones(n,1)*w \approx [X;A]
  %
  % that is:
  %
  %    [Z 1] * [R; w] \approx [X; A]
  %
  % We can find R and w by solving
  %
  %    [Z(ind_a,:) 1] * [R; w] \approx A
  %
  Za = Z((n-ra+1):n,:);
  Y  = [Za ones(ra,1)];  
  Rw = pinv(Y'*Y)*(Y'*A);
  R  = Rw(1:d,:);
  w  = Rw(d+1,:);

  % ensure that R is orthogonal
  [Ur,Sr,Vr] = svd(R);
  R = Ur*Vr';
  
  % update w
  w = mean(A-Za*R,1);

  % return the final X
  X = Z(1:rx,:)*R + repmat(w,rx,1);
