function [X] = LocalizeLowRankSDP_3D_MICRO_enhanced(rx,ra,d,nt,A,cons,X0,W0,Display,lambda)
%
% [X] = LOCALIZELOWRANKSDP_3D(rx,ra,d,nt,A,cons,X0,Display) returns
% estimated locations for sensors using low-rank SDP.
%
% Input:
%
%        rx:   size(X,1)
%
%        ra:   size(A,1)
%
%         d:   number of dimensions (= size(X,2) = size(A,2))
%
%        nt:   number of time points (= size(X,3) = size(A,3))
%
%         A:   a ra-by-d-by-nt tensor giving the locations of anchor nodes
%              at different time points
%
%      cons:   constraints.
%
%               cons.EQx:   a rx-by-rx-by-nt tensor giving equality constraints
%                           on sensor-to-sensor distance.
%                           alternatively, a rx*rx-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.EQa:   a rx-by-ra-by-nt tensor giving equality constraints
%                           on sensor-to-anchor distance.
%                           alternatively, a rx*ra-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.UBx:   a rx-by-rx-by-nt tensor giving upper bounds
%                           on sensor-to-sensor distance.
%                           alternatively, a rx*rx-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.UBa:   a rx-by-ra-by-nt tensor giving upper bounds
%                           on sensor-to-anchor distance.
%                           alternatively, a rx*ra-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.LBx:   a rx-by-rx-by-nt tensor giving lower bounds
%                           on sensor-to-sensor distance.
%                           alternatively, a rx*rx-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.LBa:   a rx-by-ra-by-nt tensor giving lower bounds
%                           on sensor-to-anchor distance.
%                           alternatively, a rx*ra-by-nt sparse matrix can be 
%                           used to save space
%
%               cons.TC:    the temporal constraint matrix (size(cons.T,2) = nt)
%                           |reshape(X,rx*d,nt)*cons.TC'|_F^2 is included in the 
%                           objective to ensure temporal smoothness
%                           You need to properly scale TC to achieve the desired
%                           weight for the smoothness penalty term
%
%               cons.LRank: location rank constraint on reshape(X,rx*d,nt)
%                           (default: 0)
%
%               cons.gamma: weight for location rank constraint.  That is,
%                           we include a penalty term 
%                           gamma*|reshape(X,rx*d,nt) - U*V'|_F^2 
%                           in the objective. (default: gamma = 1e-2)
%
%               cons.t_inc: Factor of fitting error/temporal error by which 
%                           we want to increase the tweight.
%
%             Note: such constraints need not be strict -- we are minimizing
%             the squared sum of constraint violations
%
%        X0:  an initial solution (which need not be feasible)
%             [ size: rx-by-d-by-nt.  default: empty <==> randomly chosen ]
%
%   Display:  level of display [ off | final | (iter) | full | excessive ]
%
%    lambda:  include a regularization term lambda*sum((x-x0).^2) in the objective
%             when calling LowRankSDP(...).  when x0 is not specified, simply
%             regularize towards the center of anchor nodes (i.e. mean(A,1))
%
% Output:
%
%         X:  the estimated sensor location (size: rx-by-d-by-nt)
%
% References:
% file:        LocalizeLowRankSDP_3D.m

  if nargin < 7, X0 = [];          end
  if nargin < 8, Display = 'iter'; end
  if nargin < 9, lambda = 0;       end
  
  EQx   = sparse(rx*rx,nt);
  EWx   = sparse(rx*rx,nt); % add by jie  
  EQa   = sparse(rx*ra,nt); 
  EWa   = sparse(rx*ra,nt); % add by jie
  UBx   = sparse(rx*rx,nt); 
  UBa   = sparse(rx*ra,nt); 
  LBx   = sparse(rx*rx,nt); 
  LBa   = sparse(rx*ra,nt); 
  TC    = sparse(0,nt);
  lrank = 0;
  gamma = 1e-2;
  if (isfield(cons,'EQx')) 
    if (ndims(cons.EQx) == 3)
      EQx = sparse(reshape(cons.EQx,rx*rx,nt));
    else
      EQx = sparse(cons.EQx);
    end
  end
  if (isfield(cons,'EWx')) % add by jie
    if (ndims(cons.EWx) == 3)
      EWx = sparse(reshape(cons.EWx,rx*rx,nt));
    else
      EWx = sparse(cons.EWx);
    end
  end
  if (isfield(cons,'EQa'))
    if (ndims(cons.EQa) == 3)
      EQa = sparse(reshape(cons.EQa,rx*ra,nt));
    else
      EQa = sparse(cons.EQa);
    end
  end
  if (isfield(cons,'EWa'))
    if (ndims(cons.EWa) == 3)
      EWa = sparse(reshape(cons.EWa,rx*ra,nt));
    else
      EWa = sparse(cons.EWa);
    end
  end
  if (isfield(cons,'UBx')) 
    if (ndims(cons.UBx) == 3)
      UBx = sparse(reshape(cons.UBx,rx*rx,nt));
    else
      UBx = sparse(cons.UBx);
    end
  end
  if (isfield(cons,'UBa'))
    if (ndims(cons.UBa) == 3)
      UBa = sparse(reshape(cons.UBa,rx*ra,nt));
    else
      UBa = sparse(cons.UBa);
    end
  end
  if (isfield(cons,'LBx')) 
    if (ndims(cons.LBx) == 3)
      LBx = sparse(reshape(cons.LBx,rx*rx,nt));
    else
      LBx = sparse(cons.LBx);
    end
  end
  if (isfield(cons,'LBa'))
    if (ndims(cons.LBa) == 3)
      LBa = sparse(reshape(cons.LBa,rx*ra,nt));
    else
      LBa = sparse(cons.LBa);
    end
  end
  if (isfield(cons,'TC'))
    TC = sparse(cons.TC);
  end
  if (isfield(cons,'LRank'))
    lrank = cons.LRank;
  end
  if (isfield(cons,'gamma'))
    gamma = cons.gamma;
  end
  if (gamma == 0)
    lrank = 0;
  end
  if (isfield(cons,'t_inc'))
    t_inc = cons.t_inc;
  else
    t_inc = 1;
  end
  
  % make EQx, UBx, LBx symmetric
  [ij t v] = find(EQx);
  [i j] = ind2sub([rx rx],ij);
  ji = sub2ind([rx rx],j,i);
  EQx = max(EQx,sparse(ji,t,v,rx*rx,nt));

  [ij t v] = find(EWx); %add by jie
  [i j] = ind2sub([rx rx],ij);
  ji = sub2ind([rx rx],j,i);
  EWx = max(EWx,sparse(ji,t,v,rx*rx,nt));
  
  [ij t v] = find(UBx);
  [i j] = ind2sub([rx rx],ij);
  ji = sub2ind([rx rx],j,i);
  UBx = max(UBx,sparse(ji,t,v,rx*rx,nt));

  [ij t v] = find(LBx);
  [i j] = ind2sub([rx rx],ij);
  ji = sub2ind([rx rx],j,i);
  LBx = max(LBx,sparse(ji,t,v,rx*rx,nt));

  % reset the upper triagle
  [ij t v] = find(EQx);
  [i j] = ind2sub([rx rx],ij);
  ind = find(i < j);
  EQx = sparse(ij(ind),t(ind),v(ind),rx*rx,nt);
  
  [ij t v] = find(EWx); % add by jie
  [i j] = ind2sub([rx rx],ij);
  ind = find(i < j);
  EWx = sparse(ij(ind),t(ind),v(ind),rx*rx,nt);

  [ij t v] = find(UBx);
  [i j] = ind2sub([rx rx],ij);
  ind = find(i < j);
  UBx = sparse(ij(ind),t(ind),v(ind),rx*rx,nt);

  [ij t v] = find(LBx);
  [i j] = ind2sub([rx rx],ij);
  ind = find(i < j);
  LBx = sparse(ij(ind),t(ind),v(ind),rx*rx,nt);

  % length of equality constraints
  m_eqx = nnz(EQx);
  m_eqa = nnz(EQa);

  % length of upper bound constraints
  m_ubx = nnz(UBx);
  m_uba = nnz(UBa);

  % length of lower bound constraints
  m_lbx = nnz(LBx);
  m_lba = nnz(LBa);
  
  % number of block diagnonal matrices
  n = nt;
  
  % total number of constraints 
  m = m_eqx + m_eqa + m_ubx + m_uba + m_lbx + m_lba; 

  % set the size of each block matrix
  x_size = repmat([rx d],n,1);

  % length of location rank constraint
  if (lrank > 0)
    n    = n + 1;
    m_lr = rx*d*nt;
    m    = m + m_lr;
    x_size = [x_size; (rx*d+nt) lrank];
  end
  
  nX = sum(x_size(:,1).*x_size(:,2));
  rX = sum(x_size(:,1));
  cX = sum(x_size(:,2));

  % re-center X0 and A
  if (isempty(X0))
    X_avg = mean(A,1);
    x0    = randn(nX,1);
    w0    = zeros(nX,1); % add by jie
    x_reg = zeros(nX,1);
  else
    X_avg = mean([X0;A],1);
    X0    = X0 - repmat(X_avg,[rx 1 1]);
    if (lrank == 0)
      x0 = X0(:);
      w0 = W0(:); % add by jie
    else
%       [U0,S0,V0] = svds(reshape(X0,[],nt),lrank); %add by cj, Make
%       sure the U0, S0 and V0 have enough number of rows and columns
      currentX = reshape(X0,[],nt);
      [U0,S0,V0] = svds(currentX,lrank); 
      [oriS0row oriS0col] = size(S0);
      if oriS0row < lrank || oriS0col < lrank
          oriS0 = S0;
          oriU0 = U0;
          oriV0 = V0;
          S0 = zeros(lrank);
          U0 = zeros(size(currentX,1),lrank);
          V0 = zeros(size(currentX,2),lrank);
          S0(1:oriS0row,1:oriS0col) = oriS0;
          U0(:,1:size(oriU0,2)) = oriU0;
          V0(:,1:size(oriV0,2)) = oriV0;
      end
      
      S0 = sparse(S0);
      U0 = U0*sqrt(S0);
      V0 = V0*sqrt(S0);
      x0 = [X0(:); U0(:); V0(:)];
      sizeOfU0V0 = length([U0(:); V0(:)]);
      w0 = [W0(:); zeros(sizeOfU0V0,1)];
    end
    x_reg = x0;
  end
  if (isempty(X_avg))
    X_avg = zeros(1,d,nt);
  end  
  A = A - repmat(X_avg,[ra 1 1]);

  A2 = sum(A.^2,2);

  % initialize the set of constraints
  CT = zeros(m,1);
  H  = cell(1,m);
  G  = sparse(m,nX);
  b  = zeros(m,1);

  % number of constraints
  nc = 0;

  % objective:            // any (i,j) in EQx
  %   minimize [ X(i,:,t)*X(i,:,t)' + X(j,:,t)*X(j,:,t)' - 2*X(i,:,t)*X(j,:,t)' - EQx(i,j,t) ]^2
  [IJ T V] = find(EQx);
  [I J] = ind2sub([rx rx],IJ);
  for c = 1:m_eqx
    nc = nc + 1;
    i = I(c);
    j = J(c);
    t = T(c);
    v = V(c);
    w = sqrt(EWx(IJ(c),T(c))); % add by jie
    
    h_bas = (t-1)*rx;
    H{nc} = sparse([i j i j]+h_bas, [i j j i]+h_bas, [1 1 -1 -1]*w, rX, rX);
    b(nc) = v*w;
  end

  % objective:              // any (i,k) in EQa
  %   minimize [ X(i,:,t)*X(i,:,t)' + A2(k,1,t) - 2*A(k,:,t)*X(i,:,t)' - EQa(i,k,t) ]^2
  [IK T V] = find(EQa);
  [I K] = ind2sub([rx ra],IK);
  for c = 1:m_eqa
    nc = nc + 1;
    i = I(c);
    k = K(c);
    t = T(c);
    v = V(c);
    w = sqrt(EWa(IK(c),T(c)));

    h_bas = (t-1)*rx;
    g_bas = (t-1)*rx*d;
    H{nc} = sparse(i+h_bas, i+h_bas, 1*w, rX, rX);
    G(nc,sub2ind([rx d],i*ones(1,d),1:d)+g_bas) = -2*A(k,:,t)*w;
    b(nc) = (v - A2(k,1,t))*w;
  end
  
  % objective:            // any (i,j) in UBx
  %   minimize max[ 0, X(i,:,t)*X(i,:,t)' + X(j,:,t)*X(j,:,t)' - 2*X(i,:,t)*X(j,:,t)' - UBx(i,j,t) ]^2
  [IJ T V] = find(UBx);
  [I J] = ind2sub([rx rx],IJ);
  for c = 1:m_ubx
    nc = nc + 1;
    i = I(c);
    j = J(c);
    t = T(c);
    v = V(c);

    h_bas = (t-1)*rx;
    H{nc} = sparse([i j i j]+h_bas,  [i j j i]+h_bas, [1 1 -1 -1], rX, rX);
    b(nc) = v;
    CT(nc) = 1;
  end

  % objective:              // any (i,k) in UBa
  %   minimize max[ 0, X(i,:,t)*X(i,:,t)' + A2(k,1,t) - 2*A(k,:,t)*X(i,:,t)' - UBa(i,k,t) ]^2
  [IK T V] = find(UBa);
  [I K] = ind2sub([rx ra],IK);
  for c = 1:m_uba
    nc = nc + 1;
    i = I(c);
    k = K(c);
    t = T(c);
    v = V(c);

    h_bas = (t-1)*rx;
    g_bas = (t-1)*rx*d;
    H{nc} = sparse(i+h_bas, i+h_bas, 1, rX, rX);
    G(nc,sub2ind([rx d],i*ones(1,d),1:d)+g_bas) = -2*A(k,:,t);    
    b(nc) = v - A2(k,1,t);
    CT(nc) = 1;
  end
    
  % objective:            // any (i,j) in LBx
  %   minimize min[ 0, X(i,:,t)*X(i,:,t)' + X(j,:,t)*X(j,:,t)' - 2*X(i,:,t)*X(j,:,t)' - LBx(i,j,t) ]^2
  [IJ T V] = find(LBx);
  [I J] = ind2sub([rx rx],IJ);
  for c = 1:m_lbx
    nc = nc + 1;
    i = I(c);
    j = J(c);
    t = T(c);
    v = V(c);

    h_bas = (t-1)*rx;
    H{nc} = sparse([i j i j]+h_bas, [i j j i]+h_bas, [1 1 -1 -1], rX, rX);
    b(nc) = v;
    CT(nc) = -1;
  end

  % objective:              // any (i,k) in LBa
  %   minimize min[ 0, X(i,:,t)*X(i,:,t)' + A2(k,1,t) - 2*A(k,:,t)*X(i,:,t)' - LBa(i,k,t) ]^2
  [IK T V] = find(LBa);
  [I K] = ind2sub([rx ra],IK);
  for c = 1:m_lba
    nc = nc + 1;
    i = I(c);
    k = K(c);
    t = T(c);
    v = V(c);

    h_bas = (t-1)*rx;
    g_bas = (t-1)*rx*d;
    H{nc} = sparse(i+h_bas, i+h_bas, 1, rX, rX);
    G(nc,sub2ind([rx d],i*ones(1,d),1:d)+g_bas) = -2*A(k,:,t);    
    b(nc) = v - A2(k,1,t);
    CT(nc) = -1;
  end

  %
  % objective:
  %
  %   minimize gamma*|reshape(X,rx*d,nt)-U*V'|_F^2
  %
  % that is:
  %
  %   minimize sum_c |sqrt(gamma)*x(c) - sqrt(gamma)/2*U(i,:)*V(t,:)'
  %                                    - sqrt(gamma)/2*V(t,:)*U(i,:)'|_F^2 
  %   (where [i t] = ind2sub([rx*d nt], c))
  %
  if (lrank > 0)
    [I,T] = ind2sub([rx*d nt], 1:m_lr);
    % convert T row index in V into row index in [U;V]
    T     = T + rx*d;
    h_bas = nt*rx;
    gam2  = 0.5*sqrt(gamma);
    for c = 1:m_lr
      nc = nc + 1;
      i = I(c);
      t = T(c);
      H{nc}   = sparse([i t]+h_bas, [t i]+h_bas, [-gam2 -gam2], rX, rX);
      G(nc,c) = 2*gam2;
    end
  end
  
  %
  % objective:
  %    minimize |reshape(X,rx*d,nt)*TC'|_F^2
  %
  % Below we find matrix R such that |R*x|_2^2 = |reshape(X,rx*d,nt)*TC'|_F^2 
  %
  % Let xt = reshape(reshape(X,rx*d,nt)',[],1).  
  % Let R0 = blkdiag(TC_rep{:}) where TC_rep = repmat({TC},1,rx*d).
  % Let ind_xt2x = reshape(reshape(1:(rx*d*nt),nt,rx*d)',[],1)
  % Let R = R0(:,ind_xt2x)
  % We then have |R*x|_2^2 = |TC*reshape(X,rx*d,nt)'|_F^2.
  %
  for iter = 1:2
%     iter
    ind_xt2x = reshape(reshape(1:(rx*d*nt),rx*d,nt)',[],1);
    [r_tc, c_tc] = size(TC);
    [I,J,V] = find(TC);
    n_tc = length(I);
    K = reshape(repmat(0:(rx*d-1),n_tc,1),[],1);
    I = repmat(I(:),rx*d,1) + r_tc*K;
    J = repmat(J(:),rx*d,1) + c_tc*K;
    V = repmat(V(:),rx*d,1);
    R = sparse(I, ind_xt2x(J), V, r_tc*rx*d, nX);

    % form P and q
%     P = lambda*speye(nX,nX) + R'*R;
%     q = -2*lambda*x_reg;
    P = lambda*sparse(diag(w0)) + R'*R; % revised by jie
    q = -2*lambda*(x_reg.*w0);
    
    % garbage collection
    clear I J K V R;
  
    % solve the low-rank SDP
    [x,obj,flag,obj_lst] = LowRankSDP(CT,H,G,b,x_size,P,q,x0,Display);
    ratio = obj_lst(1)/obj_lst(2)*t_inc;
%     disp('Fitting Error')
%     obj_lst(1)
%     disp('Temporal Error')
%     obj_lst(2)
%     ratio
    if (isnan(ratio) || ratio <= 1 || isinf(ratio))
      break
    else
      ratio = max(min(ratio,10),1/10);
      TC = TC*ratio;
    end
  end

  % turn x back into tensor format
  X = reshape(x(1:(rx*d*nt)),[rx d nt]);
  
  % go back to the original center
  X = X + repmat(X_avg, [rx 1 1]);

