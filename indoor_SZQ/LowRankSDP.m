function [x,obj,flag,obj_lst] = LowRankSDP(CT,H,A,b,blk_size,P,q,x0,Display)
%
% [x,obj,flag] = LOWRANKSDP(H,A,b,blk_size,C,f,x0,Display) performs 
% unconstrained optimization of the following quartic objective function:
%
%   x'*P*x + q'*x + sum_{i=1:m} PosNeg( CT(i), trace(H{i}*X*X') + A(i,:)*x - b(i) )^2
%
% where
%
%     X = blkdiag(X1,X2,...,Xn) is a block diagonal matrix (size: rx-by-cx)
%     x = [X1(:); X2(:); ...; Xn(:)] is the concatenation of all Xi's (length: nx)
%
%     PosNeg(ct,v) =        v             if ct =  0
%                    max(0,v)             if ct =  1
%                    min(0,v)             if ct = -1
%
% Equivalently, 
%
%     PosNeg(ct,v) =        0             if sign(v) == -ct
%                           v             otherwise
%
% Input:
%
%        CT:   a m-by-1 vector that gives the type of each constraint.  CT(i) \in {-1, 0, 1}
%
%         H:   a 1-by-m cell of symmetric matrices H{i} specifying the quadratic 
%              component for each squared term (each has size: rx*rx)
%
%         A:   a m-by-rx matrix specifying the linear component for each
%              squared term 
%
%         b:   a m-by-1 vector specifying the constant component for each
%              squared term
%
%  blk_size:   a n-by-2 matrix, where blk_size(j,1:2) = size(Xj).  Thus 
%              nx = sum(blk_size(:,1).*blk_size(:,2))
%
%         P:   a nx-by-nx matrix specifying a quadratic term in the objective
%              (P can also be a scalar, in which case x'*P*x = sum(x.*(P*x)))
%
%         q:   a nx-by-1 vector specifying a linear term in the objective
%              (q can be a scalar, in which case q'*x = sum(q.*x))
%
%        x0:   an initial solution (default: x0 = randn(nx,1))
%
%   Display:   level of display [ off | final | (iter) | full | excessive ]
%
% Output:
%
%       x:     the solution
%
%     obj:     the final objective
%
%    flag:     flag returned by the optimizer minFunc()
%
% Reference:
%
%    S. Burer and R.D.C. Monteiro. 
%    A Nonlinear Programming Algorithm for Solving Semidefinite Programs 
%    Via Low-Rank Factorization. Mathematical Programming (series B), 
%    95(2):329-357, 2003.
%
%    S. Burer and C. Choi.
%    Computational Enhancements in Low-Rank Semidefinite Programming. 
%    Optimization Methods and Software, 21(3):493-512, 2006.
%    http://dollar.biz.uiowa.edu/~sburer/papers/013-lowrank3.pdf
%
% file:        LowRankSDP.m

  if (nargin < 6), x0 = [];          end
  if (nargin < 7), Display = 'iter'; end

  % get length information
  m = length(b);
  n = size(blk_size,1);

  % get the negation of CT
  NegCT = -sign(CT(:));
  
  % create index for converting x between vector and block diagonal format
  % we have:   X = sparse(Ix,Jx,x,rx,cx)    
  x_cell = cell(1,n);
  for i = 1:n
    % need to make x_cell{i} sparse so that blkdiag(x_cell{:}) is sparse
    x_cell{i} = sparse(ones(blk_size(i,1),blk_size(i,2)));
  end
  x_mat = blkdiag(x_cell{:});
  clear x_cell;
  [rx,cx] = size(x_mat);
  [Ix,Jx] = find(x_mat);
  nx = length(Ix);
  clear x_mat;
  
  %
  % Let xt = [X1'(:); X2'(:); ...; Xn'(:)] be a nx-by-1 column vector formed 
  % by collapsing all the block diagonal matrices in X'. Then each row of X is 
  % stored consecutively in xt.  Create indices for finding these rows in xt.
  %
  %    x_beg(i):   beginning position of non-zero elements of X(i,:) in xt
  %    x_end(i):   ending position of non-zero elements of X(i,:) in xt
  %    x_len(i):   x_end(i)-x_beg(i)+1
  %    x_blk(i):   index j for block matrix Xj that contains all non-zero 
  %                    elements in this row
  %
  x_beg = zeros(rx,1);
  x_end = zeros(rx,1);
  x_len = zeros(rx,1);
  x_blk = zeros(rx,1);
  row_bas = 0;   % base row index
  xt_bas  = 0;   % base index for elements in xt
  for i = 1:n
    ri = blk_size(i,1);
    ci = blk_size(i,2);
    x_beg((row_bas+1):(row_bas+ri)) = xt_bas+(1:ci:(ci*ri));
    x_end((row_bas+1):(row_bas+ri)) = xt_bas+(ci:ci:(ci*ri));
    x_len((row_bas+1):(row_bas+ri)) = ci;
    x_blk((row_bas+1):(row_bas+ri)) = i;
    row_bas = row_bas + ri;
    xt_bas = xt_bas + ri*ci;
  end

  % set up lookup table for converting index in xt to x
  ind_xt2x = nonzeros(sparse(Jx,Ix,1:nx,cx,rx));

  %
  % Next: we need to do some preprocessing so that we can 
  % efficiently compute H*X inside calcObjGrad and exactLineSearch
  %

  % 
  % The following code is functionally equivalent to:
  %
  %    Hcat = vertcat(H{:});
  %    [Ih0,Jh,Vh] = find(Hcat);
  %
  % However, vertcat() a large number of large sparse matrices
  % can be very expensive and much slower than the following code
  % (despite its use of slow matlab for-loops)
  %
  Ih0_cell = cell(1,m);
  Jh_cell  = cell(1,m);
  Vh_cell  = cell(1,m);
  for i = 1:m
    [Ih0_i,Jh_cell{i},Vh_cell{i}] = find(H{i});
    Ih0_cell{i} = Ih0_i + rx*(i-1);
  end
  Ih0 = vertcat(Ih0_cell{:});
  Jh  = vertcat(Jh_cell{:});
  Vh  = vertcat(Vh_cell{:});
  clear Ih0_cell Jh_cell Vh_cell;

  %
  % map Ih0 back to 
  %   (1) Ch: constraint index (i.e. i if an element comes from H{i})
  %   (2) Ih: original row index in an al H{i}
  %
  Ih = mod(Ih0-1,rx)+1;
  Ch = (Ih0-Ih)/rx + 1;

  %
  % Next, for each (i,j,v) in vertcat(H{:}), find the indices for all elements
  % in x(:) that needs to be multiplied with (i,j,v).  Specifically, in order to
  % compute Hcat*X, we should multiply element (i,j,v) with j-th row of X, whose
  % elements have indices ind_xt2x(x_beg(j):x_end(j)) = ind_xt2x(jk_vec).  The 
  % product should be added to i-th row of X, whose elements have indices
  % ind_xt2x(x_beg(i):x_end(i)) = ind_xt2x(ik_vec)
  %
  % since we are computing trace(H{i}*X*X'), and X is block diagonal, we only
  % need to consider those i,j s.t. x_blk(i) == x_blk(j), which is given by
  % ind = find(x_blk(Ih) == x_blk(Jh))';
  %
  % After everything is done, we can compute Hcat*X and convert the result into
  % a m-by-nx matrix simply as:
  %
  %    HX = sparse(H_st.C, H_st.IK, H_st.V .* x(H_st.JK), m, nx);
  %
  % NOTE: sparse(...) function in matlab ensures that when there are multiple
  % values associated with the same entry, all these values are summed up, which
  % is precisely what we want.
  %
  ind = find(x_blk(Ih) == x_blk(Jh))';  % xxx: ind needs to be a row vector
  vec_len = sum(x_len(Ih(ind)));
  h_vec = zeros(vec_len,1);
  jk_vec = zeros(vec_len,1);
  ik_vec = zeros(vec_len,1);
  bas = 0;
  for h = ind
    i = Ih(h);
    j = Jh(h);
    len = x_len(j);
    h_vec((bas+1):(bas+len)) = h;
    ik_vec((bas+1):(bas+len)) = x_beg(i):x_end(i);
    jk_vec((bas+1):(bas+len)) = x_beg(j):x_end(j);
    bas = bas + len;
  end
  
  %
  % NOTE: [C IK] may have duplicate rows.  So sparse(C,IK,1) may have fewer 
  % non-zero entries than length(C).  To handle this inconvenience, we
  % use H_st.HX0 to provide non-zero positions of sparse(C,IK,1), and use
  % H_st.HX0_idx(k) to indicate which non-zero of H_st.HX0 entry 
  % [C(k) IK(k)] gets mapped into.
  %
  C  = Ch(h_vec);
  IK = ind_xt2x(ik_vec);
  % XXX: must use [IK C] instead of [C IK] => sorted by column first
  [IKC_u,I,J] = unique([IK C], 'rows'); 
  H_st.HX0 = sparse(IKC_u(:,2), IKC_u(:,1), 1, m, nx);
  H_st.HX0_nz = nnz(H_st.HX0);
  H_st.HX0_idx = uint32(J);
  H_st.V = Vh(h_vec);   
  H_st.JK = uint32(ind_xt2x(jk_vec));
  
  % garbage collection
  clear h_vec ik_vec jk_vec Ih0 Ih Jh Ch Vh ind_xt2x C IK IKC_u I J ind;
  clear x_beg x_end x_len x_blk;
  
  % 
  % We only use calcObjGrad for initial computation of obj + grad.
  % Later, we use exactLineSearch to return the new x and new g.
  %
  funObj = @(x) calcObjGrad(x,m,NegCT,H_st,A,b,nx,P,q);

  % set the initial solution x0 if user didn't provide one
  if (isempty(x0))
    x0 = randn(nx,1);
  end

  % make sure x0 is full
  if (issparse(x0))
    x0 = full(x0);
  end
  
  [x,obj,flag] = FMin(funObj,x0);
  [obj,grad,obj_lst] = funObj(x);
  
%
% Objective: minimize
%
%  x'*P*x + q'*x + sum_{i=1:m} PosNeg( -NegCT(i), trace(H{i}*X*X') + A(i,:)*x - b(i) )^2
%
function [obj,grad,obj_lst] = calcObjGrad(x,m,NegCT,H_st,A,b,nx,P,q)

  hx = accumvec(H_st.HX0_idx, H_st.V.*distrvec(H_st.JK,x), H_st.HX0_nz);
  HX = spsetnz(H_st.HX0, hx);
  HXxAxb = HX*x + A*x - b;
  HXxAxb(sign(HXxAxb)==NegCT) = 0;

  Px   = P*x;
  obj  = sum(HXxAxb.^2) + sum(Px.*x) + sum(q.*x);
  grad = 4*(HX'*HXxAxb) + 2*(A'*HXxAxb) + 2*Px + q;
  
  if (nargout == 3)
    obj_lst(1) = sum(HXxAxb.^2);
    obj_lst(2) = sum(Px.*x);
    obj_lst(3) = sum(q.*x);
  end


