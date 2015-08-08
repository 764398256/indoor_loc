function [A_hat E_hat B_hat iter] = inexact_alm_rpca_stable_fix_rank(D, xi, lambda, tol, maxIter,r)

% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA.
%
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% Initialize A,E,Y,u
% while ~converged 
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end
%
% Minming Chen, October 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

% addpath PROPACK;

[m n] = size(D);

%%%%%%%%%%add by jie, generate stable matrix
temp1 = diag(ones(1,n));
temp2 = diag(-2*ones(1,n-1),1);
temp3 = diag(ones(1,n-2),2);
stable = temp1 + temp2 + temp3;
stable(end-1:end,:) = 0;
%%%%%%%%%%add by jie, generate stable matrix

xi = 1;
p1left = sqrt(xi) * stable';

if nargin < 2
    xi = 1;
end

if nargin < 3
    lambda = 1 / sqrt(m);
end

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;
W = zeros( m,n);

A_hat = zeros( m, n);
E_hat = zeros( m, n);
B_hat = zeros( m, n);
mu = 1.25/norm_two; % this one can be tuned
nu = mu;
mu_bar = mu * 1e7;
nu_bar = mu_bar;
rho = 1.5;          % this one can be tuned
rhonu = 1.1;
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = r; % revised by jie, let sv starting from r
while ~converged       
    iter = iter + 1;
    
    temp_T = D - A_hat + (1/mu)*Y;
    E_hat = max(temp_T - lambda/mu, 0);
    E_hat = E_hat+min(temp_T + lambda/mu, 0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%update A%%%%%%%%%%%%%%
    R = 1/(mu+nu) * (mu*(D - E_hat + 1/mu*Y)+nu*(B_hat - 1/nu*W));
    
    if choosvd(n, sv) == 1
        [U S V] = lansvd(R, sv, 'L');
    else
        [U S V] = svd(R, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/(mu+nu)));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    if svp > r  % add by jie
        sv = r; 
        svp = r;    
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/(mu+nu)) * V(:, 1:svp)';
    
    p1 = [p1left sqrt(nu) * eye(n)];
    p2 = [zeros(m,n) sqrt(nu)*(A_hat + 1/nu * W)];
    B_hat = p2 * p1' * (p1*p1')^(-1);

    total_svd = total_svd + 1;
    
    Z = D - A_hat - E_hat;
    
    Y = Y + mu*Z;
    
    Znonnegative = A_hat - B_hat;
    W = W + nu * Znonnegative;
    
    mu = min(mu*rho, mu_bar);
    nu = min(nu*rhonu, nu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    nonnegativeStopCriterion = norm(Znonnegative,'fro') / d_norm;
    if (stopCriterion < tol) && (nonnegativeStopCriterion < tol)
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
%         disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
%             ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
%             ' stopCriterion ' num2str(stopCriterion)]);
    end    
    
    if ~converged && iter >= maxIter
%         disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end
