function [x,obj,flag] = FMin(funObj,x0,solver)
%
% [x,obj,flag] = FMIN(funObj,x0,solver,options) is a wrapper around solvers
% like minFunc and fmincon. 
%
% file:        FMin.m
% directory:   /u/yzhang/Anomaly/Matlab/
% created:     Wed Nov 18 2009 
% author:      Yin Zhang 
% email:       yzhang@cs.utexas.edu
%

  SOLVER_MINFUNC = 1;
  SOLVER_FMINCON = 2;
  
  if nargin < 3 
    if (~isempty(which('minFunc')))
      solver = SOLVER_MINFUNC;
    else
      solver = SOLVER_FMINCON;
    end
  end

  switch(solver)
   
   case SOLVER_MINFUNC
    options.Display = 'off';     % display level
    options.Method = 'lbfgs';    % method
    options.Corr = 5;            % number of past gradients in L-BFGS
    options.MaxIter = 1000;      % max iterations
    options.MaxFunEvals = inf;   % max function evaluations
    options.TolX = 1e-10;        % tolerance on function/parameter change
    options.TolFun = 1e-6;       % tolerance on first-order optimality
    [x,obj,flag] = minFunc(funObj,x0,options);
   
   case SOLVER_FMINCON 
    options = ...
        optimset('Display', 'off', ...
                 'Algorithm', 'interior-point', ...
                 'GradObj', 'on', ...
                 'Hessian', {'lbfgs',5}, ...
                 'MaxIter', 1000, ...
                 'MaxFunEvals', inf, ...
                 'TolX', 1e-10, ...
                 'TolFun', 1e-6);
    [x,obj,flag] = fmincon(funObj,x0,[],[],[],[],-inf,[],[],options);  
  
  end
