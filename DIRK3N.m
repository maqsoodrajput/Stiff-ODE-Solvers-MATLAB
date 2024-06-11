function [t,y] = DIRK3N(f,J,t,y0)
   %Two Stage Order 3 Diagonally Implicit Runge Kutta (DIRK) Scheme
   %Scheme: 
   %y1=   un + (3+sqrt(3))/6*dt*f(y1);
   %y2=   un - sqrt(3)/3*dt*f(y1) +  (3+sqrt(3))/6*dt*f(y2);
   %unp1= un + .5*dt*f(y1) +.5*dt*f(y2)
   %Tutorial:
   %First Stage implicit Solve (top of Page 106 with gammas=2 p=3)
   %F(x) x-(un+(3+sqrt(3))/6*dt*f(x))
   %Second Stage implicit Solve
   %F(x)   x-(un-sqrt(3)/3*dt*f(y1) +  (3+sqrt(3))/6*dt*f(x)))
   %Then, Update time step
   a = zeros(2,2);
   a(1,1) = (3 + sqrt(3))/6;
   a(1,2) = 0.0;
   a(2,1) = -sqrt(3)/3;
   a(2,2) = (3 + sqrt(3))/6;
   b = zeros(2,1);
   b(1,1) = 1/2;
   b(2,1) = 1/2;
   tol = 1e-13;  
   %options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
   options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
   y = zeros(length(y0),length(t));
   y(:,1) = y0;
   dt = t(2) - t(1);
   %p=3;
   tn = t(1);
   un.p = y0;
   yn=un.p;
   F1 = @(k1)  k1 - f(yn + a(1,1)*dt.*k1)
   F2 = @(k2)  k2 - f(yn + a(2,1)*dt.*k1 + a(2,2)*dt.*k2);
   gnparams = struct('maxit',1000,'toler',1.0e-4,'lsmethod','chol')
   for i =2:length(t)
	%[y1,tmp1] = newton( @(x)  x - ( un + a(1,1)*dt*f(x)),un,J);
   	%[y2,tmp2] = newton( @(x)  x - ( un - a(2,1)*dt*f(y1) +  a(2,2)*dt*f(x)),un,J);
        [inform,Y1] = GaussN(f,J,F1,un.p,gnparams);
        [inform,Y2] = GaussN(f,J,F2,un.p,gnparams);
   	y(:,i)= un.p + b(1,1)*dt*f(Y1) + b(2,1)*dt*f(Y2);
	un = y(:,i);            %Updating the initial guess	
   end
end


function [inform, x] = GaussN(fun,Jac, resid, x, gnparams)
%  Implements the Gauss-Newton algorithm for solving non-linear least squares
%  problems.
%
%  Input:
%
%    fun      - a pointer to a function
%    resid    - a pointer to residual function
%    x        - the following structure:
%               * x.p - the starting point values
%               * x.f - the function value of x.p
%               * x.g - the gradient value of x.p
%               * x.h - the Hessian value of x.p
%               with only x.p guaranteed to be set.
%    gnparams - the following structure, as an example:
%
%                 gnparams = struct('maxit',1000,'toler',1.0e-4,
%                                  'lsmethod','chol');
%              
%              Note that 'lsmethod' can be 'chol', 'qr', or 'svd'.
%
%  Output:
%
%    inform - structure containing two fields:
%      * inform.status - 1 if gradient tolerance was
%                        achieved;
%                        0 if not.
%      * inform.iter   - the number of steps taken
%    x      - the solution structure, with the solution point along with
%             function, gradient, and Hessian evaluations thereof.

%  Number of function, gradient, and Hessian evaluations.
global numf numg numh
numf = 0;
numg = 0;
numh = 0;

%  Populate local caching of gnparams parameters.
toler = gnparams.toler;  % Set gradient tolerance.
maxit = gnparams.maxit;  % Set maximum number of allowed iterations.
lsmethod = gnparams.lsmethod;  % Set method to find search direction.

%  Initialize parameter structure for StepSizeSW function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
                'stpmax', 1e20, 'maxfev', 10000);

alfa = 1;  % Initial value of alpha.
xc.p = x;  % Set the current point to the initial point, x.p.

for i = 1:maxit
    %  Compute function and gradient at current point.
    %xc.f = feval(fun, xc.p, 1);
    %xc.g = feval(fun, xc.p, 2);
    xc.f = feval(fun,xc.p);
    xc.g = feval(Jac,xc.p);
    %  Check for termination condition: norm of gradient less than toler.
    if norm(xc.g) < toler
        inform.status = 1;  % Indicates success.
        inform.iter = i;  % Number of iterations.
        x.p = xc.p;
        x.f = xc.f;
        x.g = xc.g;
        return;
    end

    %  Calculate residual function value and Jacobian at current point.
    %r = feval(resid, xc.p, 1);
    %J = feval(resid, xc.g, 2);
    r = feval(resid, xc.p);
    J = feval(resid, xc.g);


    %  Determine which method should be used to find the search direction, p.
    switch lsmethod
        case 'chol'  % Use Cholesky factorization of J'*J to get p.
            %  If J'*J is positive-definite, R upper triangular matrix such
            %  that R'*R = J'*J.
            %  Else, Cholesky failed and flag is a positive integer.
            [R, flag] = chol(J'*J);
            %  If the Cholesky factorization failed, return with fail status.
            if flag ~= 0
                inform.status = 0;  % Update status to failure indicator, 0.
                inform.iter = i;  % Number of iterations.
                x.p = xc.p;
                x.f = xc.f;
                return;
            end
            %  J'Jp = -J'r and R'R = -J'r, hence p = -R\(R'\(J'*r)).
            p = -R \ (R' \ (J'*r));
        case 'qr'  % Use QR factorization of J to get p.
            %  P permutation matrix, Q unitay matrix, and R upper triangular
            %  matrix with diagonal arranged in absolute decreasing order.
            [Q, R, P] = qr(J);
            n = size(J, 2);
            Q1 = Q(1:end, 1:n);
            R = R(1:n, 1:end);
            %  p = argmin ||J*p+r||^2 = solution of R*P'*p + Q1'*r = 0, hence
            %  p = -P' \ (R \ (Q1'*r)) = -P * (R \ (Q1' * r)).
            p = -P * (R \ (Q1'*r));
        case 'svd'  % Use SVD factorization of J to get p.
            %  U and V unitary matrices, S diagonal matrix.
            [U, S, V] = svd(full(J));
            n = size(J, 2);
            U1 = U(1:end, 1:n);
            S = S(1:n, 1:n);
            %  Since the Moore-Penrose inverse is pinv(J) = V*inv(S)*U1' and 
            %  p = -pinv(J)*r, we have that p = -V*inv(S)*U1'*r.
            p = -V * inv(S) * U1' * r;
    end
    
    %  Get step size that satisfies strong Wolfe conditions.
    [alfa, x] = StepSizeSW(fun, xc, p, alfa, params);
    %  Update current point in p-direction with step size alpha.
    xc.p = xc.p + alfa * p;
end
%  If reached, method failed.
inform.status = 0;  % Update status to failure indicator, 0.
inform.iter = maxit;  % Number of iterations i = maxit at this point.
x.p = xc.p;
%x.f = feval(fun, x.p, 1);
%x.g = feval(fun, x.p, 2);
x.f = feval(fun, x.p);
x.g = feval(J, x.p);

return;  % Return inform and final point x
end





function [x,xhist] = newton(fun,x0,Jac,varargin)
% Solves a system of nonlinear equations  fun(x) = 0 by Newton's method 
% (a.k.a. the Newton-Raphson method) Only real-valued solutions accepted.
% Usage:
%  x = newton(fun,x0)
%    fun is a handle to a function returning a column vector f of function 
%    values and optionally the Jacobian matrix J, where J(i,j) = df(i)/dx(j)
%       f = fun(x) or [f,J] = fun(x) 
%    x0 is the column vector of starting values for the search.
%    x is the solution, if one is found. Otherwise, newton issues a 
%       warning; "No convergence" and returns a vector of NaNs
%    f, x0 and x are all column vectors of length n. J is n by n.
%  x = newton(fun,x0,name,value, ...)
%    allows the user to control the iteration.  
%    The following names are allowed:
%      bounds: n by 2 matrix where bounds(i,1) and bounds(i,2) are 
%        lower and upper bounds, respectively, for variable i.   
%        Use -Inf and/or Inf to indicate no bound.  Default: No bounds.
%      maxiter: Maximum number of iterations,  Default: 50
%      tolx:    Minimum length of last iteration step. Default: 1e-4
%      tolfun:  Minimum value for norm(f). Default: 1e-5
%    Example: x = newton(fun,x0,'maxiter',100,'tolfun',1e-3) 
%  [x,xhist] = newton(fun,x0)
%     xhist contains the search history. x(k,i) is x(i) at iteration k
%   See also: - broyden 
%               (https://se.mathworks.com/matlabcentral/fileexchange/54667)
%             - fsolve (optimization toolbox)
% Are Mjaavatten, November 2020
% Version 2.0: January 2021
%   Option for numeric Jacobian
%   Option for bounds and for modifying iteration parameters
  
  % Parse any name-value pairs that will override the defaults:
  p = inputParser;
  addParameter(p,'tolfun',1e-5);  % Last argument is the default value
  addParameter(p,'tolx',1e-4);
  addParameter(p,'maxiter',50);
  addParameter(p,'bounds',[]);
  parse(p,varargin{:})
  opt = p.Results;    
  % Fields of opt take the default values unless overriden by a 
  % name/value pair
  
  bnds = ~isempty(opt.bounds);  % Bounds on x
  if any(diff(opt.bounds')<0)
    error('newton:bounds',...
      'bounds(i,2)-bounds(i,1) must be > 0 for all i')
  end  
  if bnds && ~all(isreal(opt.bounds))
    error('newton:complexbounds','Bounds cannot be complex numbers')
  end
  %no_jacobian = nargout(fun) < 2;
  x0 = x0(:);  % Make sure x0 is a column vector
  %if no_jacobian
  %  f = fun(x0);
  %  J = jacobi(fun,x0);
  %else
  %  [f,J] = fun(x0);
  %end
  f = fun(x0);
  J = Jac(x0);
  x = x0(:);  % Make sure x0 is a column vector
  if ~(size(x) == size(f))
    error('newton:dimension',...
      'fun must return a column vector of the same size as x0')
  end  
  xhist = zeros(opt.maxiter,length(x));
  for i = 1:opt.maxiter
    xhist(i,:) = x';  % Search history
    %if any(isnan(J(:))) || rcond(J) < 1e-15
    %  warning('newton:singular',...
    %    'Singular jacobian at iteration %d.\n Iteration stopped.',i);
    %  x = NaN*x0;
    %  return
    %end     
    dx = -J\f(:);
    if bnds
      xnew = real(x + dx);
      xnew = min(max(xnew,opt.bounds(:,1)),opt.bounds(:,2));
      dx = xnew-x;
    end
    x = x + dx; 
    if norm(f(:))<opt.tolfun && norm(dx)<opt.tolx
      xhist = xhist(1:i,:);
      if norm(fun(real(x)))> opt.tolfun  % Complex part not negligible
        warning('newton:complex','Converged to complex solution.')
      else
        x = real(x);
      end
      return
    end    
    %if no_jacobian
    %  f = fun(x);
    %  J = jacobi(fun,x);
    %else
    %  [f,J] = fun(x);
    %end
    f = fun(x);
    J = Jac(x);
  end
  if any(abs(xhist(end-1,:)-xhist(end,:))<= 1e-8) 
    s1 = 'Search stuck. ';
    s2 = 'Possibly because Newton step points out of bounded region.';
    if bnds
      warning('newton:stuck',[s1,s2,'\nNo convergence']);
    else
      warning('newton:stuck',[s1,'\nNo convergence']);
    end
    fprintf('Last x:\n')
    for k = 1:numel(x)
      fprintf('%f\n',x(k))
    end
    x = x*NaN;
  else
    warning('newton:noconv','No convergence.')
    x = x*NaN;
  end
end


function [x, xhist] = broyden(fun,x0,Jac,varargin)
% Solves system of nonlinear equations fun(X) = 0 by Broyden's method
% Broyden's method is a quasi-Newton method which updates an approximate 
% Jacobian at each new Newton step. Only real-valued solutions accepted.
%  
% Usage:
%  x = broyden(fun,x0)
%    fun is a handle to a function returning a column vector f
%       of function values
%    x0 is the column vector of starting values for the search.
%    x is the solution, if one is found. Otherwise, broyden issues a 
%       warning; "No convergence" and returns a vector of NaNs
%    f, x0 and x are all column vectors of length n. J is n by n.
%  x = broyden(fun,x0,name,value, ...)
%    allows the user to control the iteration.  
%    The following names are allowed:
%      bounds: n by 2 matrix where bounds(i,1) and bounds(i,2) are 
%        lower and upper bounds, respectively, for variable i.   
%        Use -Inf and/or Inf to indicate no bound.  Default: No bounds.
%      maxiter: Maximum number of iterations,  Default: 50
%      tolx:    Minimum length of last iteration step. Default: 1e-4
%      tolfun:  Minimum value for norm(f). Default: 1e-5
%    Example: x = broyden(fun,x0,'maxiter',100,'tolfun',1e-3) 
%  [x, xhist] = f(fun,x0,...)
%     xhist contains the search history. x(k,i) is x(i) at iteration k
%
% Example:
%  fun = @(x) [x(1)^2 + x(2)^2 - 4; exp(x(1)) + x(2) - 1];
%  x = broyden(fun,[1;1])
%  A second root has x(1) > 0, x(2) < 0.  This may be found by selecting  
%  another initial value x0, or by using bounds.  E.g.:
%  [x,xhist] = broyden(fun,[1;1],'bounds',[0,2;-3,-1]);
%
%  See also: newton
%             (https://se.mathworks.com/matlabcentral/fileexchange/82320)
%             fsolve (in optimization toolbox)
%  Author: Are Mjaavatten, 
% 
%  The update formula for the Jacobian J is taken from p399 in 
%  Matthias Heinkenschloss: Lecture Notes - Numerical Analysis II
%  https://bpb-us-e1.wpmucdn.com/blogs.rice.edu/dist/8/4754/files/2019/01/CAAM-454-554-1lvazxx.pdf
%
% Version 1:   2015-12-29
% Version 1.1: 2016-11-29 Added iteration history
% Version 1.2: 2020-11-04 New handling of convergence parameters
% Version 1.3: 2020-12-19 New handling of bounds
% Version 2.0: January 2021 Improved robustness. Better error handling
% Version 2.01: 2023-05-18; Replaced broken link to Matthias Heinkenschloss
  % Parse any name-value pairs that will override the defaults:
  p = inputParser;
  addParameter(p,'tolfun',1e-5);  % Last argument is the default value
  addParameter(p,'tolx',1e-4);
  addParameter(p,'maxiter',50);
  addParameter(p,'bounds',[]);
  parse(p,varargin{:})
  opt = p.Results;    
  % Fields of opt take the default values unless overriden by a 
  % name/value pair
  
  bnds = ~isempty(opt.bounds);  % Bounds on x
  if any(diff(opt.bounds')<0)
    error('newton:bounds',...
      'bounds(i,2)-bounds(i,1) must be > 0 for all i')
  end  
  if bnds && ~all(isreal(opt.bounds))
    warning('newton:complexbounds','Bounds cannot be complex numbers')
  end							 
  x = x0(:);
  f = fun(x);
  J = Jac(x0);
  if ~(size(x) == size(f))
    error('broyden:dimension',...
      'fun must return a column vector of the same size as x0')
  end  
  %J = jacobi(fun,x);  % Initial Jacobian matrix
  xhist = zeros(opt.maxiter,length(x));
  for i = 1:opt.maxiter
    xhist(i,:) = x';  % Search history			  bnds
    %if any(isnan(J(:))) || rcond(J) < 1e-15
    %  J = jacobi(fun,x); % Try with a full Jacobian
      %if isnan(J) || rcond(J) < 1e-15
      %  warning('broyden:singular',...'
	%	        'Singular jacobian at iteration %d.\n Iteration stopped.',i);
        %x = NaN*x0;
        %return
      %end
    %end
    dx = -J\f(:);
    if bnds
      xnew = real(x + dx);
      xnew = min(max(xnew,opt.bounds(:,1)),opt.bounds(:,2));
      dx = xnew-x;
    end
    x  = x+dx;
    if norm(f(:))<opt.tolfun && norm(dx)<opt.tolx
      xhist = xhist(1:i,:);
      if norm(fun(real(x)))> opt.tolfun  % Complex part not negligible
        warning('Converged to complex solution')
      else
        x = real(x);
      end
      return
    end 
    f  = fun(x);
    J  = J + f*dx'/(dx'*dx);  
  end
  if any(abs(xhist(end-1,:)-xhist(end,:))<= 1e-8) 
    s1 = 'Search stuck. ';
    s2 = 'Possibly because Newton step points out of bounded region.';
    if bnds
      warning('broyden:stuck',[s1,s2,'\nNo convergence']);
    else
      warning('broyden:stuck',[s1,'\nNo convergence']);
    end
    fprintf('Last x:\n')
    for k = 1:numel(x)
      fprintf('%f\n',x(k))
    end
	x = x*NaN;
  else
    warning('broyden:noconv','No convergence.')
  end
end
function J = jacobi(f,x,delta)
% Simple Jacobian matrix estimate using central differences
%  f: function handle
%  x: column vector of independet variables
%  delta: Optional step leghth.  Default: 1e-7*sqrt(norm(x))
%
% See also: John D'Errico's jacobianest.m in
%   https://www.mathworks.com/matlabcentral/fileexchange/13490
%   (Robust and high accuracy)
  if nargin < 3
      delta = 1e-7*sqrt(norm(x));
  end
  y0 = feval(f,x);
  n = length(y0);
  m = length(x);
  J = zeros(n,m);
  X = repmat(x,[1,m]);
  d = delta/2*eye(m);
  for i = 1:m
    J(:,i) = (f(X(:,i)+d(:,i))-f(X(:,i)-d(:,i)))/delta;
  end
end
