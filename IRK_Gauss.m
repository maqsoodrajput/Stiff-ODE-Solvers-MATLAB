function [tvals,y] = IRK_Gauss(fcn,Jac,tvals,y0)

% fcn is the RHS in y'=f(t,y)
% Jac is the Jacobian evaluated at the internal stages Y_i 
% $J_i = \partial f(Y_i)/ \partial y$

N = length(tvals);
y = zeros(length(y0),N);
y(:,1) = y0;

t = tvals(1);
y_guess = y0;
h = tvals(2) - tvals(1);


Tend = tvals(end);
% iterate over time steps
for i = 2:length(tvals)
	 t = tvals(i);
	 y(:,i+1) = IRKsolver(fcn,Jac,y_guess,h,t);
	 y_guess = y(:,i+1);
end
end


function y = IRKsolver(fcn,Jac,y0,h,t)
	v = @(x) double(x);
%	A = [5/36,2/9-sqrt(15)/15,5/36-sqrt(15)/30; ...
%      	5/36+sqrt(15)/24,2/9,5/36-sqrt(15)/24; ...
%      5/36+sqrt(15)/30, 2/9+sqrt(15)/15,5/36];
%   c = [0.5-sqrt(15)/10; 0.5; 0.5+sqrt(15)/10];
%   b = [5/18, 8/18, 5/18];

   A = [v(5)/v(36), v(2)/v(9)-sqrt(v(15))/v(15), v(5)/v(36)-sqrt(v(15))/v(30); ...
      v(5)/v(36)+sqrt(v(15))/v(24), v(2)/v(9), v(5)/v(36)-sqrt(v(15))/v(24); ...
      v(5)/v(36)+sqrt(v(15))/v(30), v(2)/v(9)+sqrt(v(15))/v(15), v(5)/v(36)];
   c = [v(0.5)-sqrt(v(15))/v(10); v(0.5); v(0.5)+sqrt(v(15))/v(10)];
   b = [v(5)/v(18), v(8)/v(18), v(5)/v(18)];


   % q = 6; order of accuracy
   s = 3;
   m = length(y0);
   
   newt_maxit = 100;           % max number of Newton iterations
   newt_ftol  = 1e-13;        % Newton solver residual tolerance
   newt_stol  = 1e-13;        % Newton solver solution tolerance

   % set Newton initial guesses as previous step solution
   z = zeros(s*m,1);
   for i = 0:s-1
      z(i*m+1:(i+1)*m) = y0;
   end
   %Jrhs = Jh(Jac,z,t,h,A,b,c);
   %Frhs = F_IRK(fcn,y0,z,t,h,z,A,b,c);
   %newton's method
   Z = newton(fcn,Jac,y0,z,t,h,A,b,c,newt_ftol, newt_stol, newt_maxit);
   y = Y_IRK(fcn,y0,Z,t,h,A,b,c,s);
end

function y = newton(fcn,Jac,y0,z,t,h,A,b,c,ftol,stol,maxit)
% initialize result, increment vector, residual, statistics
%y = y0;
%s = ones(size(y));

% perform iterations
for i=1:maxit

   % compute residual at current guess
   F = F_IRK(fcn,y0,z,t,h,A,b,c);

   % check residual and increment for stopping
   if ((norm(s,inf) < stol) | (norm(F,inf) < ftol))
      ierr = 0;
      return
   end

   % compute Jacobian
   A = Jh(Jac,z,t,h,A,b,c);

   % perform Newton update
   s = A\F;
   y = y - s;

end

% if we've made it to this point, the Newton iteration did not converge
% ierr = 1;
%fprintf('\nnewton warning: nonconvergence after %i iterations (|F| = %g)\n',maxit,norm(F,inf));

end


function Amat = Jh(Jac,z,t,h,A,b,c)
% Inputs:  z = current guesses for [z1, ..., zs]
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage IRK method, ODE Jacobian function.

	s = length(b);
	zlen = length(z);
	nvar = floor(zlen/s);
	z = reshape(z,nvar,s);
	J = cell(s);
	for is=1:s
		t = t + h*c(is);
		J{is} = Jac(t,z(:,is)); 
	end
	
	Amat = zeros(nvar*s);
	
	for j=1:s
        	for i=1:s
         		Amat(nvar*(i-1)+1:nvar*i,nvar*(j-1)+1:nvar*j) = A(i,j)*J{j};
      		end
   	end
		
	Amat = eye(nvar*s) - h*Amat;
end
		
function F = F_IRK(fcn,y0,z,t,h,A,b,c)
% Inputs:  z = current guesses for [z1, ..., zs]
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for each intermediate
% stage solution, through calling the user-supplied ODE
% right-hand side function.
	
	s = length(b);
	zlen = length(z);
	nvar = floor(zlen/s);
	z = reshape(z,nvar,s);
	% call f at our guesses. This f is from y' =f(t,y)
   	f = zeros(nvar,s);
   	for is=1:s
      		t = t + h*c(is);
      		f(:,is) = fcn(t, z(:,is));
  	end

   	% form the IRK residuals
  	 %    Fs = zs - y_n - h*sum(a(s,j)*fj)
   	F = zeros(nvar,s);
   	for is=1:s
      		F(:,is) = z(:,is) - y0;
      		for j=1:s
         		F(:,is) = F(:,is) - h*A(is,j)*f(:,j);
      		end
   	end

   % reshape our output
   F = reshape(F, nvar*s, 1);

% end of function
end
	

function y = Y_IRK(fcn,y0,z,t,h,A,b,c,s)
% this function is used after newton's method
	zlen = length(z);
   	nvar = floor(zlen/s);

   	% reshape our z arguments into separate vectors for each stage
   	z = reshape(z,nvar,s);

   	% call f at our stages
   	f = zeros(nvar,s);
   	for is=1:s
      		t = t + h*c(is);
      		f(:,is) = fcn(t, z(:,is));
   	end

   % form the solution
   %    ynew = yold + h*sum(b(j)*fj)
   %disp(size(f))
   %disp(size(b))
   y = y0 + h*dot(f,b);
end




