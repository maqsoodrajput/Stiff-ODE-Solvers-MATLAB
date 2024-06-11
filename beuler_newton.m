function [t,y,iter] = beuler_newton(f,J,t,y0)
dt = t(2)-t(1); %time step
%tol = 1e-13
%options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
%options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
y = zeros(length(y0),length(t));
% Midpoint for intial guess
y_guess = y0;
y(:,1) = y0;
tic
iter=0;
for i = 2:length(t)
	%ymid = y(:,i) + dt*f(t(i+1) + dt/2,y(:,i)+ f(t(i),y(:,i))*dt/2); 
	%fun = @(yn) y(:,i-1) + dt*f(t(i),yn)-yn;  %set up function x-f(x)= g(x) = 0
	%y(:,i) = fsolve(fun,y_guess,options); %use fsolve to find next point
        [y(:,i),isConverged,itr]=newtonSDIRK(f,J,y(:,i-1),y_guess,dt);
	y_guess = y(:,i);
        iter=iter+itr;
end
%fprintf('\nThe Backward Euler method at t = 3000 gave y1,y2 = %f %f \n',y(:,end))
toc
end


%https://sites.pitt.edu/~kimwong/lab03/index.html
function [Y,isConverged,iter]=newtonSDIRK(f,fPartial,yk,Y,h)
% [Y,isConverged]=newton4euler(f_ode,xkp1,yk,Y,h)
% special function to evaluate Newton's method for back_euler

% your name and the date
iter=0;
TOL = 1.e-1;
MAXITS = 500;
 
isConverged= (1==0);  % starts out FALSE
for n=1:MAXITS
  iter=iter+1;
  %[fValue fPartial] = f_ode( xkp1, Y);
  fValue = f(Y);
  PartialVal = fPartial(Y);
  %F = yk + h * fValue - Y;
  %J = h * PartialVal - eye(numel(Y));
  %increment=J\F;
  %Y = Y - increment;
  F = -yk - h * fValue + Y;
  J = -h * PartialVal+eye(numel(Y));
  increment=J\F;
  Y0=Y;
  Y = Y - increment;

  %if norm(increment,2) < TOL*norm(Y,2)
  if norm(abs(Y0-Y)) < TOL
    isConverged= (1==1);  % turns TRUE here
    return
  end
end
end
