function [t,y] = beuler(f,t,y0)
dt = t(2)-t(1); %time step
tol = 1e-13
options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
%options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
y = zeros(length(y0),length(t));
% Midpoint for intial guess
y_guess = y0;
y(:,1) = y0;
tic
for i = 2:length(t)
	%ymid = y(:,i) + dt*f(t(i+1) + dt/2,y(:,i)+ f(t(i),y(:,i))*dt/2); 
	fun = @(yn) y(:,i-1) + dt*f(t(i),yn)-yn;  %set up function x-f(x)= g(x) = 0
	y(:,i) = fsolve(fun,y_guess,options); %use fsolve to find next point
	y_guess = y(:,i);
end
%fprintf('\nThe Backward Euler method at t = 3000 gave y1,y2 = %f %f \n',y(:,end))
toc
end
