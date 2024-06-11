function [tvals,y] = RadauIIA(fcn,tvals,y0)

% fcn is the RHS in y'=f(t,y)
% Jac is the Jacobian evaluated at the internal stages Y_i 
% $J_i = \partial f(Y_i)/ \partial y$

y = zeros(length(y0),length(tvals));
y(:,1) = y0;
t = tvals(1);
y_guess = y0;
dt = tvals(2) - tvals(1);

% iterate over time steps
for i = 2:length(tvals)
	 t = tvals(i-1);
	 y(:,i) = solver(fcn,y_guess,dt,t);
	 y_guess = y(:,i);
end
end

function y = solver(f,yn,dt,tn)
   a = zeros(3,3);
   a(1,1) = (88-7*sqrt(6))/360;
   a(1,2) = (296 - 169*sqrt(6))/1800;
   a(1,3) = (-2+3*sqrt(6))/225;

   a(2,1) = (296 + 169*sqrt(6))/1800;
   a(2,2) = (88 + 7*sqrt(6))/360;
   a(2,3) = (-2-3*sqrt(6))/225;
   
   a(3,1) = (16-sqrt(6))/36;
   a(3,2) = (16+sqrt(6))/36;
   a(3,3) = (1)/9;

   c = zeros(3,1);
   c(1,1) = (4 - sqrt(6))/10;
   c(2,1) = (4 + sqrt(6))/10;
   c(3,1) = 1.0;

   b = zeros(3,1); 
   b(1,1) = (16-sqrt(6))/36;
   b(2,1) = (16+sqrt(6))/36;
   b(3,1) = 1/9;
   tol = 1e-13;  
   
   %options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
   options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
   %options = optimset('disp','on','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);

   
   %k1 = f(tn + c(1,2)*dt,yn+a(1,1)*dt*k1 + a(1,2)*dt*k2);
   %k2 = f(tn + c(2,1)*dt,yn+a(2,1)*dt*k1 + a(2,2)*dt*k2);
   %y = yn + b(1,2)*dt*k1 + b(2,1)*dt*k2;
   
    
   %F = @(k)  [k(1) - f(tn + c(1,2)*dt,yn + a(1,1)*dt.*k(1) + a(1,2)*dt.*k(2));
   %           k(2) - f(tn + c(2,1)*dt,yn + a(2,1)*dt.*k(1) + a(2,2)*dt.*k(2))];  
   F = @(g)  [g(1) - yn - a(1,1)*dt.*f(g(1)) - a(1,2)*dt.*f(g(2)) - a(1,3)*dt.*f(g(3));
              g(2) - yn - a(2,1)*dt.*f(g(1)) - a(2,2)*dt.*f(g(2)) - a(2,3)*dt.*f(g(3));
	      g(3) - yn - a(3,1)*dt.*f(g(1)) - a(3,2)*dt.*f(g(2)) - a(3,3)*dt.*f(g(3))];  
   
 
   Y0 = [yn,yn,yn]; % initiall guess
   G = fsolve(F,Y0,options);
   y = yn + b(1,1)*dt.*f(G(1)) + b(2,1)*dt.*f(G(2)) + b(3,1)*dt.*f(G(3));
end


function y = test_solver(f,yn,dt,tn)
   a = zeros(3,3);
   a(1,1) = (88-7*sqrt(6))/360;
   a(1,2) = (296 - 169*sqrt(6))/1800;
   a(1,3) = (-2+3*sqrt(6))/225;

   a(2,1) = (296 + 169*sqrt(6))/1800;
   a(2,2) = (88 + 7*sqrt(6))/360;
   a(2,3) = (-2-3*sqrt(6))/225;
   
   a(3,1) = (16-sqrt(6))/36;
   a(3,2) = (16+sqrt(6))/36;
   a(3,3) = (1)/9;

   c = zeros(3,1);
   c(1,1) = (4 - sqrt(6))/10;
   c(2,1) = (4 + sqrt(6))/10;
   c(3,1) = 1.0;

   b = zeros(3,1); 
   b(1,1) = (16-sqrt(6))/36;
   b(2,1) = (16+sqrt(6))/36;
   b(3,1) = 1/9;
   tol = 1e-13;  
   
   %options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
   options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
   %options = optimset('disp','on','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);

   
   %k1 = f(tn + c(1,2)*dt,yn+a(1,1)*dt*k1 + a(1,2)*dt*k2);
   %k2 = f(tn + c(2,1)*dt,yn+a(2,1)*dt*k1 + a(2,2)*dt*k2);
   %y = yn + b(1,2)*dt*k1 + b(2,1)*dt*k2;
   
    
   %F = @(k)  [k(1) - f(tn + c(1,2)*dt,yn + a(1,1)*dt.*k(1) + a(1,2)*dt.*k(2));
   %           k(2) - f(tn + c(2,1)*dt,yn + a(2,1)*dt.*k(1) + a(2,2)*dt.*k(2))];  
   F = @(k)  [k(1) - f(yn + a(1,1)*dt.*k(1) + a(1,2)*dt.*k(2) + a(1,3)*dt.*k(3));
              k(2) - f(yn + a(2,1)*dt.*k(1) + a(2,2)*dt.*k(2) + a(2,3)*dt.*k(3));
	      k(3) - f(yn + a(3,1)*dt.*k(1) + a(3,2)*dt.*k(2) + a(3,3)*dt.*k(3))];  
   
 
   K0 = [yn,yn,yn]; % initiall guess
   K = fsolve(F,K0,options);
   y = yn + b(1,1)*dt.*K(1) + b(2,1)*dt.*K(2) + b(3,1)*dt.*K(3);
end






