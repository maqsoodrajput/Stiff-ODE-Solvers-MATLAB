function [tvals,y] = EDIRK2(fcn,tvals,y0)

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
   a = zeros(2,2);
   a(1,1) = 1/3;
   a(1,2) = 0.0;

   a(2,1) = 3/4;
   a(2,2) = 1/4;

   c = zeros(2,1);
   c(1,1) = 1/3;
   c(2,1) = 1.0;

   b = zeros(2,1); 
   b(1,1) = 3/4;
   b(2,1) = 1/4;

   tol = 1e-13;  
   
   %options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
   options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
   %options = optimset('disp','on','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);

   
   %k1 = f(tn + c(1,2)*dt,yn+a(1,1)*dt*k1 + a(1,2)*dt*k2);
   %k2 = f(tn + c(2,1)*dt,yn+a(2,1)*dt*k1 + a(2,2)*dt*k2);
   %y = yn + b(1,2)*dt*k1 + b(2,1)*dt*k2;
   
    
   %F = @(k)  [k(1) - f(tn + c(1,2)*dt,yn + a(1,1)*dt.*k(1) + a(1,2)*dt.*k(2));
   %           k(2) - f(tn + c(2,1)*dt,yn + a(2,1)*dt.*k(1) + a(2,2)*dt.*k(2))];  
   F = @(k)  [k(1) - f(yn + a(1,1)*dt.*k(1) + a(1,2)*dt.*k(2));
              k(2) - f(yn + a(2,1)*dt.*k(1) + a(2,2)*dt.*k(2))];
   
   K0 = [yn,yn]; % initiall guess
   K = fsolve(F,K0,options);
   y = yn + b(1,1)*dt.*K(1) + b(2,1)*dt.*K(2);
end
