function [tspan, y] = TRBDF2(f,tspan,y0)

N = length(y0);
y = zeros(length(y0),length(tspan));
y(:,1) = y0;
dt = tspan(2) - tspan(1);
Gamma = 2 - sqrt(2);

tol = 1e-13;
options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
%options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');

for i = 2:length(tspan)
	%% TR-stage
	u_g = fsolve( @(u) u - y0 - Gamma*dt*0.5*f(u) - Gamma*dt*0.5*f(y0),y0,options);
	%% BDF2-stage
	y(:,i) = fsolve(@(u) u - ((1-Gamma)/(2-Gamma))*dt*f(u) - (1/(Gamma*(2-Gamma)))*u_g + ((1-Gamma)^2/(Gamma*(2-Gamma)))*y0,y0,options);
	y0 = y(:,i);
end
