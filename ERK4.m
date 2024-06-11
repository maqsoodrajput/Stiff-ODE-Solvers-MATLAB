function [t,y] = ERK4(fcn,t,y0)

dt = t(2) - t(1);
y = zeros(length(y0),length(t));
y(:,1) = y0;
tic
for i = 2:length(t)
	tn = t(i-1);
	yn = y(:,i-1);
	k1 = fcn(tn,yn);
	k2 = fcn(tn+dt/2,yn+dt*k1/2);
	k3 = fcn(tn+dt/2,yn+dt*k2/2);
	k4 = fcn(tn + dt,yn+ dt*k3);
	y(:,i) = yn + (k1+2*k2 + 2*k3 + k4)*dt/6;
end	
toc
end
