function [tspan,y] = ROS2(f,J,tspan,y0)

N = length(y0);
y = zeros(N,length(tspan));
y(:,1) = y0;
h = tspan(2) - tspan(1);
gamma = 1 + 1/sqrt(2);
gamma21 = -2*gamma;
alpha21 = 1;
b1 = 0.5;
b2 = 0.5;

for i = 2:length(tspan)
	Jac0 = J(y0);
	fcn0 = f(y0);
        rhs = (eye(N)-h*gamma*Jac0);
	K1 = rhs\fcn0;
	K2 = rhs\(f(y0 + h*alpha21*K1) + h*Jac0*gamma21*K1);
	y(:,i) = y0 + h*b1*K1 + h*b2*K2;
        y0 = y(:,i);
end
end
