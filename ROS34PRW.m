function [tspan,y] = ROS34PRW(f,J,tspan,y0)

N = length(y0);
y = zeros(N,length(tspan));
y(:,1) = y0;
h = tspan(2) - tspan(1);
Gamma = 0.435866521508459;
alpha = zeros(4,4);
alpha(2,1) = 0.87173304301691801;
alpha(3,1) = 1.4722022879435914;
alpha(3,2) = -0.31840250568090289;
alpha(4,1) = 0.81505192016694938;
alpha(4,2) = 0.5;
alpha(4,3) = -0.31505192016694938;
gamma = zeros(4,4);
gamma(2,1) = -alpha(2,1);
gamma(3,1) = -1.2855347382089872;
gamma(3,2) =  0.50507005541550687;
gamma(4,1) =  -0.50507005541550687;
gamma(4,2) = 0.21793326075422950;
gamma(4,3) = -0.17178529043404503;


b = zeros(4,1);
b(1) = 0.33303742833830591;
b(2) = 0.71793326075422947;
b(3) = -0.48683721060099439;
b(4) = Gamma;

for i = 2:length(tspan)
	Jac0 = J(y0);
	fcn0 = f(y0);
        rhs = (eye(N)-h*Gamma*Jac0);
	K1 = rhs\fcn0;
	K2 = rhs\(f(y0 + h*alpha(2,1)*K1) + h*gamma(2,1).*Jac0*K1);
	K3 = rhs\(f(y0 + h*alpha(3,1)*K1  + h*alpha(3,2)*K2) +  h*gamma(3,1).*Jac0*K1 + h*gamma(3,2)*Jac0*K2);
	K4 = rhs\(f(y0 + h*alpha(4,1)*K1  + h*alpha(4,2)*K2 + h*alpha(4,3)*K3) +  h*gamma(4,1).*Jac0*K1 + h*gamma(4,2)*Jac0*K2 + h*gamma(4,3)*Jac0*K3);
	y(:,i) = y0 + h*b(1).*K1 + h*b(2).*K2 + h*b(3).*K3 + h*b(4).*K4;
        y0 = y(:,i);
end
end
