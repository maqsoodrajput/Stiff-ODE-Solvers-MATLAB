function [tspan,y] = ROSW(f,J,tspan,y0)

N = length(y0);
y = zeros(N,length(tspan));
y(:,1) = y0;
h = tspan(2) - tspan(1);
Gamma = 0.25;
alpha = zeros(6,6);
alpha(2,1) = 0.1453095851778752E+00;
alpha(3,1) = ;
alpha(3,2) = 0.0;
gamma = zeros(3,3);
gamma(2,1) = -alpha(2,1);
gamma(3,1) = -0.67075317547305480;
gamma(3,2) = -0.17075317547305482;

b = zeros(3,1);
b(1) = 0.10566243270259355;
b(2) = 0.049038105676657971;
b(3) = 0.84529946162074843;

for i = 2:length(tspan)
	Jac0 = J(y0);
	fcn0 = f(y0);
        rhs = (eye(N)-h*Gamma*Jac0);
	K1 = rhs\fcn0;
	K2 = rhs\(f(y0 + h*alpha(2,1)*K1) + h*gamma(2,1).*Jac0*K1);
	K3 = rhs\(f(y0 + h*alpha(3,1)*K1  + h*alpha(3,2)*K2) +  h*gamma(3,1).*Jac0*K1 + h*gamma(3,2)*Jac0*K2);
	y(:,i) = y0 + h*b(1).*K1 + h*b(2).*K2 + h*b(3).*K3;
        y0 = y(:,i);
end
end
