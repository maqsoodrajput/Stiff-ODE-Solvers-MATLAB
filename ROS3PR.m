function [tspan,y] = ROS3PR(f,J,tspan,y0)


a = zeros(3,3);
a(2,1) = 2.3660254037844388;
a(3,1) = 0.0000000000000000;
a(3,2) = 1.0000000000000000;

b = zeros(3,1);
b(1,1) = 2.9266384402395124e-01;
b(2,1) = -8.1338978618764143e-02;
b(3,1) = 7.8867513459481287e-01;


%b(1,1) = (7/18) - (1/18)*sqrt(3);
%b(2,1) = (1/9) - (1/9)*sqrt(3);
%b(3,1) = (1/2) + (1/6)*sqrt(3);

%b(1,1) = (13/36) + (1/12)*sqrt(3);
%b(2,1) = (1/3)   - (7/27)*sqrt(3);
%b(3,1) = (11/36) - (19/108)*sqrt(3);

gamma = zeros(3,3);

gamma(2,1) = -2.3660254037844388e+00;
gamma(3,1) = -2.8468642516567449e-01;
gamma(3,2) = -1.0813389786187642e+00;

h = tspan(2) - tspan(1);

Gamma = 7.8867513459481287e-01;
%Gamma = (3+sqrt(3))/6;

y = zeros(length(y0),length(tspan));
y(:,1) = y0;
N = length(y0);
for i = 2:length(tspan)
	y(:,i) = solver_test(N,a,b,gamma,Gamma,f,J,y0,h);
	y0 = y(:,i);	
end

end 
 

function y = solver(N,a,b,gamma,Gamma,f,J,y0,h)
fcn0 = f(y0);
Jac0 = J(y0);
lhs = eye(N) - h*Gamma*Jac0;
%1st stage
rhs = h*fcn0;
k1 = lhs\rhs;
%2nd stage
fcn = f(y0);
rhs = fcn + h.*Jac0*(gamma(2,1).*k1);
k2 = lhs\rhs;
%3rd stage
fcn = f(y0);
rhs = fcn + h.*Jac0*(gamma(3,1).*k1 + gamma(3,2).*k2);
k3 = lhs\rhs;

y = y0 + b(1,1).*k1 + b(2,1).*k2 + b(3,1).*k3;
end




function y = solver_test(N,a,b,gamma,Gamma,f,J,y0,h)
%1st stage
fcn1 = f(y0);
Jac1 = J(y0);
lhs = eye(N) - h*Gamma*Jac1;
rhs = h*fcn1;
k1 = lhs\rhs;
%2nd stage
g2 = y0 + a(2,1).*k1;
fcn2 = f(g2);
Jac2 = J(g2);
rhs = fcn2 + h.*Jac2*(gamma(2,1).*k1);
lhs = eye(N) - h*Gamma*Jac2;
k2 = lhs\rhs;
%3rd stage
g3 = y0 + a(3,1).*k1 + a(3,2).*k2;
fcn3 = f(g3);
Jac3 = J(g3);
rhs = fcn3 + h.*Jac3*(gamma(3,1).*k1 + gamma(3,2).*k2);
lhs = eye(N) -h*Gamma*Jac3;
k3 = lhs\rhs;

y = y0 + b(1,1).*k1 + b(2,1).*k2 + b(3,1).*k3;
end




