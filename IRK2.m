function [t, y] = IRK2(f, tspan, y0)

% f: function (RHS in y' = f(t,x))
% tspan: [t0, tf], where t0 is the initial time and tf is the final time
% y0: initial value of the solution
% t: vector of time points
% y: vector of solution values

N = length(tspan)-1;
% Determine time step size
h = (tspan(2) - tspan(1))/(N);

% Initialize solution vectors
t = linspace(tspan(1), tspan(2), N+1);
y = zeros(length(y0),N+1);
y(:,1) = y0;

% IRK2 coefficients
a = 1/2 - sqrt(3)/6;
b = 1/2 + sqrt(3)/6;
c = 1/2;

% Fixed-point iteration parameters
maxiter = 10000; % maximum number of iterations
tol = 1e-3; % tolerance level

for n = 1:N
    % Compute y_{n+1/2}
    yhalf = (y(:,n) + y(:,n+1))/2;
    
    % Initialize fixed-point iteration
    yk = y(:,n+1);
    err = 1;
    iter = 0;
    
    % Fixed-point iteration
    while (err > tol) && (iter < maxiter)
        yk1 = y(:,n) + h*a*f(t(n+1), yhalf) + h*b*f(t(n+1), yk);
        %err = norm(yk1 - yk);
	err = l2norm(yk1-yk);
        yk = yk1;
        iter = iter + 1;
    end
    
    % Update solution vector
    y(:,n+1) = yk;
end

end

function y = l2norm(y)
y = sum(y.^2)/length(y);
end
