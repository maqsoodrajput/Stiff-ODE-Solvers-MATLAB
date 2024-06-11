function [t,y] = feuler(f,t,y0)
dt = t(2)-t(1); %time step
y = zeros(length(y0),length(t));
y(:,1) = y0;
tic
for i = 2:length(t)
	y(:,i) = y(:,i-1) + dt*f(t(i-1),y(:,i-1)); 
end
if t(i)~=t(end)
disp('error')
end
toc
end
