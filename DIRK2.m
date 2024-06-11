function [t,y] = DIRK2(f,t,y0)
   %Two Stage Order 2 DIRK Method
   %Scheme: 
   %y1 = un + (2-sqrt(2))/2*dt*f(y1)
   %y2 = un + sqrt(2)/2*dt*f(y1) + (2-sqrt(2))/2*dt*f(y2)
   %unp1 = un+ sqrt(2)/2*dt*f(y1) + (2-sqrt(2))/2*dt*f(y2)
                    
   %p = 2
   %y1=fsolve( @(x) x- ( un + (2-sqrt(2))/2*dt*f(x)) ,un);
   %y2= fsolve(@(x) x-( un + sqrt(2)/2*dt*f(y1) + (2-sqrt(2))/2*dt*f(x)),un);
   %unp1 = un+ sqrt(2)/2*dt*f(y1) + (2-sqrt(2))/2*dt*f(y2)
   
   tol = 1e-13;  
   %options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter');
   options = optimset('disp','off','TolFun',tol,'MaxFunEvals',1e5,'Maxiter',1e5);
   y = zeros(length(y0),length(t));
   y(:,1) = y0;
   dt = t(2) - t(1);
   p=3;
   tn = t(1);
   un = y0;
   a = zeros(2,2);
   a(1,1) = (2 - sqrt(2))/2;
   a(1,2) = 0.0;
   a(2,1) = sqrt(2)/2;
   a(2,2) = (2 - sqrt(2))/2;
   b = zeros(2,1);
   b(1,1) = sqrt(2)/2;
   b(2,1) = (2 - sqrt(2))/2;
   for i = 2:length(t)
   	y1 = fsolve( @(x) x - ( un + a(1,1)*dt*f(x)) ,un,options);
   	y2 = fsolve(@(x)  x - ( un + a(2,1)*dt*f(y1) +  a(2,2)*dt*f(x)),un,options);
   	y(:,i)= un + b(1,1)*dt*f(y1) + b(2,1)*dt*f(y2);
	un = y(:,i);            %Updating the initial guess	
   end
end
