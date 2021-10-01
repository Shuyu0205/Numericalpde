%Q3, Poisson equation, SOR method with change of coordinates
%dx=dy=h!
clear
q=10;
s=10;
r=8;
h=0.05;
L=((s+q)/h);
N=(r/h);
phi= zeros(N+1,L+1);
%i=2:L;
%j=2:N;
residual=[];
w=1.5;
%Different from the Jacobi method we just change all phi_new to phi will do
%the G-S method work.
residual_new=1;
n=0;%count iteration
tolerance=-7;
while residual_new>=tolerance;
    phi_old=phi;
    %process four corner points
    phi(1,1)=0.25*w*(2*phi(2,1)+2*phi(1,2))+(1-w)*phi_old(1,1);
    phi(1,L+1)=0.25*w*(2*phi(2,L+1)+2*phi(1,L))+(1-w)*phi_old(1,L+1);
    phi(N+1,L+1)=0.25*w*(2*phi(N,L+1)+2*phi(N+1,L))+(1-w)*phi_old(N+1,L+1);
    phi(N+1,1)=0.25*w*(2*phi(N,1)+2*phi(N+1,2))+(1-w)*phi_old(N+1,1);
    for j=2:N
        
        phi(j,1)=0.25*w*(phi(j+1,1)+phi(j-1,1)+2*phi(j,2))+(1-w)*phi_old(j,1);%left barrier    
    end
    for i=2:L
        
        phi(N+1,i)=0.25*w*(phi(N+1,i+1)+phi(N+1,i-1)+2*phi(N,i))+(1-w)*phi_old(N+1,i);%Upper barrier    
    end
    for j=2:N
        
        phi(j,L+1)=0.25*w*(phi(j+1,L+1)+phi(j-1,L+1)+2*phi(j,L))+(1-w)*phi_old(j,L+1);%right barrier
    end
    for i2=2:1:L  
         %Special lower boundary identification
         %free (-q,0] & (1,s)
         if (-q+((i2-1)*h))<=0 | (-q+((i2-1)*h))>1
             phi(1,i2)=0.25*w*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2))+(1-w)*phi_old(1,i2);
                %boundary layer (1,0)
         else
             x=-q+((i2-1)*h);
             phi(1,i2)=(w/(2+2*d_y(x)+2))*(phi(1,i2+1)+phi(1,i2-1)+(d_y(x)+1)*(phi(2,i2)+(phi(2,i2)-(d_y(x)/(1+(d_y(x))^2))*(2*h+phi(1,i2+1)-phi(1,i2-1))))...
             +(d_y(x)/2)*(phi(2,i2+1)-(phi(2,i2+1)-(d_y(x)/(1+(d_y(x))^2))*(2*h+phi(1,i2+2)-phi(1,i2)))-phi(2,i2-1)+(phi(2,i2-1)-(d_y(x)/(1+(d_y(x))^2))*(2*h+phi(1,i2)-phi(1,i2-2)))))...
             +(1-w)*phi_old(1,i2);
         end
    end
    for j=2:N
        for i=2:L
            %interior points
            phi(j,i)=0.25*w*(phi(j+1,i)+phi(j,i+1)+phi(j,i-1)+phi(j-1,i))+(1-w)*phi_old(j,i);
        end
    end
    %calaulate the log-10 error
    residual_new=log10(max(max(abs(phi-phi_old))));
    residual=[residual;residual_new];
    n=n+1;%count the iteration
end
%uncomment to plot the error
plot(residual);
%hold on
%waterfall(-q:h:s,0:h:r,phi)
fprintf('Number of iteration to reach residual %6.0f',tolerance)
n
%uncomment the following line to verify the solution
%Proj2_part_1of_Ascertainfunction_01935446(phi,h,q,s)
%hold on

function y=d_y(x)
tau=0.01;%can be changed, the thickness of boundary layer
A=0.298222773;
B=0.127125232;
C=0.357907906;
D=0.291984971;
E=0.105174606;
y=tau.*(A*(1./(2.*sqrt(x)))-B-2.*C.*x+3.*D.*x.^2-4.*E.*x.^3);
end