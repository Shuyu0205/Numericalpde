%Q1, Poisson equation, Jacobi method 
%dx=dy=h!
clear
q=10;
s=10;
r=8;
h=0.05;
L=((s+q)/h);
N=(r/h);
phi= zeros(N+1,L+1);
phi_new=phi;
i=2:L;
j=2:N;
residual=[];
residual_new=1;
n=0;%count iteration
tolerance=-7;
%while n<=3000;
while residual_new>=tolerance;
    %process four corner points
    phi_new(1,1)=0.25*(2*phi(2,1)+2*phi(1,2));
    phi_new(1,L+1)=0.25*(2*phi(2,L+1)+2*phi(1,L));
    phi_new(N+1,L+1)=0.25*(2*phi(N,L+1)+2*phi(N+1,L));
    phi_new(N+1,1)=0.25*(2*phi(N,1)+2*phi(N+1,2));
    %interior points
    phi_new(j,i)=0.25*(phi(j+1,i)+phi(j,i+1)+phi(j,i-1)+phi(j-1,i));
    %boundary conditions
    phi_new(j,1)=0.25*(phi(j+1,1)+phi(j-1,1)+2*phi(j,2));%left barrier
    phi_new(j,L+1)=0.25*(phi(j+1,L+1)+phi(j-1,L+1)+2*phi(j,L));%right barrier
    phi_new(N+1,i)=0.25*(phi(N+1,i+1)+phi(N+1,i-1)+2*phi(N,i));%Upper barrier
    %Special boundary identification
    for i2=2:1:L        
        %free (-q,0] & (1,s)
        if (-q+((i2-1)*h))<=0 | (-q+((i2-1)*h))>1
           phi_new(1,i2)=0.25*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2));
        %boundary layer (0,1]
        else
           phi_new(1,i2)=0.25*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2)-((2*h+phi(1,i2+1)-phi(1,i2-1))*(d_y((-q+((i2-1)*h))))));
        end;
    end
    %calaulate the log-10 error
    residual_new=log10(max(max(abs(phi_new-phi))));
    phi=phi_new;
    residual=[residual;residual_new];
    n=n+1;%count the iteration
end 
%uncomment to plot the error
plot(residual);
hold on %remain this and run G_S.m and SOR.m to place everything in one plot
%waterfall(-q:h:s,0:h:r,phi)%plot the solution
fprintf('Number of iteration to reach residual %6.0f',tolerance)
n
%uncomment the following line to verify the solution
%Proj2_part_1of_Ascertainfunction_01935446(phi,h,q,s)
function y=d_y(x)
tau=0.05;
A=0.298222773;
B=0.127125232;
C=0.357907906;
D=0.291984971;
E=0.105174606;
y=tau.*(A*(1./(2.*sqrt(x)))-B-2.*C.*x+3.*D.*x.^2-4.*E.*x.^3);
end

