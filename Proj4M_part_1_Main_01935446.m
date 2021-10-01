clear
q=5;
s=5;
r=8;
h=0.02;
N=((s+q)/h);
M=(r/h);
Ome=1;
Ma=0.01;%Mach number
phi=zeros(M+1,N+1);
rho=zeros(M+1,N+1);
%i=2:L;
%j=2:N;
residual=[];
%Different from the Jacobi method we just change all phi_new to phi will do
%the G-S method work.
residual_new=1;
k=0;%count iteration
tolerance=-6;
%while n<=3000;
tic
while residual_new>=tolerance;
    phi_old=phi;
    %process four corner points
    phi(1,1)=0.25*(2*phi(2,1)+2*phi(1,2));
    phi(1,N+1)=0.25*(2*phi(2,N+1)+2*phi(1,N));
    phi(M+1,N+1)=0.25*(2*phi(M,N+1)+2*phi(M+1,N));
    phi(M+1,1)=0.25*(2*phi(M,1)+2*phi(M+1,2));
    %boundary conditions, should obey left, upper, right, lower to make all
    %points updates
    for m=2:M
        
        phi(m,1)=0.25*(phi(m+1,1)+phi(m-1,1)+2*phi(m,2));%left barrier  
    end
    for n=2:N
        
        phi(M+1,n)=0.25*(phi(M+1,n+1)+phi(M+1,n-1)+2*phi(M,n));%Upper barrier
    end
    for m=2:M
        phi(m,N+1)=0.25*(phi(m+1,N+1)+phi(m-1,N+1)+2*phi(m,N));%right barrier
    end
    %Special lower boundary identification
    for i2=2:1:N        
        %free (-q,0] & (1,s)
        if (-q+((i2-1)*h))<=0 | (-q+((i2-1)*h))>1
           phi(1,i2)=0.25*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2));
        %aerofoil (1,0)
        else
           phi(1,i2)=0.25*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2)-((2*h+phi(1,i2+1)-phi(1,i2-1))*(d_y((-q+((i2-1)*h))))));
        end;
    end
    for m=2:M
        for n=2:N          
    %interior points
            a=(1/(2*h))*(phi(m,n+1)-phi(m,n-1));%central difference
            b=(1/(2*h))*(rho(m,n+1)-rho(m,n-1));
            c=(1/(2*h))*(phi(m+1,n)-phi(m-1,n));
            d=(1/(2*h))*(rho(m+1,n)-rho(m-1,n));
            phi(m,n)=0.25*(phi(m+1,n)+phi(m,n+1)+phi(m,n-1)+phi(m-1,n))+((h^2)/(4*(1+rho(m,n))))*(((1+a)*b)+c*d);
        end
    end
    for m=2:M
        for n=2:N
            a=(1/(2*h))*(phi(m,n+1)-phi(m,n-1));
            c=(1/(2*h))*(phi(m+1,n)-phi(m-1,n));
            rho(m,n) = (1-0.2*Ma^2*(2*a+a^2+c^2))^(2.5)-1;
        end
    end
    
    %calaulate the log-10 residual
    n_p=(M+1)*(N+1);
    sum=0;
    for m=1:1:M+1
        for n=1:1:N+1
            det=(abs(phi_old(m,n))+abs(phi(m,n)))*0.5;
            sum=sum+(abs(abs(phi_old(m,n))-abs(phi(m,n))))/det;
        end
    end
    residual_new=log10(max(max(abs(phi-phi_old))));
    residual=[residual;residual_new];
    k=k+1;%count the iteration
end
toc
%uncomment to plot the error
%plot(residual);
%hold on
surf(-q:h:s,0:h:r,phi)
shading flat
fprintf('Number of iteration to reach residual %6.0f',tolerance)
k
%uncomment the following line to verify the solution
%Proj4M_part_2_perturbation_01935446(phi,rho,h,q,s)
%hold on
function y=d_y(x)
tau=0.05;
A=0.298222773;
B=0.127125232;
C=0.357907906;
D=0.291984971;
E=0.105174606;
y=tau.*(A*(1./(2.*sqrt(x)))-B-2.*C.*x+3.*D.*x.^2-4.*E.*x.^3);
end