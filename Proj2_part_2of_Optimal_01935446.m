%Q2, SOR method ,finding w
%dx=dy=h!
clear
q=10;
s=10;
r=8;
h=0.05;
L=((s+q)/h);
N=(r/h);
n_1=[];
%i=2:L;
%j=2:N;
%Different from the Jacobi method we just change all phi_new to phi will do
%the G-S method work.
tolerance=-6;
for w=1.65:0.01:1.7
    residual_new=1;
    n=0;%count iteration
    residual=[];
    phi= zeros(N+1,L+1);
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
             phi(1,i2)=0.25*w*(phi(1,i2+1)+phi(1,i2-1)+2*phi(2,i2)-((2*h+phi(1,i2+1)-phi(1,i2-1))*(d_y((-q+((i2-1)*h))))))+(1-w)*phi_old(1,i2);
         end;
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
n_1=[n_1,n];
end
n_1
plot(1.65:0.01:1.7,n_1)
function y=d_y(x)
tau=0.05;
A=0.298222773;
B=0.127125232;
C=0.357907906;
D=0.291984971;
E=0.105174606;
y=tau.*(A*(1./(2.*sqrt(x)))-B-2.*C.*x+3.*D.*x.^2-4.*E.*x.^3);
end