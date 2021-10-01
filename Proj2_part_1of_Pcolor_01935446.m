%take a matrix phi and then plot dphi/dx and dphi/dy at y=0
function u=cplot(phi,L,N,h)
    u=zeros(N+1,L+1);
    v=zeros(N+1,L+1);
    u(1,1)=0;
    u(1,L+1)=0;
    u(N+1,L+1)=0;
    u(N+1,1)=0;
    v(1,1)=0;
    v(1,L+1)=0;
    v(N+1,L+1)=0;
    v(N+1,1)=0;
    for j=2:N
        for i=2:L
            u(j,i)=(phi(j,i+1)-phi(j,i-1))/(2*h);
            v(j,i)=(phi(j+1,i)-phi(j-1,i))/(2*h);
        end
    end
    subplot(1,2,1)
    waterfall(-10:h:10,0:h:10,u);
    subplot(1,2,2)
    waterfall(-10:h:10,0:h:10,v);
    %subplot(2,2,3)
    %waterfall(-10:h:10,0:h:10,u);
    %subplot(2,2,4)
    %waterfall(-10:h:10,0:h:10,v);
            
end