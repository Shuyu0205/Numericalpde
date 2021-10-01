%take a matrix phi and then plot dphi/dx and dphi/dy at y=0
function u=Proj4_part_2_Pcolor_01935446(phi,q,s,r,N,M,h)
    u=zeros(M+1,N+1);
    v=zeros(M+1,N+1);
    u(1,1)=0;
    u(1,N+1)=0;
    u(M+1,N+1)=0;
    u(M+1,1)=0;
    v(1,1)=0;
    v(1,N+1)=0;
    v(M+1,N+1)=0;
    v(M+1,1)=0;
    for j=2:M
        for i=2:N
            u(j,i)=(phi(j,i+1)-phi(j,i-1))/(2*h);
            v(j,i)=(phi(j+1,i)-phi(j-1,i))/(2*h);
        end
    end
    subplot(1,2,1)
    surf(-q:h:s,0:h:r,u);
    
    subplot(1,2,2)
    surf(-q:h:s,0:h:r,v);
    shading falt;
    %subplot(2,2,3)
    %waterfall(-10:h:10,0:h:10,u);
    %subplot(2,2,4)
    %waterfall(-10:h:10,0:h:10,v);
            
end