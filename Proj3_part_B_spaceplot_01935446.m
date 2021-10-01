%plot in the three plane y=-0.2, y=0.2 and x=0
T=25;
c=1;
h=0.02;
k=0.008;
y_0=0.5;
q=(c*k)/(h);%here we let q=0.4
N=(4/h)+1;
M=(4/h)+1;
J=(T/k)+1;
delta=0.2;%wave source scale
d=0.5;%the gap width
%build up the empty solution matrix;
U= zeros(M,N,J);
%Start setting initial condition to derive the first page of the solution
for m=1:1:M %y
    for n=1:1:N%x
        r=sqrt((-2+((n-1)*h))^2+(((-2+(m-1)*h))-y_0)^2);
        if (r>=-delta) && (r<=delta)%locate the wave source
            U(m,n,1)=cos((pi*r)/(2*delta));
        end
    end
end
%Second page of the solution
for m=2:1:M-1%y
    for n=2:1:N-1%x
        U(m,n,2)=0.5*((q^2)*((U(m,n+1,1)-2*U(m,n,1)+U(m,n-1,1))+(U(m+1,n,1)-2*U(m,n,1)+U(m-1,n,1)))+2*U(m,n,1));%first time layer
    end
end
%Bottom solid wall 
for j=1:1:J
    for n=1:1:N
        U(1,n,j)=0;
    end
end
%Left solid wall
for j=1:1:J
    for m=1:1:M
        if (-2+(m-1)*h)<=0
            U(m,1,j)=0;
        end
    end
end
%right solid wall boundary
for j=1:1:J
    for m=1:1:M
        if (-2+(m-1)*h)<=0
            U(m,N,j)=0;
            m_location=m;%save it for upper solid wall
        end
    end
end
%Upper solid wall boundary
for j=1:1:J
    for n=1:1:N
        if (-2+(n-1)*h)<=-d |  (-2+(n-1)*h)>=d
            U(m_location,n,j)=0;
        end
    end
end
%main calculation
for j=3:1:J
    for m=2:1:M-1
        for n=2:1:N-1
            if m==m_location
                if (-2+(n-1)*h)<=d && (-2+(n-1)*h)>=-d
                    U(m,n,j)=(q^2)*((U(m,n+1,j-1)-2*U(m,n,j-1)+U(m,n-1,j-1))+(U(m+1,n,j-1)-2*U(m,n,j-1)+U(m-1,n,j-1)))+2*U(m,n,j-1)-U(m,n,j-2);
                end
            else
                U(m,n,j)=(q^2)*((U(m,n+1,j-1)-2*U(m,n,j-1)+U(m,n-1,j-1))+(U(m+1,n,j-1)-2*U(m,n,j-1)+U(m-1,n,j-1)))+2*U(m,n,j-1)-U(m,n,j-2);
            end                    
        end
    end
    %non-reflection boundary
    %left and right
    for m=m_location+1:1:M-1
        U(m,1,j)=(1/(1+(q/2)))*((q/2)*(U(m,2,j)-U(m,2,j-2)+U(m,1,j-2))+((q^2)/2)*(U(m+1,1,j-1)-2*U(m,1,j-1)+U(m-1,1,j-1))+2*U(m,1,j-1)-U(m,1,j-2));
        U(m,N,j)=(1/(1+(q/2)))*(-(q/2)*(-U(m,N-1,j)-U(m,N,j-2)+U(m,N-1,j-2))+((q^2)/2)*(U(m+1,N,j-1)-2*U(m,N,j-1)+U(m-1,N,j-1))+2*U(m,N,j-1)-U(m,N,j-2));
    end
    %Top
    for n=2:1:N-1
        U(M,n,j)=(1/(1+(q/2)))*(-(q/2)*(-U(M-1,n,j)-U(M,n,j-2)+U(M-1,n,j-2))+((q^2)/2)*(U(M,n+1,j-1)-2*U(M,n,j-1)+U(M,n-1,j-1))+2*U(M,n,j-1)-U(M,n,j-2));
    end
    %Two corner
    U(M,1,j)=(1/(-3-2*q))*(q*((-U(M-1,1,j)-U(M,1,j-2)+U(M-1,1,j-2))-(U(M,2,j)-U(M,2,j-2)+U(M,1,j-2)))-6*U(M,1,j-1)+3*U(M,1,j-2));%(-2,2) point
    U(M,N,j)=(1/(-3-2*q))*(q*((-U(M-1,N,j)-U(M,N,j-2)+U(M-1,N,j-2))+(-U(M,N-1,j)-U(M,N,j-2)+U(M,N-1,j-2)))-6*U(M,N,j-1)+3*U(M,N,j-2));%(2,2) point
end
U_x_1=zeros(J,N);
U_x_2=zeros(J,N);
U_y_1=zeros(M,J);
%at y=0.2
m_1=(2.2/h)+1;
for j=1:1:J
    for n=1:1:N
    U_x_1(j,n)=U(m_1,n,j);
    end
end
%at y=-0.2
m_2=(1.8/h)+1;
for j=1:1:J
    for n=1:1:N
    U_x_2(j,n)=U(m_2,n,j);
    end
end
%at x=0
n_1=(2/h)+1;
for j=1:1:J
    for m=1:1:M
    U_y_1(m,j)=U(m,n_1,j);
    end
end
subplot(1,3,1)
surf(-2:h:2,0:k:T,U_x_1)
shading flat
subplot(1,3,2)
surf(-2:h:2,0:k:T,U_x_2)
shading flat
subplot(1,3,3)
surf(0:k:T,-2:h:2,U_y_1)
shading flat

