%This script calculate the error by comparing the result from doubling the
%points
T=25;
c=1;
h_1=0.016;
k_1=0.008;
y_0=0.5;
q=(c*k_1)/(h_1);%here we let q=0.5
N_1=(4/h_1)+1;
M_1=(4/h_1)+1;
J_1=(T/k_1)+1;
delta=0.2;%wave source scale
d_1=0.2;%the width of slit
%build up the empty solution matrix;
U_1= zeros(M_1,N_1,J_1);
%Start setting initial condition to derive the first page of the solution
for m=1:1:M_1 %y
    for n=1:1:N_1%x
        r=sqrt((-2+((n-1)*h_1))^2+(((-2+(m-1)*h_1))-y_0)^2);
        if (r>=-delta) && (r<=delta)%locate the wave source
            U_1(m,n,1)=cos((pi*r)/(2*delta));
        end
    end
end
%Second page of the solution
for m=2:1:M_1-1%y
    for n=2:1:N_1-1%x
        U_1(m,n,2)=0.5*((q^2)*((U_1(m,n+1,1)-2*U_1(m,n,1)+U_1(m,n-1,1))+(U_1(m+1,n,1)-2*U_1(m,n,1)+U_1(m-1,n,1)))+2*U_1(m,n,1));%first time layer
    end
end
%Bottom solid wall 
for j=1:1:J_1
    for n=1:1:N_1
        U_1(1,n,j)=0;
    end
end
%Left solid wall
for j=1:1:J_1
    for m=1:1:M_1
        if (-2+(m-1)*h_1)<=0
            U_1(m,1,j)=0;
        end
    end
end
%right solid wall boundary
for j=1:1:J_1
    for m=1:1:M_1
        if (-2+(m-1)*h_1)<=0
            U_1(m,N_1,j)=0;
            m_location=m;%save it for upper solid wall
        end
    end
end
%Upper solid wall boundary
for j=1:1:J_1
    for n=1:1:N_1
        if (-2+(n-1)*h_1)<=-d_1 |  (-2+(n-1)*h_1)>=d_1
            U_1(m_location,n,j)=0;
        end
    end
end
%main calculation
for j=3:1:J_1
    for m=2:1:M_1-1
        for n=2:1:N_1-1
            if m==m_location
                if (-2+(n-1)*h_1)<=d_1 && (-2+(n-1)*h_1)>=-d_1
                    U_1(m,n,j)=(q^2)*((U_1(m,n+1,j-1)-2*U_1(m,n,j-1)+U_1(m,n-1,j-1))+(U_1(m+1,n,j-1)-2*U_1(m,n,j-1)+U_1(m-1,n,j-1)))+2*U_1(m,n,j-1)-U_1(m,n,j-2);
                end
            else
                U_1(m,n,j)=(q^2)*((U_1(m,n+1,j-1)-2*U_1(m,n,j-1)+U_1(m,n-1,j-1))+(U_1(m+1,n,j-1)-2*U_1(m,n,j-1)+U_1(m-1,n,j-1)))+2*U_1(m,n,j-1)-U_1(m,n,j-2);
            end                    
        end
    end
    %non-reflection boundary
    %left and right
    for m=m_location+1:1:M_1-1
        U_1(m,1,j)=(1/(1+(q/2)))*((q/2)*(U_1(m,2,j)-U_1(m,2,j-2)+U_1(m,1,j-2))+((q^2)/2)*(U_1(m+1,1,j-1)-2*U_1(m,1,j-1)+U_1(m-1,1,j-1))+2*U_1(m,1,j-1)-U_1(m,1,j-2));
        U_1(m,N_1,j)=(1/(1+(q/2)))*(-(q/2)*(-U_1(m,N_1-1,j)-U_1(m,N_1,j-2)+U_1(m,N_1-1,j-2))+((q^2)/2)*(U_1(m+1,N_1,j-1)-2*U_1(m,N_1,j-1)+U_1(m-1,N_1,j-1))+2*U_1(m,N_1,j-1)-U_1(m,N_1,j-2));
    end
    %Top
    for n=2:1:N_1-1
        U_1(M_1,n,j)=(1/(1+(q/2)))*(-(q/2)*(-U_1(M_1-1,n,j)-U_1(M_1,n,j-2)+U_1(M_1-1,n,j-2))+((q^2)/2)*(U_1(M_1,n+1,j-1)-2*U_1(M_1,n,j-1)+U_1(M_1,n-1,j-1))+2*U_1(M_1,n,j-1)-U_1(M_1,n,j-2));
    end
    %Two corner
    U_1(M_1,1,j)=(1/(-3-2*q))*(q*((-U_1(M_1-1,1,j)-U_1(M_1,1,j-2)+U_1(M_1-1,1,j-2))-(U_1(M_1,2,j)-U_1(M_1,2,j-2)+U_1(M_1,1,j-2)))-6*U_1(M_1,1,j-1)+3*U_1(M_1,1,j-2));%(-2,2) point
    U_1(M_1,N_1,j)=(1/(-3-2*q))*(q*((-U_1(M_1-1,N_1,j)-U_1(M_1,N_1,j-2)+U_1(M_1-1,N_1,j-2))+(-U_1(M_1,N_1-1,j)-U_1(M_1,N_1,j-2)+U_1(M_1,N_1-1,j-2)))-6*U_1(M_1,N_1,j-1)+3*U_1(M_1,N_1,j-2));%(2,2) point
end
c=1;
h_1=0.008;
k_1=0.004;
y_0=0.5;
q=(c*k_1)/(h_1);%here we let q=0.5
N_2=(4/h_1)+1;
M_2=(4/h_1)+1;
J_2=(T/k_1)+1;
delta=0.2;%wave source scale
d_2=0.2;%the width of slit
%build up the empty solution matrix;
U_2= zeros(M_2,N_2,J_2);
%Start setting initial condition to derive the first page of the solution
for m=1:1:M_2 %y
    for n=1:1:N_2%x
        r=sqrt((-2+((n-1)*h_1))^2+(((-2+(m-1)*h_1))-y_0)^2);
        if (r>=-delta) && (r<=delta)%locate the wave source
            U_2(m,n,1)=cos((pi*r)/(2*delta));
        end
    end
end
%Second page of the solution
for m=2:1:M_2-1%y
    for n=2:1:N_2-1%x
        U_2(m,n,2)=0.5*((q^2)*((U_2(m,n+1,1)-2*U_2(m,n,1)+U_2(m,n-1,1))+(U_2(m+1,n,1)-2*U_2(m,n,1)+U_2(m-1,n,1)))+2*U_2(m,n,1));%first time layer
    end
end
%Bottom solid wall 
for j=1:1:J_2
    for n=1:1:N_2
        U_2(1,n,j)=0;
    end
end
%Left solid wall
for j=1:1:J_2
    for m=1:1:M_2
        if (-2+(m-1)*h_1)<=0
            U_2(m,1,j)=0;
        end
    end
end
%right solid wall boundary
for j=1:1:J_2
    for m=1:1:M_2
        if (-2+(m-1)*h_1)<=0
            U_2(m,N_2,j)=0;
            m_location=m;%save it for upper solid wall
        end
    end
end
%Upper solid wall boundary
for j=1:1:J_2
    for n=1:1:N_2
        if (-2+(n-1)*h_1)<=-d_2 |  (-2+(n-1)*h_1)>=d_2
            U_2(m_location,n,j)=0;
        end
    end
end
%main calculation
for j=3:1:J_2
    for m=2:1:M_2-1
        for n=2:1:N_2-1
            if m==m_location
                if (-2+(n-1)*h_1)<=d_2 && (-2+(n-1)*h_1)>=-d_2
                    U_2(m,n,j)=(q^2)*((U_2(m,n+1,j-1)-2*U_2(m,n,j-1)+U_2(m,n-1,j-1))+(U_2(m+1,n,j-1)-2*U_2(m,n,j-1)+U_2(m-1,n,j-1)))+2*U_2(m,n,j-1)-U_2(m,n,j-2);
                end
            else
                U_2(m,n,j)=(q^2)*((U_2(m,n+1,j-1)-2*U_2(m,n,j-1)+U_2(m,n-1,j-1))+(U_2(m+1,n,j-1)-2*U_2(m,n,j-1)+U_2(m-1,n,j-1)))+2*U_2(m,n,j-1)-U_2(m,n,j-2);
            end                    
        end
    end
    %non-reflection boundary
    %left and right
    for m=m_location+1:1:M_2-1
        U_2(m,1,j)=(1/(1+(q/2)))*((q/2)*(U_2(m,2,j)-U_2(m,2,j-2)+U_2(m,1,j-2))+((q^2)/2)*(U_2(m+1,1,j-1)-2*U_2(m,1,j-1)+U_2(m-1,1,j-1))+2*U_2(m,1,j-1)-U_2(m,1,j-2));
        U_2(m,N_2,j)=(1/(1+(q/2)))*(-(q/2)*(-U_2(m,N_2-1,j)-U_2(m,N_2,j-2)+U_2(m,N_2-1,j-2))+((q^2)/2)*(U_2(m+1,N_2,j-1)-2*U_2(m,N_2,j-1)+U_2(m-1,N_2,j-1))+2*U_2(m,N_2,j-1)-U_2(m,N_2,j-2));
    end
    %Top
    for n=2:1:N_2-1
        U_2(M_2,n,j)=(1/(1+(q/2)))*(-(q/2)*(-U_2(M_2-1,n,j)-U_2(M_2,n,j-2)+U_2(M_2-1,n,j-2))+((q^2)/2)*(U_2(M_2,n+1,j-1)-2*U_2(M_2,n,j-1)+U_2(M_2,n-1,j-1))+2*U_2(M_2,n,j-1)-U_2(M_2,n,j-2));
    end
    %Two corner
    U_2(M_2,1,j)=(1/(-3-2*q))*(q*((-U_2(M_2-1,1,j)-U_2(M_2,1,j-2)+U_2(M_2-1,1,j-2))-(U_2(M_2,2,j)-U_2(M_2,2,j-2)+U_2(M_2,1,j-2)))-6*U_2(M_2,1,j-1)+3*U_2(M_2,1,j-2));%(-2,2) point
    U_2(M_2,N_2,j)=(1/(-3-2*q))*(q*((-U_2(M_2-1,N_2,j)-U_2(M_2,N_2,j-2)+U_2(M_2-1,N_2,j-2))+(-U_2(M_2,N_2-1,j)-U_2(M_2,N_2,j-2)+U_2(M_2,N_2-1,j-2)))-6*U_2(M_2,N_2,j-1)+3*U_2(M_2,N_2,j-2));%(2,2) point
end
max(max(U_1(:,:,J_1)))
max(max(U_2(:,:,J_2)))
error=abs(max(max(U_1(:,:,J_1)))-max(max(U_2(:,:,J_2))))