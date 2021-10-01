%part A solid wall
clear
T=5;
c=4;
h=0.004;
k=0.001;
q=(c*k)/(h);%here we let q=1
N=(4/h)+1;
L=(T/k)+1;
delta=0.1;%wave source length
%build up the empty solution matrix;
U= zeros(L,N);
%Start setting initial condition
for n=1:1:N
    if ((-2+((n-1)*h))>=-delta) && ((-2+((n-1)*h))<=delta)
        U(1,n)=cos((pi*(-2+(n-1)*h))/(2*delta));
    end
end

for n=2:1:N-1  
    U(2,n)=0.5*((q^2)*(U(1,n+1)-2*U(1,n)+U(1,n-1))+2*U(1,n));%first time layer
end
for j=3:1:L
    U(j,1)=0;%Boundary condition
    U(j,N)=0;%Boundary condition
    for n=2:1:N-1
            U(j,n)=(q^2)*(U(j-1,n+1)-2*U(j-1,n)+U(j-1,n-1))+2*U(j-1,n)-U(j-2,n);%interior points
    end
end

surf(-2:h:2,0:k:T,U)
shading flat
