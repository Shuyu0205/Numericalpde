%Solve the problem with 1000*1000 points
clear
T=1;
c=4;
h=0.004;
k=0.001;
q=(c*k)/(h);%here we let q=1
N=(4/h)+1;
L=(T/k)+1;
delta=0.5;%wave source length
%build up the empty solution matrix;
U_1= zeros(L,N);
%Start setting initial condition
for n=1:1:N
    if ((-2+((n-1)*h))>=-delta) && ((-2+((n-1)*h))<=delta)
        U_1(1,n)=cos((pi*(-2+(n-1)*h))/(2*delta));
    end
end

for n=2:1:N-1  
    U_1(2,n)=0.5*((q^2)*(U_1(1,n+1)-2*U_1(1,n)+U_1(1,n-1))+2*U_1(1,n));%first time layer
end
for j=3:1:L
    U_1(j,1)=(1/(1+q))*((q^2)*(U_1(j-1,2)-2*U_1(j-1,1)+U_1(j-1,2)+(1/q)*U_1(j-2,1))+2*U_1(j-1,1)-U_1(j-2,1));%left boundary use Forward scheme 
    U_1(j,N)=(1/(1+q))*((q^2)*(U_1(j-1,N-1)-2*U_1(j-1,N)+U_1(j-1,N-1)+(1/q)*U_1(j-2,N))+2*U_1(j-1,N)-U_1(j-2,N));%right boundary use backward scheme
    for n=2:1:N-1
            U_1(j,n)=(q^2)*(U_1(j-1,n+1)-2*U_1(j-1,n)+U_1(j-1,n-1))+2*U_1(j-1,n)-U_1(j-2,n);%interior points
    end
end
%double the number of points
h_2=0.002;%We change the resoluation
k_2=0.0005;
q=(c*k_2)/(h_2);%here we let q=1
N=(4/h_2)+1;
L=(T/k_2)+1;
delta=0.5;%wave source length
%build up the empty solution matrix;
U_2= zeros(L,N);
%Start setting initial condition
for n=1:1:N
    if ((-2+((n-1)*h_2))>=-delta) && ((-2+((n-1)*h_2))<=delta)
        U_2(1,n)=cos((pi*(-2+(n-1)*h_2))/(2*delta));
    end
end

for n=2:1:N-1  
    U_2(2,n)=0.5*((q^2)*(U_2(1,n+1)-2*U_2(1,n)+U_2(1,n-1))+2*U_2(1,n));%first time layer
end
for j=3:1:L
    U_2(j,1)=(1/(1+q))*((q^2)*(U_2(j-1,2)-2*U_2(j-1,1)+U_2(j-1,2)+(1/q)*U_2(j-2,1))+2*U_2(j-1,1)-U_2(j-2,1));%left boundary use Forward scheme 
    U_2(j,N)=(1/(1+q))*((q^2)*(U_2(j-1,N-1)-2*U_2(j-1,N)+U_2(j-1,N-1)+(1/q)*U_2(j-2,N))+2*U_2(j-1,N)-U_2(j-2,N));%right boundary use backward scheme
    for n=2:1:N-1
            U_2(j,n)=(q^2)*(U_2(j-1,n+1)-2*U_2(j-1,n)+U_2(j-1,n-1))+2*U_2(j-1,n)-U_2(j-2,n);%interior points
    end
end
%double the number of points
h_3=0.001;%We change the resoluation
k_3=0.00025;
q=(c*k_3)/(h_3);%here we let q=1
N=(4/h_3)+1;
L=(T/k_3)+1;
delta=0.5;%wave source length
%build up the empty solution matrix;
U_3= zeros(L,N);
%Start setting initial condition
for n=1:1:N
    if ((-2+((n-1)*h_3))>=-delta) && ((-2+((n-1)*h_3))<=delta)
        U_3(1,n)=cos((pi*(-2+(n-1)*h_3))/(2*delta));
    end
end

for n=2:1:N-1  
    U_3(2,n)=0.5*((q^2)*(U_3(1,n+1)-2*U_3(1,n)+U_3(1,n-1))+2*U_3(1,n));%first time layer
end
for j=3:1:L
    U_3(j,1)=(1/(1+q))*((q^2)*(U_3(j-1,2)-2*U_3(j-1,1)+U_3(j-1,2)+(1/q)*U_3(j-2,1))+2*U_3(j-1,1)-U_3(j-2,1));%left boundary use Forward scheme 
    U_3(j,N)=(1/(1+q))*((q^2)*(U_3(j-1,N-1)-2*U_3(j-1,N)+U_3(j-1,N-1)+(1/q)*U_3(j-2,N))+2*U_3(j-1,N)-U_3(j-2,N));%right boundary use backward scheme
    for n=2:1:N-1
            U_3(j,n)=(q^2)*(U_3(j-1,n+1)-2*U_3(j-1,n)+U_3(j-1,n-1))+2*U_3(j-1,n)-U_3(j-2,n);%interior points
    end
end
%Both solutions at x=-2
max(U_1(:,1))
max(U_2(:,1))
error_1=abs(max(U_1(:,1))-max(U_2(:,1)))
%Both solutions at x=-1.5
max(U_1(:,(0.5/h)+1))
max(U_2(:,(0.5/h_2)+1))
error_2=abs(max(U_1(:,(0.5/h)+1))-max(U_2(:,(0.5/h_2)+1)))
%Both solutions at x=0
max(U_1(:,(2/h)+1))
max(U_2(:,(2/h_2)+1))
error_3=abs(max(U_1(:,(2/h)+1))-max(U_2(:,(2/h_2)+1)))
%Both solutions at x=1.5
max(U_1(:,(3.5/h)+1))
max(U_2(:,(3.5/h_2)+1))
error_4=abs(max(U_1(:,(3.5/h)+1))-max(U_2(:,(3.5/h_2)+1)))
%Both solutions at x=2
max(U_1(:,(4/h)+1))
max(U_2(:,(4/h_2)+1))
error_5=abs(max(U_1(:,(4/h)+1))-max(U_2(:,(4/h_2)+1)))
plot(0:k:T,U_1(:,(3.5/h)+1))
hold on
plot(0:k_2:T,U_2(:,(3.5/h_2)+1))
hold on
plot(0:k_3:T,U_3(:,(3.5/h_3)+1))

