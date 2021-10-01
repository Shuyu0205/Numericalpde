%part C implicit
clear
N=600;
L=6;
T=0.1;
h=L/(N);
D=0.005;%can be changed
U_constant=0.8;%can be changed
k=0.001;
miu=k/h;
r=k/(h)^2;
A=[];%Stiffness matrix A
U=[];%Solution matrix U
%build up the initial condition vector;
U_0=[];
U_exact=[];
U_exact_0=[];
x_0=1;
%Start of calculation
for n=1:1:(L/h)-1
    U_0=[U_0;exp((-(n*h-x_0)^2)/D)];
end
U_boundary=[exp((-(0-x_0)^2)/D);U_0;exp((-(L-x_0)^2)/D)];%set boundary
%Set the starting solution
U_0(1)=U_0(1)+((((U_constant*miu)/2)+D*r)*(1/sqrt(4*k+1))*exp(-(0-x_0-U_constant*k)^2/(D*(4*k+1))));%boundary process
U_0(end)=U_0(end)-(((U_constant*miu)/2)-D*r)*(1/sqrt(4*k+1))*exp(-(L-x_0-U_constant*k)^2/(D*(4*k+1)));%allocate boundary condition
U=U_0;


%Assembly the stiffness matrix A
for n=1:1:(L/h)-1
    C=[];%a row vector
    for j=1:1:(L/h)-1
        if n==j-1%assign element
            C=[C,-(((U_constant*miu)/2)+D*r)];
        elseif n==j%assign element
            C=[C,((2*D*r)+1)];
        elseif n==j+1%assign element
            C=[C,(((U_constant*miu)/2)-D*r)];
        else%assign element
            C=[C,0];
        end
    end%assign row
    A=[A;C];
end

%Solve the matrix system
for j=1:1:(T/k)
    %Calculate the solution
    U_temp=tridiag(A,(U(:,j)))';
    U_temp_boundary=[(1/sqrt(4*j*k+1))*exp(-(0-x_0-U_constant*j*k)^2/(D*(4*j*k+1))); U_temp;(1/sqrt(4*j*k+1))*exp(-(L-x_0-U_constant*j*k)^2/(D*(4*j*k+1)))];%add boundary condition
    U_temp(1)=U_temp(1)+(((U_constant*miu)/2)+D*r)*(1/sqrt(4*(j+1)*k+1))*exp(-((0-x_0-U_constant*(j+1)*k)^2)/(D*(4*(j+1)*k+1)));%boundary process
    U_temp(end)=U_temp(end)-(((U_constant*miu)/2)-D*r)*(1/sqrt(4*(j+1)*k+1))*exp(-(L-x_0-U_constant*(j+1)*k)^2/(D*(4*(j+1)*k+1)));%boundary process
    %Use Euler method to approximate the boundary.
    U=[U,U_temp];
    U_boundary=[U_boundary, U_temp_boundary];%The matrix after assembly the boundary conditions
end

%x=L/2 error
C_numerical=[];
n_2=(L/2)/h;
C_numerical=U_boundary(n_2+1,end);
error=log10(abs(C_numerical-(1/sqrt(4*T+1))*exp(-((L/2)-x_0-U_constant*T)^2/(D*(4*T+1)))))

%Uncomment and comment the waterfall to plot x-C, x-error.
%C_exact=(1/sqrt(4*T+1))*exp(-((0:h:L)-x_0-U_constant*T).^2/(D*(4*T+1)));
%error_C=U_boundary(:,end)-C_exact';
%subplot(1,2,1)
%plot(0:h:L,error_C)
%subplot(1,2,2)
%plot(0:h:L,U_boundary(:,end))

%plot solution
waterfall(0:k:T,0:h:L,U_boundary);