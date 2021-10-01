%part A
clear
N=100;
L=1;
T=5;
h=L/(N);
k=0.0002;
r=k/(h)^2;
A=[];%Stiffness matrix A
U=[];%Solution matrix U
%build up the initial condition vector;
U_0=[];

%Start of calculation
for n=1:1:(L/h)-1
    U_0=[U_0;sin(pi*n*h)];
end
%Set the starting solution
U=[0;U_0;0];%set boundary
tic
for j=1:1:(T/k)
    %Calculate the solution
    U_temp=[];
    for n=1:1:(L/h)
        if n==1
            U_temp_n=0;
        else
            U_temp_n=((r/(pi^2))*U(n+1,j))+(1-((2*r)/(pi^2)))*U(n,j)+(r/(pi^2))*U(n-1,j);
        end
        U_temp=[U_temp;U_temp_n];      
    end
    U_temp=[U_temp;0];
    U=[U,U_temp];
end
toc
%x=0.1 error
C_numerical=[];
n_1=0.1/h;
C_numerical=U(n_1+1,end);
C_exact_1=exp(-(0:k:T))*sin(pi*0.1);
error_1=log10(abs(C_numerical-(exp(-5)*sin(pi*0.1))))
%x=0.5 error
C_numerical=[];
n_2=0.5/h;
C_numerical=U(n_2+1,end);
C_exact_2=exp(-(0:k:T))*sin(pi*0.5);
error_2=log10(abs(C_numerical-(exp(-5)*sin(pi*0.5))))


waterfall(0:k:T,0:h:L,U);


