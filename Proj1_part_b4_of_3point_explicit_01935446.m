%part B modified 
clear
N=100;
L=1;
T=5;
h=L/(N);
k=0.002;
r=(k)/(h)^2;
A=[];%Stiffness matrix A
U=[];%Solution matrix U
%build up the initial condition vector;
U_0=[];
U_exact=[];
U_exact_0=[];

%Start of calculation
for n=1:1:(L/h)-1
    U_0=[U_0;sin(pi*n*h)];
end
%Set the starting solution
U=U_0;
U_boundary=[0;U_0;0];%set boundary

%use explicit method once to approximate row no.2, ie t=k
U_temp=[];
for n=2:1:L/h
    U_temp_n=((r/(pi^2))*U_boundary(n+1,1))+(1-((2*r)/(pi^2)))*U_boundary(n,1)+(r/(pi^2))*U_boundary(n-1,1);
    U_temp=[U_temp;U_temp_n];
end
U=[U,U_temp];
U_temp_boundary=[0;U_temp;0];
U_boundary=[U_boundary,U_temp_boundary];
%rewrite matrix A
A=[];
%Assembly the stiffness matrix A
for n=1:1:(L/h)-1
    C=[];%a row vector
    for j=1:1:(L/h)-1
        if n==j-1%assign element
            C=[C,-(2*r/(pi^2))];
        elseif n==j%assign element
            C=[C,(3+((4*r)/(pi^2)))];
        elseif n==j+1%assign element
            C=[C,-(2*r/(pi^2))];
        else%assign element
            C=[C,0];
        end
    end%assign row
    A=[A;C];
end

%Solve the matrix system
for j=2:1:(T/k)
    %Calculate the solution
    U_temp=tridiag(A,(4*U(:,j)-U(:,j-1)))';
    %Use Euler method to approximate the boundary.
    U_temp_boundary=[0;U_temp;0];%add boundary condition
    U=[U,U_temp];
    U_boundary=[U_boundary, U_temp_boundary];%The matrix after assembly the boundary conditions
end


%x=0.5 error
C_numerical=[];
n_2=0.5/h;
C_numerical=U_boundary(n_2+1,end);
error_2=log10(abs(C_numerical-(exp(-5)*sin(pi*0.5))))


C_exact=(exp(-(5))*sin((0:h:L)*pi))
error_explicit=U_boundary(:,end)-C_exact'
%plot(0:h:L,error_explicit')
%hold on

waterfall(0:k:T,0:h:L,U_boundary);
    