%part B semiexplicit
clear
N=180;
L=1;
T=5;
h=L/(N);
k=0.001;
r=k/(h)^2;
A=[];%Stiffness matrix A
B=[];%RHS matrix B
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

%Assembly the stiffness matrix A
for n=1:1:(L/h)-1
    C=[];%a row vector
    for j=1:1:(L/h)-1
        if n==j-1%assign element
            C=[C,-(r/(pi^2))];
        elseif n==j%assign element
            C=[C,(2+((2*r)/(pi^2)))];
        elseif n==j+1%assign element
            C=[C,-(r/(pi^2))];
        else%assign element
            C=[C,0];
        end
    end%assign row
    A=[A;C];
end

%Assembly the Right hand side matrix B
for n=1:1:(L/h)-1
    C=[];%a row vector
    for j=1:1:(L/h)-1
        if n==j-1%assign element
            C=[C,(r/(pi^2))];
        elseif n==j%assign element
            C=[C,(2-((2*r)/(pi^2)))];
        elseif n==j+1%assign element
            C=[C,(r/(pi^2))];
        else%assign element
            C=[C,0];
        end
    end%assign row
    B=[B;C];
end

tic
%Solve the matrix system
for j=1:1:(T/k)
    %Calculate the solution
    RHS=(B*(U(:,j)))';%produce RHS formed by C-N method
    U_temp=tridiag(A,RHS)';
    %Use Euler method to approximate the boundary.
    U_temp_boundary=[0;U_temp;0];%add boundary condition
    U=[U,U_temp];
    U_boundary=[U_boundary, U_temp_boundary];%The matrix after assembly the boundary conditions
end
toc

%x=0.1 error
C_numerical=[];
n_1=0.1/h;
C_numerical=U_boundary(n_1+1,end);
error_1=log10(abs(C_numerical-(exp(-5)*sin(pi*0.1))))
%x=0.5 error
C_numerical=[];
n_2=0.5/h;
C_numerical=U_boundary(n_2+1,end);
error_2=log10(abs(C_numerical-(exp(-5)*sin(pi*0.5))))

%plot
%waterfall(0:h:L,0:k:T,U_exact_boundary);
waterfall(0:k:T,0:h:L,U_boundary);
    