% Solve 1D Laplace equation -uxx=f(x) in [a,b]
clear all
clc
close all
ax=0.0;
bx=1.0;
N=4;% Number of control volume

M=6;% number of iteration when refine mesh
norml2=zeros(M,1); % norm l2;
normh1=zeros(M,1); % norrm h1

% lambda0=1;
% lambda1=0;
lambda0=1/2;
lambda1=11/2;
% lambda0=pi/2;
% lambda1=pi/2;

ll=zeros(M,1);


for jj=1:M
    dx=(bx-ax);
    
    % Create the mesh point
    x=zeros(N+1,1);
    for i_iter=1:N+1
        x(i_iter)=ax+(1-cos((i_iter-1)*pi/(2*N)))*dx;
    end
    
    % create control point
    
    x_cp=zeros(N+2,1);
    for i_iter=1:N+2
        if(i_iter==1)
            x_cp(i_iter)=x(i_iter);
        else
            if(i_iter==N+2)
                x_cp(i_iter)=x(i_iter-1);
            else
                x_cp(i_iter)=2/3*x(i_iter-1)+1/3*x(i_iter);
            end
        end
    end

    % Creare the Matrix
    A=zeros(N,N);
    for i_iter=1:N
        a1=-1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+1)-x_cp(i_iter)));
        b1=-1/((x(i_iter+1)-x(i_iter))*(x_cp(i_iter+2)-x_cp(i_iter+1)));
        if(i_iter==1)
            A(i_iter,i_iter+1)=b1;
            A(i_iter,i_iter)=a1/(1+lambda0*(x_cp(i_iter+1)-x_cp(i_iter)))-(a1+b1);
        else
            if(i_iter==N)
                A(i_iter,i_iter-1)=a1;
                A(i_iter,i_iter)=b1/(1+lambda1*(x_cp(i_iter+2)-x_cp(i_iter+1)))-(a1+b1);
            else
                A(i_iter,i_iter-1)=a1;
                A(i_iter,i_iter+1)=b1;
                A(i_iter,i_iter)=-(a1+b1);
            end
        end
    end
    
    % Create vector b
    b=zeros(N,1);
    for i_iter=1:N
%        b(i_iter)=f((x(i_iter)+x(i_iter+1))/2.0);% Midpoint rule
        b(i_iter)=(f(x(i_iter))+f(x(i_iter+1)))/2.0; % Trepozoidal rule          
    end
    
    u=zeros(N,1);
    u=A\b;
    
    u_ex=zeros(N+2,1);
    for i_iter=1:N+2
        u_ex(i_iter)=u_exact(x_cp(i_iter));
    end
    
    u_dis=zeros(N+2,1);
    u_dis(1)=u_ex(1);
    u_dis(N+2)=u_ex(N+2);
    for i_iter=1:N
        u_dis(i_iter+1)=u(i_iter);
    end
    
    figure
    plot(x_cp,u_dis,'red',x_cp,u_ex);
    
    for i_iter=1:N
        norml2(jj)=norml2(jj)+(u_dis(i_iter+1)-u_ex(i_iter+1))^2*(x(i_iter+1)-x(i_iter));
    end
    norml2(jj)=sqrt(norml2(jj));
    
    for i_iter=1:N+1
        normh1(jj)=normh1(jj)+((u_dis(i_iter+1)-u_ex(i_iter+1))-(u_dis(i_iter)-u_ex(i_iter)))^2/(x_cp(i_iter+1)-x_cp(i_iter));
    end
    
    normh1(jj)=sqrt(normh1(jj));
    
    ll(jj)=N;
    
    N=2*N;
end

figure
plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll),1.5*log(ll)+2, 'black', log(ll), 2*log(ll)+1.5,'green');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')