%%                                                                       %%
%%     Solve 1D Laplace equation u_t-u_xx=f(x,t) in [a,b]x[0,T]          %%
%%                                                                       %%

%%
clear 
clc
close all
format long
tic
%% Input
a=0.0;
b=1.0;
time=0.5;
kappa = 1;   % Choose a number of kappa
% kappa = 1/16;   % Choose a number of kappa

N=6;         % Number of control volume
T=300;       % Number of time iterior
M=4;         % Number of iteration when refine mesh
ll=zeros(M,1);
normh1=zeros(M,1);      % norrm h1;
norml2_x=zeros(M,1);    % norm l2;
normh1_x=zeros(M,1);    % norrm h1;

%%
for jj=1:M
    %% Create the mesh point
    dx=(b-a)/N;
    k=time/T;
    x=zeros(N+1,1);   % x_(i+1/2)
    t=0:k:time;
    for i_iterx=1:N+1
        x(i_iterx)=a+(i_iterx-1)*dx;
    end
    
    %% create control point
    x_cp=zeros(N+2,1);  % x(i)
    for i_iterx=1:N+2
        if(i_iterx==1)
            x_cp(i_iterx)=x(i_iterx);
        else
            if(i_iterx==N+2)
                x_cp(i_iterx)=x(i_iterx-1);
            else
                x_cp(i_iterx)=(x(i_iterx-1)+x(i_iterx))/2.0;
            end
        end
    end
    
    %% Create form for computing
    A = zeros(N,N);  % Matrix A
    u = zeros(N,T);  % Matrix value of u
    u0 = zeros(N,1); % u(x,0)
    phi_a = zeros(1,T);  %u(a,t)
    phi_b = zeros(1,T);  %u(b,t)
    
    %% Input iterior
    for i_iterx=1:N
        u0(i_iterx) = u_ex(x_cp(i_iterx+1),0);
    end
    
    for i_itert = 1:T
        phi_a(i_itert) = u_ex(a,t(i_itert));
    end
    
    for i_itert = 1:T
        phi_b(1,i_itert) = u_ex(b,t(i_itert));
    end
    
    %% Matrix of function f
    F = zeros(N,T);
    for i_itert =1:T
        for i_iterx =1:N
            if(i_iterx ==1)
                alpha = kappa/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+1)-x_cp(i_iterx)));
                F(i_iterx,i_itert) = (fun(x(i_iterx+1),t(i_itert))+fun(x(i_iterx),t(i_itert)))/2 + alpha*phi_a(1,i_itert); 
            elseif (i_iterx ==N)
                gamma = kappa/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+2)-x_cp(i_iterx+1)));
                F(i_iterx,i_itert) = (fun(x(i_iterx+1),t(i_itert))+fun(x(i_iterx),t(i_itert)))/2 + gamma*phi_b(1,i_itert); 
            else
                F(i_iterx,i_itert) = (fun(x(i_iterx+1),t(i_itert))+fun(x(i_iterx),t(i_itert)))/2;

            end
        end
    end
    
    %% Matrix A
    for i_iterx = 1:N
        alpha = kappa/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+1)-x_cp(i_iterx)));
        gamma = kappa/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+2)-x_cp(i_iterx+1)));
        if(i_iterx ==1)
            A(i_iterx,i_iterx) = -alpha - gamma;
            A(i_iterx,i_iterx+1) = gamma;
        elseif(i_iterx == N)
            A(i_iterx,i_iterx) = -alpha - gamma;
            A(i_iterx,i_iterx-1) = alpha;
        else
            A(i_iterx,i_iterx) = -alpha - gamma;
            A(i_iterx,i_iterx-1) = alpha;
            A(i_iterx,i_iterx+1) = gamma;
        end
    end
    
    
    %% Forward Euler method
    u(:,1) = u0(:,1);
    for i_itert = 2:T
        u(:,i_itert) = (eye(N,N)+k*A)*u(:,i_itert-1) + k*F(:,i_itert-1);        
    end        
    %% Create a matrix for exact solution

    u_exact = zeros(N+2,T);
    t=0:k:time;
    for i_itert=1:T
        for i_iterx = 1:N+2
            u_exact(i_iterx,i_itert) = u_ex(x_cp(i_iterx),t(i_itert));
        end
    end
    
    %% Create a matrix for approximate solution
    u_dis = zeros(N+2,T);
    t=0:k:time;
    for i_itert=1:T
        for i_iterx = 1:N+2
            if(i_iterx == 1)
                u_dis(i_iterx,i_itert) = phi_a(1,i_itert);
            elseif(i_iterx == N+2)
                u_dis(i_iterx,i_itert) = phi_b(1,i_itert);
            else
                u_dis(i_iterx,i_itert) = u(i_iterx - 1, i_itert);
            end
        end
    end
    figure
    plot(x_cp,u_dis(:,T),'b-*',x_cp,u_exact(:,T),'r');
    title(['Compare exact and approximate solution at h = ',num2str(dx), ',k =' ,num2str(k)]);
    legend('Approximate solution', 'Exact solution')


     %% Calculate error for all period
    normh0 = 0;
    for i_itert =1:T
       for i_iterx = 1:N
           normh0 = normh0 +(u_dis(i_iterx+1,i_itert)-u_exact(i_iterx+1,i_itert)-(u_dis(i_iterx,i_itert)-u_exact(i_iterx,i_itert)))^2/(x_cp(i_iterx+1)-x_cp(i_iterx));
       end
    end
    normh1(jj) = normh1(jj) + sqrt(normh0)*k;
    normh1(jj) = sqrt(normh1(jj));
       
    %% Calculate error independent on t
    for i_iterx = 1:N
        norml2_x(jj) = norml2_x(jj) +(u_dis(i_iterx+1,i_itert)-u_exact(i_iterx+1,i_itert))^2*(x(i_iterx+1)-x(i_iterx));
    end
    norml2_x(jj) = sqrt(norml2_x(jj));
    
    for i_iterx = 1:N
        normh1_x(jj) = normh1_x(jj) +(u_dis(i_iterx+1,i_itert)-u_exact(i_iterx+1,i_itert)-(u_dis(i_iterx,i_itert)-u_exact(i_iterx,i_itert)))^2/(x_cp(i_iterx+1)-x_cp(i_iterx));
    end
    normh1_x(jj) = sqrt(normh1_x(jj));
    
    %% Refine mesh
    ll(jj) = N*T;
    N=N*2;
    T=T*2;
end

%% Plot error for all period
figure
plot(log(ll.^(1/2)),-log(normh1),'b',log(ll.^(1/2)),log(ll.^(1/2)),'r');
title('Error in L2 and H1 norms');
legend('H1 Norm','x' )

%% Plot error independent on t
figure
plot(log(ll.^(1/2)),-log(norml2_x),'r',log(ll.^(1/2)),2*log(ll.^(1/2)),'black',log(ll.^(1/2)),-log(normh1_x),'b',log(ll.^(1/2)),1.5*log(ll.^(1/2)),'c');
title('Error in L2 and H1 norms independent on t');
legend('L2 Norm', '2x', 'H1 Norm','3x/2' )

toc












