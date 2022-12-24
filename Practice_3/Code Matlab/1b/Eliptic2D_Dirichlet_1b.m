%%                                                                       %%
%   Code duoi day se giai quyet bai toan vv -u_xx(x,y)-u_yy(x,y)=f(x,y) 
%   bang Phuong phap the tich huu han voi bien Dirichlet khong thuan nhat
%%                                                                       %%
clear all
close all
clc
format long
%% Thong tin dau vao, Input.
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;
N = 8; % So cach chia tren truc Ox
M = 9; % So cach chia tren truc Oy
iteration=4;
ll = zeros(iteration,1);
norml2=zeros(iteration,1);
normh1=zeros(iteration,1);

%% Xay dung ma tran, qua trinh duoc lap lai 4 lan.
for jj=1:iteration
    
    dx = (bx-ax)/N;
    dy = (by-ay)/M;
    
    
    %% Tao luoi  
    x=zeros(N+1,1);
    for i=1:N+1
        x(i)=(i-1)*dx;
    end
    
    x_cp=zeros(N+2,1);
    for i=1:N+2
        if(i==1)
            x_cp(i)=x(i);
        else
            if(i==N+2)
                x_cp(i)=x(i-1);
            else
                x_cp(i)=(x(i-1)+x(i))/2.0;
            end
        end
    end
    
    y = zeros(M+1,1);
    for j=1:M+1
        y(j)=(j-1)*dy;
    end
    
    y_cp=zeros(M+2,1);
    for j=1:M+2
        if(j==1)
            y_cp(j)=y(j);
        else
            if(j==M+2)
                y_cp(j)=y(j-1);
            else
                y_cp(j)=(y(j-1)+y(j))/2.0;
            end
        end
    end
    
    %% Tao ma tran A
    A=zeros(N*M,N*M);
    B=sparse(N,N);
    C=sparse(N,N);
    D=sparse(N,N); 
    for j=1:M
        c1 = -1/((y(j+1) - y(j))*(y_cp(j+1) - y_cp(j)));
        for i=1:N
            C(i,i)=c1;
        end
        d1 = -1/((y(j+1) - y(j))*(y_cp(j+2) - y_cp(j+1)));
        for i=1:N
            D(i,i)=d1;
        end
        for i=1:N
            a1 = -1/((x(i+1) - x(i))*(x_cp(i+1) - x_cp(i)));
            b1 = -1/((x(i+1) - x(i))*(x_cp(i+2) - x_cp(i+1)));
            if(i==1)
                B(i,i)=-(a1+ b1 + c1 + d1);
                B(i,i+1)=b1;
            else
                if(i==N)
                    B(i,i)=-(a1+ b1 + c1 + d1);
                    B(i,i-1)=a1;
                else
                    B(i,i)=-(a1+ b1 + c1 + d1);
                    B(i,i-1)=a1;
                    B(i,i+1)=b1;
                end
            end
        end
        if(j==1)
            A((j-1)*N+1:j*N,(j-1)*N+1:j*N)=B;
            A((j-1)*N+1:j*N,j*N+1:(j+1)*N)=D;
        else
            if(j==M)
                A((j-1)*N+1:j*N,(j-1)*N+1:j*N)=B;
                A((j-1)*N+1:j*N,(j-2)*N+1:(j-1)*N)=C;
            else
                A((j-1)*N+1:j*N,(j-1)*N+1:j*N)=B;
                A((j-1)*N+1:j*N,j*N+1:(j+1)*N)=D;
                A((j-1)*N+1:j*N,(j-2)*N+1:(j-1)*N)=C;
            end
        end
    end
    
    %% Tao vector ve phai F
    F=zeros(N*M,1);
    for j=1:M
        c0 = -1/((y(j+1) - y(j))*(y_cp(j+1) - y_cp(j)));
        d0 = -1/((y(j+1) - y(j))*(y_cp(j+2) - y_cp(j+1)));       
        if(j==1)
            for i =1:N
                a0 = -1/((x(i+1) - x(i))*(x_cp(i+1) - x_cp(i)));
                b0 = -1/((x(i+1) - x(i))*(x_cp(i+2) - x_cp(i+1)));
                if(i==1)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - c0*u_exact(x_cp(i+1),0) - a0*u_exact(0,y_cp(j+1));
                elseif(i==N)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - c0*u_exact(x_cp(i+1),0) - b0*u_exact(1,y_cp(j+1));
                else
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - c0*u_exact(x_cp(i+1),0);
                end
            end
        elseif(j==M)
            for i =1:N
                a0 = -1/((x(i+1) - x(i))*(x_cp(i+1) - x_cp(i)));
                b0 = -1/((x(i+1) - x(i))*(x_cp(i+2) - x_cp(i+1)));
                if(i==1)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - d0*u_exact(x_cp(i+1),1) - a0*u_exact(0,y_cp(j+1));
                elseif(i==N)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - d0*u_exact(x_cp(i+1),1) - b0*u_exact(1,y_cp(j+1));
                else
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - d0*u_exact(x_cp(i+1),1) ;
                end
            end
        else
            for i =1:N
                a0 = -1/((x(i+1) - x(i))*(x_cp(i+1) - x_cp(i)));
                b0 = -1/((x(i+1) - x(i))*(x_cp(i+2) - x_cp(i+1)));
                if(i==1)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - a0*u_exact(0,y_cp(j+1));
                elseif(i==N)
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4 - b0*u_exact(1,y_cp(j+1));
                else
                    F((j-1)*N+i)=(f(x(i),y(j)) + f(x(i + 1),y(j)) + f(x(i),y(j+1)) + f(x(i+1),y(j+1)))/4;
                end
            end
        end
    end
   
    %% Tim nghiem xap xi u
    u=A\F;
    
    
     %% Tao ma tran cho ket qua chinh xac
    u_ex=zeros(M+2,N+2);
    for j=1:M+2
        for i=1:N+2
            u_ex(j,i)=u_exact(x_cp(i),y_cp(j));
        end
    end
    
    %% Tao ma tran chua ket qua xap xi va bien
    u_dis=u_ex;
    for j=1:M
        for i=1:N
            u_dis(j+1,i+1)=u((j-1)*N+i);
        end
    end
    %% Ve nghiem xap xi
     figure
     surf(x_cp,y_cp,u_dis)

    %% Tinh sai so

    for i=1:N
        for j = 1:M
            norml2(jj)=norml2(jj)+(u_dis(j+1,i+1)-u_ex(j+1,i+1))^2*(x(i+1)-x(i))*(y(j+1)-y(j));
        end
    end
    norml2(jj)=sqrt(norml2(jj));
    for i=1:N
        for j=1:M
            normh1(jj)=normh1(jj)+((u_dis(j+1,i)-u_ex(j+1,i))-(u_dis(j,i)-u_ex(j,i)))^2*(x(i+1)-x(i))/(y_cp(j+1)-y_cp(j)) + ((u_dis(j,i+1)-u_ex(j,i+1))-(u_dis(j,i)-u_ex(j,i)))^2*(y_cp(j+1)-y_cp(j))/(x(i+1)-x(i));
        end
    end
    
    normh1(jj)=sqrt(normh1(jj));
    ll(jj) = N*M;
    %% Refine mesh
    N=2*N;
    M =2*M;
end

%% Ve bac sai so va so sanh voi duong thang 2x, 3/2x
figure
plot(log(ll.^(1/2)),-log(norml2),'r', log(ll.^(1/2)), -log(normh1),'blue', log(ll.^(1/2)),1.5*log(ll.^(1/2))+2, 'black', log(ll.^(1/2)), 2*log(ll.^(1/2))+2,'green');
title('Error');
legend('L^2 Norm', 'H^1 norm', '3/2x', '2x')