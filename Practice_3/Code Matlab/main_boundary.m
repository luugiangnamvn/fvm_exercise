% Solve 1D Laplace equation -uxx=f(x) in [a,b]
clear all
clc
close all
ax=0.0;
bx=1.0;
ay=0.0;
by=1.0;

N=12;% Number of control volume
K=10;
M=1;% Number of iteration when refine mesh
% norml2=zeros(M,1); % norm l2;
% normh1=zeros(M,1); % norrm h1

ll=zeros(M,1);
for jj=1:M
    dx=(bx-ax)/N;
    dy=(by-ay)/K;
    % Create the mesh point
    x=zeros(N+1,1);
    y=zeros(K+1,1);
    for i_iterx=1:N+1
        x(i_iterx)=ax+(i_iterx-1)*dx;
    end
    for i_itery=1:K+1
        y(i_itery)=ay+(i_itery-1)*dy;
    end
    
    % create control point
    
    x_cp=zeros(N+2,1);
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
    y_cp=zeros(K+2,1);
    for i_itery=1:K+2
        if(i_itery==1)
            y_cp(i_itery)=y(i_itery);
        else
            if(i_itery==K+2)
                y_cp(i_itery)=y(i_itery-1);
            else
                y_cp(i_itery)=(y(i_itery-1)+y(i_itery))/2.0;
            end
        end
    end

    % Creare the Matrix
    B=zeros(N,N);
    C=zeros(N,N);
    D=zeros(N,N);
    ai=zeros(1,N);
    bi=zeros(1,N);
    cj=zeros(1,K);
    dj=zeros(1,K);
    for i_iterx=1:N
        a1=-1/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+1)-x_cp(i_iterx)));
        ai(1,i_iterx)=a1;
        b1=-1/((x(i_iterx+1)-x(i_iterx))*(x_cp(i_iterx+2)-x_cp(i_iterx+1)));
        bi(1,i_iterx)=b1;
        for i_itery=1:K            
            c1=-1/((y(i_itery+1)-y(i_itery))*(y_cp(i_itery+1)-y_cp(i_itery)));
            cj(1,i_itery)=c1;
            d1=-1/((y(i_itery+1)-y(i_itery))*(y_cp(i_itery+2)-y_cp(i_itery+1)));       
            dj(1,i_itery)=d1;
        end
    end
end
A=zeros(N*K,N*K);
for j=1:K
    if(j==1)
        B=zeros(N,N);
        for i=1:N
            if(i==1)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i+1)=-bi(1,i);
            elseif (i==N)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
            else
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
                B(i,i+1)=-bi(1,i);
            end
        end
        D=zeros(N,N);
        for i=1:N
            D(i,i)=-dj(1,j);
        end
        A((j-1)*N+1:j*N ,(j-1)*N+1:j*N)=B;
        A((j-1)*N+1:j*N ,j*N+1:(j+1)*N)=D;
    elseif(j==K)
        B=zeros(N,N);
        for i=1:N
            if(i==1)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i+1)=-bi(1,i);
            elseif (i==N)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
            else
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
                B(i,i+1)=-bi(1,i);
            end
        end
        C=zeros(N,N);
        for i=1:N
            C(i,i)=-cj(1,j);
        end
        A((j-1)*N+1:j*N ,(j-1)*N+1:j*N)=B;
        A((j-1)*N+1:j*N ,(j-2)*N+1:(j-1)*N)=C;
    else
        B=zeros(N,N);
        for i=1:N
            if(i==1)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i+1)=-bi(1,i);
            elseif (i==N)
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
            else
                B(i,i)=ai(1,i)+bi(1,i)+cj(1,j)+dj(1,j);
                B(i,i-1)=-ai(1,i);
                B(i,i+1)=-bi(1,i);
            end
        end
        C=zeros(N,N);
        for i=1:N
            C(i,i)=-cj(1,j);
        end
        D=zeros(N,N);
        for i=1:N
            D(i,i)=-dj(1,j);
        end
        A((j-1)*N+1:j*N ,(j-1)*N+1:j*N)=B;
        A((j-1)*N+1:j*N ,(j-2)*N+1:(j-1)*N)=C;
        A((j-1)*N+1:j*N ,j*N+1:(j+1)*N)=D;
    end
end

    % Create vector b
    b=zeros(N*K,1);
    for j=1:K
        for i=1:N
        b((j-1)*N+i)=(f(x(i),y(j))+f(x(i),y(j+1))+f(x(i+1),y(j))+f(x(i+1),y(j+1)))/4.0; % Trepozoidal rule  
        end
    end
%     a0=-1/((x(2)-x(1))*(x_cp(2)-x_cp(1)));
%     bN=-1/((x(N+1)-x(N))*(x_cp(N+2)-x_cp(N+1)));
%     b(1) = (f(x(1))+f(x(2)))/2.0 - a0*u_exact(ax);
%     b(N) = (f(x(N))+f(x(N+1)))/2.0 - bN*u_exact(bx);
% 
    u=zeros(N*K,1);
    u=A\b;
%     
    u_ex=zeros((N+2)*(K+2),1);
    for j=1:K+2
        for i=1:N+2
            u_ex((j-1)*(N+2)+i)=u_exact(x_cp(i),y_cp(j));
        end
    end
    
    u_dis=zeros((N+2),K+2);
%     u_dis(1,:)=u_ex(1,:);
%     u_dis(N+2,:)=u_ex((N+2),:);
%     u_dis(1,:)=u_ex(1,:);
%     u_dis(N+2,:)=u_ex((N+2),:);
    for j=1:K+2
        if(j==1)
            for i=1:N+2
                if(i==1)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                elseif (i==N+2)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                end
            end
        elseif (j==K+2)
            for i=1:N+2
                if(i==1)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                elseif (i==N+2)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                end
            end
        else
            for i=1:N+2
                if(i==1)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                elseif (i==N+2)
                    u_dis(i,j) = u_exact(x_cp(i),y_cp(j));
                else
                    u_dis(i,j) = u((j-2)*N+i-1);
                end
            end
        end
    end
    [X,Y]=meshgrid(x_cp,y_cp);  
    figure
%     surf(x_cp,y_cp,u_dis,'red',x_cp,y_cp,u_ex);
    surf(X,Y,u_dis');
    legend('Discrese solution','Exact solution')
%     norm12=zeros(N*K,1);
%     for j=1:K
%         for i=1:N
%             norml2((j-1)*N+i)=norml2((j-1)*N+i)+(u((j-1)*N+i)-u_ex((j-1)*N+i)^2*(x(i_iter+1)-x(i_iter));
%         end
%     end
%     norml2(jj)=sqrt(norml2(jj));
%     
%     for i_iter=1:N+1
%         normh1(jj)=normh1(jj)+((u_dis(i_iter+1)-u_ex(i_iter+1))-(u_dis(i_iter)-u_ex(i_iter)))^2/(x_cp(i_iter+1)-x_cp(i_iter));
%     end
%     
%     normh1(jj)=sqrt(normh1(jj));
%     
%     ll(jj)=N;
%     
%     N=2*N;
% 
% figure
% plot(log(ll),-log(norml2),'r', log(ll), -log(normh1),'blue', log(ll), 2*log(ll)+1.5,'green');
% title('Error');
% legend('L^2 Norm', 'H^1 norm', '2x')