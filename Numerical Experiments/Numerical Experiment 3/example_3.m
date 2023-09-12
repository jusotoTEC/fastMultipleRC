function example_general_1()
    clc; clear; close all
    
    intervalo=20:20:400;
    time1=zeros(length(10:10:50),1); time2=zeros(length(10:10:50),1);
    e1=zeros(length(10:10:50),1); e2=zeros(length(10:10:50),1);
    k=0;    
    for m=intervalo
        m
        k=k+1;
        r=floor(m/4); p=2; 
        A=rand(m); B=rand(m,m/2,p); C=rand(m/2,m,p); 
        Aux=norm(A,'fro');
        tic
        [~,error1]=multiple_low_rank_opt1(A,B,C,r);
        t1=toc;
        time1(k)=t1;
        e1(k)=error1/Aux;
        tic
        [~,error2]=multiple_low_rank_opt2(A,B,C,r);
        t2=toc;   
        time2(k)=t2;
        e2(k)=error2/Aux;
    end
    
    figure
    hold on
    area(intervalo, time1)
    area(intervalo, time2)
    
    figure
    hold on
    plot(intervalo, e1,'r')
    plot(intervalo, e2,'b')
    
end

function Y=pinv_new(X)
    [m,n]=size(X);
    if m>n
        Y=(X'*X)\X';
    else
        Y=X'/(X*X');
    end
end

function Y=GBRP(X,r)
   n=size(X,2);
   Y2=randn(n,r);
   for j=1:3
       Y1=X*Y2;
       Y2=X'*Y1;
   end
   [Qr,~]=qr(Y2,0);   
   Y=X*(Qr*Qr');
end

function Y=low_rank_opt1(A,B,C,r)
    Bp=pinv(B); Cp=pinv(C); 
    T=B*Bp*A*Cp*C;
    [U,S,V]=svd(T);
    Y=Bp*U(:,1:r)*S(1:r,1:r)*(V(:,1:r))'*Cp;
end

function Y=low_rank_opt2(A,B,C,r)
    Bp=pinv_new(B); Cp=pinv_new(C); 
    T=B*Bp*A*Cp*C;
    Tr=GBRP(T,r);
    Y=Bp*Tr*Cp;
end

function y=obj_func(A,B,C,X)
    p=size(B,3);    
    Aux=A;
    for i=1:p
        Aux=Aux-B(:,:,i)*X(:,:,i)*C(:,:,i);
    end
    y=norm(Aux,'fro')^2;
end

function [X,error1]=multiple_low_rank_opt1(A,B,C,r)
    iterMax=1000;    
    m=size(A,1); p=size(B,3);
    
    X=zeros(m/2,m/2,p);
    for i=1:p
        X(:,:,i)=rand(m/2,r)*rand(r,m/2);
    end
    
    e_old=obj_func(A,B,C,X);    
    
    for k=1:iterMax
        for i=1:p
            At=A;
            for j=1:p
                if j~=i
                    At=At-B(:,:,j)*X(:,:,j)*C(:,:,j);
                end
            end
            X(:,:,i)=low_rank_opt1(At,B(:,:,i),C(:,:,i),r);
        end
        e_new=obj_func(A,B,C,X);
        error1=abs(e_new-e_old);
        if error1<10^-5
            break
        end       
        e_old=e_new;
    end        
end

function [X,error1]=multiple_low_rank_opt2(A,B,C,r)
    iterMax=1000;    
    m=size(A,1); p=size(B,3);
    
    X=zeros(m/2,m/2,p);
    for i=1:p
        X(:,:,i)=rand(m/2,r)*rand(r,m/2);
    end
    
    e_old=obj_func(A,B,C,X);    
    
    for k=1:iterMax
        for i=1:p
            At=A;
            for j=1:p
                if j~=i
                    At=At-B(:,:,j)*X(:,:,j)*C(:,:,j);
                end
            end
            X(:,:,i)=low_rank_opt2(At,B(:,:,i),C(:,:,i),r);
        end
        e_new=obj_func(A,B,C,X);
        error1=abs(e_new-e_old);
        if error1<10^-5
            break
        end       
        e_old=e_new;
    end        
end