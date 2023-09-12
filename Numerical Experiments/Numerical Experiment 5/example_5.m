function example_distributed_1()
    clc; clear; close all
    
    tam=50:25:200; p=2; s=1000; 
    
    h=length(tam);
    
    time1=zeros(h,1); time2=zeros(h,1);
    e1=zeros(h,1); e2=zeros(h,1);
    
    k=0;
    
    for m=tam
    
        X=randn(m,s);
        Exx=(1/s)*(X*X');

        Enn=zeros(m,m,p); A=zeros(m,m,p); 
        A_d=[]; Enn_d=[]; Exz=[];

        for i=1:p
            sigma=rand(1)/4;
            Enn(:,:,i)=sigma^2*eye(m);
            A(:,:,i)=rand(m);
            A_d=blkdiag(A_d,A(:,:,i));
            Enn_d=blkdiag(Enn_d,Enn(:,:,i));        
            Exz=[Exz Exx*(A(:,:,i))'];
        end

        Ezz=A_d*kron(ones(p),Exx)*A_d'+Enn_d;

        M=Exz*pinv(sqrtm(Ezz)); 

        Ezz_sqrtm=sqrtm(Ezz);

        C=zeros(m,m*p,p);

        B=zeros(m,m,p);

        for i=1:p
           C(:,:,i)=Ezz_sqrtm((i-1)*m+1:m*i,:);
           B(:,:,i)=eye(m);
        end

        k=k+1; r=round(3*m/4);

        tic
        [~,error1]=multiple_low_rank_opt1(M,B,C,r);
        t1=toc;
        time1(k)=t1;
        e1(k)=error1;
        tic
        [~,error2]=multiple_low_rank_opt2(M,B,C,r);
        t2=toc;   
        time2(k)=t2;
        e2(k)=error2;
    end
    
    figure
    hold on
    area(tam, time1)
    area(tam, time2)
    title('Time')
    
    figure
    hold on
    plot(tam, e1,'r')
    plot(tam, e2,'b')
    title('Error')
    
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
    
    X=zeros(m,m,p);
    for i=1:p
        X(:,:,i)=rand(m,r)*rand(r,m);
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
    
    X=zeros(m,m,p);
    for i=1:p
        X(:,:,i)=rand(m,r)*rand(r,m);
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