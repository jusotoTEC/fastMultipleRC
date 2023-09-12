% Example 2: Execution time in secons versus matrix dimension

clc; clear; close all

inter=200:200:4000;

time1=zeros(1,length(inter)); time2=zeros(1,length(inter)); 
error1=zeros(1,length(inter)); error2=zeros(1,length(inter));

k=0;
for m=inter
   k=k+1;
   n=2*m;
   K=randn(m,n);
   r=m/8;
   %SVD Method
   tic
   [U1,S1,V1]=svd(K);
   K1=U1(:,1:r)*S1(1:r,1:r)*(V1(:,1:r))';
   t1=toc;
   time1(k)=t1;
   
   %Proposed Method
   tic
   Y2=randn(n,r);
   for j=1:3
       Y1=K*Y2;
       Y2=K'*Y1;
   end
   [Qr,~]=qr(Y2,0);   
   K2=K*(Qr*Qr');
   t2=toc;
   time2(k)=t2;
   
   error1(k)=norm(K-K1,'fro');
   error2(k)=norm(K-K2,'fro');
end

%Diagrams 
figure
hold on
area(inter,time1)
area(inter,time2)
grid on
xlabel('Dimension (m)')
ylabel('Execution Time (s)')

figure
hold on
plot(inter,error1,'r')
plot(inter,error2,'b')
grid on
xlabel('Dimension (m)')
ylabel('Error')
