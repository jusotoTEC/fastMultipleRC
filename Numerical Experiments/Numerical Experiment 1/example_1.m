% Example 1: Execution time in secons versus matrix dimension

clc; clear; close all

dimension=[100 500:500:5000];

time1=zeros(1,length(dimension)); time2=zeros(1,length(dimension)); error=zeros(1,length(dimension));

k=0;
for m=dimension
   k=k+1;
   r=m/4;
   D=rand(m,r)*rand(r,m);   
   %SVD Method
   tic
   P1=pinv(D);
   t1=toc;
   time1(k)=t1;
   
   %Proposed Method
   tic
   alpha=0.05;
   P2=(D'*D+alpha*eye(m))\D';
   t2=toc;
   time2(k)=t2;
   
   error(k)=norm(P1-P2,'fro');
end

%Diagrams 
hold on
area(dimension,time1)
area(dimension,time2)
grid on
xlabel('Dimension (m)')
ylabel('Execution Time (s)')

figure
area(dimension,error)
grid on
xlabel('Dimension (m)')
ylabel('Error')
