clc; clear; close all

folder='imagen_nasa_saturno';

A1=imread([folder '\saturno (1).jpg']);
[m,n]=size(A1);
a=dir([folder '\*.jpg']);
num_img=size(a,1);
O=zeros(m,n,num_img); N=zeros(m,n,num_img);
A=zeros(m*n,num_img); C=zeros(m*n,num_img);

sigma=0.2;
for i=1:num_img
    B=im2double(imread([folder '\saturno (' num2str(i) ').jpg']));
    O(:,:,i)=B; A(:,i)=reshape(O(:,:,i),[m*n 1]);
    N(:,:,i)=B+sigma*randn(m,n); C(:,i)=reshape(N(:,:,i),[m*n 1]);
end

m_img=m*n; n_img=num_img;


porc=0.75; B=eye(m*n);
r=round(min([m_img n_img])*porc);

display('Start Methods')

tic
F1=low_rank_M2(A,B,C,r);
time1=toc;
tic
F2=low_rank_M1(A,B,C,r);
time2=toc;

display('End Methods')

e1=norm(A-F1*C,'fro');
e2=norm(A-F2*C,'fro');

num_noise=809;
op=2;
if op==1
    img_orig=A(:,num_noise);
    img_noise=img_orig+sigma*randn(m*n,1);
elseif op==2
    img_orig=A(:,num_noise);
    img_noise=C(:,num_noise)+0.22*randn(m*n,1);
else
    Aux1=im2double(imread('img_new.jpg'));
    img_orig=double(reshape(Aux1,[m*n 1]));
    img_noise=img_orig+sigma^2*randn(m*n,1);
end

X1=F1*img_noise;
X2=F2*img_noise;


img_o=reshape(img_orig,[m,n]);
figure
imshow(img_o)
title('Original')

img_n=reshape(img_noise,[m,n]);
figure
imshow(img_n)
title('With Noise')

img1=reshape(X1,[m,n]);
figure
imshow(img1)
title('Reconstruction with SVD')
xlabel(['MSE = ' num2str(norm(img_o-img1,'fro')) ', time = ' num2str(time1)])

img2=reshape(X2,[m,n]);
figure
imshow(img2)
title('Reconstruction with New Method')
xlabel(['MSE = ' num2str(norm(img_o-img2,'fro')) ', time = ' num2str(time2)])

%%%Speedup, Prcent Difference, Error and SSIM
speedup=time2./time1;
percent_difference=100*((time2-time1)./(time2));
ssim_glrma=ssim(img_o, img1);
ssim_fast_glrma=ssim(img_o, img2);
mse_glrma=norm(img_o-img1,'fro')^2;
mse_fast_glrma=norm(img_o-img1,'fro')^2;
save('example_1.mat','speedup','percent_difference','ssim_glrma','ssim_fast_glrma','mse_glrma','mse_fast_glrma')
