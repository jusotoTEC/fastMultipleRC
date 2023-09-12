function [NN] = LQQT(L,r,q)
    [~,n]=size(L);   
        Y2=randn(n,r);
    for i=1:q
        Y1=L*Y2;
        Y2=L'*Y1;
    end
    [Q1,~]=qr(Y2,0);
    NN=mmtimes(L,Q1,Q1');
    end