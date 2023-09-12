function X=low_rank_M1(A,B,C,r)
    Bp=pinv(B);
    Cp=pinv(C);
    M=B*Bp*A*Cp*C;
    [U,S,V]=svd(M);    
    X=Bp*U(:,1:r)*S(1:r,1:r)*(V(:,1:r))'*Cp;
end