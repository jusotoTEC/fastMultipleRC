function X=low_rank_M2(A,B,C,r)
    Bp=TPM(B);
    Cp=TPM(C);
    M=mmtimes(B,Bp,A,Cp,C);    
    Mr=LQQT(M,r,3);
    X=mmtimes(Bp,Mr,Cp);    
end