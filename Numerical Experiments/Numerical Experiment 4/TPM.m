function T = TPM(Z)
    [m,n] = size(Z);   
    alpha=5*10^-4;
    if m>=n
        T=(Z'*Z+alpha*eye(n))\Z';
    else
        T=((Z*Z'+alpha*eye(m))\Z)';
    end
end
