function P = ProdEngine(i,j,s,Matrices)
% Perform matrix product from the optimal order, recursive engine
if i==j
    P = Matrices{i};
else
    k = s(j,i);
    P = ProdEngine(i,k-1,s,Matrices)*ProdEngine(k,j,s,Matrices);
end
end
