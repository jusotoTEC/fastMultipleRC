function expr = Prodexpr(i,j,s,vnames)
% Return the string expression of the optimal order 
if i==j
    expr = vnames{i};
else
    k = s(j,i);
    expr = ['(' Prodexpr(i,k-1,s,vnames) '*' Prodexpr(k,j,s,vnames) ')'];
end
end % Prodexpr