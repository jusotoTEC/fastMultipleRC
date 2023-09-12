function [s, qmin] = MatrixChainOrder(szmats)
% Find the best ordered chain-product, the best splitting index
% of M(i)*...*M(j) is stored in s(j,i) of the array s (only the lower
% part is filled)
% Top-down dynamic programming, complexity O(n^3)
n = length(szmats)-1;
s = zeros(n);
pk = szmats(2:n);
ij = (0:n-1)*(n+1)+1;
left = zeros(1,n-1);
right = zeros(1,n-1);
L = 1;
while true % off-diagonal offset
    q = zeros(size(pk));
    for j=1:n-L % this is faster and BSXFUN or product with DIAGONAL matrix
        q(:,j) = (szmats(j)*szmats(j+L+1))*pk(:,j);
    end
    q = q + left + right;
    [qmin, loc] = min(q, [], 1);
    s(ij(1:end-L)+L) = (1:n-L)+loc;
    
    if L<n-1
        pk = [pk(:,1:end-1);
              pk(end,2:end)];
        left = [left(:,1:end-1);
                qmin(1:end-1)];
        right = [qmin(2:end);
                 right(:,2:end)];
        L = L+1;
    else
        break
    end % if
end % while-loop
end % MatrixChainOrder
