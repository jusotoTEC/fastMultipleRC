function [P, orderstruct] = mmtimes(varargin)
% P = mmtimes(M1, M2, ... Mn)
%   return a chain matrix product P = M1*M2* ... *Mn
%
% {Mi} are matrices with compatible dimension: size(Mi,2) = size(Mi+1,1)
% 
% Because the matrix multiplication is associative; the chain product can
% be carried out with different order, leading to the same result (up to
% round-off error). MMTIMES uses "optimal" order of binary product to
% reduce the computational effort (probably the accuracy is also improved).
%
% The function assumes the cost of the product of (m x n) with (n x p)
% matrices is (m*n*p). This assumption is typically true for full matrix.
%
% Notes:
%   Scalar matrix are groupped together, and the rest will be
%   multiplied with optimal order.
%
%   To get the the structure that stores the best order, call with the
%   second outputs:
%   >> [P orderstruct] = mmtimes(M1, M2, ... Mn);
%   % This structure can be used later if the input matrices have the
%   % same sizes as those in the first call (but with different contents)
%   >> P = mmtimes(M1, M2, ... Mn, orderstruct);
%
% See also: mtimes
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Orginal: 19-Jun-2010
%          20-Jun-2010: quicker top-down algorithm
%          23-Jun-2010: treat the case of scalars
%          16-Aug-2010: passing optimal order as output/input argument
Matrices = varargin;
buildexpr = false;
if ~isempty(Matrices) && isstruct(Matrices{end})
    orderstruct = Matrices{end};
    Matrices(end) = [];
else
    % Detect scalars
    iscst = cellfun('length',Matrices) == 1;
    if any(iscst)
        % scalars are multiplied apart
        cst = prod([Matrices{iscst}]);
        Matrices = Matrices(~iscst);
    else
        cst = 1;
    end
    % Size of matrices
    szmats = [cellfun('size',Matrices,1) size(Matrices{end},2)];
    s = MatrixChainOrder(szmats);
    orderstruct = struct('cst', cst, ...
                         's', s, ...
                         'szmats', szmats);
                     
    if nargout>=2
        % Prepare to build the string expression
        vnames = arrayfun(@inputname, 1:nargin, 'UniformOutput', false);
        % Default names, e.g., M1, M2, ..., for inputs that is not single variable 
        noname = cellfun('isempty', vnames);
        vnames(noname) = arrayfun(@(i) sprintf('M%d', i), find(noname), 'UniformOutput', false);
        if any(iscst)
            % String '(M1*M2*...)' for constants
            cstexpr = strcat(vnames(iscst),'*');
            cstexpr = strcat(cstexpr{:});
            cstexpr = ['(' cstexpr(1:end-1) ')'];
        else
            cstexpr = '';
        end
        vnames = vnames(~iscst);
        buildexpr = true;
    end
end
if ~isempty(Matrices)
    P = ProdEngine(1,length(Matrices),orderstruct.s,Matrices);    
    if orderstruct.cst~=1
        P = orderstruct.cst*P;
    end
    if buildexpr
        expr = Prodexpr(1,length(Matrices),orderstruct.s,vnames);
        if ~isempty(cstexpr)
            % Concatenate the constant expression in front
            expr = [cstexpr '*' expr];
        end
        orderstruct.expr = expr;
    end
else
    P = orderstruct.cst;
    if nargout>=2
        orderstruct.expr = cstexpr;       
    end
end
end % mmtimes
