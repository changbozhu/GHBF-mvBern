function  y = sgm(x, varargin)
% Multivariate Sigmoid function
% INPUT:
%   x       -- a column vector
%   zeta    -- an optional scale parameter
%
if nargin > 1
    alpha = varargin{1};
    [m,n] = size(alpha);
    if m~=n && (  (m==1)||( n==1)  )
        alpha = diag(alpha);
    end    
else
    alpha = eye(length(x));
end
if size(x, 1) < size(x, 2)
    x=x';
end
y = 1 ./(1 + exp(-alpha*x));
return;
end

