function A = wish(h,n)
% Command:  s=wish(h,n)
% Purpose:  Draws an m x m matrix from a wishart distribution
%           with scale matrix h and degrees of freedom nu = n.
%           This procedure uses Bartlett's decomposition.
% Inputs:   h     -- m x m scale matrix.
%           n     -- scalar degrees of freedom.
% Outputs:  s     -- m x m matrix draw from the wishart
%                    distribution.
% Note: Parameterized so that mean is n*h

%[L U]=lu(h);
%A = randn(size(h,1),n);
%A = U'*A*A'*U;
%A = chol(h)'*randn(size(h,1),n);
A=zeros(size(h,1),n);
for i=1:n
A(:,i) = mvnrnd(zeros(size(h,1),1),h,1);
end
A = A*A';