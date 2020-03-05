function xM = henon(n,x0V)
% xM = henon(n,x0V)
% Generates the henon map for the parameter values that
% give rise to chaos.
% INPUT
% n   : length for the generated Henon trajectory
% x0V : vector of initial values. If omitted a valid random vector
%       is chosen.
% OUTPUT
% xM  : matrix of dimension n x 2, first variable in the first 
%       column, second variable in the second column.

if nargin == 1 
    x0limV = [-1 1 -0.3 0.3]'; 
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) (x0limV(4)-rand)/(x0limV(4)-x0limV(3))]';
else
    x0V = x0V(:);
end
aV = [0.3 -1.4]'; % Henon parameters for chaos
ntrans = 100;
        
xM = NaN*ones(ntrans+n,2);
xM(1,:) = x0V';
for i=2:n+ntrans
    xM(i,1) = 1+aV(2)*xM(i-1,1)^2+xM(i-1,2);
    xM(i,2) = aV(1)*xM(i-1,1);
end
xM = xM(ntrans+1:n+ntrans,:);
