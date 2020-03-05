function xM = ikeda(n,x0V)
% xM = ikeda(n,x0V)
% Generates the ikeda map for the parameter values that
% give rise to chaos. Note that the two variable map is 
% actually a map on a complex variable. 
% INPUT
% n   : length for the generated Ikeda trajectory
% x0V : vector of initial values. If omitted a valid random vector
%       is chosen. (x0V has two components, the first is the real
%       part and the second the uimaginary part of a complex number).
% OUTPUT
% xM  : matrix of dimension n x 2, the real part of the complex variable
%       in the first column, the imaginary part in the second column.
if nargin == 1 
    x0limV = [-1 2 -0.3 0.3]'; 
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) (x0limV(4)-rand)/(x0limV(4)-x0limV(3))]';
else
    x0V = x0V(:);
end
ntrans = 100;
        
zV = NaN*ones(ntrans+n,1);
zV(1) = x0V(1) + i*x0V(2);
for k=2:n+ntrans
    zV(k) = 1+0.9*zV(k-1)*exp(0.4*i-6*i/(1+abs(zV(k-1))^2));
end
xM = [real(zV(ntrans+1:n+ntrans))'; imag(zV(ntrans+1:n+ntrans))']';
