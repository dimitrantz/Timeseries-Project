function [xV,phiallV,thetaallV]= generateSARMAts(phiV,thetaV,phisV,thetasV,s,n,sdnoise)
% xV = generateSARMAts(phiV,thetaV,phisV,thetasV,s,n,sdnoise)
% Generate a SARMA(p,q)x(ps,qs)_s time series of length 'n' with Gaussian
% input noise. 
% INPUTS
% - phiV    : the coefficients of the AR part, 
%             phiV = [phi(0) phi(1) ... phi(p)]' and phi(0) is the constant
%             term. The coefficients define the AR delay polynomial
%             phi(B) = 1-phi(1)*B-phi(2)*B^2-...-phi(p)*B^p
%             and the mean of the time series is
%             mx=phi(0)/(1-phi(1)-phi(2)-...-phi(p))
% - thetaV  : the coefficients of the MA part, [theta(1) ... theta(q)]'. 
%             The coefficients define the MA delay polynomial
%             theta(B) = 1-theta(1)*B-theta(2)*B^2-...-theta(q)*B^q
% - phisV   : the coefficients of the seasonal AR part, 
%             phisV = [phis(1) ... phis(ps)]'. The coefficients define the 
%             seasonal AR delay polynomial
%             phis(B_s) = 1-phis(1)*B_s-phis(2)*B_s^2-...-phi(ps)*B_s^ps
% - thetasV : the coefficients of the seasonal MA part, 
%             thetasV = [thetas(1) ... thetas(qs)]'. The coefficients
%             define the seasonal MA delay polynomial
%             thetas(B_s) = 1-thetas(1)*B_s-thetas(2)*B_s^2-...-theta(qs)*B_s^qs
% - s       : the period (seasonality).
% - n       : the length of the time series to be generated
% - sdnoise : the SD of the input noise (if left out then sdnoise=1).
% The generating SARMA(p,q)x(ps,qs)_s process reads in delay polynomial
% form
% phi(B)*phis(B_s)*x(t) = theta(B)*thetas(B_s)*z(t), z(t) ~ WN(0,sdnoise^2)

if nargin==6
    sdnoise=1;
end
phiV = phiV(:);
phisV = phisV(:);
thetaV = thetaV(:);
thetasV = thetasV(:);
p = length(phiV)-1;
q = length(thetaV);
ps = length(phisV);
qs = length(thetasV);
pbig = 1+ps*s+p;
pbig = 1+ps*s+p;

% Form the AR delay polynomial of SARMA from the convolution of the two AR
% delay polynomials for delay one and delay s.
if p==0 && ps==0
    phiallV = [];
else
    phi1V = [1;-phiV(2:end)];
    phi2V = zeros(ps*s+1,1);
    phi2V(1) = 1;
    phi2V(1+s:s:ps*s+1)=-phisV;
    phi3V = conv(phi1V,phi2V);
    phiallV = -phi3V(2:end);
    pall = length(phiallV);
    phiallV = [phiV(1);phiallV];  
    rootarV = roots([1;-phiallV(2:end)]);
    if any(abs(rootarV)>=1)
        fprintf('The AR(%d) part of the SARMA process is not stationary.\n',pall);
    end
end
% The same for MA.
if q==0 && qs==0
    thetaallV = [];
else
    theta1V = [1;-thetaV];
    theta2V = zeros(qs*s+1,1);
    theta2V(1) = 1;
    theta2V(1+s:s:qs*s+1)=-thetasV;
    theta3V = conv(theta1V,theta2V);
    thetaallV = -theta3V(2:end);
    qall = length(thetaallV);
    rootmaV = roots([1;-thetaallV]);
    if any(abs(rootmaV)>=1)
        fprintf('The MA(%d) part of the SARMA model is not reversible.\n',qall);
    end
end
pq = max(pall,qall);
ntrans = 100+pq;
x0V = sdnoise*randn(pq,1);
zV = randn(n+ntrans,1) * sdnoise;
xV = NaN*ones(n+ntrans,1);
xV(1:pq) = x0V;
if p==0 && ps==0
    for i=pq+1:n+ntrans
        xV(i)=phiallV(1)+zV(i)-thetaallV'*flipud(zV(i-qall:i-1));
    end
elseif q==0 && qs==0    
    for i=pq+1:n+ntrans
        xV(i)=phiallV(1)+phiallV(2:end)'*flipud(xV(i-pall:i-1))+zV(i);
    end
else
    for i=pq+1:n+ntrans
        xV(i)=phiallV(1)+phiallV(2:end)'*flipud(xV(i-pall:i-1))+zV(i)-thetaallV'*flipud(zV(i-qall:i-1));
    end
end
xV = xV(ntrans+1:n+ntrans);
