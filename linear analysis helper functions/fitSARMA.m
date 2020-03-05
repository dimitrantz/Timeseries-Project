function [nrmseV,phiV,thetaV,SDz,aicS,fpeS,sarmamodel]=fitSARMA(xV,p,q,ps,qs,s,Tmax)
% [nrmseV,phiV,thetaV,phisV,thetasV,SDz,aicS,fpeS,sarmamodel]=fitSARMA(xV,p,q,ps,qs,s,
% Tmax)
% FITSARMA fits a seasonal autoregressive moving average (ARMA) model and
% computes the fitting error (normalized root mean square error) for a
% given number of steps ahead. 
% The SARMA(p,q)x(ps,qs)_s model has the form
% x(t) = phi(0) + phi(1)*x(t-1) + ... + phi(p)*x(t-p) + 
%        +phi(s)*x(t-s) + ... + phi(s+p-1)*x(t-p-s+1) + ... 
%        +phi(ps*s)*x(t-ps*s) + ... + phi(ps*s+p-1)*x(t-p-ps*s+1) + ... 
%        +z(t) - theta(1)*z(t-1) + ... + theta(q)*z(t-p) - 
%        -theta(s)*x(t-s) - ... - theta(s+p-1)*x(t-p-s+1) - ... 
%        -theta(qs*s)*x(t-qs*s) - ... - phi(qs*s+p-1)*x(t-p-qs*s+1) - ... 
% z(t) ~ WN(0,sdnoise^2).
% INPUTS:
%  xV      : vector of the scalar time series
%  p       : the order of the AR part of the model.
%  q       : the order of the MA part of the model.
%  ps      : the order of the seasonal AR part of the model.
%  qs      : the order of the seasonal MA part of the model.
%  s       : the period of the seasonal component.    
%  Tmax    : the prediction horizon, the fit error is computed for
%            T=1...Tmax steps ahead.
% OUTPUT: 
%  nrmseV  : vector of length Tmax, the nrmse of the fit for T-mappings,
%            T=1...Tmax. 
%  phiV    : the coefficients of the estimated AR part (of length
%            (ps*s+p+1) with phi(0) as first component. Note that these are
%            the coefficients of both AR and seasonal AR.
%  thetaV  : the coefficients of the estimated MA part (of length qs*s+q). 
%            Note that these are the coefficients of both MA and seasonal MA.
%  SDz     : the standard deviation of the noise term.
%  aicS    : the AIC value for the model.
%  fpeS    : the FPE value for the model.
%  sarmamodel : the model structure (contains all the above apart from
%               nrmseV)

if nargin==6
    Tmax = 1;
elseif nargin==5
    Tmax = 1;
    s=1;
elseif nargin==4
    Tmax = 1;
    s=1;
    qs=0;
elseif nargin==3
    Tmax = 1;
    s=1;
    qs=0;
    ps=0;
elseif nargin==2
    Tmax = 1;
    s=1;
    qs=0;
    ps=0;
    q=0;
end
if isempty(p)
    p = 0;
end
if isempty(q)
    q = 0;
end
if isempty(ps)
    ps = 0;
end
if isempty(qs)
    qs = 0;
end

if p>s && ps>0
    error('The AR order is larger than the period s. Use an ARMA model');
end
if q>s && qs>0
    error('The MA order is larger than the period s. Use an ARMA model');
end

xV = xV(:);
n = length(xV);
mx = mean(xV);
xxV = xV-mx;

sarmamodel = sarma(xxV,p,q,ps,qs,s);
if p==0 && ps==0
    phiV = [];
else
    phi0 = (1+sum(sarmamodel.a(2:1+ps*s+p)))*mx;
    phiV = [phi0 -sarmamodel.a(2:1+ps*s+p)];
    rootarV = roots(sarmamodel.a);
    if any(abs(rootarV)>=1)
        fprintf('The estimated AR(%d) part of the SARMA model is not stationary.\n',p);
    end
end
if q==0 && qs==0
    thetaV = [];
else
    thetaV = -sarmamodel.c(2:1+qs*s+q);
    rootmaV = roots(sarmamodel.c);
    if any(abs(rootmaV)>=1)
        fprintf('The estimated MA(%d) part of the SARMA model is not reversible.\n',q);
    end
end
SDz = sqrt(sarmamodel.NoiseVariance);
aicS = aic(sarmamodel);
fpeS = sarmamodel.EstimationInfo.FPE;
nrmseV = NaN*ones(Tmax,1);
for T=1:Tmax
    tmpS = predict(sarmamodel,xxV,T);
    if myversion==0
        xpreV = tmpS+mx;   
    elseif myversion==1
        xpreV = tmpS.OutputData+mx;   
    else
        xpreV = tmpS{1}+mx;
    end
    nrmseV(T) = nrmse(xV(q+1:n),xpreV(q+1:n));
end

