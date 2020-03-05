function [nrmseV,preM,phiV,thetaV] = predictSARMAnrmse(xV,p,q,ps,qs,s,Tmax,nlast,tittxt)
% [nrmseV,preM,phiV,thetaV] = predictSARMAnrmse(xV,p,q,ps,qs,s,Tmax,nlast,tittxt)
% PREDICTSARMANRMSE makes predictions with a SARMA model on a last part
% of a given time series and computes the prediction error (NRMSE measure)
% for T-step ahead predictions. The SARMA(p,q)x(ps,qs)_s model is
% x(t) = phi(0) + phi(1)*x(t-1) + ... + phi(p)*x(t-p) + 
%        +phi(s)*x(t-s) + ... + phi(s+p-1)*x(t-p-s+1) + ... 
%        +phi(ps*s)*x(t-ps*s) + ... + phi(ps*s+p-1)*x(t-p-ps*s+1) + ... 
%        +z(t) - theta(1)*z(t-1) + ... + theta(q)*z(t-p) - 
%        -theta(s)*x(t-s) - ... - theta(s+p-1)*x(t-p-s+1) - ... 
%        -theta(qs*s)*x(t-qs*s) - ... - phi(qs*s+p-1)*x(t-p-qs*s+1) - ... 
% z(t) ~ WN(0,sdnoise^2).
% Note that if ps=0 and qs=0, the model reduces to ARMA(p,q).
% If also q=0 (and  ps=0 and qs=0) it reduces to AR(p) (autoregressive
% model of order p), and if p=0 (and  ps=0 and qs=0) it reduces to MA(q)
% (moving average model of order q).
% INPUTS:
%  xV      : vector of the scalar time series
%  p       : the order of AR part of the model.
%  q       : the order of MA part of the model.
%  ps      : the order of the seasonal AR part of the model.
%  qs      : the order of the seasonal MA part of the model.
%  s       : the period of the seasonal component.    
%  Tmax    : the predictions in the test set are repeated for each of the 
%            prediction steps T=1...Tmax
%  nlast   : the size of the test set to compute the prediction error on.
%          : If not specified, it is half the length of the time series
%  tittxt  : string to be displayed in the title of the figure.
%            If not specified, no plot is made
% OUTPUT: 
%  nrmseV  : vector of length Tmax, the nrmse for the predictions for time
%            steps T=1...Tmax, on the test set.
%  preM    : matrix of nlast columns and Tmax rows, having the T-ahead 
%            predictions at column T, T=1...Tmax (the first T-1 components
%            are NaN).
%  phiV    : the coefficients of the estimated AR part (of length
%            (ps*s+p+1) with phi(0) as first component. Note that these are
%            the coefficients of both AR and seasonal AR.
%  thetaV  : the coefficients of the estimated MA part (of length qs*s+q). 
%            Note that these are the coefficients of both MA and seasonal MA.

n = length(xV);
xV = xV(:);
if nargin==8
    tittxt = [];
elseif nargin==7
    tittxt = [];
    nlast = round(n/2);
end
if isempty(nlast)
    nlast = round(n/2);
end
if nlast>=n-2*q,
    error('test set is too large for the given time series!')
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

n1 = n-nlast;  % size of training set
x1V = xV(1:n1); 
mx1 = mean(x1V(1:n1)); 
xx1V = x1V(1:n1)-mx1; % set mean of the training set to zero.
sarmamodel = sarma(xx1V,p,q,ps,qs,s);
pall = p+ps*s;
qall = q+qs*s;
if pall==0
    phiV = [];
else
    phi0 = (1+sum(sarmamodel.a(2:1+pall)))*mx1;
    phiV = [phi0;-sarmamodel.a(2:pall+1)']; % Note that the AR coefficients are for the centered time series.
    rootarV = roots(sarmamodel.a);
    if any(abs(rootarV)>=1)
        fprintf('The estimated AR(%d) part of the SARMA model is not stationary.\n',pall);
    end
end
if qall==0
    thetaV = [];
else
    thetaV = -sarmamodel.c(2:qall+1)'; % Note that the MA coefficients are for the centered time series.
    rootmaV = roots(sarmamodel.c);
    if any(abs(rootmaV)>=1)
        fprintf('The estimated MA(%d) part of the SARMA model is not reversible.\n',qall);
    end
end
preM = NaN*ones(n+Tmax-1,Tmax); % for simplicity use the indices for the whole
                                % time series, the first n1 will be ignored
preM = NaN*ones(n,Tmax);
xxV = xV-mx1;
for T=1:Tmax
    tmpS = predict(sarmamodel,xxV,T);
    if myversion==0
        preM(:,T) = tmpS+mx1;   
    elseif myversion==1
        preM(:,T) = tmpS.OutputData+mx1;   
    else
        preM(:,T) = tmpS{1}+mx1;
    end
end                                
nrmseV = ones(Tmax,1);
for T=1:Tmax
    nrmseV(T) = nrmse(xV(n1+T:n),preM(n1+T:n,T));
end
preM = preM(n1+1:n,:);

if ~isempty(tittxt)
	figno = gcf;
	figure(figno)
	clf
	plot([1:Tmax]',nrmseV,'.-k')
	hold on
	plot([1 Tmax],[1 1],'y')
	xlabel('prediction time T')
	ylabel('NRMSE(T)')
    if qall==0
    	title(sprintf('%s, NRMSE(T) for SAR(%d)x(%d)_%d  prediction, n=%d, nlast=%d',...
            tittxt,p,ps,s,n,nlast))
    elseif p==0
        title(sprintf('%s, NRMSE(T) for SMA(%d)x(%d)_%d prediction, n=%d, nlast=%d',...
           tittxt,q,qs,s,n,nlast))
    else
        title(sprintf('%s, NRMSE(T) for SARMA(%d,%d)x(%d,%d)_%d prediction, n=%d, nlast=%d',...
            tittxt,p,q,ps,qs,s,n,nlast))
    end
end
