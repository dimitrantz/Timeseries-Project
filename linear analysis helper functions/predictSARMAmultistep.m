function [preV,phiV,thetaV] = predictSARMAmultistep(xV,n1,p,q,ps,qs,s,Tmax,tittxt)
% [preV,phiV,thetaV] = predictSARMAmultistep(xV,n1,p,q,ps,qs,s,Tmax,tittxt)
% PREDICTSARMAMULTISTEP makes multi-step ahead predictions using the 
% SARMA(p,q)x(ps,qs)_s model
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
%  n1      : size of training set, number of samples of the segment of xV 
%            (x(1),x(2),...,x(n1)) to which the SARMA model is fitted.
%            If n1 is empty, then n1 = length(xV).
%  p       : the order of the AR part of the model.
%  q       : the order of the MA part of the model
%  ps      : the order of the seasonal AR part of the model.
%  qs      : the order of the seasonal MA part of the model.
%  s       : the period of the seasonal component.    
%  Tmax    : the prediction horizon, predictions are made for T=1...Tmax
%            steps ahead, i.e. for times n1+1,n1+2,...,n1+Tmax
%  tittxt  : string to be displayed in the title of the figure (if not
%            given no figure is displayed) 
% OUTPUT: 
%  preV    : vector of length Tmax of the predicted values 
%            x(n1+1),x(n1+2),...,x(n1+Tmax)
%  phiV    : the coefficients of the estimated AR part (of length
%            (ps*s+p+1) with phi(0) as first component. Note that these are
%            the coefficients of both AR and seasonal AR.
%  thetaV  : the coefficients of the estimated MA part (of length qs*s+q). 
%            Note that these are the coefficients of both MA and seasonal MA.
% The actual and predicted values are plotted (provided that the length of
% 'xV' is at least 'n1+Tmax').

if nargin==8
    tittxt = [];
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
n = length(xV);
if isempty(n1)
    n1=n;
end
if n1>n
    error('The given size of training set exceeds the length of time series.');
end
xV = xV(:);
mx = mean(xV);
xxV = xV-mx;

sarmamodel = sarma(xxV(1:n1),p,q,ps,qs,s);
pall = p+ps*s;
qall = q+qs*s;
if pall==0
    phiV = [];
else
    phi0 = (1+sum(sarmamodel.a(2:1+pall)))*mx;
    phiV = [phi0; -sarmamodel.a(2:1+pall)'];
    rootarV = roots(sarmamodel.a);
    if any(abs(rootarV)>=1)
        fprintf('The estimated AR(%d) part of the SARMA model is not stationary.\n',pall);
    end
end
if qall==0
    thetaV = [];
else
    thetaV = -sarmamodel.c(2:1+qall)';
    rootmaV = roots(sarmamodel.c);
    if any(abs(rootmaV)>=1)
        fprintf('The estimated MA(%d) part of the SARMA model is not reversible.\n',qall);
    end
end
pq = max(pall,qall);
preV = NaN*ones(pq+Tmax,1);
if qall>0
    tmpS = predict(sarmamodel,xxV(1:n1),1);
    if myversion==0
        xpreV = tmpS+mx;   
    elseif myversion==1
        xpreV = tmpS.OutputData+mx;   
    else
        xpreV = tmpS{1}+mx;
    end
    zV = xV(1:n1) - xpreV;
    zpreV = zeros(pq+Tmax,1);
    zpreV(pq-qall+1:pq)=zV(n1-qall+1:n1);
    % Now we can make predictions based also on the MA fit residuals, if
    % necessary
end
if p>0
    preV(pq-pall+1:pq)=xxV(n1-pall+1:n1);
end
if qall==0
    for T=pq+1:pq+Tmax
        preV(T)=mx + phiV(2:pall+1)'*preV(T-1:-1:T-pall);
    end
elseif pall==0
    for T=pq+1:pq+Tmax
        preV(T)=mx -thetaV'*zpreV(T-1:-1:T-qall);
    end
else
    for T=pq+1:pq+Tmax
        preV(T)=mx + phiV(2:pall+1)'*preV(T-1:-1:T-pall)-thetaV'*zpreV(T-1:-1:T-qall);
    end
end
preV = preV(pq+1:pq+Tmax);
if ~isempty(tittxt)
    iV = [n1+1:n1+Tmax]';
    if length(xV)<n1+Tmax
        i2V = [n1+1:length(xV)]';
        oriV = xV(i2V);
    elseif length(xV)==n1
        i2V = [];
        oriV = [];
    else
        i2V = iV;
        oriV = xV(i2V);
    end
    figure(gcf)
    clf
    plot(iV,preV,'.-r')
    hold on
    xlabel('T')
    ylabel('x(t+T)')
    if qall==0
        title(sprintf('%s, multi-step SAR(%d)x(%d)_%d prediction',tittxt,p,ps,s))
    elseif pall==0
        title(sprintf('%s, multi-step SMA(%d)x(%d)_%d prediction',tittxt,q,qs,s))
    else
        title(sprintf('%s, multi-step SARMA(%d,%d)x(%d,%d)_%d prediction',tittxt,p,q,ps,qs,s))
    end
    if isempty(oriV)
        legend('predicted',0)
    else
        plot(i2V,oriV,'.-k')
        legend('real','predicted',0)
    end
end
