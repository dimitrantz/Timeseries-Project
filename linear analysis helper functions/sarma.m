function sarmamodel = sarma(xV,p,q,ps,qs,s)
% sarmamodel = sarma(xV,p,q,ps,qs,s)
% Estimates a SARMA model on the given time series and gives out an idpoly
% model. It actually prepares the delay polynomial and then calls 'armax'.
% INPUTS:
%  xV      : vector of the scalar time series
%  p       : the order of the AR part of the model.
%  q       : the order of the MA part of the model.
%  ps      : the order of the seasonal AR part of the model.
%  qs      : the order of the seasonal MA part of the model.
%  s       : the period of the seasonal component.    
% OUTPUT:
%  sarmamodel: the idpoly model.

% Prepare the delay polynomials to be set as input for the ARMA model of
% extended AR and MA orders to include the seasonal polynomials.
if s==0
    s=1; % if for any reason s=0 is given, turn it to 1 to avoid stack of the program
end
phiV=ones(p+1,1);       
phisV = zeros(ps*s+1,1);
phisV(1:s:ps*s+1)=ones(ps+1,1);
thetaV=ones(q+1,1);       
thetasV = zeros(qs*s+1,1);
thetasV(1:s:qs*s+1)=ones(qs+1,1);
% The final delay polynomials are the product of time-incremental and
% seasonal delay polynomials 
phiallV=conv(phiV,phisV);
thetaallV=conv(thetaV,thetasV);

% The zeros in the extended delay polynomials indicate that the respective
% coefficients should be fixed (set to 0) and not to be estimated.
settozeroV = [];
indV = find(phiallV==0);
if ~isempty(indV)
    settozeroV = [settozeroV; indV-1]; % The first is always 1 (coefficient of X(t)) 
end
indV = find(thetaallV==0);
if ~isempty(indV)
    settozeroV = [settozeroV; indV-1+length(phiallV)-1]; 
    % The first is always 1 (coefficient of Z(t)), to this add the length of phi-vector 
end

% Pass the delay polynomials in the appropriate format using idpoly to be 
% input in armax. Also update the property for fixed parameters.
phiinputV=zeros(1,length(phiallV));
phiinputV(1)=1;
thetainputV=zeros(1,length(thetaallV));
thetainputV(1)=1;
sarmapol=idpoly(phiinputV,[],thetainputV,[]); 
set(sarmapol,'FixedParameter',settozeroV);   
sarmamodel=armax(xV,sarmapol); 