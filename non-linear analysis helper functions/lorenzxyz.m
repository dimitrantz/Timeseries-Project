function xM = lorenzxyz(n,taus,parV,x0V,h)
% xM = lorenzxyz(n,taus,parV,x0V,h)
% Input:
% - n   : trajectory length
% - taus    : (optional) sampling time
% - parV    : (optional) parameters for Lorenz system
% - x0V     : (optional) initial condition for the generated trajectory
% - h       : (optional) discretization step for the numerical integration

x0limV = [-10 10 -20 20 5 30]';
ntrans = 1000;
dim=3; % dummy parameter called in the integration function
if nargin==4
    h=0.01;
elseif nargin==3
    h=0.01;
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) ...
            (x0limV(4)-rand)/(x0limV(4)-x0limV(3)) ...
            (x0limV(6)-rand)/(x0limV(6)-x0limV(5))];
elseif nargin==2
    h=0.01;
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) ...
            (x0limV(4)-rand)/(x0limV(4)-x0limV(3)) ...
            (x0limV(6)-rand)/(x0limV(6)-x0limV(5))];
    parV = [10 28 8/3]; % chaos
    % [ -8/3 0 0; 0 -10 10; 0 28 -1 ]; % chaos
    % parV = [10 13.962 8/3]; % stable point
    % parV = [10 20 8/3]; % stable point
elseif nargin==1
    h=0.01;
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) ...
            (x0limV(4)-rand)/(x0limV(4)-x0limV(3)) ...
            (x0limV(6)-rand)/(x0limV(6)-x0limV(5))];
    parV = [10 28 8/3]; % chaos
    taus=h;
end
if isempty(parV)
    parV = [10 28 8/3]; % chaos
end    
if isempty(x0V)
    x0V = [(x0limV(2)-rand)/(x0limV(2)-x0limV(1)) ...
            (x0limV(4)-rand)/(x0limV(4)-x0limV(3)) ...
            (x0limV(6)-rand)/(x0limV(6)-x0limV(5))];
end    
% Integrate the equations over length steps to remove any transient behavior
for i=1:ntrans
    xoutV=step_it('dLorenz63',x0V,parV,h,dim);
    x0V=xoutV;
end
% Integrate the equations over length steps and record all states
resample = taus / h;
if resample - round(resample) ~= 0
    error('The sampling time should be multiple of the integration step!')
end
xM = NaN*ones(n,3);
if resample>1
    for i=1:n
        for j=1:resample
            xoutV=step_it('dLorenz63',x0V,parV,h,dim);
            x0V=xoutV;
        end
        xM(i,:)=xoutV;
    end
else
    for i=1:n
        xoutV=step_it('dLorenz63',x0V,parV,h,dim);
        xM(i,:)=xoutV;
        x0V=xoutV;
    end
end
