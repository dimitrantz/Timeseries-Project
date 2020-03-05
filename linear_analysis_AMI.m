% Project Timeseries 2017-2018
% Team 29
% Diamanti Maria 8133
% Ntzioni Dimitra 8209

%% Timeseries
[matrix, filtered] = extremes(VarName1, 2);

AMA = [];
for i = 1:length(matrix)
    if matrix(i,3) == 1;
        AMA = [AMA,matrix(i,2)];
    end
end

figure(1)
plot(AMA)
title('AMA')

AMI = [];
for i = 1:length(matrix)
    if matrix(i,3) == -1;
        AMI = [AMI,matrix(i,2)];
    end
end

figure(2)
plot(AMI)
title('AMI')

AMD = [];
for i = 1:length(AMI)
    AMD(i) = abs(AMI(i) - AMA(i+1));
   
end

figure(3)
plot(AMD)
title('AMD')

TMA = [];
for i = 1:(length(matrix)/2)
    TMA(i) = abs(matrix(2*i,1) - matrix(2*i+1,1));
end

figure(4)
plot(TMA)
title('TMA')

TMI = [];
TMI(1) = abs(matrix(1,1) - matrix(2,1));
for i = 2:(length(matrix)/2)
    TMI(i) = abs(matrix(2*i-1,1) - matrix(2*i,1));
end

figure(5)
plot(TMI)
title('TMI')

TBP = [];
for i = 1:(length(matrix)/2)
    TBP(i) = abs(matrix(2*i-1,1) - matrix(2*i+1,1));
end

figure(6)
plot(TBP)
title('TBP')

%% plot autocorrelation AMI

n = length(AMI);
sdnoise = 10;
mux2= 20;
perseason = 12;
maxtau = 100;
alpha = 0.05;

zalpha = norminv(1-alpha/2);

acxM = autocorrelation(AMI, maxtau);
autlim = zalpha/sqrt(n);
figure(9)
clf
hold on
for ii=1:maxtau
    plot(acxM(ii+1,1)*[1 1],[0 acxM(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('autocorrelation AMI'))

%% Remove trend using first differences

x1 = ones(length(AMI),1);
for i=2:1:length(AMI) % first differences
    x1(i) = AMI(i) - AMI(i-1);
end
x1(1) = AMI(1);

%% Autocorrelation of the detrended time series by first differences

rXt1 = autocorrelation(x1,maxtau);
autlim = zalpha/sqrt(n);
figure(10)
clf
hold on
for ii=1:maxtau
    plot(rXt1(ii+1,1)*[1 1],[0 rXt1(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title(sprintf('detrended time series by first differences, autocorrelation'))


%% Partial autocorrelation second differences

parRX2 = parautocor(x1,maxtau);
    
figure(12)
clf
hold on
for i=1:(maxtau-1)
  plot(i*[1 1],[0 parRX2(i)],'b','linewidth',2)
end
plot([0 maxtau+1],[0 0],'k')
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau\tau}')
title('Partial autocorrelation, first differences')

%% Define the optimal model from AIC index

figure(13)
AIC = []; %check AR order
for j = 1:10
   for i = 1:10
        [nrmseV,phiV,thetaV,SDz,aicS,fpeS,armamodel] = fitARMA(x1,i,j,5); %%fit with ARMA(p,q) AIC is minimum
        AIC(i,j) = aicS;
   end
   hold on
   plot(AIC(:,j))
end

[p,q] = find(AIC == min(min(AIC)))

% chosen: ARMA(8,7) for data29v1
% chosen: ARMA(2,5) for data29v2

%% Prediction for 2 steps ahead

[nrmseV, preM, phiV, thetaV] = predictARMAnrmse(x1,p,q,2);
nrmseV

%% NRMSE


premAMI = ones(length(preM),2);
for i=2:1:length(preM) % first differences
    premAMI(i,:) = preM(i,:) + preM(i-1,:);
end
premAMI(1,:) = preM(1,:);
premAMI1 = premAMI(:,1)';
premAMI2 = premAMI(:,2)';

start = floor(length(AMI)/2)+1;


nrmseAMI1 = nrmse(AMI(start:end),premAMI1);
nrmseAMI1
nrmseAMI2 = nrmse(AMI(start:end),premAMI2);
nrmseAMI2

figure(30)
hold on
plot(AMI(start:end))
plot(premAMI1)
hold off
title('AMI - Prediction one step ahead')

figure(31)
hold on
plot(AMI(start:end))
plot(premAMI2)
hold off
title('AMI - Prediction two steps ahead')
