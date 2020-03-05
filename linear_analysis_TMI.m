% Project Timeseries 2017-2018
% Team 29
% Diamanti Maria 8133
% Ntzioni Dimitra 8209

%% Timeseries
[matrix, filtered] = extremes(VarName2, 2);

TMI = [];
TMI(1) = abs(matrix(1,1) - matrix(2,1));
for i = 2:(length(matrix)/2)
    TMI(i) = abs(matrix(2*i-1,1) - matrix(2*i,1));
end

figure(5)
plot(TMI)
title('TMI')

%% plot autocorrelation TMI

n = length(TMI);
sdnoise = 10;
mux2= 20;
perseason = 12;
maxtau = 100;
alpha = 0.05;

zalpha = norminv(1-alpha/2);

acxM = autocorrelation(TMI, maxtau);
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
title(sprintf('autocorrelation TMI'))

%% Remove trend using first differences

x1 = ones(length(TMI),1);
for i=2:1:length(TMI) % first differences
    x1(i) = TMI(i) - TMI(i-1);
end
x1(1) = TMI(1);

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
TMI = TMI';
parRX2 = parautocor(TMI,maxtau);
    
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
title('Partial autocorrelation, TMI')

%% Define the optimal model from AIC index

figure(13)
AIC = []; %check AR order
for j = 1:10
   for i = 1:10
        [nrmseV,phiV,thetaV,SDz,aicS,fpeS,armamodel] = fitARMA(TMI,i,j,5); %%fit with ARMA(p,q) AIC is minimum
        AIC(i,j) = aicS;
   end
   hold on
   plot(AIC(:,j))
end

[p,q] = find(AIC == min(min(AIC)))


%% Prediction for 2 steps ahead

[nrmseV, preM, phiV, thetaV] = predictARMAnrmse(TMI,p,q,2);
nrmseV

%% NRMSE

premTMI1 = preM(:,1)';
premTMI2 = preM(:,2)';

start = floor(length(TMI)/2)+1;


TMI = TMI';
nrmseTMI1 = nrmse(TMI(start:end),premTMI1);
nrmseTMI1
nrmseTMI2 = nrmse(TMI(start:end),premTMI2);
nrmseTMI2

figure(30)
hold on
plot(TMI(start:end))
plot(premTMI1)
hold off
title('TMI - Prediction one step ahead')

figure(31)
hold on
plot(TMI(start:end))
plot(premTMI2)
hold off
title('TMI - Prediction two steps ahead')
