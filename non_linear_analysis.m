% Project Timeseries 2017-2018
% Team 29
% Diamanti Maria 8133
% Ntzioni Dimitra 8209


%% a

data = VarName1;

%2-d diagram

m = 2;
tau = 100;

[xM2] = embeddelays(data, m, tau);
figure (1)
plotd2d3(xM2,'Data29v1 - 2d');

%3-d diagram

m = 3;

[xM3] = embeddelays(data, m, tau);
figure (2)
plotd2d3(xM3,'Data29v2 - 3d');

%% b

% Find ô from the autocorrelation diagram

n = length(data);
sdnoise = 10;
mux2= 20;
perseason = 12;
maxtau = 100;
alpha = 0.05;

zalpha = norminv(1-alpha/2);

acxM = autocorrelation(data, maxtau);
autlim = zalpha/sqrt(n);
figure(4)
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
title(sprintf('autocorrelation'))


%Diagram of the mutual Information

[mutM] = mutualinformation(data, tau);
[locextM,xV] = extremes(mutM(:,2),0,1,0,0,1);
chosen_t = locextM(1,1);

% Find m

mmax = 10;
[fnnM,mdistV,sddistsV] = falsenearest(data,chosen_t,mmax,10,0,'m');
chosen_m = 5; %% chosen from the diagram of the false nearest neighboors

%% c


% Diagrams to specify the value of m and to find the v(m)
[rcM,cM,rdM,dM,nuM] = correlationdimension(data,chosen_t,mmax);


%% d

% Prediction

chosen_t = 1;  % the chosen value of ô
chosen_m = 5;  % the chosen value of m
n1 = length(data)/2;
Tmax = 10;
q = 0;
start = length(data)/2;
stop =  length(data)/2 + 9;

nrmseValue = [];
for i = 1:20
    [preV] = localpredictmultistep(data,n1,chosen_t,chosen_m,Tmax,i,q);
    nrmseV = nrmse(data(start:stop),preV);
    nrmseValue = [nrmseValue, nrmseV];
end

nrmseValue
index = find(nrmseValue == min(nrmseValue))

[preV] = localpredictmultistep(data,n1,chosen_t,chosen_m,Tmax,index,q,'Prediction');

%

plot(nrmseValue);
