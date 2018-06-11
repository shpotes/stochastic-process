% Finite Difference Method

%% Clears Workspace
clear all
clc
close all

%% Arbitrary Parameters
r = 0.05;
sigma = 0.7;
K = 100;
S_0 = 100;
F = 2;
% dS = 0.05;
T = 1/252;
F_max = 2;

%figure;
FF = [];

Call = blsprice(S_0, K, r, T, sigma);
for TS=1:10
    [F, P, v, t] = dif_fin(r, sigma, K, S_0, F_max, TS, T);
    FF = [FF F];
end
plot(FF*0 + Call, 'r')
hold on
plot(FF, 'ok')

%%
% 
proceso = xlsread('datos_sin_fecha.csv'); 
plot(proceso);
sample_data = xlsread('retornos.csv'); 
plot(sample_data);

% Right
figure
parcorr(proceso)

% Wrong
figure
parcorr(sample_data)
