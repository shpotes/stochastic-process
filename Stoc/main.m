% Finite Difference Method

%% Clears Workspace
clear all
clc
close all

%% Read data
data = csvread('../data/matlab.csv');
stock = data(:,1);
ret = data(:, 2);


%% Nosotros Arriba
r = 0.05
dT = 1/252
sigma = sqrt(var(ret)/dT)
S_0 = stock(end)
K = S_0
TS = 3;
F_max = 2
T = 1/4

%% Juan
%S_0 = 9.71;
%sigma = 0.4473;
%K = S_0;
%TS = 100;

% Sergio Abajo
%S_0 = 284.73;
%K = S_0;
%sigma = 0.4628;
%S_0 = round(S_0);

FF = [];
tic
for TS=1:10
  TS
  call = blsprice(S_0, K, r, T, sigma);
  F = dif_fin(S_0, K, r, T, sigma, F_max, TS);
  error = (F - call)/call;
  FF = [[TS, F, call, error]; FF];
  %A = [T, call, F, P, v, t];
  %AA = [AA; A];
end
toc

FF
'Nosotros'
csvwrite('../data/error.csv', FF)
