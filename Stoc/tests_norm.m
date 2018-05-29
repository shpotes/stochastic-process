% Tests de Normalidad
%
clc 
clear all
close all

% Parameters
ALPHA = 0.05;                           % Significance Level

% Data
sample_data = xlsread('retornos.csv');  % Instantaneous Returns

% Shapiro Wilks
% H = 0 => Do not reject the null hypothesis at significance level ALPHA.
% H = 1 => Reject the null hypothesis at significance level ALPHA.
% pValue : p_Value > ALPHA
% SWstatistic : Test statistic 
[H_sw, pValue, SWstatistic] = swtest(sample_data, ALPHA);    % Not rejected

% Jarque-Bera
H_jb = jbtest(sample_data,ALPHA)                             % Not rejected

% Kolmogorov - Smirnov
% Lilliefors

% Histogram
histogram(sample_data)
title('Histograma de Retornos Instantaneos');
xlabel('Conteo');
ylabel('Retornos Instantaneos');

