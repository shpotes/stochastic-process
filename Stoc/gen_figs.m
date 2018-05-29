clc
clear all
close all

A=[];
indices = 0 : 4;
for iter = indices
r = 0.05;                       % Tasa de interes libre de riesgo
sigma = 0.2;                    % Volatilidad
K = 100;                        % Precio de ejercicio
S_0 = 123.456;                  % Valor de accion en tiempo actual
S_max = ceil(2*S_0);            % Cota superior de la accion
dS = 2^(-iter);                 % Tamaño de intervalos del precio  
M = S_max/dS;                   % Intervalos del precio 
T = 1/12;                       % Tiempo de muestreo
dT_temp = (dS/(sigma*S_max))^2;   %   
N_temp = T/dT_temp;             %    
N = ceil(N_temp);               % Intervalos del tiempo
dT = T/N;                       % Tamaño de intervalos del tiempo

% Black Scholes
d_1 = (log(S_0/K) + (r+(1/2)*sigma^2)*T)/(sigma* sqrt(T));
d_2 = d_1 - sigma*sqrt(T);

f_call = S_0 * normcdf(d_1) - K*exp(-r*T)*normcdf(d_2);
%[Call,Put] = blsprice(S_0, K, r, T, sigma);

% Ponderación
a = zeros(M-1, 1);
b = zeros(M-1, 1);
c = zeros(M-1, 1);

for j = 1 : M-1
    % Ponderados
    a(j) = (dT/(1+r*dT)) * ((sigma^2*j^2)/2 - (r*j/2));
    b(j) = (dT/(1+r*dT)) * ((1/dT) - sigma^2*j^2);
    c(j) = (dT/(1+r*dT)) * ((sigma^2*j^2)/2 + (r*j/2));
    sum_j(j) = a(j) + b(j) + c(j);
end
% % 
% plot(sum_j);
% title('a_j + b_j + c_j');
% %
% plot(c)
% title('c_j');

% Vector de incrementos de s
vect_S = 0 : dS : S_max;

% Crear Malla F
F = zeros(N+1, M+1);

% Condiciones de frontera
F(:,1) = max(0-K, 0);                   % Cond 1

for j = 1 : M+1
    F(N+1, j) = max((j-1)*dS - K, 0);   % Cond 2
end

F(:, M+1) = max(S_max - K, 0);          % Cond 3

% Llenar Malla 
for i = N : -1 : 1
    for j = 2 : M
        F(i,j) = a(j-1) * F(i+1, j-1) + b(j-1) * F(i+1, j) + c(j-1) * F(i+1, j+1);
    end
end

% numerical vs exact
P = round(S_0/dS) + 1;
f_num = F(1,P);
A = [A f_num];

% % Proceso - datos
% figure
% datos = xlsread('datos_sin_fecha.csv'); 
% plot(datos);
% title('Serie de Tiempo');
% xlabel('Tiempo');
% ylabel('Precio del Activo');
% 
% % Retornos
% figure
% retornos = xlsread('retornos.csv'); 
% plot(retornos);
% title('Retornos Instantaneos');
% xlabel('Tiempo');
% ylabel('Retorno');
% 
% % PACF
% figure;
% parcorr(datos);
% 
% figure;
% parcorr(retornos);
end
figure
hold on

plot(indices, A, 'ok')
plot([indices(1), indices(end)], [f_call, f_call], 'r')

title('Soluciones númericas vs Black Scholes');
xlabel('dS = 2^{(-num)}');
ylabel('Precio de la opción');



