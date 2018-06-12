function F=dif_fin(S_0, K, r, T, sigma, F_max, TS)
 % Input
 % r     : Tasa de interes libre de riesgo  (0.05)
 % sigma : Volatilidad                      (estimarla)
 % K     : Precio de ejercicio              (Cualquiera)
 % S0    : Valor de accion en tiempo actual (Ultimo valor de la serie)
 % FMax  : Smax = S0 * FMAX
 % TS    : 1/Tamaños del inter de muestreo
 % T     : Termino de venicimineto
 %
 % Output
 % F     : Solucion de la ecuacion de BS

  dS = 1/TS; % Tamaño del intervalo
  S_max = S_0*F_max; % Precio de activo lo suficientemente grande
  M = round(S_max/dS);
  dT_temp = (dS/(sigma*S_max))^2; % Se escoge un dt correspondiente
  N_temp = T/dT_temp;
  N = ceil(N_temp);
  dt = T/N; % Se redonde el valor correspondiente de dt

  % En lugar de generar la matriz de (M+1)(N+1), se trabajará con un vector de M+1
  % El vector J sirve para facilitar la vectorizacion de operaciones
  J = 2:M; 
  
  a = @(j) dt/(1 + r*dt).*(sigma^2 .* j.^2 - r .*j)./2;
  b = @(j) dt/(1 + r*dt).*(1/dt - sigma^2 .*j.^2);
  c = @(j) dt/(1 + r*dt).*(sigma^2 .* j.^2 + r .*j)./2;

  tmp = J*dS - K; % St - K
  tmp = (tmp).*(tmp > 0) % max (St - k, 0)
  v = [0, tmp, S_max - K]; % Condiciones iniciales

  for i=1:(N+1)
    % Ejecucion del metodo, respetando las condiciones de frontera
    v = [0, a(J).*v(J-1) + b(J).*v(J) + c(J).*v(J+1), S_max-K];
  end

  P = round(S_0*TS);
  F = v(P);
end
