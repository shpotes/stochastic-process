function F=dif_fin(S_0, K, r, T, sigma, F_max, TS)
 % Input
 % r     : Tasa de interes libre de riesgo  (0.05)
 % sigma : Volatilidad                      (estimarla)
 % K     : Precio de ejercicio              (Cualquiera)
 % S0    : Valor de accion en tiempo actual (Ultimo valor de la serie)
 % FMax  : Smax = S0 * FMAX
 % dS    : 1/Tamaños del inter de muestreo
 % T     : tiempo de muestreo
 %
 % Output
 % F     : 
  
  dS = 1/TS;
  S_max = S_0*F_max;              
  M = round(S_max/dS);
  dT_temp = (dS/(sigma*S_max))^2;
  N_temp = T/dT_temp;
  N = ceil(N_temp);
  dt = T/N;

  J = 2:M-1;
  
  a = @(j) dt/(1 + r*dt).*(sigma^2 .* j.^2 - r .*j)./2;
  b = @(j) dt/(1 + r*dt).*(1/dt - sigma^2 .*j.^2);
  c = @(j) dt/(1 + r*dt).*(sigma^2 .* j.^2 + r .*j)./2;
  
  tmp = J*dS - K;
  v = [0, (tmp).*(tmp > 0), S_max - K];

  for i=1:N
    v = [0, a(J).*v(J-1) + b(J).*v(J) + c(J).*v(J+1), S_max-K];
  end
  P = round(S_0*TS);
  F = v(P);
end
