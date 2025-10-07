d = readmatrix("data_motor.csv");
t = d(:,2);
u = d(:,3);
y = d(:,4);

%se elimina los valores NaN
valid_idx = ~isnan(t) & ~isnan(u) & ~isnan(y); 
t = t(valid_idx);
u = u(valid_idx);
y = y(valid_idx);

yi = y(1); %Primer valor de la respuesta del proceso
yf = y(end); %Ultimo valor de la respeusta del proceso
ui = u(1); % Primer valor de la se;al escalon
uf = u(end); % Ultimo valor de la se;al escalon
y_inicio = find(y ~= 0, 1);%buscamos el primer valor diferente de cero
t0 = t(y_inicio);  % tiempo donde inicia el experimento


% Método de Ziegler-Nichols
v_prom = diff(y)./diff(t);
m = max(v_prom);
v_max = find(v_prom == m, 1);

t1 = (t(v_max) + t(v_max+1))/2;
y1 = (y(v_max) + y(v_max+1))/2;

Rt = y1 + m.*(t - t1);

k = (yf - yi)/(uf - ui); % Ganancia
theta = t1 + (yi - y1)/m - t0; %tiempo mueto
tau = t1 + (yf - y1)/m; %Tehta + tau
tau1 = tau - theta; %constante de tiempo

G1 = tf(k, [tau1 1], 'InputDelay', theta);
fm = lsim(G1, u,t);
error1 = mean((y-fm).^2); % se calclula el error cuadratico medio

%Metodo dce miller, el valor de k y theta son los mismo que 
k2 = k; 
theta2 = theta;
y2 = yi + 0.632*(yf-yi)-t0; %y2 es el valor donde inica el 63.2%
[a, ind] = min(abs(y - y2)); %indice del valor mas cercano a 63.27%
t2 = t(ind); %t2, valor del tiempo del 63.27% de y
tau2 = t2 - theta2;

G2 = tf(k,[tau2 1],'InputDelay',theta2);
fm2 = lsim(G2 ,u,t);

error2 = mean((y-fm2).^2); %error cuadratico medio del segundo metodo
%metodo analitico
y3 = y2 ; % el 63.27% desde que inicia el expereimento
y3_i = yi + 0.284*(yf-yi) - t0; % el 28.4% desde que incia el experimento

[b, ind2] = min(abs(y - y3_i)); % indice del valor del 28.4%
t3 = t(ind2); % tiempo en el que inica el experimento del 28.4%
tau3 = 3/2 * (t2 - t3);
theta3 = t2 - tau3 ;

G3 = tf(k,[tau3 1],'InputDelay',theta3);
fm3 = lsim(G3 ,u,t);
error3 = mean((y-fm3).^2); %Error cuadratico medio
figure
hold on
plot(t, u, 'Color', [0 0.447 0.741], 'LineWidth', 1.5)
plot(t, fm, 'Color', [0.929 0.694 0.125], 'LineWidth', 1.5)
plot(t, fm2, 'Color', [0.850 0.325 0.098], 'LineWidth', 1.5)
plot(t, fm3, 'Color', [0.301 0.745 0.933], 'LineWidth', 1.5)
plot(t, y, 'Color', [0.466 0.674 0.188], 'LineWidth', 1.5)
plot(t, Rt, '--', 'Color', [0.494 0.184 0.556])
yline(yf, 'g--')
yline(yi, 'r--')
plot(theta + t0, yi, 'ro', 'MarkerFaceColor', 'r')
plot(tau1 + t0, yf, 'wo', 'MarkerFaceColor', 'w')
xlabel('Tiempo (s)')
ylabel('Amplitud')
title('Integración Grafica')
legend('Señal escalón','Metodo de Ziegler y Nichols','Metodo de Miller','Metodo Analitico',...
    'Salida medida','Recta tangente','Línea 100','Línea 0')
ylim([-0.1 1.6])
grid on

