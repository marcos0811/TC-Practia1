% Se realiza la lectura de los datos
d = readmatrix("data_motor.csv");
t = d(:,2); % Tiempo
u = d(:,3); % Escalón 
y = d(:,4); % Salida del sistema

%% Definición de variables base
y0 = y(2);           % valor inicial 
yf = y(end);         % valor final 
u0 = u(2);           % valor inicial de la entrada
uf = u(end);         % valor final de la entrada
U  = uf - u0;        % amplitud del escalón
Dy = yf - y0;        % variación de la salida

% Ganancia del sistema
k = Dy / U;          % ganancia del proceso

%% Cálculo de pendiente (diferencia finita)
v_prom = [];
for i = 2:length(y)-1
    prom = (y(i+1)-y(i))/(t(i+1)-t(i));
    v_prom = [v_prom, prom];
end

m = max(v_prom); % pendiente máxima (punto de inflexión)
v_max = find(v_prom == m,1); 
t1 = (t(v_max) + t(v_max+1))/2;
y1 = (y(v_max) + y(v_max+1))/2;

% Recta tangente
e = y1 + m*(t - t1);

%% Método de Ziegler-Nichols
t_theta = t1 + (y0 - y1)/m;   % tiempo muerto
t_tau   = t1 + (yf - y1)/m;   % instante cuando alcanza yf
tau     = t_tau - t_theta;    % constante de tiempo aparente

y_model1 = zeros(size(t));
for i = 1:length(t)
    if t(i) < t_theta
        y_model1(i) = y0;
    else
        y_model1(i) = y0 + k*U*(1 - exp(-(t(i)-t_theta)/tau));
    end
end

%% Método de Miller
y_632 = y0 + 0.632*Dy;
[~, idx_63] = min(abs(y - y_632));
t_63 = t(idx_63);

tau2 = t_63 - t_theta;

y_model2 = zeros(size(t));
for i = 1:length(t)
    if t(i) < t_theta
        y_model2(i) = y0;
    else
        y_model2(i) = y0 + k*U*(1 - exp(-(t(i)-t_theta)/tau2));
    end
end

%% Método Analítico
y_284 = y0 + 0.284*Dy;
[~, idx_284] = min(abs(y - y_284));
t_284 = t(idx_284);

tau3 = 3/2 * (t_63 - t_284);
theta3 = t_63 - tau3;

y_model3 = zeros(size(t));
for i = 1:length(t)
    if t(i) < theta3
        y_model3(i) = y0;
    else
        y_model3(i) = y0 + k*U*(1 - exp(-(t(i)-theta3)/tau3));
    end
end

%% Gráficas
figure;
plot(t, y, 'b', 'LineWidth',1.5);           % datos reales
hold on;
plot(t, y_model1, 'r', 'LineWidth',1.5); % modelo zegler nichols
plot(t, y_model2, 'g', 'LineWidth',1.5); % modelo miller
plot(t, y_model3, 'm', 'LineWidth',1.5); % modelo analítico
plot(t, e, 'k:', 'LineWidth',1.5);         % recta tangente
plot(t,u, 'g', 'LineWidth',1.5); % señal escalon


yline(yf, '--', 'LineWidth',2)             % línea 100%
yline(y0, '--','LineWidth',2, 'Color','yellow')    % línea base
plot(t_theta, y0, 'ro', 'MarkerFaceColor','r'); % corte con base
plot(t_tau, yf, 'go', 'MarkerFaceColor','g');   % corte con 100%

grid on;
legend('Respuesta al proceso','Ziegler y Nichols','Miller','Analítico','Tangente','Señal escalon','Linea 100%','Linea Base','\theta','\theta+\tau');
xlabel('Tiempo [s]');
ylabel('Salida');
ylim([-0.1 1.6]);
title('Identificación Grafica');
