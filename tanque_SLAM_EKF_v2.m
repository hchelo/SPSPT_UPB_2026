%% tanque_sim_EKF_SLAM_rayos.m
% EKF-SLAM con el robot oruga Hiwonder siguiendo la ruta por WAYPOINTS
% (cuadrado + diagonal de regreso) del script de Python, escalada al
% tamaño del robot.
%   - Predicción: odometría de encoders (con deslizamiento y ruido)
%   - Sensor    : LiDAR 2D que solo impacta beacons (sin paredes)
%   - Mapa      : 40 beacons aleatorios, DESCONOCIDOS para el filtro
%   - Velocidad constante; el control solo ajusta el giro hacia el waypoint
%   - La diagonal final vuelve a la zona inicial -> cierre de lazo
%
% Estado:  X = [x; y; phi; B_a_x; B_a_y; B_b_x; B_b_y; ...]
% Asociación de datos conocida (cada beacon tiene ID), como en Python.
%
% Se dibujan los rayos del LiDAR: rojos los que impactan un beacon,
% claros los que no encuentran nada dentro del alcance. Las líneas
% verde-azuladas punteadas marcan los beacons que usa el EKF.

clear; close all; clc;

%% === ESCALA RESPECTO AL SCRIPT DE PYTHON ===
% Python usa ±75 m y alcance 30 m. Con escala 1/15 queda en ±5 m y el
% LiDAR en 2 m, proporciones idénticas pero a tamaño del tanque.
escala = 1/15;

%% === VEHÍCULO (Hiwonder) ===
B_ext = 0.19343;            % ancho exterior (de oruga a oruga) [m]
B_int = 0.102;              % ancho interior [m]
B = (B_ext + B_int)/2;      % separación efectiva orugas [m]
L = 0.306;                  % largo chasis [m]

%% === SIMULACIÓN ===
dt    = 0.02;
Tsim  = 600;                % tiempo máximo; se corta al terminar la ruta
t     = 0:dt:Tsim;
N     = length(t);
V_cte = 0.35;               % velocidad constante [m/s]
w_max = 1.5;                % [rad/s]
k_w   = 2.5;                % ganancia de orientación (como en Python)
plot_every = 10;            % refresco de la animación (pasos)

%% === ENCODERS Y DESLIZAMIENTO ===
sigma_encoder = 0.02;
io_mean = 0.08; ii_mean = 0.10; sigma_slip = 0.05;

%% === WAYPOINTS (mismos que Python, escalados) ===
WP = escala * [  0    0;
                75    0;
                75   75;
               -75   75;
               -75  -75;
               -75  -75;
                75  -75;
               -70   70];
umbral_wp = 5*escala;       % distancia para dar por alcanzado un waypoint

%% === BEACONS REALES (aleatorios, como Python) ===
rng(7);
n_beacons     = 40;
beacons       = escala * (-95 + 190*rand(n_beacons,2));
beacon_radius = 0.05;
rng(1);                     % semilla para los ruidos de la simulación

%% === LIDAR ===
lidar_rate_hz   = 10;
lidar_every     = max(1, round(1/(lidar_rate_hz*dt)));
lidar_max_range = 30*escala;           % 2 m
lidar_n_rays    = 180;
lidar_angles    = linspace(-pi, pi, lidar_n_rays+1); lidar_angles(end) = [];
sigma_lidar     = 0.02;
rayos_libres_cada = 2;          % dibujar 1 de cada 2 rayos sin impacto (1 = todos)

sigma_b_range   = 0.03;
sigma_b_bearing = deg2rad(1.0);
R_beacon  = diag([sigma_b_range^2, sigma_b_bearing^2]);
gate_chi2 = 9.21;

%% === EKF-SLAM: ESTADO INICIAL ===
X = [0; 0; 0];
P = diag([0.01, 0.01, deg2rad(1)]).^2;
Q_enc = diag([0.02, 0.02, deg2rad(0.3)]).^2;
Q_u   = diag([sigma_encoder^2, (sigma_encoder/(B/2))^2]);

lm_idx   = zeros(n_beacons,1);  % fila de B_j en X (0 = aún no en el mapa)
lm_order = [];                  % IDs en el orden en que entraron al estado

%% === REGISTROS ===
x_real = zeros(1,N); y_real = zeros(1,N); phi_real = zeros(1,N);
odom_x = zeros(1,N); odom_y = zeros(1,N); odom_phi = zeros(1,N);
X_hist = zeros(3,N);
err_mec = zeros(1,N); err_fus = zeros(1,N); sig3 = zeros(1,N);
n_obs_log = zeros(1,N); n_map_log = zeros(1,N); err_map_log = nan(1,N);
cnt_gate = 0;

scan_x = nan; scan_y = nan; last_obs = [];
wp_idx = 1; kf = N;
lim = 110*escala;

%% === FIGURA 1: ANIMACIÓN ===
fig1 = figure(1); set(fig1,'Position',[40,40,1100,850]); hold on; grid on; axis equal;
plot(WP(:,1), WP(:,2), '--', 'Color',[0.5 0 0.5], 'HandleVisibility','off');
plot(WP(:,1), WP(:,2), 's', 'Color',[0.5 0 0.5], 'MarkerSize',9, ...
     'MarkerFaceColor',[0.5 0 0.5], 'DisplayName','Waypoints');
for i = 1:size(WP,1)
    text(WP(i,1)+0.1, WP(i,2)+0.15, sprintf('WP%d',i-1), 'FontSize',7, 'Color',[0.5 0 0.5]);
end
plot(beacons(:,1), beacons(:,2), 'k+', 'MarkerSize',8, 'LineWidth',1.2, 'DisplayName','Beacons reales');
xlim([-lim lim]); ylim([-lim lim]);
xlabel('X [m]'); ylabel('Y [m]');
ax_map = gca;

hRayMiss = plot(nan,nan,'-', 'Color',[1 0.87 0.78], 'LineWidth',0.5, 'DisplayName','Rayo sin impacto');
hRayHit  = plot(nan,nan,'-', 'Color',[1 0.35 0.35], 'LineWidth',0.8, 'DisplayName','Rayo con impacto');
hScan  = plot(nan,nan,'.', 'Color',[0.85 0 0], 'MarkerSize',9, 'DisplayName','Impactos LiDAR');
hRange = plot(nan,nan,':', 'Color',[1 0.5 0], 'HandleVisibility','off');
hRays  = plot(nan,nan,':', 'Color',[0 0.6 0.6], 'LineWidth',1.2, 'DisplayName','Beacons usados por el EKF');
hReal  = plot(nan,nan,'k-', 'LineWidth',1.6, 'DisplayName','Real');
hMec   = plot(nan,nan,'b--','LineWidth',1.2, 'DisplayName','Odom. mecánica');
hEst   = plot(nan,nan,'g-', 'LineWidth',1.5, 'DisplayName','EKF-SLAM');
hCovR  = plot(nan,nan,'Color',[0 0.5 0], 'DisplayName','3\sigma robot');
hLmPt  = plot(nan,nan,'r.', 'MarkerSize',12, 'DisplayName','Beacons estimados');
hLmEl  = gobjects(n_beacons,1);
for j = 1:n_beacons
    hLmEl(j) = plot(nan,nan,'r-','LineWidth',0.9,'HandleVisibility','off');
end
plot(nan,nan,'r-','DisplayName','3\sigma beacon');
hTank = crear_tanque(ax_map);          % tanque (orugas + chasis + frente)
legend('Location','northeastoutside');

% Vista ampliada que sigue al tanque
ax_zoom = axes('Parent',fig1,'Position',[0.76 0.08 0.21 0.32]);
hold(ax_zoom,'on'); box(ax_zoom,'on'); grid(ax_zoom,'on'); axis(ax_zoom,'equal');
plot(ax_zoom, beacons(:,1), beacons(:,2), 'k+', 'MarkerSize',8, 'LineWidth',1.2);
hZReal = plot(ax_zoom, nan,nan,'k-','LineWidth',1.2);
hZEst  = plot(ax_zoom, nan,nan,'g-','LineWidth',1.2);
hZMiss = plot(ax_zoom, nan,nan,'-','Color',[1 0.87 0.78],'LineWidth',0.5);
hZHit  = plot(ax_zoom, nan,nan,'-','Color',[1 0.35 0.35],'LineWidth',0.8);
hZScan = plot(ax_zoom, nan,nan,'.','Color',[0.85 0 0],'MarkerSize',12);
hZLm   = plot(ax_zoom, nan,nan,'r.','MarkerSize',14);
hTankZ = crear_tanque(ax_zoom);
title(ax_zoom, 'Vista del tanque', 'FontSize',8);
zoom_w = 0.9;                          % semiancho de la vista ampliada [m]

%% === FIGURA 2: MÉTRICAS ===
fig2 = figure(2); set(fig2,'Position',[252,200,1032,600]);
ax1 = subplot(3,1,1); hold on; grid on; box on;
hErrM = plot(nan,nan,'b--','LineWidth',1.4,'DisplayName','Error mecánica');
hErrF = plot(nan,nan,'g-', 'LineWidth',1.6,'DisplayName','Error EKF-SLAM');
hSig  = plot(nan,nan,'k:', 'LineWidth',1.2,'DisplayName','Cota 3\sigma');
ylabel('Error robot [m]'); legend('Location','northwest');
ax2 = subplot(3,1,2); hold on; grid on; box on;
hNobs = stairs(nan,nan,'Color',[1 0.5 0],'LineWidth',1.4,'DisplayName','Observados ahora');
hNmap = stairs(nan,nan,'r-','LineWidth',1.4,'DisplayName','En el mapa');
ylabel('# beacons'); ylim([0 n_beacons+1]); legend('Location','northwest');
ax3 = subplot(3,1,3); hold on; grid on; box on;
hEmap = plot(nan,nan,'r-','LineWidth',1.4);
ylabel('Error medio mapa [m]'); xlabel('Tiempo [s]');

%% === SIMULACIÓN ===
for k = 1:N-1
    %% --- Control por waypoints (velocidad constante, pose real como en Python) ---
    [V_cmd, omega, wp_idx, fin_ruta] = control_waypoints(x_real(k), y_real(k), phi_real(k), ...
                                         WP, wp_idx, V_cte, umbral_wp, k_w, w_max);
    if fin_ruta, kf = k; break; end

    %% --- Deslizamiento + ground truth ---
    if t(k) >= 15 && t(k) <= 17
        io = 0.15+0.05*randn; ii = 0.18+0.05*randn;
    else
        io = io_mean+sigma_slip*randn; ii = ii_mean+sigma_slip*randn;
    end
    Vo = (V_cmd + (B/2)*omega)*(1-io);
    Vi = (V_cmd - (B/2)*omega)*(1-ii);
    x_real(k+1)   = x_real(k) + dt*(Vo+Vi)/2*cos(phi_real(k));
    y_real(k+1)   = y_real(k) + dt*(Vo+Vi)/2*sin(phi_real(k));
    phi_real(k+1) = wrapToPi(phi_real(k) + dt*(Vo-Vi)/B);

    %% --- Encoders y odometría mecánica ---
    Vom = Vo + sigma_encoder*randn; Vim = Vi + sigma_encoder*randn;
    Venc = (Vom+Vim)/2; Wenc = (Vom-Vim)/B;
    odom_phi(k+1) = wrapToPi(odom_phi(k) + dt*Wenc);
    odom_x(k+1)   = odom_x(k) + dt*Venc*cos(odom_phi(k));
    odom_y(k+1)   = odom_y(k) + dt*Venc*sin(odom_phi(k));

    %% --- EKF-SLAM: PREDICCIÓN (solo se mueve el robot) ---
    ph = X(3);
    F  = [1 0 -dt*Venc*sin(ph); 0 1 dt*Venc*cos(ph); 0 0 1];
    Wu = [dt*cos(ph) 0; dt*sin(ph) 0; 0 dt];
    X(1:3) = [X(1)+dt*Venc*cos(ph); X(2)+dt*Venc*sin(ph); wrapToPi(ph+dt*Wenc)];
    P(1:3,1:3) = F*P(1:3,1:3)*F' + Wu*Q_u*Wu' + Q_enc;
    if numel(X) > 3
        P(1:3,4:end) = F*P(1:3,4:end);      % correlación robot-mapa
        P(4:end,1:3) = P(1:3,4:end)';
    end

    %% --- LiDAR + EKF-SLAM: CORRECCIÓN / AUMENTO DEL ESTADO ---
    if mod(k, lidar_every) == 0
        xr = x_real(k+1); yr = y_real(k+1); pr = phi_real(k+1);

        rng_scan = lidar_scan(xr, yr, pr, lidar_angles, lidar_max_range, beacons, beacon_radius);
        hit = rng_scan < lidar_max_range;
        rng_scan(hit) = rng_scan(hit) + sigma_lidar*randn(1, nnz(hit));
        scan_x = xr + rng_scan(hit).*cos(pr + lidar_angles(hit));
        scan_y = yr + rng_scan(hit).*sin(pr + lidar_angles(hit));

        last_obs = [];
        for j = 1:n_beacons
            dxr = beacons(j,1)-xr; dyr = beacons(j,2)-yr;
            rr = hypot(dxr,dyr);
            if rr > lidar_max_range, continue; end
            z = [rr + sigma_b_range*randn;
                 wrapToPi(atan2(dyr,dxr) - pr + sigma_b_bearing*randn)];
            last_obs(end+1) = j; %#ok<AGROW>

            if lm_idx(j) == 0
                % ===== Beacon nuevo: aumentar estado y covarianza =====
                a  = z(2) + X(3);
                lx = X(1) + z(1)*cos(a);
                ly = X(2) + z(1)*sin(a);
                Gr = [1 0 -z(1)*sin(a); 0 1 z(1)*cos(a)];       % d(lm)/d(robot)
                Gz = [cos(a) -z(1)*sin(a); sin(a) z(1)*cos(a)];  % d(lm)/d(z)
                n  = numel(X);
                Pll = Gr*P(1:3,1:3)*Gr' + Gz*R_beacon*Gz';
                Plx = Gr*P(1:3,:);
                X = [X; lx; ly]; %#ok<AGROW>
                P = [P, Plx'; Plx, Pll];
                lm_idx(j) = n+1;
                lm_order(end+1) = j; %#ok<AGROW>
            else
                % ===== Beacon conocido: actualización EKF =====
                i  = lm_idx(j);
                dx = X(i)-X(1); dy = X(i+1)-X(2);
                q  = dx^2+dy^2; r = sqrt(q);
                z_hat = [r; wrapToPi(atan2(dy,dx) - X(3))];
                n = numel(X);
                H = zeros(2,n);
                H(:,1:3)   = [-dx/r, -dy/r,  0;  dy/q, -dx/q, -1];
                H(:,i:i+1) = [ dx/r,  dy/r;     -dy/q,  dx/q];
                nu = z - z_hat; nu(2) = wrapToPi(nu(2));
                S  = H*P*H' + R_beacon;
                if nu'/S*nu > gate_chi2
                    cnt_gate = cnt_gate + 1; continue;
                end
                K = P*H'/S;
                X = X + K*nu; X(3) = wrapToPi(X(3));
                IKH = eye(n) - K*H;
                P = IKH*P*IKH' + K*R_beacon*K';
                P = (P+P')/2;
            end
        end
    end

    %% --- Registros ---
    X_hist(:,k+1) = X(1:3);
    err_mec(k+1)  = hypot(odom_x(k+1)-x_real(k+1), odom_y(k+1)-y_real(k+1));
    err_fus(k+1)  = hypot(X(1)-x_real(k+1), X(2)-y_real(k+1));
    sig3(k+1)     = 3*sqrt(max(eig(P(1:2,1:2))));
    n_obs_log(k+1) = numel(last_obs);
    n_map_log(k+1) = numel(lm_order);
    if ~isempty(lm_order)
        ii_ = lm_idx(lm_order);
        err_map_log(k+1) = mean(hypot(X(ii_)-beacons(lm_order,1), X(ii_+1)-beacons(lm_order,2)));
    end
    kf = k+1;

    %% --- Gráficos ---
    if mod(k,plot_every) == 0
        kk = k+1;
        set(hReal,'XData',x_real(1:kk),'YData',y_real(1:kk));
        set(hMec, 'XData',odom_x(1:kk),'YData',odom_y(1:kk));
        set(hEst, 'XData',X_hist(1,1:kk),'YData',X_hist(2,1:kk));
        % Rayos LiDAR en vivo desde la pose real (visualización)
        pr_v = [x_real(kk); y_real(kk); phi_real(kk)];
        rg_v = lidar_scan(pr_v(1), pr_v(2), pr_v(3), lidar_angles, lidar_max_range, beacons, beacon_radius);
        [hx, hy, mx, my, px, py] = rayos(pr_v, rg_v, lidar_angles, lidar_max_range, rayos_libres_cada);
        set(hRayHit, 'XData',hx,'YData',hy);
        set(hRayMiss,'XData',mx,'YData',my);
        set(hScan,   'XData',px,'YData',py);
        set(hZHit,   'XData',hx,'YData',hy);
        set(hZMiss,  'XData',mx,'YData',my);
        th = linspace(0,2*pi,60);
        set(hRange,'XData',x_real(kk)+lidar_max_range*cos(th), ...
                   'YData',y_real(kk)+lidar_max_range*sin(th));
        if isempty(last_obs)
            set(hRays,'XData',nan,'YData',nan);
        else
            nd = numel(last_obs);
            lx_ = [repmat(x_real(kk),1,nd); beacons(last_obs,1)'; nan(1,nd)];
            ly_ = [repmat(y_real(kk),1,nd); beacons(last_obs,2)'; nan(1,nd)];
            set(hRays,'XData',lx_(:),'YData',ly_(:));
        end
        actualizar_tanque(hTank,  x_real(kk), y_real(kk), phi_real(kk), L, B_ext);
        actualizar_tanque(hTankZ, x_real(kk), y_real(kk), phi_real(kk), L, B_ext);
        set(hZReal,'XData',x_real(1:kk),'YData',y_real(1:kk));
        set(hZEst, 'XData',X_hist(1,1:kk),'YData',X_hist(2,1:kk));
        set(hZScan,'XData',px,'YData',py);
        if ~isempty(lm_order)
            iz = lm_idx(lm_order);
            set(hZLm,'XData',X(iz),'YData',X(iz+1));
        end
        xlim(ax_zoom, x_real(kk) + [-zoom_w zoom_w]);
        ylim(ax_zoom, y_real(kk) + [-zoom_w zoom_w]);

        [ex,ey] = error_ellipse(P(1:2,1:2), X(1:2)', 3);
        set(hCovR,'XData',ex,'YData',ey);
        if ~isempty(lm_order)
            ii_ = lm_idx(lm_order);
            set(hLmPt,'XData',X(ii_),'YData',X(ii_+1));
            for j = lm_order
                i = lm_idx(j);
                [ex,ey] = error_ellipse(P(i:i+1,i:i+1), X(i:i+1)', 3);
                set(hLmEl(j),'XData',ex,'YData',ey);
            end
        end
        title(ax_map, sprintf('EKF-SLAM  |  t = %.1f s  |  WP objetivo: WP%d  |  beacons en el mapa: %d/%d', ...
              t(kk), min(wp_idx,size(WP,1))-1, numel(lm_order), n_beacons));

        set(hErrM,'XData',t(1:kk),'YData',err_mec(1:kk));
        set(hErrF,'XData',t(1:kk),'YData',err_fus(1:kk));
        set(hSig, 'XData',t(1:kk),'YData',sig3(1:kk));
        set(hNobs,'XData',t(1:kk),'YData',n_obs_log(1:kk));
        set(hNmap,'XData',t(1:kk),'YData',n_map_log(1:kk));
        set(hEmap,'XData',t(1:kk),'YData',err_map_log(1:kk));
        title(ax1, sprintf('Error medio robot EKF-SLAM = %.4f m', mean(err_fus(2:kk))));
        drawnow limitrate;
    end
end

%% === RESULTADOS ===
idx = 2:kf;
fprintf('\n=== EKF-SLAM (ruta cuadrada, tanque Hiwonder) ===\n');
fprintf('Tiempo de recorrido     : %.1f s\n', t(kf));
fprintf('Error medio mecánica    : %.4f m\n', mean(err_mec(idx)));
fprintf('Error final mecánica    : %.4f m\n', err_mec(kf));
fprintf('Error medio EKF-SLAM    : %.4f m\n', mean(err_fus(idx)));
fprintf('Error final EKF-SLAM    : %.4f m\n', err_fus(kf));
fprintf('Error dentro de 3sigma  : %.1f %% del tiempo\n', 100*mean(err_fus(idx) <= sig3(idx)));
fprintf('Beacons en el mapa      : %d/%d\n', numel(lm_order), n_beacons);
fprintf('Rechazos por gating     : %d\n', cnt_gate);
fprintf('\nBeacon |  x real  y real |  x est   y est |  error [m] | 3sigma máx [m]\n');
for j = lm_order
    i = lm_idx(j);
    fprintf('  B%-3d | %6.2f %6.2f | %6.2f %6.2f |   %.4f   |   %.4f\n', j, ...
        beacons(j,1), beacons(j,2), X(i), X(i+1), ...
        hypot(X(i)-beacons(j,1), X(i+1)-beacons(j,2)), ...
        3*sqrt(max(eig(P(i:i+1,i:i+1)))));
end

%% === FIGURA 3: MAPA FINAL ===
fig3 = figure(3); set(fig3,'Position',[100,60,900,850]); hold on; grid on; axis equal;
plot(beacons(:,1), beacons(:,2), 'k+', 'MarkerSize',10, 'LineWidth',1.6, 'DisplayName','Beacons reales');
plot(x_real(1:kf), y_real(1:kf), 'k-', 'LineWidth',1.2, 'DisplayName','Real');
plot(odom_x(1:kf), odom_y(1:kf), 'b--', 'LineWidth',1.0, 'DisplayName','Odom. mecánica');
plot(X_hist(1,1:kf), X_hist(2,1:kf), 'g-', 'LineWidth',1.2, 'DisplayName','EKF-SLAM');
if ~isempty(lm_order)
    ii_ = lm_idx(lm_order);
    plot(X(ii_), X(ii_+1), 'r.', 'MarkerSize',14, 'DisplayName','Beacons estimados + 3\sigma');
    for j = lm_order
        i = lm_idx(j);
        [ex,ey] = error_ellipse(P(i:i+1,i:i+1), X(i:i+1)', 3);
        plot(ex, ey, 'r-', 'LineWidth',1.1, 'HandleVisibility','off');
    end
end
xlim([-lim lim]); ylim([-lim lim]);
xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('EKF-SLAM — Mapa final (%d/%d beacons detectados)', numel(lm_order), n_beacons));
legend('Location','northeastoutside');
saveas(fig3, 'ekf_slam_mapa.png');

%% === FIGURA 4: MATRIZ DE COVARIANZA ===
n   = numel(X);
nL  = n - 3;
lbl = {'x_r','y_r','\theta_r'};
for j = lm_order
    lbl = [lbl, {sprintf('B%d_x',j), sprintf('B%d_y',j)}]; %#ok<AGROW>
end
tk = 1:max(1, floor(n/20)):n;     % hasta ~20 etiquetas

fig4 = figure(4); set(fig4,'Position',[60,60,1400,620]);
sgtitle(sprintf('Matriz de covarianza EKF-SLAM — estado %d\\times%d (%d beacons)', n, n, nL/2));

subplot(1,2,1);
imagesc(log10(abs(P) + 1e-12)); axis square; colormap(gca, parula);
cb = colorbar; cb.Label.String = 'log_{10}|P_{ij}|';
hold on;
rectangle('Position',[0.5 0.5 3 3], 'EdgeColor','c', 'LineWidth',2);
if nL > 0
    rectangle('Position',[3.5 3.5 nL nL], 'EdgeColor','g', 'LineWidth',1.5);
    rectangle('Position',[3.5 0.5 nL 3],  'EdgeColor',[1 0.6 0], 'LineWidth',1.5);
    rectangle('Position',[0.5 3.5 3 nL],  'EdgeColor',[1 0.6 0], 'LineWidth',1.5);
end
set(gca,'XTick',tk,'XTickLabel',lbl(tk),'YTick',tk,'YTickLabel',lbl(tk), ...
        'XTickLabelRotation',90,'FontSize',7);
title('Covarianza (log_{10}) — cian: robot, verde: mapa, naranja: robot-mapa');

subplot(1,2,2);
d = sqrt(diag(P));
C = P ./ (d*d');
imagesc(C, [-1 1]); axis square;
cm = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1);
      ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
colormap(gca, cm);
cb = colorbar; cb.Label.String = '\rho_{ij}';
set(gca,'XTick',tk,'XTickLabel',lbl(tk),'YTick',tk,'YTickLabel',lbl(tk), ...
        'XTickLabelRotation',90,'FontSize',7);
title('Correlación normalizada \rho_{ij} = P_{ij}/(\sigma_i\sigma_j)');
saveas(fig4, 'ekf_slam_cov.png');

%% === FIGURA 5: MAPA DE BEACONS (reales vs estimados, sin rutas) ===
fig5 = figure(5); set(fig5,'Position',[120,80,850,850]); hold on; grid on; box on; axis equal;
plot(beacons(:,1), beacons(:,2), 'ko', 'MarkerSize',7, 'MarkerFaceColor','k', ...
     'DisplayName','Beacons reales');
if ~isempty(lm_order)
    ii_ = lm_idx(lm_order);
    plot(X(ii_), X(ii_+1), 'ro', 'MarkerSize',6, 'MarkerFaceColor','r', ...
         'DisplayName','Centros estimados (EKF-SLAM)');
end
xlim([-lim lim]); ylim([-lim lim]);
xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Mapa de beacons — reales (negro) vs estimados (rojo)  |  %d/%d detectados', ...
      numel(lm_order), n_beacons));
legend('Location','northeastoutside');
saveas(fig5, 'ekf_slam_beacons.png');

fprintf('\n[ok] Figuras guardadas en: %s\n', pwd);

%% === FUNCIONES ===
function [V, w, wp_idx, fin] = control_waypoints(x, y, phi, WP, wp_idx, V_cte, umbral, k_w, w_max)
    % Go-to-waypoint con velocidad constante (equivale a get_control_waypoints
    % de Python, pero sin frenar al acercarse)
    fin = false;
    while wp_idx <= size(WP,1) && hypot(WP(wp_idx,1)-x, WP(wp_idx,2)-y) < umbral
        wp_idx = wp_idx + 1;
    end
    if wp_idx > size(WP,1)
        V = 0; w = 0; fin = true; return;
    end
    e = wrapToPi(atan2(WP(wp_idx,2)-y, WP(wp_idx,1)-x) - phi);
    V = V_cte;
    w = max(min(k_w*e, w_max), -w_max);
end

function r = lidar_scan(xs, ys, phis, angles, rmax, beacons, rb)
    a = phis + angles; dxr = cos(a); dyr = sin(a);
    r = rmax*ones(size(a));
    for b = 1:size(beacons,1)
        mx = xs-beacons(b,1); my = ys-beacons(b,2);
        bb = mx*dxr + my*dyr;
        disc = bb.^2 - (mx^2 + my^2 - rb^2);
        tt = -bb - sqrt(max(disc,0));
        ok = disc >= 0 & tt > 0 & tt < r;
        r(ok) = tt(ok);
    end
end

function [hx, hy, mx, my, px, py] = rayos(pose, rg, angles, rmax, cada)
    % Líneas de los rayos (separadas por NaN) para dibujar con un solo objeto
    a = pose(3) + angles;
    hit = rg < rmax;
    ex = pose(1) + rg.*cos(a); ey = pose(2) + rg.*sin(a);
    n = nnz(hit);
    hx = [repmat(pose(1),1,n); ex(hit); nan(1,n)]; hy = [repmat(pose(2),1,n); ey(hit); nan(1,n)];
    hx = hx(:); hy = hy(:);
    miss = find(~hit); miss = miss(1:cada:end); m = numel(miss);
    mx = [repmat(pose(1),1,m); ex(miss); nan(1,m)]; my = [repmat(pose(2),1,m); ey(miss); nan(1,m)];
    mx = mx(:); my = my(:);
    px = ex(hit); py = ey(hit);
    if n == 0, hx = nan; hy = nan; px = nan; py = nan; end
    if m == 0, mx = nan; my = nan; end
end

function a = wrapToPi(a)
    a = mod(a+pi, 2*pi) - pi;
end

function h = crear_tanque(ax)
    % Orugas, chasis e indicador de frente del Hiwonder
    h.trackL = patch(ax,0,0,[0.15 0.15 0.15],'EdgeColor','k','HandleVisibility','off');
    h.trackR = patch(ax,0,0,[0.15 0.15 0.15],'EdgeColor','k','HandleVisibility','off');
    h.treadL = plot(ax,nan,nan,'-','Color',[0.55 0.55 0.55],'HandleVisibility','off');
    h.treadR = plot(ax,nan,nan,'-','Color',[0.55 0.55 0.55],'HandleVisibility','off');
    h.body   = patch(ax,0,0,[0.2 0.4 0.8],'FaceAlpha',0.9,'EdgeColor','k','HandleVisibility','off');
    h.front  = patch(ax,0,0,[1 0.85 0.1],'EdgeColor','k','HandleVisibility','off');
end

function actualizar_tanque(h, x, y, phi, L, W)
    % L: largo, W: ancho exterior (de oruga a oruga)
    Rm = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    tw = 0.045;                              % ancho de cada oruga [m]
    yc = W/2 - tw/2;                         % centro lateral de cada oruga
    trk = [-L/2 L/2 L/2 -L/2; -tw/2 -tw/2 tw/2 tw/2];
    tL = Rm*(trk + [0;  yc]);
    tR = Rm*(trk + [0; -yc]);
    wi = W/2 - tw;                           % semiancho del chasis
    bd = Rm*[-0.42*L 0.42*L 0.42*L -0.42*L; -wi -wi wi wi];
    fr = Rm*[0.40*L 0.18*L 0.18*L; 0 0.6*wi -0.6*wi];
    % líneas de las zapatas de la oruga
    xs = linspace(-L/2, L/2, 9);
    tx = [xs; xs; nan(1,9)]; ty = [-tw/2*ones(1,9); tw/2*ones(1,9); nan(1,9)];
    zL = Rm*[tx(:)'; ty(:)' + yc];
    zR = Rm*[tx(:)'; ty(:)' - yc];
    set(h.trackL,'XData',tL(1,:)+x,'YData',tL(2,:)+y);
    set(h.trackR,'XData',tR(1,:)+x,'YData',tR(2,:)+y);
    set(h.treadL,'XData',zL(1,:)+x,'YData',zL(2,:)+y);
    set(h.treadR,'XData',zR(1,:)+x,'YData',zR(2,:)+y);
    set(h.body,  'XData',bd(1,:)+x,'YData',bd(2,:)+y);
    set(h.front, 'XData',fr(1,:)+x,'YData',fr(2,:)+y);
end

function [x, y] = error_ellipse(C, c, ns)
    C = (C+C')/2;
    [V, D] = eig(C); [d, i] = sort(max(diag(D),1e-12),'descend'); V = V(:,i);
    t = linspace(0,2*pi,100);
    e = V*[ns*sqrt(d(1))*cos(t); ns*sqrt(d(2))*sin(t)];
    x = e(1,:)+c(1); y = e(2,:)+c(2);
end