%% tanque_sim_EKF_beacons_lidar_v2.m
% Robot oruga (Hiwonder) en pista sinusoidal de 12 m con:
%   - Odometría mecánica (encoders con deslizamiento y ruido)  -> PREDICCIÓN EKF
%   - LiDAR 2D simulado (ray casting SOLO contra beacons)      -> visualización
%   - Beacons de posición CONOCIDA medidos en range-bearing     -> CORRECCIÓN EKF
%
% Cambios respecto a la versión anterior:
%   * Sin paredes: el entorno es abierto, el LiDAR solo impacta beacons.
%   * Velocidad comandada constante desde t = 0 (sin rampa de aceleración).
%   * La simulación se detiene cuando el robot llega al final de la pista.
%
% Supuestos:
%   * Asociación de datos conocida: cada beacon tiene ID.
%   * El LiDAR "ve" un beacon si está dentro de su alcance.
%   * Entre x ~ 5 m y x ~ 7 m no hay beacons al alcance: el EKF queda solo
%     con encoders y se ve crecer la elipse de incertidumbre.

clear; close all; clc;
rng(1);   % reproducible (comentar para resultados aleatorios)

%% === PARÁMETROS DEL VEHÍCULO (Hiwonder) ===
B_ext = 0.19343;            % Distancia entre orugas (exterior) [m]
B_int = 0.102;              % Distancia entre orugas (interior) [m]
B     = (B_ext + B_int)/2;  % Separación efectiva [m]
L     = 0.306;              % Longitud del chasis [m]

%% === PARÁMETROS DE SIMULACIÓN ===
dt      = 0.02;             % Paso de tiempo [s]
Tsim    = 75;               % Tiempo máximo [s]
t       = 0:dt:Tsim;
N       = length(t);
V_cte   = 0.35;             % Velocidad CONSTANTE del tanque [m/s]
w_max   = 1.0;              % [rad/s]

%% === CONTROLADOR STANLEY ===
k1 = 0.3; k2 = 0.3;

%% === DESLIZAMIENTO Y RUIDO DE ENCODERS ===
sigma_encoder = 0.02;       % [m/s]
io_mean = 0.08; ii_mean = 0.10; sigma_slip = 0.05;

%% === TRAYECTORIA SINUSOIDAL 12 m ===
dist_total = 12;
X_path = linspace(0, dist_total, 1200);
A = 1.2; waves = 3; f = waves/dist_total;
Y_path = A*sin(2*pi*f*X_path);
phi_traj = atan2(gradient(Y_path), gradient(X_path));

%% === ENTORNO: SOLO BEACONS (sin paredes) ===
beacons = [-0.5  2.0;               % posiciones conocidas [x y]
            1.0 -2.0;
            2.5  2.0;
            9.5 -2.0;
           11.0  2.0;
           12.5 -2.0];
n_beacons     = size(beacons,1);
beacon_radius = 0.05;               % poste cilíndrico [m]

%% === LIDAR 2D ===
lidar_rate_hz     = 10;
lidar_every       = max(1, round(1/(lidar_rate_hz*dt)));
lidar_max_range   = 3.0;            % alcance [m] (bajarlo agranda la zona ciega)
lidar_n_rays      = 180;
lidar_angles      = linspace(-pi, pi, lidar_n_rays+1); lidar_angles(end) = [];
sigma_lidar_range = 0.02;           % ruido por rayo [m]

% Medición a beacons (range-bearing)
sigma_b_range   = 0.03;             % [m]
sigma_b_bearing = deg2rad(1.0);     % [rad]
R_beacon  = diag([sigma_b_range^2, sigma_b_bearing^2]);
gate_chi2 = 9.21;                   % gating Mahalanobis (2 GDL, 99%)

%% === ESTADOS ===
x_real = zeros(1,N); y_real = zeros(1,N); phi_real = zeros(1,N);
odom_mec_x = zeros(1,N); odom_mec_y = zeros(1,N); odom_mec_phi = zeros(1,N);

X_est = zeros(3,N);                 % [x; y; phi]
P_est = zeros(3,3,N);
P_est(:,:,1) = diag([0.01, 0.01, deg2rad(1)]).^2;
Q_enc = diag([0.02, 0.02, deg2rad(0.3)]).^2;
Q_u   = diag([sigma_encoder^2, (sigma_encoder/(B/2))^2]);

err_mec   = zeros(1,N);
err_fus   = zeros(1,N);
sig3_fus  = zeros(1,N);
n_vis_log = zeros(1,N);

scan_x = nan; scan_y = nan; last_det = []; last_n_vis = 0;
kf = N;                             % último índice simulado

%% === FIGURA 1: MAPA ===
figure(1); set(gcf,'Position',[58,36,1400,480]); hold on; grid on; axis equal;
plot(X_path, Y_path, 'k--', 'LineWidth',1.0, 'DisplayName','Trayectoria deseada');
plot(beacons(:,1), beacons(:,2), 'k^', 'MarkerSize',10, ...
     'MarkerFaceColor',[1 0.85 0], 'DisplayName','Beacons');
for j = 1:n_beacons
    text(beacons(j,1)+0.12, beacons(j,2), sprintf('B%d',j), 'FontSize',8);
end
xlim([-1.3, 13.3]); ylim([-2.8, 2.8]);
xlabel('X [m]'); ylabel('Y [m]');
title('Robot oruga – EKF con encoders + LiDAR/beacons');

hScan   = plot(nan,nan,'.', 'Color',[1 0.5 0], 'MarkerSize',8, 'DisplayName','Scan LiDAR');
hRange  = plot(nan,nan,':', 'Color',[1 0.5 0], 'HandleVisibility','off');
hRays   = plot(nan,nan,'-', 'Color',[0 0.7 0.7], 'LineWidth',1.0, 'DisplayName','Beacon detectado');
hBody   = patch(0,0,'b','FaceAlpha',0.5,'EdgeColor','k','HandleVisibility','off');
hTrackL = patch(0,0,[0.2 0.2 0.2],'FaceAlpha',0.7,'EdgeColor','none','HandleVisibility','off');
hTrackR = patch(0,0,[0.2 0.2 0.2],'FaceAlpha',0.7,'EdgeColor','none','HandleVisibility','off');
hPathReal = plot(nan,nan,'g-', 'LineWidth',1.6, 'DisplayName','Real');
hPathMec  = plot(nan,nan,'r--','LineWidth',1.2, 'DisplayName','Odom. mecánica');
hPathEst  = plot(nan,nan,'b-', 'LineWidth',1.5, 'DisplayName','EKF (enc + beacons)');
hCov      = plot(nan,nan,'Color',[0.5 0 0],'LineWidth',1.0,'DisplayName','Elipse 3\sigma EKF');
legend('Location','northeastoutside');

%% === FIGURA 2: ERRORES Y BEACONS VISIBLES ===
figure(2); set(gcf,'Position',[252,380,1032,420]);
ax1 = subplot(2,1,1); hold on; grid on; box on;
hErrM = plot(nan,nan,'r--','LineWidth',1.4,'DisplayName','Error mecánica');
hErrF = plot(nan,nan,'b-', 'LineWidth',1.6,'DisplayName','Error EKF');
hSig  = plot(nan,nan,'k:', 'LineWidth',1.2,'DisplayName','Cota 3\sigma EKF');
ylabel('Error posición [m]'); legend('Location','northwest');
ax2 = subplot(2,1,2); hold on; grid on; box on;
hNvis = stairs(nan,nan,'Color',[0 0.6 0.6],'LineWidth',1.4);
xlabel('Tiempo [s]'); ylabel('# beacons'); ylim([0 n_beacons]);

%% === SIMULACIÓN ===
for k = 1:N-1
    %% --- Control (velocidad constante, Stanley solo corrige el giro) ---
    V_cmd = V_cte;
    [e_lat, theta_err] = lateral_error_and_heading(x_real(k), y_real(k), phi_real(k), X_path, Y_path, phi_traj);
    delta = -(k1*theta_err + atan(k2*e_lat/V_cmd));
    delta = max(min(delta, pi/3), -pi/3);
    if abs(tan(delta)) < 1e-6
        omega_yaw = 0;
    else
        omega_yaw = V_cmd / (L/tan(delta));
    end
    omega_yaw = max(min(omega_yaw, w_max), -w_max);

    %% --- Deslizamiento ---
    if t(k) >= 15 && t(k) <= 17
        io = 0.15 + 0.05*randn;  ii = 0.18 + 0.05*randn;
    else
        io = io_mean + sigma_slip*randn;  ii = ii_mean + sigma_slip*randn;
    end
    Vo_true = (V_cmd + (B/2)*omega_yaw)*(1-io);
    Vi_true = (V_cmd - (B/2)*omega_yaw)*(1-ii);

    %% --- Ground truth (k -> k+1) ---
    Vcg_real = (Vo_true + Vi_true)/2;
    w_real   = (Vo_true - Vi_true)/B;
    x_real(k+1)   = x_real(k) + dt*Vcg_real*cos(phi_real(k));
    y_real(k+1)   = y_real(k) + dt*Vcg_real*sin(phi_real(k));
    phi_real(k+1) = wrapToPi(phi_real(k) + dt*w_real);

    %% --- Encoders ---
    Vo_meas = Vo_true + sigma_encoder*randn;
    Vi_meas = Vi_true + sigma_encoder*randn;
    Vcg_enc = (Vo_meas + Vi_meas)/2;
    w_enc   = (Vo_meas - Vi_meas)/B;

    odom_mec_phi(k+1) = wrapToPi(odom_mec_phi(k) + dt*w_enc);
    odom_mec_x(k+1)   = odom_mec_x(k) + dt*Vcg_enc*cos(odom_mec_phi(k));
    odom_mec_y(k+1)   = odom_mec_y(k) + dt*Vcg_enc*sin(odom_mec_phi(k));

    %% --- Predicción EKF ---
    xp = X_est(:,k); Pp = P_est(:,:,k); ph = xp(3);
    F  = [1 0 -dt*Vcg_enc*sin(ph);
          0 1  dt*Vcg_enc*cos(ph);
          0 0  1];
    Wu = [dt*cos(ph) 0; dt*sin(ph) 0; 0 dt];
    x_upd = [xp(1) + dt*Vcg_enc*cos(ph);
             xp(2) + dt*Vcg_enc*sin(ph);
             wrapToPi(ph + dt*w_enc)];
    P_upd = F*Pp*F' + Wu*Q_u*Wu' + Q_enc;

    %% --- LiDAR + corrección con beacons ---
    if mod(k, lidar_every) == 0
        xr = x_real(k+1); yr = y_real(k+1); pr = phi_real(k+1);

        % Scan (solo impacta beacons; rayos sin impacto = rango máximo)
        ranges = lidar_scan(xr, yr, pr, lidar_angles, lidar_max_range, beacons, beacon_radius);
        hit = ranges < lidar_max_range;
        ranges(hit) = ranges(hit) + sigma_lidar_range*randn(1, nnz(hit));
        scan_x = xr + ranges(hit).*cos(pr + lidar_angles(hit));
        scan_y = yr + ranges(hit).*sin(pr + lidar_angles(hit));

        % Detección de beacons y actualización secuencial del EKF
        last_det = [];
        for j = 1:n_beacons
            dxr = beacons(j,1) - xr; dyr = beacons(j,2) - yr;
            rr  = hypot(dxr, dyr);
            if rr > lidar_max_range, continue; end

            z = [rr + sigma_b_range*randn;
                 wrapToPi(atan2(dyr,dxr) - pr + sigma_b_bearing*randn)];

            dxe = beacons(j,1) - x_upd(1); dye = beacons(j,2) - x_upd(2);
            q = dxe^2 + dye^2; re = sqrt(q);
            z_hat = [re; wrapToPi(atan2(dye,dxe) - x_upd(3))];
            Hb = [-dxe/re, -dye/re,  0;
                   dye/q,  -dxe/q,  -1];

            innov = z - z_hat; innov(2) = wrapToPi(innov(2));
            S = Hb*P_upd*Hb' + R_beacon;
            if innov'/S*innov > gate_chi2, continue; end   % outlier

            K = P_upd*Hb'/S;
            x_upd = x_upd + K*innov;
            x_upd(3) = wrapToPi(x_upd(3));
            IKH = eye(3) - K*Hb;
            P_upd = IKH*P_upd*IKH' + K*R_beacon*K';         % forma de Joseph
            last_det(end+1) = j; %#ok<AGROW>
        end
        last_n_vis = numel(last_det);
    end
    X_est(:,k+1)   = x_upd;
    P_est(:,:,k+1) = (P_upd + P_upd')/2;
    n_vis_log(k+1) = last_n_vis;

    %% --- Errores ---
    err_mec(k+1)  = hypot(odom_mec_x(k+1)-x_real(k+1), odom_mec_y(k+1)-y_real(k+1));
    err_fus(k+1)  = hypot(X_est(1,k+1)-x_real(k+1),   X_est(2,k+1)-y_real(k+1));
    sig3_fus(k+1) = 3*sqrt(max(eig(P_est(1:2,1:2,k+1))));

    fin_pista = x_real(k+1) >= dist_total;

    %% --- Gráficos ---
    if mod(k,5) == 0 || fin_pista
        kk = k+1;
        set(hPathReal,'XData',x_real(1:kk),     'YData',y_real(1:kk));
        set(hPathMec, 'XData',odom_mec_x(1:kk), 'YData',odom_mec_y(1:kk));
        set(hPathEst, 'XData',X_est(1,1:kk),    'YData',X_est(2,1:kk));
        set(hScan, 'XData',scan_x, 'YData',scan_y);

        th = linspace(0,2*pi,60);
        set(hRange,'XData',x_real(kk)+lidar_max_range*cos(th), ...
                   'YData',y_real(kk)+lidar_max_range*sin(th));

        if isempty(last_det)
            set(hRays,'XData',nan,'YData',nan);
        else
            nd = numel(last_det);
            lx = [repmat(x_real(kk),1,nd); beacons(last_det,1)'; nan(1,nd)];
            ly = [repmat(y_real(kk),1,nd); beacons(last_det,2)'; nan(1,nd)];
            set(hRays,'XData',lx(:),'YData',ly(:));
        end

        [bx,by] = dibujar_rect(x_real(kk), y_real(kk), phi_real(kk), L, B*0.8);
        set(hBody,'XData',bx,'YData',by);
        [txL,tyL] = createOval(L*0.9, min(0.06,0.4*B), 0, -B/2);
        [txR,tyR] = createOval(L*0.9, min(0.06,0.4*B), 0,  B/2);
        Rm = [cos(phi_real(kk)) -sin(phi_real(kk)); sin(phi_real(kk)) cos(phi_real(kk))];
        tL = Rm*[txL;tyL]; tR = Rm*[txR;tyR];
        set(hTrackL,'XData',tL(1,:)+x_real(kk),'YData',tL(2,:)+y_real(kk));
        set(hTrackR,'XData',tR(1,:)+x_real(kk),'YData',tR(2,:)+y_real(kk));

        cxy = P_est(1:2,1:2,kk);
        if all(isfinite(cxy(:))) && all(eig(cxy) > 0)
            [ex,ey] = error_ellipse(cxy, X_est(1:2,kk)', 3);
            set(hCov,'XData',ex,'YData',ey);
        end

        set(hErrM,'XData',t(1:kk),'YData',err_mec(1:kk));
        set(hErrF,'XData',t(1:kk),'YData',err_fus(1:kk));
        set(hSig, 'XData',t(1:kk),'YData',sig3_fus(1:kk));
        set(hNvis,'XData',t(1:kk),'YData',n_vis_log(1:kk));
        title(ax1, sprintf('Error medio EKF = %.4f m', mean(err_fus(1:kk))));
        drawnow limitrate;
    end

    if fin_pista
        kf = k+1;
        break;
    end
end

%% === RESULTADOS ===
idx = 2:kf;
fprintf('\n=== RESULTADOS (pista 12 m, EKF encoders + beacons) ===\n');
fprintf('Tiempo de recorrido  : %.1f s\n', t(kf));
fprintf('Error medio mecánica : %.4f m\n', mean(err_mec(idx)));
fprintf('Error medio EKF      : %.4f m\n', mean(err_fus(idx)));
fprintf('Error máximo EKF     : %.4f m\n', max(err_fus(idx)));
fprintf('Tiempo sin beacons   : %.1f %%\n', 100*mean(n_vis_log(idx)==0));
fprintf('Error EKF dentro de 3sigma: %.1f %% del tiempo\n\n', ...
        100*mean(err_fus(idx) <= sig3_fus(idx)));

%% === FUNCIONES AUXILIARES ===
function r = lidar_scan(xs, ys, phis, angles, rmax, beacons, rb)
    % Ray casting 2D solo contra círculos (beacons); sin impacto -> rmax
    a = phis + angles; dxr = cos(a); dyr = sin(a);
    r = rmax*ones(size(a));
    for b = 1:size(beacons,1)
        mx = xs - beacons(b,1); my = ys - beacons(b,2);
        bb = mx*dxr + my*dyr;
        disc = bb.^2 - (mx^2 + my^2 - rb^2);
        tt = -bb - sqrt(max(disc,0));
        ok = disc >= 0 & tt > 0 & tt < r;
        r(ok) = tt(ok);
    end
end

function [e_lat, theta_err] = lateral_error_and_heading(xv, yv, phi_v, Xp, Yp, phi_p)
    d2 = (Xp - xv).^2 + (Yp - yv).^2;
    [~, idx] = min(d2);
    vx = xv - Xp(idx); vy = yv - Yp(idx);
    cross_z = cos(phi_p(idx))*vy - sin(phi_p(idx))*vx;
    e_lat = sign(cross_z)*sqrt(d2(idx));
    theta_err = wrapToPi(phi_v - phi_p(idx));
end

function angle = wrapToPi(angle)
    angle = mod(angle + pi, 2*pi) - pi;
end

function [X, Y] = dibujar_rect(cx, cy, theta, L, W)
    rect = [-L/2 -W/2; L/2 -W/2; L/2 W/2; -L/2 W/2]';
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    pts = R*rect;
    X = pts(1,:) + cx; Y = pts(2,:) + cy;
end

function [x, y] = createOval(a, b, cx, cy)
    a = a/4*3; b = b/2;
    x = [a, a, -a, -a, a] + cx;
    y = [b, -b, -b, b, b] + cy;
end

function [x, y] = error_ellipse(cov, center, nsigma)
    [V, D] = eig(cov);
    [Dv, idx] = sort(diag(D),'descend'); V = V(:,idx);
    th = linspace(0, 2*pi, 100);
    el = V * [nsigma*sqrt(Dv(1))*cos(th); nsigma*sqrt(Dv(2))*sin(th)];
    x = el(1,:) + center(1); y = el(2,:) + center(2);
end