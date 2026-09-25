%% tanque_sim_EKFSLAM_sin_paredes.m
% EKF-SLAM con el robot oruga Hiwonder sobre el MISMO entorno sin paredes
% que la versión de SLAM de grafos. Sin toolboxes.
%
% Idea del método:
%   1) ESTADO    : mu = [x y phi | m1x m1y | m2x m2y | ...] y su covarianza P.
%                  El robot y TODOS los landmarks se estiman juntos.
%   2) PREDICCIÓN: en cada dt se propaga la pose con la odometría de los
%                  encoders (modelo diferencial de orugas).
%   3) LANDMARKS : como no hay beacons, se extraen del LiDAR las ESQUINAS de
%                  las cajas (split-and-merge + intersección de rectas).
%                  Una esquina es un punto fijo del mundo: no depende del
%                  punto de vista (a diferencia del centroide de un cluster).
%   4) ASOCIACIÓN: vecino más cercano con distancia de Mahalanobis.
%                  Si ningún landmark es compatible -> landmark nuevo.
%   5) CORRECCIÓN: actualización EKF secuencial con modelo rango-rumbo.
%                  Re-observar un landmark antiguo = "cierre de lazo" implícito.
%   6) MAPA      : grilla de ocupación con los scans integrados en la pose
%                  estimada (cada 0.3 m o 15°, igual que los nodos del grafo).
%
% Diferencia clave con la versión de grafos: aquí NO se corrigen poses
% pasadas; la grilla se construye con la pose estimada en ese instante.

clear; close all; clc;
rng(1);

%% === VEHÍCULO (Hiwonder) ===
B_ext = 0.19343; B_int = 0.102;
B = (B_ext + B_int)/2;            % separación efectiva orugas [m]
L = 0.306;                        % largo chasis [m]

%% === SIMULACIÓN ===
dt    = 0.02;
Tsim  = 600;                      % se corta al terminar la ruta
t     = 0:dt:Tsim;
N     = length(t);
V_cte = 0.35;                     % velocidad constante [m/s]
w_max = 1.5; k_w = 2.5;
plot_every = 10;                  % un cuadro cada 10*dt (siempre igual)
vel_anim   = 1;                   % 1 = tiempo real, 2 = x2 ... (0 = lo más rápido posible)

%% === ENCODERS Y DESLIZAMIENTO ===
sigma_encoder = 0.02;
io_mean = 0.08; ii_mean = 0.10; sigma_slip = 0.05;

%% === WAYPOINTS (ruta cuadrada del script de Python, escalada) ===
escala = 1/15;
WP = escala * [0 0; 75 0; 75 75; -75 75; -75 -75; -75 -75; 75 -75; -70 70];
umbral_wp = 5*escala;

%% === ENTORNO: espacio abierto con obstáculos (SIN paredes) ===
room = 7;
segs = zeros(0,4);
obst = [ % cx    cy    ancho alto
        -2.0   6.2   0.6   0.6;
         2.5   6.2   0.8   0.5;
         6.2   2.0   0.5   0.8;
         6.2  -3.0   0.6   0.6;
        -6.2  -1.0   0.5   1.0;
        -6.2   3.5   0.6   0.6;
         1.0  -6.2   1.0   0.5;
        -3.0  -6.2   0.6   0.6;
         2.5   2.5   1.0   1.0;
        -2.5  -2.5   1.0   1.0;
         3.5  -2.0   0.4   0.4;
        -3.5   1.5   0.4   0.4;
         0.0   3.5   0.5   0.5;
        -1.5   3.8   0.6   0.3;
         1.5  -3.5   0.5   0.5;
        -1.0  -1.8   0.4   0.4;
         4.0   1.5   0.4   0.4;
         6.5   6.5   0.8   0.8;
        -6.5   6.5   0.8   0.8;
         6.5  -6.5   0.8   0.8;
        -6.5  -6.5   0.8   0.8];
C_true = zeros(2,0);                           % esquinas reales (solo para evaluar)
for i = 1:size(obst,1)
    segs = [segs; rect_segs(obst(i,:))]; %#ok<AGROW>
    o = obst(i,:);
    x1 = o(1)-o(3)/2; x2 = o(1)+o(3)/2; y1 = o(2)-o(4)/2; y2 = o(2)+o(4)/2;
    C_true = [C_true, [x1 x2 x2 x1; y1 y1 y2 y2]]; %#ok<AGROW>
end

%% === LIDAR (tipo RPLIDAR) ===
lidar_max_range = 5.0;
lidar_n_rays    = 360;                         % 1°
lidar_angles    = linspace(-pi, pi, lidar_n_rays+1); lidar_angles(end) = [];
sigma_lidar     = 0.02;
lidar_every     = 5;                           % scan cada 5*dt = 10 Hz
rayos_libres_cada = 3;

%% === EXTRACTOR DE ESQUINAS ===
ext.d_clust   = 0.20;          % dist. máx. entre puntos consecutivos de un cluster [m]
ext.thr_split = 0.06;          % umbral split-and-merge [m] (~3 sigma_lidar)
ext.min_pts   = 4;             % puntos mínimos por cara
ext.min_len   = 0.12;          % largo mínimo de cada cara [m]
ext.max_cos   = 0.30;          % |cos| entre caras -> ángulo entre ~72° y ~108°
ext.max_gap   = 0.10;          % distancia máx. intersección <-> punto de quiebre [m]
ext.r_max     = 4.5;           % solo esquinas a menos de 4.5 m

%% === PARÁMETROS DEL EKF ===
sig_V = 1.5*sigma_encoder/sqrt(2);             % ruido de V de encoders (inflado x1.5)
sig_W = 1.5*sqrt(2)*sigma_encoder/B;           % ruido de omega de encoders
Qu = diag([sig_V^2, sig_W^2]);
R  = diag([0.05^2, deg2rad(2)^2]);             % ruido medición esquina [rango rumbo]
chi2_asoc   = 9.21;            % gate de asociación (chi2, 2 gdl, 99%)
chi2_nuevo  = 25;              % por encima de esto -> candidato a landmark nuevo
d_min_nuevo = 0.25;            % no crear landmark a menos de 0.25 m de otro
r_cand      = 1.0;             % solo se prueban landmarks a < 1 m de la predicción
t_relazo    = 30;              % re-observar tras >30 s cuenta como cierre de lazo
Nlm_max     = 200;

mu = zeros(3,1);
P  = diag([1e-6 1e-6 1e-6]);
lm_visto = zeros(1,0);         % último instante en que se vio cada landmark
n_asoc = 0; n_nuevos = 0; n_rech = 0; n_relazo = 0; n_esq = 0; t_lc = [];
obs_ids = [];

%% === MAPA (keyframes como los nodos del grafo) ===
nodo_dist = 0.30; nodo_ang = deg2rad(15);
pose_kf = []; n_kf = 0;

%% === GRILLA DE OCUPACIÓN ===
map.res  = 0.05;
map.xmin = -room-0.5; map.ymin = -room-0.5;
map.nx   = round(2*(room+0.5)/map.res); map.ny = map.nx;
l_occ = 0.85; l_free = -0.4; l_max = 5;
Lmap = zeros(map.ny, map.nx);
xs_map = map.xmin + ((1:map.nx)-0.5)*map.res;
ys_map = map.ymin + ((1:map.ny)-0.5)*map.res;

%% === REGISTROS ===
x_real = zeros(1,N); y_real = zeros(1,N); phi_real = zeros(1,N);
odom_x = zeros(1,N); odom_y = zeros(1,N); odom_phi = zeros(1,N);
x_slam = zeros(3,N);
err_odo = zeros(1,N); err_slam = zeros(1,N); sig3 = zeros(1,N);
wp_idx = 1; kf = N;

%% === FIGURA 1 ===
fig1 = figure(1); set(fig1,'Position',[30,60,1500,720]);
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal'); box(ax1,'on');
hRayMiss = plot(ax1, nan,nan,'-', 'Color',[1 0.87 0.78], 'LineWidth',0.5, 'DisplayName','Rayo sin impacto');
hRayHit  = plot(ax1, nan,nan,'-', 'Color',[1 0.35 0.35], 'LineWidth',0.6, 'DisplayName','Rayo con impacto');
for i = 1:size(obst,1)
    dibujar_caja(ax1, obst(i,:), [0.35 0.35 0.35]);
end
plot(ax1, C_true(1,:), C_true(2,:), '+', 'Color',[0.6 0.6 0.6], 'MarkerSize',5, 'DisplayName','Esquinas reales');
plot(ax1, WP(:,1), WP(:,2), 's--', 'Color',[0.5 0 0.5], 'DisplayName','Waypoints');
hScan  = plot(ax1, nan,nan,'.', 'Color',[0.85 0 0], 'MarkerSize',7, 'DisplayName','Impactos LiDAR');
hReal  = plot(ax1, nan,nan,'k-', 'LineWidth',1.5, 'DisplayName','Real');
hOdo   = plot(ax1, nan,nan,'b--','LineWidth',1.2, 'DisplayName','Odometría');
hSlam  = plot(ax1, nan,nan,'g-', 'LineWidth',1.5, 'DisplayName','EKF-SLAM');
hObs   = plot(ax1, nan,nan,'c-', 'LineWidth',0.8, 'DisplayName','Asociaciones');
hEll   = plot(ax1, nan,nan,'m-', 'LineWidth',0.7, 'HandleVisibility','off');
hLm    = plot(ax1, nan,nan,'md', 'MarkerFaceColor','m', 'MarkerSize',4, 'DisplayName','Landmarks estimados (3\sigma)');
hEllR  = plot(ax1, nan,nan,'-', 'Color',[1 0.5 0], 'LineWidth',1.2, 'DisplayName','Pose 3\sigma');
hTank1 = crear_tanque(ax1);
xlim(ax1, [-room-0.5 room+0.5]); ylim(ax1, [-room-0.5 room+0.5]);
xlabel(ax1,'X [m]'); ylabel(ax1,'Y [m]'); title(ax1,'Mundo real, rayos LiDAR y landmarks (esquinas)');
legend(ax1,'Location','southoutside','Orientation','horizontal','NumColumns',4);

ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); box(ax2,'on');
hMap = imagesc(ax2, xs_map, ys_map, 0.5*ones(map.ny,map.nx));
set(ax2,'YDir','normal'); colormap(ax2, gray); caxis(ax2,[0 1]);
hSlam2 = plot(ax2, nan,nan,'g-', 'LineWidth',1.3);
hScan2 = plot(ax2, nan,nan,'.', 'Color',[1 0.3 0.3], 'MarkerSize',6);
hEll2  = plot(ax2, nan,nan,'m-', 'LineWidth',0.7);
hLm2   = plot(ax2, nan,nan,'md', 'MarkerFaceColor','m', 'MarkerSize',4);
hTank2 = crear_tanque(ax2);
xlim(ax2, [-room-0.5 room+0.5]); ylim(ax2, [-room-0.5 room+0.5]);
xlabel(ax2,'X [m]'); ylabel(ax2,'Y [m]');
title(ax2,'Mapa construido (grilla de ocupación)');

%% === FIGURA 2: ERRORES ===
fig2 = figure(2); set(fig2,'Position',[250,300,1000,380]);
ax3 = axes(fig2); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
hEo = plot(ax3, nan,nan,'b--','LineWidth',1.4,'DisplayName','Error odometría');
hEs = plot(ax3, nan,nan,'g-', 'LineWidth',1.6,'DisplayName','Error EKF-SLAM');
hSg = plot(ax3, nan,nan,'g:', 'LineWidth',1.2,'DisplayName','Cota 3\sigma EKF');
hEl = plot(ax3, nan,nan,'rv', 'MarkerFaceColor','r','DisplayName','Re-observación (cierre de lazo)');
xlabel(ax3,'Tiempo [s]'); ylabel(ax3,'Error de posición [m]');
legend(ax3,'Location','northwest');

%% === SIMULACIÓN ===
mapa_pendiente = false;
t_reloj = tic;
for k = 1:N-1
    %% --- Control por waypoints ---
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
    % Velocidad lineal REAL constante (= V_cte): el deslizamiento solo
    % perturba la velocidad de giro, no el avance.
    W_real = (Vo - Vi)/B;
    Vo = V_cmd + (B/2)*W_real;
    Vi = V_cmd - (B/2)*W_real;
    x_real(k+1)   = x_real(k) + dt*(Vo+Vi)/2*cos(phi_real(k));
    y_real(k+1)   = y_real(k) + dt*(Vo+Vi)/2*sin(phi_real(k));
    phi_real(k+1) = wrapToPi(phi_real(k) + dt*(Vo-Vi)/B);

    %% --- Encoders -> odometría ---
    Vom = Vo + sigma_encoder*randn; Vim = Vi + sigma_encoder*randn;
    Venc = (Vom+Vim)/2; Wenc = (Vom-Vim)/B;
    odom_phi(k+1) = wrapToPi(odom_phi(k) + dt*Wenc);
    odom_x(k+1)   = odom_x(k) + dt*Venc*cos(odom_phi(k));
    odom_y(k+1)   = odom_y(k) + dt*Venc*sin(odom_phi(k));

    %% ===== EKF: PREDICCIÓN (solo afecta al bloque del robot) =====
    phi = mu(3);
    mu(1) = mu(1) + dt*Venc*cos(phi);
    mu(2) = mu(2) + dt*Venc*sin(phi);
    mu(3) = wrapToPi(mu(3) + dt*Wenc);
    Fx = [1 0 -dt*Venc*sin(phi); 0 1 dt*Venc*cos(phi); 0 0 1];
    Gu = [dt*cos(phi) 0; dt*sin(phi) 0; 0 dt];
    P(1:3,:) = Fx*P(1:3,:);
    P(:,1:3) = P(:,1:3)*Fx';
    P(1:3,1:3) = P(1:3,1:3) + Gu*Qu*Gu';

    %% ===== EKF: CORRECCIÓN con esquinas del LiDAR =====
    if mod(k, lidar_every) == 0
        pr = [x_real(k+1); y_real(k+1); phi_real(k+1)];
        rg = lidar_scan(pr, lidar_angles, lidar_max_range, segs);
        hit = rg < lidar_max_range;
        rg(hit) = rg(hit) + sigma_lidar*randn(1, nnz(hit));

        Z = extraer_esquinas(rg, lidar_angles, lidar_max_range, ext);   % [r; b]
        n_esq = n_esq + size(Z,2);
        obs_ids = [];
        hubo_relazo = false;

        for q = 1:size(Z,2)
            z = Z(:,q);
            nlm = (numel(mu)-3)/2;
            a  = mu(3) + z(2);
            mz = mu(1:2) + z(1)*[cos(a); sin(a)];      % esquina en el mundo (predicha)

            best = 0; d2best = inf; deu_min = inf;
            if nlm > 0
                LM  = reshape(mu(4:end), 2, []);
                deu = hypot(LM(1,:)-mz(1), LM(2,:)-mz(2));
                deu_min = min(deu);
                for i = find(deu < r_cand)
                    id = [1 2 3 2+2*i 3+2*i];
                    [zh, Hs] = modelo_obs(mu(id));
                    S  = Hs*P(id,id)*Hs' + R;
                    nu = z - zh; nu(2) = wrapToPi(nu(2));
                    d2 = nu'/S*nu;
                    if d2 < d2best, d2best = d2; best = i; end
                end
            end

            if d2best < chi2_asoc
                [mu, P] = ekf_update(mu, P, best, z, R);
                n_asoc = n_asoc + 1;
                obs_ids(end+1) = best; %#ok<AGROW>
                if t(k+1) - lm_visto(best) > t_relazo
                    n_relazo = n_relazo + 1; hubo_relazo = true;
                end
                lm_visto(best) = t(k+1);
            elseif d2best > chi2_nuevo && deu_min > d_min_nuevo && nlm < Nlm_max
                [mu, P] = ekf_nuevo(mu, P, z, R);
                lm_visto(end+1) = t(k+1); %#ok<AGROW>
                n_nuevos = n_nuevos + 1;
                obs_ids(end+1) = nlm+1; %#ok<AGROW>
            else
                n_rech = n_rech + 1;                 % ambigua: se descarta
            end
        end
        if hubo_relazo, t_lc(end+1) = t(k+1); end %#ok<AGROW>

        %% --- Grilla: integrar el scan en keyframes (pose estimada) ---
        if isempty(pose_kf) || norm(mu(1:2)-pose_kf(1:2)) > nodo_dist || ...
           abs(wrapToPi(mu(3)-pose_kf(3))) > nodo_ang
            Lmap = integrar_scan(Lmap, map, mu(1:3), rg, lidar_angles, ...
                                 lidar_max_range, l_occ, l_free, l_max);
            pose_kf = mu(1:3); n_kf = n_kf + 1;
            mapa_pendiente = true;
        end
    end

    %% --- Registros ---
    x_slam(:,k+1) = mu(1:3);
    err_odo(k+1)  = hypot(odom_x(k+1)-x_real(k+1), odom_y(k+1)-y_real(k+1));
    err_slam(k+1) = hypot(mu(1)-x_real(k+1), mu(2)-y_real(k+1));
    sig3(k+1)     = 3*sqrt(P(1,1)+P(2,2));
    kf = k+1;

    %% --- Gráficos ---
    if mod(k, plot_every) == 0
        kk = k+1;
        set(hReal,'XData',x_real(1:kk),'YData',y_real(1:kk));
        set(hOdo, 'XData',odom_x(1:kk),'YData',odom_y(1:kk));
        set(hSlam,'XData',x_slam(1,1:kk),'YData',x_slam(2,1:kk));
        pr_v = [x_real(kk); y_real(kk); phi_real(kk)];
        rg_v = lidar_scan(pr_v, lidar_angles, lidar_max_range, segs);
        [hx, hy, mx, my, px, py] = rayos(pr_v, rg_v, lidar_angles, lidar_max_range, rayos_libres_cada);
        set(hRayHit, 'XData',hx,'YData',hy);
        set(hRayMiss,'XData',mx,'YData',my);
        set(hScan,   'XData',px,'YData',py);
        hit_v = rg_v < lidar_max_range;
        pe = transformar(mu(1:3), [rg_v(hit_v).*cos(lidar_angles(hit_v)); ...
                                   rg_v(hit_v).*sin(lidar_angles(hit_v))]);
        set(hScan2,'XData',pe(1,:),'YData',pe(2,:));

        nlm = (numel(mu)-3)/2;
        if nlm > 0
            LM = reshape(mu(4:end), 2, []);
            [ex, ey] = elipses_lm(mu, P, nlm, 3);
            set(hLm, 'XData',LM(1,:),'YData',LM(2,:)); set(hLm2,'XData',LM(1,:),'YData',LM(2,:));
            set(hEll,'XData',ex,'YData',ey);           set(hEll2,'XData',ex,'YData',ey);
            if ~isempty(obs_ids)
                m = numel(obs_ids);
                ox = [repmat(mu(1),1,m); LM(1,obs_ids); nan(1,m)];
                oy = [repmat(mu(2),1,m); LM(2,obs_ids); nan(1,m)];
                set(hObs,'XData',ox(:),'YData',oy(:));
            else
                set(hObs,'XData',nan,'YData',nan);
            end
        end
        [erx, ery] = elipse(mu(1:2), P(1:2,1:2), 3);
        set(hEllR,'XData',erx,'YData',ery);
        actualizar_tanque(hTank1, x_real(kk), y_real(kk), phi_real(kk), L, B_ext);

        if mapa_pendiente
            set(hMap,'CData', 1 - 1./(1+exp(-Lmap))); mapa_pendiente = false;
        end
        set(hSlam2,'XData',x_slam(1,1:kk),'YData',x_slam(2,1:kk));
        actualizar_tanque(hTank2, mu(1), mu(2), mu(3), L, B_ext);
        title(ax2, sprintf('Mapa construido  |  landmarks: %d  |  re-observaciones (lazo): %d', ...
              nlm, n_relazo));

        set(hEo,'XData',t(1:kk),'YData',err_odo(1:kk));
        set(hEs,'XData',t(1:kk),'YData',err_slam(1:kk));
        set(hSg,'XData',t(1:kk),'YData',sig3(1:kk));
        if ~isempty(t_lc)
            set(hEl,'XData',t_lc,'YData',interp1(t(1:kk), err_slam(1:kk), t_lc));
        end
        % Ritmo constante en pantalla: esperar hasta que el reloj alcance t(kk)
        if vel_anim > 0
            espera = t(kk)/vel_anim - toc(t_reloj);
            if espera > 0, pause(espera); end
        end
        drawnow limitrate;
    end
end

%% === RESULTADOS ===
idx = 2:kf;
nlm = (numel(mu)-3)/2;
LM  = reshape(mu(4:end), 2, []);
e_lm = zeros(1,nlm);
for i = 1:nlm
    e_lm(i) = min(hypot(C_true(1,:)-LM(1,i), C_true(2,:)-LM(2,i)));
end
fprintf('\n=== EKF-SLAM (landmarks = esquinas) ===\n');
fprintf('Tiempo de recorrido        : %.1f s\n', t(kf));
v_real = hypot(diff(x_real(1:kf)), diff(y_real(1:kf)))/dt;
fprintf('Velocidad real             : %.4f m/s  (desv. est. %.2e m/s)\n', mean(v_real), std(v_real));
fprintf('Esquinas detectadas        : %d\n', n_esq);
fprintf('  asociadas / nuevas / desc: %d / %d / %d\n', n_asoc, n_nuevos, n_rech);
fprintf('Landmarks en el estado     : %d  (esquinas reales: %d)\n', nlm, size(C_true,2));
fprintf('Tamaño del estado          : %d  (P de %dx%d)\n', numel(mu), numel(mu), numel(mu));
fprintf('Re-observaciones (lazo)    : %d\n', n_relazo);
fprintf('Keyframes en la grilla     : %d\n', n_kf);
fprintf('Error medio odometría      : %.4f m   (final %.4f m)\n', mean(err_odo(idx)), err_odo(kf));
fprintf('Error medio EKF-SLAM       : %.4f m   (final %.4f m)\n', mean(err_slam(idx)), err_slam(kf));
fprintf('%% tiempo dentro de 3sigma  : %.1f %%\n', 100*mean(err_slam(idx) <= sig3(idx)));
if nlm > 0
    fprintf('Error medio landmarks      : %.4f m   (espurios >0.3 m: %d)\n', ...
            mean(e_lm), nnz(e_lm > 0.3));
end

%% === FIGURA 3: MAPA FINAL vs ENTORNO REAL ===
fig3 = figure(3); set(fig3,'Position',[100,60,900,850]); hold on; axis equal; box on;
imagesc(xs_map, ys_map, 1 - 1./(1+exp(-Lmap)));
set(gca,'YDir','normal'); colormap(gray); caxis([0 1]);
for s = 1:size(segs,1)
    plot(segs(s,[1 3]), segs(s,[2 4]), 'r-', 'LineWidth',0.8, 'HandleVisibility','off');
end
plot(nan,nan,'r-','DisplayName','Obstáculos reales');
plot(x_real(1:kf), y_real(1:kf), 'k-', 'LineWidth',1.2, 'DisplayName','Trayectoria real');
plot(x_slam(1,1:kf), x_slam(2,1:kf), 'g-', 'LineWidth',1.2, 'DisplayName','EKF-SLAM');
if nlm > 0
    [ex, ey] = elipses_lm(mu, P, nlm, 3);
    plot(ex, ey, 'm-', 'LineWidth',0.7, 'HandleVisibility','off');
    plot(LM(1,:), LM(2,:), 'md', 'MarkerFaceColor','m', 'MarkerSize',5, 'DisplayName','Landmarks (3\sigma)');
end
xlim([-room-0.5 room+0.5]); ylim([-room-0.5 room+0.5]);
xlabel('X [m]'); ylabel('Y [m]');
title('Mapa final EKF-SLAM (grilla + landmarks) sobre el entorno real');
legend('Location','southoutside','Orientation','horizontal');
saveas(fig3, 'ekfslam_sin_paredes_mapa.png');
fprintf('\n[ok] Mapa guardado en: %s\n', fullfile(pwd,'ekfslam_sin_paredes_mapa.png'));

%% ======================= FUNCIONES =======================

function [zh, H] = modelo_obs(v)
    % v = [x y phi mx my] -> z = [rango; rumbo] y jacobiano 2x5
    dx = v(4)-v(1); dy = v(5)-v(2);
    q = dx^2 + dy^2; r = sqrt(q);
    zh = [r; wrapToPi(atan2(dy,dx) - v(3))];
    H  = [-dx/r, -dy/r,  0,  dx/r, dy/r;
           dy/q, -dx/q, -1, -dy/q, dx/q];
end

function [mu, P] = ekf_update(mu, P, i, z, R)
    % Corrección EKF con el landmark i (se usa solo el bloque disperso de H)
    id = [1 2 3 2+2*i 3+2*i];
    [zh, Hs] = modelo_obs(mu(id));
    PHt = P(:,id)*Hs';
    S   = Hs*PHt(id,:) + R;
    K   = PHt/S;
    nu  = z - zh; nu(2) = wrapToPi(nu(2));
    mu  = mu + K*nu; mu(3) = wrapToPi(mu(3));
    P   = P - K*S*K';
    P   = (P + P')/2;
end

function [mu, P] = ekf_nuevo(mu, P, z, R)
    % Inicializa un landmark a partir de [r; b] y aumenta el estado
    a  = mu(3) + z(2);
    m  = [mu(1) + z(1)*cos(a); mu(2) + z(1)*sin(a)];
    Gx = [1 0 -z(1)*sin(a); 0 1 z(1)*cos(a)];
    Gz = [cos(a) -z(1)*sin(a); sin(a) z(1)*cos(a)];
    Plx = Gx*P(1:3,:);
    Pll = Gx*P(1:3,1:3)*Gx' + Gz*R*Gz';
    mu = [mu; m];
    P  = [P, Plx'; Plx, Pll];
end

function Z = extraer_esquinas(rg, ang, rmax, ext)
    % Clusters -> split-and-merge -> rectas por mínimos cuadrados totales
    % -> intersección de caras adyacentes casi perpendiculares = esquina.
    Z = zeros(2,0);
    n = numel(rg);
    hit = rg < rmax;
    if ~any(hit), return; end
    i0 = find(~hit, 1); if isempty(i0), i0 = 1; end
    ord = [i0:n, 1:i0-1];                    % empezar en un rayo sin impacto
    rg = rg(ord); ang = ang(ord); hit = hit(ord);
    px = rg.*cos(ang); py = rg.*sin(ang);
    k = 1;
    while k <= n
        if ~hit(k), k = k+1; continue; end
        s = k;
        while k < n && hit(k+1) && hypot(px(k+1)-px(k), py(k+1)-py(k)) < ext.d_clust
            k = k+1;
        end
        e = k; k = k+1;
        if e-s+1 < 2*ext.min_pts+1, continue; end
        Pc = [px(s:e); py(s:e)];
        br = dividir(Pc, 1, size(Pc,2), ext.thr_split);
        nb = numel(br);
        for q = 2:nb-1
            a = br(q-1); m = br(q); b = br(q+1);
            A  = Pc(:, a+(q>2) : m-1);           % se excluyen los puntos de quiebre
            Bp = Pc(:, m+1 : b-(q<nb-1));
            if size(A,2) < ext.min_pts || size(Bp,2) < ext.min_pts, continue; end
            if norm(A(:,end)-A(:,1)) < ext.min_len || norm(Bp(:,end)-Bp(:,1)) < ext.min_len
                continue;
            end
            [c1, d1] = ajustar_recta(A); [c2, d2] = ajustar_recta(Bp);
            if abs(d1'*d2) > ext.max_cos, continue; end
            ts = [d1, -d2] \ (c2 - c1);
            pc = c1 + ts(1)*d1;
            if norm(pc - Pc(:,m)) > ext.max_gap, continue; end
            r = norm(pc);
            if r > ext.r_max, continue; end
            Z(:,end+1) = [r; atan2(pc(2), pc(1))]; %#ok<AGROW>
        end
    end
end

function br = dividir(P, i, j, thr)
    % Split recursivo (iterative end-point fit). Devuelve índices de quiebre.
    if j - i < 2, br = [i j]; return; end
    d = P(:,j) - P(:,i); Ld = norm(d);
    if Ld < 1e-9, br = [i j]; return; end
    nrm = [-d(2); d(1)]/Ld;
    dist = abs(nrm'*(P(:,i+1:j-1) - P(:,i)));
    [dm, kk] = max(dist);
    if dm > thr
        m = i + kk;
        b1 = dividir(P, i, m, thr); b2 = dividir(P, m, j, thr);
        br = [b1, b2(2:end)];
    else
        br = [i j];
    end
end

function [c, d] = ajustar_recta(Q)
    c = mean(Q, 2);
    [V, D] = eig(cov(Q'));
    [~, im] = max(diag(D));
    d = V(:, im);
end

function [ex, ey] = elipse(m, C, ns)
    th = linspace(0, 2*pi, 24);
    [V, D] = eig((C + C')/2);
    c = V*sqrt(max(D,0))*ns*[cos(th); sin(th)];
    ex = [m(1) + c(1,:), nan]'; ey = [m(2) + c(2,:), nan]';
end

function [ex, ey] = elipses_lm(mu, P, nlm, ns)
    ex = nan(25, nlm); ey = ex;
    for i = 1:nlm
        id = [2+2*i, 3+2*i];
        [ex(:,i), ey(:,i)] = elipse(mu(id), P(id,id), ns);
    end
    ex = ex(:); ey = ey(:);
end

function [V, w, wp_idx, fin] = control_waypoints(x, y, phi, WP, wp_idx, V_cte, umbral, k_w, w_max)
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

function s = rect_segs(o)
    x1 = o(1)-o(3)/2; x2 = o(1)+o(3)/2; y1 = o(2)-o(4)/2; y2 = o(2)+o(4)/2;
    s = [x1 y1 x2 y1; x2 y1 x2 y2; x2 y2 x1 y2; x1 y2 x1 y1];
end

function r = lidar_scan(pose, angles, rmax, segs)
    a = pose(3) + angles; dxr = cos(a); dyr = sin(a);
    r = rmax*ones(size(a));
    for w = 1:size(segs,1)
        ex = segs(w,3)-segs(w,1); ey = segs(w,4)-segs(w,2);
        vx = segs(w,1)-pose(1);   vy = segs(w,2)-pose(2);
        den = dxr*ey - dyr*ex;
        tt = (vx*ey - vy*ex)./den;
        uu = (vx*dyr - vy*dxr)./den;
        ok = abs(den) > 1e-12 & tt > 0 & uu >= 0 & uu <= 1 & tt < r;
        r(ok) = tt(ok);
    end
end

function dibujar_caja(ax, o, color)
    x1 = o(1)-o(3)/2; x2 = o(1)+o(3)/2; y1 = o(2)-o(4)/2; y2 = o(2)+o(4)/2;
    patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], color, 'EdgeColor','k', 'HandleVisibility','off');
end

function [hx, hy, mx, my, px, py] = rayos(pose, rg, angles, rmax, cada)
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

function R = rot(th)
    R = [cos(th) -sin(th); sin(th) cos(th)];
end

function p = transformar(pose, pts)
    p = rot(pose(3))*pts + pose(1:2);
end

function Lm = integrar_scan(Lm, map, pose, ranges, angles, rmax, l_occ, l_free, l_max)
    a = pose(3) + angles;
    s = (0:ceil(rmax/map.res)-1)' * map.res;
    libre = s < (ranges - map.res);
    Xs = pose(1) + s*cos(a); Ys = pose(2) + s*sin(a);
    lin = celdas(Xs(libre), Ys(libre), map);
    Lm(lin) = Lm(lin) + l_free;
    hit = ranges < rmax;
    lin = celdas(pose(1) + ranges(hit).*cos(a(hit)), pose(2) + ranges(hit).*sin(a(hit)), map);
    Lm(lin) = Lm(lin) + l_occ;
    Lm = min(max(Lm, -l_max), l_max);
end

function lin = celdas(x, y, map)
    ix = floor((x(:) - map.xmin)/map.res) + 1;
    iy = floor((y(:) - map.ymin)/map.res) + 1;
    ok = ix >= 1 & ix <= map.nx & iy >= 1 & iy <= map.ny;
    lin = unique(sub2ind([map.ny map.nx], iy(ok), ix(ok)));
end

function a = wrapToPi(a)
    a = mod(a+pi, 2*pi) - pi;
end

function h = crear_tanque(ax)
    h.trackL = patch(ax,0,0,[0.15 0.15 0.15],'EdgeColor','k','HandleVisibility','off');
    h.trackR = patch(ax,0,0,[0.15 0.15 0.15],'EdgeColor','k','HandleVisibility','off');
    h.treadL = plot(ax,nan,nan,'-','Color',[0.55 0.55 0.55],'HandleVisibility','off');
    h.treadR = plot(ax,nan,nan,'-','Color',[0.55 0.55 0.55],'HandleVisibility','off');
    h.body   = patch(ax,0,0,[0.2 0.4 0.8],'FaceAlpha',0.9,'EdgeColor','k','HandleVisibility','off');
    h.front  = patch(ax,0,0,[1 0.85 0.1],'EdgeColor','k','HandleVisibility','off');
end

function actualizar_tanque(h, x, y, phi, L, W)
    Rm = rot(phi);
    tw = 0.045; yc = W/2 - tw/2;
    trk = [-L/2 L/2 L/2 -L/2; -tw/2 -tw/2 tw/2 tw/2];
    tL = Rm*(trk + [0;  yc]); tR = Rm*(trk + [0; -yc]);
    wi = W/2 - tw;
    bd = Rm*[-0.42*L 0.42*L 0.42*L -0.42*L; -wi -wi wi wi];
    fr = Rm*[0.40*L 0.18*L 0.18*L; 0 0.6*wi -0.6*wi];
    xs = linspace(-L/2, L/2, 9);
    tx = [xs; xs; nan(1,9)]; ty = [-tw/2*ones(1,9); tw/2*ones(1,9); nan(1,9)];
    zL = Rm*[tx(:)'; ty(:)' + yc]; zR = Rm*[tx(:)'; ty(:)' - yc];
    set(h.trackL,'XData',tL(1,:)+x,'YData',tL(2,:)+y);
    set(h.trackR,'XData',tR(1,:)+x,'YData',tR(2,:)+y);
    set(h.treadL,'XData',zL(1,:)+x,'YData',zL(2,:)+y);
    set(h.treadR,'XData',zR(1,:)+x,'YData',zR(2,:)+y);
    set(h.body,  'XData',bd(1,:)+x,'YData',bd(2,:)+y);
    set(h.front, 'XData',fr(1,:)+x,'YData',fr(2,:)+y);
end