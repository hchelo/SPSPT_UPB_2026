%% tanque_sim_graphSLAM_sin_paredes.m
% SLAM basado en grafos "estilo Cartographer" con el robot oruga Hiwonder.
% Sin toolboxes: todo el método está implementado aquí para poder explicarlo.
%
% Idea del método (igual que Cartographer / slam_toolbox, simplificado):
%   1) NODOS     : cada vez que el robot avanza 0.3 m o gira 15° se crea un
%                  nodo del grafo con su pose y su scan LiDAR.
%   2) FRONT-END : el scan nuevo se alinea (ICP) contra un SUBMAPA formado
%                  por los scans de los últimos nodos -> arista secuencial.
%   3) CIERRE DE LAZO: si el robot pasa cerca de un nodo antiguo, se alinea
%                  el scan actual con el de ese nodo -> arista de lazo.
%   4) BACK-END  : al aparecer un lazo se optimiza el grafo completo
%                  (Gauss-Newton sobre todas las poses).
%   5) MAPA      : grilla de ocupación construida con los scans en las poses
%                  optimizadas (se reconstruye después de cada optimización).
%
% A diferencia de EKF-SLAM: no hay beacons ni asociación de datos; se usa
% la geometría completa del entorno.
%
% VERSIÓN SIN PAREDES: espacio abierto con obstáculos (cajas y columnas).
% Se dibujan los rayos del LiDAR: rojos los que impactan un obstáculo,
% claros los que no encuentran nada dentro del alcance.

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
plot_every = 10;

%% === ENCODERS Y DESLIZAMIENTO ===
sigma_encoder = 0.02;
io_mean = 0.08; ii_mean = 0.10; sigma_slip = 0.05;

%% === WAYPOINTS (ruta cuadrada del script de Python, escalada) ===
escala = 1/15;
WP = escala * [0 0; 75 0; 75 75; -75 75; -75 -75; -75 -75; 75 -75; -70 70];
umbral_wp = 5*escala;

%% === ENTORNO: espacio abierto con obstáculos (SIN paredes) ===
room = 7;                                      % solo define el área dibujada
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
         % obstáculos extra: sin paredes el ICP necesita más geometría
         0.0   3.5   0.5   0.5;
        -1.5   3.8   0.6   0.3;
         1.5  -3.5   0.5   0.5;
        -1.0  -1.8   0.4   0.4;
         4.0   1.5   0.4   0.4;
         6.5   6.5   0.8   0.8;
        -6.5   6.5   0.8   0.8;
         6.5  -6.5   0.8   0.8;
        -6.5  -6.5   0.8   0.8];
for i = 1:size(obst,1)
    segs = [segs; rect_segs(obst(i,:))]; %#ok<AGROW>
end

%% === LIDAR (tipo RPLIDAR) ===
lidar_max_range = 5.0;
lidar_n_rays    = 360;                         % 1°
lidar_angles    = linspace(-pi, pi, lidar_n_rays+1); lidar_angles(end) = [];
sigma_lidar     = 0.02;
rayos_libres_cada = 3;          % dibujar 1 de cada 3 rayos sin impacto

%% === PARÁMETROS DEL SLAM DE GRAFOS ===
nodo_dist  = 0.30;             % crear nodo cada 0.3 m ...
nodo_ang   = deg2rad(15);      % ... o cada 15°
K_submapa  = 10;               % nodos que forman el submapa local (1 = scan-to-scan)
icp_iter   = 30;
icp_dmax   = 0.30;             % distancia máx. de correspondencia (front-end) [m]
lc_dmax    = 0.50;             % idem para cierre de lazo [m]
lc_radio   = 1.5;              % buscar lazos con nodos a menos de 1.5 m
lc_min_sep = 40;               % ... que sean al menos 40 nodos más antiguos
lc_cada    = 2;                % intentar cierre de lazo cada 2 nodos
lc_rmse    = 0.06;             % calidad mínima del ICP para aceptar un lazo
lc_frac    = 0.70;
Om_scan = diag([1/0.03^2, 1/0.03^2, 1/deg2rad(1)^2]);   % arista scan matching
Om_odo  = diag([1/0.10^2, 1/0.10^2, 1/deg2rad(5)^2]);   % si falla el ICP
Om_lc   = diag([1/0.05^2, 1/0.05^2, 1/deg2rad(2)^2]);   % arista de lazo

%% === GRILLA DE OCUPACIÓN ===
map.res  = 0.05;
map.xmin = -room-0.5; map.ymin = -room-0.5;
map.nx   = round(2*(room+0.5)/map.res); map.ny = map.nx;
l_occ = 0.85; l_free = -0.4; l_max = 5;
Lmap = zeros(map.ny, map.nx);
xs_map = map.xmin + ((1:map.nx)-0.5)*map.res;
ys_map = map.ymin + ((1:map.ny)-0.5)*map.res;

%% === ALMACENAMIENTO DEL GRAFO ===
Nmax = 3000;
Xn = zeros(3,Nmax);            % poses de los nodos (estimadas)
Xn_true = zeros(3,Nmax);       % poses reales (solo para evaluar)
odo_n = zeros(3,Nmax);         % odometría en cada nodo
scans = cell(1,Nmax);          % puntos del scan en marco del robot
ranges_n = cell(1,Nmax);       % rangos crudos (para el mapa)
Ei = []; Ej = []; Ez = zeros(3,0); EOm = zeros(3,3,0); Etipo = [];
nn = 0; n_lc = 0; n_fallback = 0; n_opt = 0; t_lc = [];

%% === REGISTROS ===
x_real = zeros(1,N); y_real = zeros(1,N); phi_real = zeros(1,N);
odom_x = zeros(1,N); odom_y = zeros(1,N); odom_phi = zeros(1,N);
x_slam = zeros(3,N);
err_odo = zeros(1,N); err_slam = zeros(1,N);
scan_w = [nan; nan];
wp_idx = 1; kf = N;

%% === FIGURA 1 ===
fig1 = figure(1); set(fig1,'Position',[30,60,1500,720]);
% --- Izquierda: mundo real ---
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal'); box(ax1,'on');
hRayMiss = plot(ax1, nan,nan,'-', 'Color',[1 0.87 0.78], 'LineWidth',0.5, 'DisplayName','Rayo sin impacto');
hRayHit  = plot(ax1, nan,nan,'-', 'Color',[1 0.35 0.35], 'LineWidth',0.6, 'DisplayName','Rayo con impacto');
for i = 1:size(obst,1)
    dibujar_caja(ax1, obst(i,:), [0.35 0.35 0.35]);
end
plot(ax1, WP(:,1), WP(:,2), 's--', 'Color',[0.5 0 0.5], 'DisplayName','Waypoints');
hScan  = plot(ax1, nan,nan,'.', 'Color',[0.85 0 0], 'MarkerSize',7, 'DisplayName','Impactos LiDAR');
hReal  = plot(ax1, nan,nan,'k-', 'LineWidth',1.5, 'DisplayName','Real');
hOdo   = plot(ax1, nan,nan,'b--','LineWidth',1.2, 'DisplayName','Odometría');
hSlam  = plot(ax1, nan,nan,'g-', 'LineWidth',1.5, 'DisplayName','SLAM (grafo)');
hNodes = plot(ax1, nan,nan,'g.', 'MarkerSize',9, 'HandleVisibility','off');
hSub   = plot(ax1, nan,nan,'o', 'Color',[1 0.5 0], 'MarkerSize',5, 'DisplayName','Submapa activo');
hLC    = plot(ax1, nan,nan,'r-', 'LineWidth',1.2, 'DisplayName','Aristas de lazo');
hTank1 = crear_tanque(ax1);
xlim(ax1, [-room-0.5 room+0.5]); ylim(ax1, [-room-0.5 room+0.5]);
xlabel(ax1,'X [m]'); ylabel(ax1,'Y [m]'); title(ax1,'Mundo real, rayos LiDAR y grafo de poses');
legend(ax1,'Location','southoutside','Orientation','horizontal','NumColumns',4);

% --- Derecha: mapa construido ---
ax2 = subplot(1,2,2); hold(ax2,'on'); axis(ax2,'equal'); box(ax2,'on');
hMap = imagesc(ax2, xs_map, ys_map, 0.5*ones(map.ny,map.nx));
set(ax2,'YDir','normal'); colormap(ax2, gray); caxis(ax2,[0 1]);
hSlam2  = plot(ax2, nan,nan,'g-', 'LineWidth',1.3);
hNodes2 = plot(ax2, nan,nan,'g.', 'MarkerSize',8);
hLC2    = plot(ax2, nan,nan,'r-', 'LineWidth',1.2);
hScan2  = plot(ax2, nan,nan,'.', 'Color',[1 0.3 0.3], 'MarkerSize',6);  % scan en la pose estimada
hTank2  = crear_tanque(ax2);
xlim(ax2, [-room-0.5 room+0.5]); ylim(ax2, [-room-0.5 room+0.5]);
xlabel(ax2,'X [m]'); ylabel(ax2,'Y [m]');
title(ax2,'Mapa construido (grilla de ocupación)');

%% === FIGURA 2: ERRORES ===
fig2 = figure(2); set(fig2,'Position',[250,300,1000,380]);
ax3 = axes(fig2); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
hEo = plot(ax3, nan,nan,'b--','LineWidth',1.4,'DisplayName','Error odometría');
hEs = plot(ax3, nan,nan,'g-', 'LineWidth',1.6,'DisplayName','Error SLAM (grafo)');
hEl = plot(ax3, nan,nan,'rv', 'MarkerFaceColor','r','DisplayName','Cierre de lazo');
xlabel(ax3,'Tiempo [s]'); ylabel(ax3,'Error de posición [m]');
legend(ax3,'Location','northwest');

%% === SIMULACIÓN ===
for k = 1:N-1
    %% --- Control por waypoints (velocidad constante) ---
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

    %% --- Encoders -> odometría ---
    Vom = Vo + sigma_encoder*randn; Vim = Vi + sigma_encoder*randn;
    Venc = (Vom+Vim)/2; Wenc = (Vom-Vim)/B;
    odom_phi(k+1) = wrapToPi(odom_phi(k) + dt*Wenc);
    odom_x(k+1)   = odom_x(k) + dt*Venc*cos(odom_phi(k));
    odom_y(k+1)   = odom_y(k) + dt*Venc*sin(odom_phi(k));
    odo_now = [odom_x(k+1); odom_y(k+1); odom_phi(k+1)];

    %% --- ¿Crear un nodo nuevo? ---
    if nn == 0
        crear_nodo = true;
    else
        mov = relativa(odo_n(:,nn), odo_now);
        crear_nodo = norm(mov(1:2)) > nodo_dist || abs(mov(3)) > nodo_ang;
    end

    mapa_cambio = false;
    if crear_nodo
        pr = [x_real(k+1); y_real(k+1); phi_real(k+1)];

        % Scan LiDAR desde la pose real
        rg = lidar_scan(pr, lidar_angles, lidar_max_range, segs);
        hit = rg < lidar_max_range;
        rg(hit) = rg(hit) + sigma_lidar*randn(1, nnz(hit));
        pts = [rg(hit).*cos(lidar_angles(hit)); rg(hit).*sin(lidar_angles(hit))];
        scan_w = transformar(pr, pts);

        nn = nn + 1;
        scans{nn} = pts; ranges_n{nn} = rg; odo_n(:,nn) = odo_now; Xn_true(:,nn) = pr;
        nueva_lc = false;

        if nn == 1
            Xn(:,1) = [0; 0; 0];
        else
            %% ===== FRONT-END: scan-to-submap =====
            guess = componer(Xn(:,nn-1), relativa(odo_n(:,nn-1), odo_now));
            js = max(1, nn-K_submapa):nn-1;
            dst = [];
            for j = js
                dst = [dst, transformar(Xn(:,j), scans{j})]; %#ok<AGROW>
            end
            dst = dst(:,1:2:end);
            [Tm, rmse, frac] = icp2d(pts, dst, guess, icp_iter, icp_dmax);
            if rmse < 0.08 && frac > 0.5
                Xn(:,nn) = Tm; Om = Om_scan;
            else
                Xn(:,nn) = guess; Om = Om_odo; n_fallback = n_fallback + 1;
            end
            Ei(end+1) = nn-1; Ej(end+1) = nn; %#ok<AGROW>
            Ez(:,end+1) = relativa(Xn(:,nn-1), Xn(:,nn));
            EOm(:,:,end+1) = Om; Etipo(end+1) = 0; %#ok<AGROW>

            %% ===== CIERRE DE LAZO =====
            if nn > lc_min_sep && mod(nn, lc_cada) == 0
                cand_rng = 1:nn-lc_min_sep;
                d = hypot(Xn(1,cand_rng)-Xn(1,nn), Xn(2,cand_rng)-Xn(2,nn));
                [ds, ord] = sort(d);
                cands = ord(ds < lc_radio);
                cands = cands(1:min(2,end));
                for j = cands
                    guess_lc = relativa(Xn(:,j), Xn(:,nn));
                    [Tl, rmse, frac] = icp2d(pts, scans{j}, guess_lc, icp_iter, lc_dmax);
                    corr = relativa(guess_lc, Tl);
                    if rmse < lc_rmse && frac > lc_frac && ...
                       norm(corr(1:2)) < 0.5 && abs(corr(3)) < deg2rad(20)
                        Ei(end+1) = j; Ej(end+1) = nn; %#ok<AGROW>
                        Ez(:,end+1) = Tl; EOm(:,:,end+1) = Om_lc; Etipo(end+1) = 1; %#ok<AGROW>
                        n_lc = n_lc + 1; nueva_lc = true;
                    end
                end
            end

            %% ===== BACK-END: optimización del grafo =====
            if nueva_lc
                Xn(:,1:nn) = optimizar_grafo(Xn(:,1:nn), Ei, Ej, Ez, EOm, 10);
                n_opt = n_opt + 1;
                t_lc(end+1) = t(k+1); %#ok<AGROW>
                % Reconstruir el mapa con las poses corregidas
                Lmap = zeros(map.ny, map.nx);
                for j = 1:nn
                    Lmap = integrar_scan(Lmap, map, Xn(:,j), ranges_n{j}, lidar_angles, ...
                                         lidar_max_range, l_occ, l_free, l_max);
                end
            end
        end
        if ~nueva_lc
            Lmap = integrar_scan(Lmap, map, Xn(:,nn), rg, lidar_angles, ...
                                 lidar_max_range, l_occ, l_free, l_max);
        end
        mapa_cambio = true;
    end

    %% --- Pose estimada actual: último nodo + odometría desde ese nodo ---
    x_slam(:,k+1) = componer(Xn(:,nn), relativa(odo_n(:,nn), odo_now));
    err_odo(k+1)  = hypot(odom_x(k+1)-x_real(k+1), odom_y(k+1)-y_real(k+1));
    err_slam(k+1) = hypot(x_slam(1,k+1)-x_real(k+1), x_slam(2,k+1)-y_real(k+1));
    kf = k+1;

    %% --- Gráficos ---
    if mod(k, plot_every) == 0 || mapa_cambio
        kk = k+1;
        set(hReal,'XData',x_real(1:kk),'YData',y_real(1:kk));
        set(hOdo, 'XData',odom_x(1:kk),'YData',odom_y(1:kk));
        set(hSlam,'XData',x_slam(1,1:kk),'YData',x_slam(2,1:kk));
        % Rayos LiDAR en vivo desde la pose real (solo visualización)
        pr_v = [x_real(kk); y_real(kk); phi_real(kk)];
        rg_v = lidar_scan(pr_v, lidar_angles, lidar_max_range, segs);
        [hx, hy, mx, my, px, py] = rayos(pr_v, rg_v, lidar_angles, lidar_max_range, rayos_libres_cada);
        set(hRayHit, 'XData',hx,'YData',hy);
        set(hRayMiss,'XData',mx,'YData',my);
        set(hScan,   'XData',px,'YData',py);
        % Los mismos impactos vistos desde la pose ESTIMADA, sobre el mapa
        hit_v = rg_v < lidar_max_range;
        pe = transformar(x_slam(:,kk), [rg_v(hit_v).*cos(lidar_angles(hit_v)); ...
                                        rg_v(hit_v).*sin(lidar_angles(hit_v))]);
        set(hScan2,'XData',pe(1,:),'YData',pe(2,:));
        set(hNodes,'XData',Xn(1,1:nn),'YData',Xn(2,1:nn));
        js = max(1,nn-K_submapa+1):nn;
        set(hSub,'XData',Xn(1,js),'YData',Xn(2,js));
        [lx, ly] = lineas_lazo(Xn, Ei, Ej, Etipo);
        set(hLC,'XData',lx,'YData',ly);
        actualizar_tanque(hTank1, x_real(kk), y_real(kk), phi_real(kk), L, B_ext);

        if mapa_cambio
            set(hMap,'CData', 1 - 1./(1+exp(-Lmap)));
        end
        set(hSlam2, 'XData',x_slam(1,1:kk),'YData',x_slam(2,1:kk));
        set(hNodes2,'XData',Xn(1,1:nn),'YData',Xn(2,1:nn));
        set(hLC2,'XData',lx,'YData',ly);
        actualizar_tanque(hTank2, x_slam(1,kk), x_slam(2,kk), x_slam(3,kk), L, B_ext);
        title(ax2, sprintf('Mapa construido  |  nodos: %d  |  lazos: %d  |  optimizaciones: %d', ...
              nn, n_lc, n_opt));

        set(hEo,'XData',t(1:kk),'YData',err_odo(1:kk));
        set(hEs,'XData',t(1:kk),'YData',err_slam(1:kk));
        if ~isempty(t_lc)
            set(hEl,'XData',t_lc,'YData',interp1(t(1:kk), err_slam(1:kk), t_lc));
        end
        drawnow limitrate;
    end
end

%% === RESULTADOS ===
idx = 2:kf;
e_nodos = hypot(Xn(1,1:nn)-Xn_true(1,1:nn), Xn(2,1:nn)-Xn_true(2,1:nn));
fprintf('\n=== SLAM DE GRAFOS (estilo Cartographer) ===\n');
fprintf('Tiempo de recorrido      : %.1f s\n', t(kf));
fprintf('Nodos del grafo          : %d\n', nn);
fprintf('Aristas secuenciales     : %d  (ICP fallido -> odometría: %d)\n', nnz(Etipo==0), n_fallback);
fprintf('Aristas de cierre de lazo: %d\n', n_lc);
fprintf('Optimizaciones del grafo : %d\n', n_opt);
fprintf('Error medio odometría    : %.4f m   (final %.4f m)\n', mean(err_odo(idx)), err_odo(kf));
fprintf('Error medio SLAM         : %.4f m   (final %.4f m)\n', mean(err_slam(idx)), err_slam(kf));
fprintf('Error medio de los nodos : %.4f m\n', mean(e_nodos));

%% === FIGURA 3: MAPA FINAL vs ENTORNO REAL ===
fig3 = figure(3); set(fig3,'Position',[100,60,900,850]); hold on; axis equal; box on;
imagesc(xs_map, ys_map, 1 - 1./(1+exp(-Lmap)));
set(gca,'YDir','normal'); colormap(gray); caxis([0 1]);
for s = 1:size(segs,1)
    plot(segs(s,[1 3]), segs(s,[2 4]), 'r-', 'LineWidth',0.8, 'HandleVisibility','off');
end
plot(nan,nan,'r-','DisplayName','Obstáculos reales');
plot(Xn_true(1,1:nn), Xn_true(2,1:nn), 'k-', 'LineWidth',1.2, 'DisplayName','Trayectoria real');
plot(Xn(1,1:nn), Xn(2,1:nn), 'g-', 'LineWidth',1.2, 'DisplayName','Grafo optimizado');
[lx, ly] = lineas_lazo(Xn, Ei, Ej, Etipo);
plot(lx, ly, 'm-', 'LineWidth',1.2, 'DisplayName','Aristas de lazo');
xlim([-room-0.5 room+0.5]); ylim([-room-0.5 room+0.5]);
xlabel('X [m]'); ylabel('Y [m]');
title('Mapa final (grilla) sobre el entorno real');
legend('Location','southoutside','Orientation','horizontal');
saveas(fig3, 'graphslam_sin_paredes_mapa.png');
fprintf('\n[ok] Mapa guardado en: %s\n', fullfile(pwd,'graphslam_sin_paredes_mapa.png'));

%% ======================= FUNCIONES =======================

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
    % Rectángulo [cx cy ancho alto] -> 4 segmentos [x1 y1 x2 y2]
    x1 = o(1)-o(3)/2; x2 = o(1)+o(3)/2; y1 = o(2)-o(4)/2; y2 = o(2)+o(4)/2;
    s = [x1 y1 x2 y1; x2 y1 x2 y2; x2 y2 x1 y2; x1 y2 x1 y1];
end

function r = lidar_scan(pose, angles, rmax, segs)
    % Ray casting vectorizado contra segmentos
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

function R = rot(th)
    R = [cos(th) -sin(th); sin(th) cos(th)];
end

function c = componer(a, b)
    % a ? b
    c = [a(1:2) + rot(a(3))*b(1:2); wrapToPi(a(3)+b(3))];
end

function r = relativa(a, b)
    % a^-1 ? b : pose de b vista desde a
    r = [rot(a(3))'*(b(1:2)-a(1:2)); wrapToPi(b(3)-a(3))];
end

function p = transformar(pose, pts)
    p = rot(pose(3))*pts + pose(1:2);
end

function [T, rmse, frac] = icp2d(src, dst, T, n_iter, d_max)
    % ICP punto a punto. T = pose que lleva src (marco local) al marco de dst.
    src = src(:,1:2:end);
    for it = 1:n_iter
        p = rot(T(3))*src + T(1:2);
        D = (p(1,:)' - dst(1,:)).^2 + (p(2,:)' - dst(2,:)).^2;
        [dmin, jn] = min(D, [], 2);
        ok = dmin < d_max^2;
        if nnz(ok) < 10, break; end
        A = p(:,ok); Bq = dst(:, jn(ok));
        ma = mean(A,2); mb = mean(Bq,2);
        [U, ~, V] = svd((A-ma)*(Bq-mb)');
        Rd = V*U';
        if det(Rd) < 0, V(:,2) = -V(:,2); Rd = V*U'; end
        td = mb - Rd*ma;
        dth = atan2(Rd(2,1), Rd(1,1));
        T = [Rd*T(1:2) + td; wrapToPi(T(3) + dth)];
        if norm(td) < 1e-5 && abs(dth) < 1e-6, break; end
    end
    p = rot(T(3))*src + T(1:2);
    dmin = min((p(1,:)' - dst(1,:)).^2 + (p(2,:)' - dst(2,:)).^2, [], 2);
    ok = dmin < d_max^2;
    frac = mean(ok);
    if any(ok), rmse = sqrt(mean(dmin(ok))); else, rmse = inf; end
end

function X = optimizar_grafo(X, Ei, Ej, Ez, EOm, n_iter)
    % Gauss-Newton sobre el grafo de poses (el primer nodo queda fijo)
    n = size(X,2); m = numel(Ei);
    for it = 1:n_iter
        I = zeros(36*m,1); J = I; Vv = I; c = 0;
        b = zeros(3*n,1);
        for e = 1:m
            i = Ei(e); j = Ej(e); z = Ez(:,e); Om = EOm(:,:,e);
            xi = X(:,i); xj = X(:,j);
            ci = cos(xi(3)); si = sin(xi(3));
            RiT  = [ci si; -si ci];
            dRiT = [-si ci; -ci -si];
            RzT  = rot(z(3))';
            dtr  = xj(1:2) - xi(1:2);
            err  = [RzT*(RiT*dtr - z(1:2)); wrapToPi(xj(3) - xi(3) - z(3))];
            A  = [-RzT*RiT, RzT*dRiT*dtr; 0 0 -1];
            Bm = [ RzT*RiT, [0;0];        0 0  1];
            ii = 3*i-2:3*i; jj = 3*j-2:3*j;
            blq = {ii,ii,A'*Om*A; ii,jj,A'*Om*Bm; jj,ii,Bm'*Om*A; jj,jj,Bm'*Om*Bm};
            for q = 1:4
                [rr, cc] = ndgrid(blq{q,1}, blq{q,2});
                I(c+1:c+9) = rr(:); J(c+1:c+9) = cc(:); Vv(c+1:c+9) = blq{q,3}(:);
                c = c + 9;
            end
            b(ii) = b(ii) + A'*Om*err;
            b(jj) = b(jj) + Bm'*Om*err;
        end
        H = sparse(I(1:c), J(1:c), Vv(1:c), 3*n, 3*n);
        H(1:3,1:3) = H(1:3,1:3) + 1e6*eye(3);          % anclar el nodo 1
        dx = -(H \ b);
        X = X + reshape(dx, 3, n);
        X(3,:) = wrapToPi(X(3,:));
        if norm(dx) < 1e-5, break; end
    end
end

function Lm = integrar_scan(Lm, map, pose, ranges, angles, rmax, l_occ, l_free, l_max)
    % Actualización log-odds: celdas libres a lo largo del rayo, ocupada al final
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

function [lx, ly] = lineas_lazo(Xn, Ei, Ej, Etipo)
    e = find(Etipo == 1);
    if isempty(e), lx = nan; ly = nan; return; end
    lx = [Xn(1,Ei(e)); Xn(1,Ej(e)); nan(1,numel(e))];
    ly = [Xn(2,Ei(e)); Xn(2,Ej(e)); nan(1,numel(e))];
    lx = lx(:); ly = ly(:);
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