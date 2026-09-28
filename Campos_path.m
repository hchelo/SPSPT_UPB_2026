%% ================================================================
%  PATH PLANNING — CAMPOS POTENCIALES + ROBOT KHEPERA (por pulsos)
%  Traducción desde Python a MATLAB
%
%  Convención del mapa (Mapa.txt, en la misma carpeta que este .m):
%    0  = libre
%    1  = obstáculo
%    2  = inicio del robot
%   -1  = meta
%
%  Requiere MATLAB R2016b o superior (funciones locales en scripts).
%  No requiere toolboxes (la dilatación se hace con conv2).
%
%  Indexación: M(x+1, y+1)  -> x = columna del archivo, y = fila.
%  Las coordenadas en celdas se mantienen en base 0 (como en Python)
%  para que las tablas de waypoints den exactamente los mismos números.
% ================================================================
clear; clc; close all;

BG = [0.05 0.05 0.05];   % fondo oscuro (#0d0d0d)

%% ================================================================
%  LEER MAPA DESDE ARCHIVO  Mapa.txt
% ================================================================
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta), carpeta = pwd; end
ruta_mapa = fullfile(carpeta, 'Mapa.txt');

M_raw = load(ruta_mapa, '-ascii');     % acepta tabulador o espacios
[FILAS, COLS] = size(M_raw);
fprintf('  Archivo leido: %s\n', ruta_mapa);
fprintf('  Filas=%d  Cols=%d\n', FILAS, COLS);

%% ================================================================
%  PARSEAR MAPA: inicio, meta y obstáculos
% ================================================================
M     = M_raw.';          % M(x+1, y+1)
ANCHO = COLS;             % eje X = columnas
ALTO  = FILAS;            % eje Y = filas

% find sobre la transpuesta -> mismo orden (fila-mayor) que np.argwhere
[ci, ri] = find(M_raw.' ==  2, 1);
[cg, rg] = find(M_raw.' == -1, 1);
if isempty(ci), error('No se encontro el inicio (valor 2) en el mapa'); end
if isempty(cg), error('No se encontro la meta (valor -1) en el mapa');  end

xi = ci - 1;  yi = ri - 1;     % coordenadas en celdas (base 0)
xg = cg - 1;  yg = rg - 1;

M(M ==  2) = 0;
M(M == -1) = 0;

fprintf('%s\n', repmat('=',1,55));
fprintf('   MAPA CARGADO\n');
fprintf('%s\n', repmat('=',1,55));
fprintf('  Tamano mapa    : %d x %d celdas (= 100x100 cm, 2cm/celda)\n', ANCHO, ALTO);
fprintf('  Inicio (xi,yi) : (%d, %d)\n', xi, yi);
fprintf('  Meta   (xg,yg) : (%d, %d)\n', xg, yg);
fprintf('  Obstaculos     : %d celdas\n', nnz(M == 1));
fprintf('%s\n', repmat('=',1,55));

%% ================================================================
%  DILATAR OBSTÁCULOS (margen de seguridad)
%  Equivale a skimage ball(2)[2] -> disco de radio 2
% ================================================================
RADIO_DILATE = 2;
[a, b]  = meshgrid(-RADIO_DILATE:RADIO_DILATE);
se      = double(a.^2 + b.^2 <= RADIO_DILATE^2);
M_dilat = double(conv2(double(M == 1), se, 'same') > 0);

% Proteger inicio y meta
M_dilat(xi+1, yi+1) = 0;
M_dilat(xg+1, yg+1) = 0;

%% ================================================================
%  CAMPOS POTENCIALES
%  UA = 1/2 * K * dist_a_meta
%  UR = 10*K + 1/dist_al_inicio   (solo en obstáculos)
% ================================================================
K = 3.0;
[XX, YY] = ndgrid(0:ANCHO-1, 0:ALTO-1);

dist_meta   = hypot(XX - xg, YY - yg);  dist_meta(dist_meta == 0)     = 1e-6;
dist_inicio = hypot(XX - xi, YY - yi);  dist_inicio(dist_inicio == 0) = 1e-6;

UA  = 0.5 * K * dist_meta;
UR  = zeros(size(UA));
obs = (M_dilat == 1);
UR(obs) = 10*K + 1 ./ dist_inicio(obs);
U = UA + UR;

%% ================================================================
%  DESCENSO DE GRADIENTE (8 vecinos)
% ================================================================
VECINOS = [-1 -1; -1 0; -1 1;
            0 -1;        0 1;
            1 -1;  1 0;  1 1];
MAX_IT = 50000;

path_x = zeros(MAX_IT+1, 1);  path_y = zeros(MAX_IT+1, 1);
path_x(1) = xi;  path_y(1) = yi;  n = 1;
cx = xi;  cy = yi;

for it = 1:MAX_IT
    if cx == xg && cy == yg, break; end
    mejor = inf;  nx = cx;  ny = cy;
    for v = 1:8
        ni = cx + VECINOS(v,1);
        nj = cy + VECINOS(v,2);
        if ni >= 0 && ni < ANCHO && nj >= 0 && nj < ALTO
            if U(ni+1, nj+1) < mejor
                mejor = U(ni+1, nj+1);
                nx = ni;  ny = nj;
            end
        end
    end
    if nx == cx && ny == cy
        fprintf('Minimo local en (%d,%d)\n', cx, cy);
        break;
    end
    cx = nx;  cy = ny;
    n = n + 1;
    path_x(n) = cx;  path_y(n) = cy;
end
path_x = path_x(1:n);
path_y = path_y(1:n);
fprintf('\nPath encontrado: %d puntos\n', n);

%% ================================================================
%  DETECTAR ESQUINAS (waypoints): cambio de dirección >= 20°
% ================================================================
UMBRAL_ANGULO = deg2rad(20);
wp_idx = 1;
for k = 2:n-1
    va = [path_x(k)   - path_x(k-1), path_y(k)   - path_y(k-1)];
    vs = [path_x(k+1) - path_x(k),   path_y(k+1) - path_y(k)];
    na = norm(va);  ns = norm(vs);
    if na < 1e-9 || ns < 1e-9, continue; end
    c = max(-1, min(1, dot(va, vs) / (na*ns)));
    if acos(c) >= UMBRAL_ANGULO
        wp_idx(end+1) = k; %#ok<SAGROW>
    end
end
wp_idx(end+1) = n;

% Filtrar waypoints muy cercanos (distancia mínima 3 celdas)
wp_f = wp_idx(1);
for idx = wp_idx(2:end)
    p = wp_f(end);
    if hypot(path_x(idx) - path_x(p), path_y(idx) - path_y(p)) >= 3
        wp_f(end+1) = idx; %#ok<SAGROW>
    end
end

WP_x = path_x(wp_f);
WP_y = path_y(wp_f);
NWP  = numel(WP_x);

% Orientación de cada waypoint = dirección hacia el siguiente
TH = zeros(NWP, 1);
for k = 1:NWP-1
    TH(k) = atan2(WP_y(k+1) - WP_y(k), WP_x(k+1) - WP_x(k));
end
TH(NWP) = atan2(WP_y(end) - WP_y(end-1), WP_x(end) - WP_x(end-1));

%% ================================================================
%  REPORTE EN CONSOLA
% ================================================================
fprintf('Waypoints detectados: %d\n\n', NWP);
fprintf('  %2s  %6s  %6s  %6s  %6s  %7s\n', '#', 'X(cel)', 'Y(cel)', 'X(cm)', 'Y(cm)', 'theta');
fprintf('  %s\n', repmat('-', 1, 45));
for k = 1:NWP
    lab = '';
    if k == 1,   lab = '  <- INICIO'; end
    if k == NWP, lab = '  <- META';   end
    fprintf('  %2d  %6d  %6d  %6.0f  %6.0f  %+6.1f deg%s\n', ...
        k-1, WP_x(k), WP_y(k), 2*WP_x(k), 2*WP_y(k), rad2deg(TH(k)), lab);
end

%% ================================================================
%  FIGURA 1 — VISUALIZACIÓN (6 paneles)
% ================================================================
fig1 = figure('Name', 'Campos Potenciales', 'Color', BG, 'Position', [40 60 1600 900]);
xs = [0.5, ANCHO-0.5];   ys = [0.5, ALTO-0.5];   % centros de celda (extent [0 W 0 H])

cmBlues   = mapaLineal([0.97 0.98 1.00], [0.03 0.19 0.42]);
cmOranges = mapaLineal([1.00 0.96 0.92], [0.50 0.15 0.01]);

% --- 1. Mapa original
ax = subplot(2,3,1); hold(ax, 'on');
imagesc(ax, xs, ys, M.');  colormap(ax, cmBlues);  caxis(ax, [0 1.5]);
h1 = plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
h2 = plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, 'Mapa Original', BG);
leyenda(ax, [h1 h2], {'Inicio', 'Meta'});

% --- 2. Mapa dilatado
ax = subplot(2,3,2); hold(ax, 'on');
imagesc(ax, xs, ys, M_dilat.');  colormap(ax, cmOranges);  caxis(ax, [0 1.5]);
plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, sprintf('Obstaculos Dilatados (r=%d celdas)', RADIO_DILATE), BG);

% --- 3. Potencial de atracción
ax = subplot(2,3,3); hold(ax, 'on');
imagesc(ax, xs, ys, UA.');  colormap(ax, parula);
cb = colorbar(ax); cb.Color = 'w';
plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, 'Potencial Atraccion U_A', BG);

% --- 4. Potencial de repulsión
ax = subplot(2,3,4); hold(ax, 'on');
imagesc(ax, xs, ys, min(UR, 50).');  colormap(ax, hot);
cb = colorbar(ax); cb.Color = 'w';
plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, 'Potencial Repulsion U_R', BG);

% --- 5. Potencial total + gradiente
ax = subplot(2,3,5); hold(ax, 'on');
imagesc(ax, xs, ys, min(U, 150).');  colormap(ax, jet);
cb = colorbar(ax); cb.Color = 'w';
[gy, gx] = gradient(-U);          % dim1 = x, dim2 = y
st = 4;  ix = 1:st:ANCHO;  iy = 1:st:ALTO;
quiver(ax, XX(ix,iy), YY(ix,iy), gx(ix,iy), gy(ix,iy), 'Color', [0.75 0.75 0.75]);
estilo(ax, 'Potencial Total + Gradiente', BG);

% --- 6. Ruta + waypoints
ax = subplot(2,3,6); hold(ax, 'on');
capa(ax, M_dilat == 1, [0.50 0.50 0.50], 0.4, 1);
capa(ax, M == 1,       [0.25 0.45 0.80], 0.8, 1);
hR = plot(ax, path_x, path_y, 'Color', [0 1 0], 'LineWidth', 1.5);
hW = scatter(ax, WP_x, WP_y, 90, 'y', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
for k = 1:NWP
    text(ax, WP_x(k)+0.5, WP_y(k)+0.5, num2str(k-1), 'Color', 'y', 'FontSize', 7);
    if k < NWP
        dx_ = WP_x(k+1) - WP_x(k);  dy_ = WP_y(k+1) - WP_y(k);
        quiver(ax, WP_x(k)+0.4*dx_, WP_y(k)+0.4*dy_, 0.2*dx_, 0.2*dy_, 0, ...
               'Color', [1 0.65 0], 'LineWidth', 1.5, 'MaxHeadSize', 4);
    end
end
hI = plot(ax, xi, yi, 'go', 'MarkerSize', 11, 'MarkerFaceColor', 'g');
hM = plot(ax, xg, yg, 'r*', 'MarkerSize', 13);
estilo(ax, 'Ruta + Waypoints (Esquinas)', BG);
xlim(ax, [0 ANCHO]);  ylim(ax, [0 ALTO]);  grid(ax, 'on');
leyenda(ax, [hR hW hI hM], {'Ruta', sprintf('%d Waypoints', NWP), 'Inicio', 'Meta'}, 'northwest');

% Tabla de waypoints al costado del panel 6
pos = ax.Position;
ax.Position = [pos(1)-0.04, pos(2), pos(3)*0.72, pos(4)];
lineas = {'Waypoints Khepera:', '', ...
          sprintf('  %2s  %4s  %4s  %5s  %5s  %7s', '#', 'X', 'Y', 'Xcm', 'Ycm', 'th'), ...
          ['  ' repmat('-', 1, 35)]};
for k = 1:NWP
    lineas{end+1} = sprintf('  %2d  %4d  %4d  %5.0f  %5.0f  %+6.1f', ...
        k-1, WP_x(k), WP_y(k), 2*WP_x(k), 2*WP_y(k), rad2deg(TH(k))); %#ok<SAGROW>
end
text(ax, 1.04, 0.99, lineas, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'FixedWidth', 'FontSize', 7, 'Color', 'w', 'Interpreter', 'none', ...
     'BackgroundColor', [0.10 0.10 0.18], 'EdgeColor', [0.27 0.27 0.67]);

sgtitle(fig1, 'Path Planning — Campos Potenciales  |  Mapa 100x100 cm  (2 cm/celda)', ...
        'Color', 'w', 'FontSize', 13);

%% ================================================================
%  FIGURA 2 — ROBOT KHEPERA NAVEGANDO POR LOS WAYPOINTS
%  Cinemática inversa por pulsos:
%    1. ORIENTAR   -> girar en sitio hasta apuntar al waypoint
%    2. DESPLAZAR  -> avanzar recto contando pulsos
%    3. REORIENTAR -> girar hasta la orientación final del waypoint
% ================================================================

% --- Parámetros Khepera
P.r    = 0.008;                 % radio rueda [m]
P.l    = 0.054;                 % distancia entre ruedas [m]
P.dt   = 0.05;                  % paso de tiempo [s]
P.PPR  = 600;                   % pulsos por revolución
P.MMP  = (2*pi*P.r) / P.PPR;    % metros por pulso
P.VLIN = 0.03;                  % m/s
P.VROT = 0.25;                  % rad/s

R_ROBOT_K = 0.0275;             % radio robot [m]  (diámetro 55 mm)
R_RUEDA_K = 0.005;              % grosor visual rueda
L_RUEDA_K = 0.010;              % largo visual rueda
CEL_A_M   = 0.02;               % 1 celda = 2 cm

% --- Waypoints en metros [x y theta]
wp_metros = [WP_x*CEL_A_M, WP_y*CEL_A_M, TH];

% --- Simular trayectoria completa
fprintf('\nSimulando Khepera en waypoints...\n');
enc        = [0 0];             % [EncL EncR]
reg_pulsos = [];
traj_k     = zeros(0, 3);
fase_k     = false(0, 1);       % true = girando, false = avanzando
x_k        = wp_metros(1, :);

for i = 1:NWP-1
    [tramo, x_k, enc, rg, fs] = ir_wp_k(x_k, wp_metros(i+1, :), i, P, enc);
    traj_k     = [traj_k; tramo];          %#ok<AGROW>
    fase_k     = [fase_k; fs];             %#ok<AGROW>
    reg_pulsos = [reg_pulsos, rg];         %#ok<AGROW>
end

COL_G = [1.00 0.843 0.000];   % amarillo = girando
COL_A = [0.00 0.902 0.463];   % verde    = avanzando
colores = double(fase_k).*COL_G + double(~fase_k).*COL_A;

% --- Reporte de pulsos
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('   REPORTE DE PULSOS — KHEPERA\n');
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  %-28s %10s %10s\n', 'Operacion', 'Delta', 'Pulsos');
fprintf('  %s\n', repmat('-', 1, 50));
tot_g = 0;  tot_a = 0;
for k = 1:numel(reg_pulsos)
    r = reg_pulsos(k);
    if r.tipo == 'G'
        fprintf('  %-28s %+9.1fdeg %9.1f\n', r.op, r.val, r.p);
        tot_g = tot_g + r.p;
    else
        fprintf('  %-28s %9.2fcm  %9.1f\n', r.op, r.val*100, r.p);
        tot_a = tot_a + r.p;
    end
end
fprintf('%s\n', repmat('=', 1, 60));
fprintf('  %-28s %10s  %9.1f\n', 'TOTAL GIROS',     '', tot_g);
fprintf('  %-28s %10s  %9.1f\n', 'TOTAL AVANCES',   '', tot_a);
fprintf('  %-28s %10s  %9.1f\n', 'TOTAL GENERAL L', '', enc(1));
fprintf('  %-28s %10s  %9.1f\n', 'TOTAL GENERAL R', '', enc(2));
fprintf('%s\n', repmat('=', 1, 60));

% --- Figura 2: escena
fig2 = figure('Name', 'Khepera', 'Color', BG, 'Position', [120 60 850 900]);
ax2  = axes(fig2);  hold(ax2, 'on');
set(ax2, 'Color', [0.03 0.03 0.03], 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
    'GridColor', 'w', 'GridAlpha', 0.15);

capa(ax2, M_dilat == 1, [0.50 0.50 0.50], 0.35, CEL_A_M);
capa(ax2, M == 1,       [0.25 0.45 0.80], 0.70, CEL_A_M);

hRuta = plot(ax2, path_x*CEL_A_M, path_y*CEL_A_M, 'Color', [0 1 0 0.35], 'LineWidth', 1);
for k = 1:NWP
    plot(ax2, wp_metros(k,1), wp_metros(k,2), 's', 'Color', 'y', 'MarkerSize', 7);
    text(ax2, wp_metros(k,1)+0.002, wp_metros(k,2)+0.002, num2str(k-1), 'Color', 'y', 'FontSize', 8);
end
hIni  = plot(ax2, wp_metros(1,1),   wp_metros(1,2),   'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
hMeta = plot(ax2, wp_metros(end,1), wp_metros(end,2), 'r*', 'MarkerSize', 14);

axis(ax2, 'equal');
xlim(ax2, [0 ANCHO*CEL_A_M]);  ylim(ax2, [0 ALTO*CEL_A_M]);
xlabel(ax2, 'x [m]');  ylabel(ax2, 'y [m]');  grid(ax2, 'on');
leyenda(ax2, [hRuta hIni hMeta], {'Ruta potencial', 'Inicio', 'Meta'}, 'northeast');
sgtitle(fig2, {'Robot Khepera — Cinematica Inversa por Pulsos', ...
               'Amarillo: Orientar   |   Verde: Desplazar   (3 pasos por waypoint)'}, ...
        'Color', 'w', 'FontSize', 11);

% --- Objetos gráficos del robot (se crean una vez y se actualizan)
hTraj  = scatter(ax2, NaN, NaN, 5, [0 0 0], 'filled');
ang    = linspace(0, 2*pi, 40);
hBody  = patch(ax2, NaN, NaN, COL_G, 'EdgeColor', 'none', 'FaceAlpha', 0.92);
hDir   = plot(ax2, NaN, NaN, 'w', 'LineWidth', 2);
hRueda = gobjects(2, 1);
for j = 1:2
    hRueda(j) = patch(ax2, NaN, NaN, [0.2 0.2 0.2], 'EdgeColor', [0.67 0.67 0.67], 'LineWidth', 0.8);
end
hTit = title(ax2, '', 'Color', 'w', 'FontSize', 8, 'Interpreter', 'none');

% --- Animación
N      = size(traj_k, 1);
step2  = max(1, floor(N/600));
frames = unique([1:step2:N, N]);
lados  = [+1, -1];

for fi = frames
    if ~isvalid(fig2), break; end
    x  = traj_k(fi, :);
    cx = x(1);  cy = x(2);  th = x(3);
    col = colores(fi, :);

    % Cuerpo
    set(hBody, 'XData', cx + R_ROBOT_K*cos(ang), 'YData', cy + R_ROBOT_K*sin(ang), 'FaceColor', col);
    % Flecha de dirección
    set(hDir, 'XData', [cx, cx + 1.7*R_ROBOT_K*cos(th)], ...
              'YData', [cy, cy + 1.7*R_ROBOT_K*sin(th)]);
    % Ruedas
    perp = [-sin(th), cos(th)];
    fwd  = [ cos(th), sin(th)];
    for j = 1:2
        wc = [cx, cy] + lados(j)*(P.l/2)*perp;
        esq = [wc + L_RUEDA_K*fwd + R_RUEDA_K*perp;
               wc - L_RUEDA_K*fwd + R_RUEDA_K*perp;
               wc - L_RUEDA_K*fwd - R_RUEDA_K*perp;
               wc + L_RUEDA_K*fwd - R_RUEDA_K*perp];
        set(hRueda(j), 'XData', esq(:,1), 'YData', esq(:,2));
    end
    % Estela coloreada por fase
    if fi > 1
        set(hTraj, 'XData', traj_k(1:fi-1,1), 'YData', traj_k(1:fi-1,2), ...
                   'CData', colores(1:fi-1,:));
    end
    % Título
    if fase_k(fi), txtF = 'Orientando'; else, txtF = 'Desplazando'; end
    hTit.String = sprintf('%s  |  t=%.1fs  |  pos=(%.1f,%.1f)cm  th=%.0f deg  |  EncL=%.0fp  EncR=%.0fp', ...
        txtF, (fi-1)*P.dt, cx*100, cy*100, rad2deg(th), enc(1), enc(2));

    drawnow;
    pause(0.03);
end

%% ================================================================
%  FUNCIONES LOCALES
% ================================================================

% --- Cinemática directa (diferencial)
function xn = cinem_dir(x, vL, vR, P)
    v  = P.r * (vR + vL) / 2;
    w  = P.r * (vR - vL) / P.l;
    xn = [x(1) + v*cos(x(3))*P.dt, ...
          x(2) + v*sin(x(3))*P.dt, ...
          x(3) + w*P.dt];
    xn(3) = atan2(sin(xn(3)), cos(xn(3)));
end

% --- Conversiones a pulsos de encoder
function p = m2p(metros, P)
    p = abs(metros) / P.MMP;
end

function p = a2p(rad, P)
    p = abs((P.l/2) * rad) / P.MMP;
end

% --- PASO 1 / 3: ORIENTAR (giro en sitio)
function [traj, x, pt, enc, reg] = girar_k(x, theta_dest, tag, P, enc)
    traj = x;  reg = [];
    delta = atan2(sin(theta_dest - x(3)), cos(theta_dest - x(3)));
    if abs(delta) < deg2rad(0.1)
        x(3) = theta_dest;  pt = 0;
        return;
    end
    pt = a2p(delta, P);  pe = 0;  s = sign(delta);
    wr = P.VROT / P.r;
    pp = (wr * P.r * P.dt) / P.MMP;
    while pe + pp <= pt
        enc = enc + [-s*pp, s*pp];  pe = pe + pp;
        x = cinem_dir(x, -s*wr, s*wr, P);
        traj(end+1, :) = x; %#ok<AGROW>
    end
    res = pt - pe;
    if res > 1e-6
        f = res / pp;
        enc = enc + [-s*res, s*res];
        x = cinem_dir(x, -s*wr*f, s*wr*f, P);
        traj(end+1, :) = x;
    end
    x(3) = theta_dest;
    traj(end+1, :) = x;
    reg = struct('op', ['GIRO ' tag], 'tipo', 'G', 'val', rad2deg(delta), 'p', pt);
end

% --- PASO 2: DESPLAZAR (avance recto)
function [traj, x, pt, enc, reg] = avanzar_k(x, dist, tag, P, enc)
    traj = x;
    th = x(3);  x0 = x(1);  y0 = x(2);
    pt = m2p(dist, P);  pe = 0;
    wr = P.VLIN / P.r;
    pp = (wr * P.r * P.dt) / P.MMP;
    while pe + pp <= pt
        enc = enc + [pp, pp];  pe = pe + pp;
        x = cinem_dir(x, wr, wr, P);  x(3) = th;
        traj(end+1, :) = x; %#ok<AGROW>
    end
    res = pt - pe;
    if res > 1e-6
        f = res / pp;
        enc = enc + [res, res];
        x = cinem_dir(x, wr*f, wr*f, P);  x(3) = th;
        traj(end+1, :) = x;
    end
    x(1) = x0 + dist*cos(th);
    x(2) = y0 + dist*sin(th);
    x(3) = th;
    traj(end+1, :) = x;
    reg = struct('op', ['AVANCE ' tag], 'tipo', 'A', 'val', dist, 'p', pt);
end

% --- Secuencia completa hacia un waypoint
function [traj, x, enc, regs, fase] = ir_wp_k(x, wp, idx, P, enc)
    gp  = wp(1:2);  thf = wp(3);
    d   = gp - x(1:2);
    dist     = hypot(d(1), d(2));
    th_hacia = atan2(d(2), d(1));
    [t1, x, ~, enc, r1] = girar_k  (x, th_hacia, sprintf('WP%d orient', idx), P, enc);
    [t2, x, ~, enc, r2] = avanzar_k(x, dist,     sprintf('WP%d %.1fcm', idx, dist*100), P, enc);
    x(1:2) = gp;
    [t3, x, ~, enc, r3] = girar_k  (x, thf,      sprintf('WP%d reorient', idx), P, enc);
    traj = [t1; t2(2:end, :); t3(2:end, :)];
    fase = [true(size(t1,1), 1); false(size(t2,1)-1, 1); true(size(t3,1)-1, 1)];
    regs = [r1, r2, r3];
end

% --- Utilidades gráficas
function cm = mapaLineal(c0, c1, n)
    if nargin < 3, n = 256; end
    t  = linspace(0, 1, n).';
    cm = (1 - t)*c0 + t*c1;
end

function capa(ax, mask, color, alfa, s)
    % Dibuja una máscara (indexada [x,y]) como capa de color semitransparente
    [W, H] = size(mask);
    C = repmat(reshape(color, 1, 1, 3), H, W);
    image(ax, [0.5, W-0.5]*s, [0.5, H-0.5]*s, C, 'AlphaData', alfa*double(mask.'));
    set(ax, 'YDir', 'normal');
end

function estilo(ax, ttl, BG)
    axis(ax, 'image');
    set(ax, 'Color', BG, 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
        'GridColor', 'w', 'GridAlpha', 0.15, 'YDir', 'normal');
    title(ax, ttl, 'Color', 'w');
    xlabel(ax, 'x [celdas]');  ylabel(ax, 'y [celdas]');
end

function leyenda(ax, h, etiquetas, ubic)
    if nargin < 4, ubic = 'best'; end
    legend(ax, h, etiquetas, 'Location', ubic, 'TextColor', 'w', ...
           'Color', [0.1 0.1 0.1], 'EdgeColor', [0.4 0.4 0.4], 'FontSize', 8);
end
