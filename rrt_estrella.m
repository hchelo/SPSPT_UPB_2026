%% ================================================================
%  PATH PLANNING — RRT*  +  ROBOT KHEPERA (cinemática inversa por pulsos)
%
%  Mismo mapa (Mapa.txt) y mismo robot que el script de campos potenciales.
%  Convención del mapa:
%    0  = libre
%    1  = obstáculo
%    2  = inicio del robot
%   -1  = meta
%
%  Requiere MATLAB R2016b o superior. Sin toolboxes.
%  Coordenadas continuas en celdas (base 0): x en [0, ANCHO-1], y en [0, ALTO-1]
%  Cada celda = 2 cm.
% ================================================================
clear; clc; close all;

BG = [0.05 0.05 0.05];

%% ================================================================
%  LEER Y PARSEAR MAPA
% ================================================================
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta), carpeta = pwd; end
ruta_mapa = fullfile(carpeta, 'Mapa.txt');

M_raw = load(ruta_mapa, '-ascii');
[FILAS, COLS] = size(M_raw);
fprintf('  Archivo leido: %s\n', ruta_mapa);
fprintf('  Filas=%d  Cols=%d\n', FILAS, COLS);

M     = M_raw.';          % M(x+1, y+1)
ANCHO = COLS;
ALTO  = FILAS;

[ci, ri] = find(M_raw.' ==  2, 1);
[cg, rg] = find(M_raw.' == -1, 1);
if isempty(ci), error('No se encontro el inicio (valor 2) en el mapa'); end
if isempty(cg), error('No se encontro la meta (valor -1) en el mapa');  end
xi = ci - 1;  yi = ri - 1;
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
%  DILATAR OBSTÁCULOS (margen de seguridad, disco r = 2 celdas)
% ================================================================
RADIO_DILATE = 2;
[a, b]  = meshgrid(-RADIO_DILATE:RADIO_DILATE);
se      = double(a.^2 + b.^2 <= RADIO_DILATE^2);
M_dilat = double(conv2(double(M == 1), se, 'same') > 0);

% Liberar un pequeño entorno alrededor de inicio y meta (solo celdas que
% no son obstáculo real), para que el árbol pueda salir del inicio y
% conectar con la meta aunque estén cerca de una pared.
[XX, YY] = ndgrid(0:ANCHO-1, 0:ALTO-1);
zona = hypot(XX - xi, YY - yi) <= RADIO_DILATE + 1 | ...
       hypot(XX - xg, YY - yg) <= RADIO_DILATE + 1;
M_dilat(zona & M == 0) = 0;

occ = (M_dilat == 1);           % mapa de ocupación que usa el planificador

if occ(xi+1, yi+1), error('El inicio cae sobre un obstaculo'); end
if occ(xg+1, yg+1), error('La meta cae sobre un obstaculo');   end

%% ================================================================
%  RRT*
% ================================================================
rng(7);                       % semilla (cámbiala para otra corrida)

MAX_ITER  = 4000;             % iteraciones (más = ruta más cercana a la óptima)
ETA       = 3.0;              % paso máximo de extensión [celdas]
GAMMA     = 60;               % constante del radio de vecindad
R_MAX     = 2*ETA;            % radio máximo de vecindad [celdas]
GOAL_BIAS = 0.05;             % probabilidad de muestrear la meta
TOL_META  = ETA;              % distancia para intentar conectar con la meta
RES_COL   = 0.25;             % resolución del chequeo de colisión [celdas]

inicio = [xi, yi];
meta   = [xg, yg];

nodes  = zeros(MAX_ITER+1, 2);
parent = zeros(MAX_ITER+1, 1);    % 0 = raíz
cost   = zeros(MAX_ITER+1, 1);
nodes(1,:) = inicio;
n = 1;

goalCand   = [];                  % nodos que pueden conectarse a la meta
hist_iter  = [];
hist_costo = [];
iter_primera = NaN;

fprintf('\nEjecutando RRT* (%d iteraciones)...\n', MAX_ITER);
tic;
for it = 1:MAX_ITER
    % 1) Muestreo
    if rand < GOAL_BIAS
        q = meta;
    else
        q = [rand*(ANCHO-1), rand*(ALTO-1)];
    end

    % 2) Nodo más cercano
    d_all = hypot(nodes(1:n,1) - q(1), nodes(1:n,2) - q(2));
    [dmin, inear] = min(d_all);
    if dmin < 1e-9, continue; end
    qnear = nodes(inear, :);

    % 3) Steer
    if dmin > ETA
        qnew = qnear + ETA * (q - qnear) / dmin;
    else
        qnew = q;
    end
    if ~segmentoLibre(qnear, qnew, occ, RES_COL), continue; end

    % 4) Vecindad
    r = min(GAMMA * sqrt(log(n+1) / (n+1)), R_MAX);
    dists = hypot(nodes(1:n,1) - qnew(1), nodes(1:n,2) - qnew(2));
    vec   = find(dists <= r);

    % 5) Elegir el mejor padre
    bestP = inear;
    bestC = cost(inear) + dists(inear);
    for j = vec.'
        c = cost(j) + dists(j);
        if c < bestC && segmentoLibre(nodes(j,:), qnew, occ, RES_COL)
            bestP = j;  bestC = c;
        end
    end

    n = n + 1;
    nodes(n,:) = qnew;
    parent(n)  = bestP;
    cost(n)    = bestC;

    % 6) Rewiring
    for j = vec.'
        if j == bestP, continue; end
        c = bestC + dists(j);
        if c < cost(j) - 1e-9 && segmentoLibre(qnew, nodes(j,:), occ, RES_COL)
            delta = c - cost(j);
            parent(j) = n;
            % Propagar la mejora a todos los descendientes de j
            cola = j;
            while ~isempty(cola)
                qq = cola(1);  cola(1) = [];
                cost(qq) = cost(qq) + delta;
                cola = [cola; find(parent(1:n) == qq)]; %#ok<AGROW>
            end
        end
    end

    % 7) ¿Conecta con la meta?
    if norm(qnew - meta) <= TOL_META && segmentoLibre(qnew, meta, occ, RES_COL)
        goalCand(end+1) = n; %#ok<SAGROW>
        if isnan(iter_primera), iter_primera = it; end
    end

    % 8) Registro de convergencia
    if ~isempty(goalCand) && mod(it, 20) == 0
        dG = hypot(nodes(goalCand,1) - xg, nodes(goalCand,2) - yg);
        hist_iter(end+1)  = it;                                  %#ok<SAGROW>
        hist_costo(end+1) = min(cost(goalCand) + dG);            %#ok<SAGROW>
    end
end
t_plan = toc;

if isempty(goalCand)
    error('RRT* no encontro ruta a la meta. Aumenta MAX_ITER o revisa el mapa.');
end

nodes  = nodes(1:n, :);
parent = parent(1:n);
cost   = cost(1:n);

% --- Extraer la mejor ruta
dG = hypot(nodes(goalCand,1) - xg, nodes(goalCand,2) - yg);
[costo_final, kbest] = min(cost(goalCand) + dG);
idx  = goalCand(kbest);
ruta = meta;
while idx ~= 0
    ruta = [nodes(idx,:); ruta]; %#ok<AGROW>
    idx  = parent(idx);
end
% Quitar puntos repetidos (si la meta fue muestreada exactamente)
keep = [true; hypot(diff(ruta(:,1)), diff(ruta(:,2))) > 1e-9];
ruta = ruta(keep, :);

% --- Poda por línea de vista (menos giros para el Khepera)
PODAR = true;
ruta_podada = podarRuta(ruta, occ, RES_COL);
if PODAR
    WP = ruta_podada;
else
    WP = ruta;
end

WP_x = WP(:,1);
WP_y = WP(:,2);
NWP  = numel(WP_x);

TH = zeros(NWP, 1);
for k = 1:NWP-1
    TH(k) = atan2(WP_y(k+1) - WP_y(k), WP_x(k+1) - WP_x(k));
end
TH(NWP) = atan2(WP_y(end) - WP_y(end-1), WP_x(end) - WP_x(end-1));

%% ================================================================
%  REPORTE EN CONSOLA
% ================================================================
L_rrt = longitud(ruta);
L_pod = longitud(ruta_podada);
fprintf('\n%s\n', repmat('=',1,55));
fprintf('   RESULTADO RRT*\n');
fprintf('%s\n', repmat('=',1,55));
fprintf('  Tiempo de planificacion : %.2f s\n', t_plan);
fprintf('  Nodos en el arbol       : %d\n', n);
fprintf('  Primera solucion en it. : %d\n', iter_primera);
fprintf('  Longitud ruta RRT*      : %.2f celdas = %.1f cm  (%d nodos)\n', L_rrt, 2*L_rrt, size(ruta,1));
fprintf('  Longitud ruta podada    : %.2f celdas = %.1f cm  (%d nodos)\n', L_pod, 2*L_pod, size(ruta_podada,1));
fprintf('%s\n', repmat('=',1,55));

fprintf('\nWaypoints para el Khepera: %d\n\n', NWP);
fprintf('  %2s  %7s  %7s  %6s  %6s  %7s\n', '#', 'X(cel)', 'Y(cel)', 'X(cm)', 'Y(cm)', 'theta');
fprintf('  %s\n', repmat('-', 1, 47));
for k = 1:NWP
    lab = '';
    if k == 1,   lab = '  <- INICIO'; end
    if k == NWP, lab = '  <- META';   end
    fprintf('  %2d  %7.2f  %7.2f  %6.1f  %6.1f  %+6.1f deg%s\n', ...
        k-1, WP_x(k), WP_y(k), 2*WP_x(k), 2*WP_y(k), rad2deg(TH(k)), lab);
end

%% ================================================================
%  FIGURA 1 — VISUALIZACIÓN (6 paneles)
% ================================================================
fig1 = figure('Name', 'RRT*', 'Color', BG, 'Position', [40 60 1600 900]);
xs = [0.5, ANCHO-0.5];   ys = [0.5, ALTO-0.5];

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

% --- 3. Árbol RRT*
ax = subplot(2,3,3); hold(ax, 'on');
capa(ax, occ,    [0.50 0.50 0.50], 0.4, 1);
capa(ax, M == 1, [0.25 0.45 0.80], 0.8, 1);
hijos = find(parent > 0);
TX = [nodes(hijos,1), nodes(parent(hijos),1), nan(numel(hijos),1)].';
TY = [nodes(hijos,2), nodes(parent(hijos),2), nan(numel(hijos),1)].';
hA = plot(ax, TX(:), TY(:), 'Color', [0.35 0.75 1 0.45], 'LineWidth', 0.5);
hR = plot(ax, ruta(:,1), ruta(:,2), 'Color', [1 0.2 0.6], 'LineWidth', 2);
plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, sprintf('Arbol RRT* (%d nodos)', n), BG);
xlim(ax, [0 ANCHO]);  ylim(ax, [0 ALTO]);
leyenda(ax, [hA hR], {'Arbol', 'Mejor ruta'}, 'northwest');

% --- 4. Convergencia del costo
ax = subplot(2,3,4); hold(ax, 'on');
plot(ax, hist_iter, 2*hist_costo, 'Color', [0 0.9 0.46], 'LineWidth', 1.8);
set(ax, 'Color', BG, 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
    'GridColor', 'w', 'GridAlpha', 0.15);
grid(ax, 'on');
title(ax, 'Convergencia del costo RRT*', 'Color', 'w');
xlabel(ax, 'Iteracion');  ylabel(ax, 'Longitud de la mejor ruta [cm]');

% --- 5. Ruta RRT* vs ruta podada
ax = subplot(2,3,5); hold(ax, 'on');
capa(ax, occ,    [0.50 0.50 0.50], 0.4, 1);
capa(ax, M == 1, [0.25 0.45 0.80], 0.8, 1);
hR1 = plot(ax, ruta(:,1), ruta(:,2), '.-', 'Color', [1 0.2 0.6], 'LineWidth', 1.2, 'MarkerSize', 10);
hR2 = plot(ax, ruta_podada(:,1), ruta_podada(:,2), '-', 'Color', [0 1 0], 'LineWidth', 2);
plot(ax, xi, yi, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(ax, xg, yg, 'r*', 'MarkerSize', 12);
estilo(ax, 'Ruta RRT* vs Ruta Podada', BG);
xlim(ax, [0 ANCHO]);  ylim(ax, [0 ALTO]);
leyenda(ax, [hR1 hR2], {sprintf('RRT* (%.1f cm)', 2*L_rrt), ...
                        sprintf('Podada (%.1f cm)', 2*L_pod)}, 'northwest');

% --- 6. Waypoints finales + tabla
ax = subplot(2,3,6); hold(ax, 'on');
capa(ax, occ,    [0.50 0.50 0.50], 0.4, 1);
capa(ax, M == 1, [0.25 0.45 0.80], 0.8, 1);
hRr = plot(ax, WP_x, WP_y, 'Color', [0 1 0], 'LineWidth', 1.5);
hW  = scatter(ax, WP_x, WP_y, 90, 'y', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.8);
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
estilo(ax, 'Waypoints Khepera', BG);
xlim(ax, [0 ANCHO]);  ylim(ax, [0 ALTO]);  grid(ax, 'on');
leyenda(ax, [hRr hW hI hM], {'Ruta', sprintf('%d Waypoints', NWP), 'Inicio', 'Meta'}, 'northwest');

pos = ax.Position;
ax.Position = [pos(1)-0.04, pos(2), pos(3)*0.72, pos(4)];
lineas = {'Waypoints Khepera:', '', ...
          sprintf('  %2s  %5s  %5s  %5s  %5s  %7s', '#', 'X', 'Y', 'Xcm', 'Ycm', 'th'), ...
          ['  ' repmat('-', 1, 37)]};
for k = 1:NWP
    lineas{end+1} = sprintf('  %2d  %5.1f  %5.1f  %5.1f  %5.1f  %+6.1f', ...
        k-1, WP_x(k), WP_y(k), 2*WP_x(k), 2*WP_y(k), rad2deg(TH(k))); %#ok<SAGROW>
end
text(ax, 1.04, 0.99, lineas, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'FixedWidth', 'FontSize', 7, 'Color', 'w', 'Interpreter', 'none', ...
     'BackgroundColor', [0.10 0.10 0.18], 'EdgeColor', [0.27 0.27 0.67]);

sgtitle(fig1, 'Path Planning — RRT*  |  Mapa 100x100 cm  (2 cm/celda)', ...
        'Color', 'w', 'FontSize', 13);

%% ================================================================
%  FIGURA 2 — ROBOT KHEPERA NAVEGANDO POR LOS WAYPOINTS
%    1. ORIENTAR   -> girar en sitio hasta apuntar al waypoint
%    2. DESPLAZAR  -> avanzar recto contando pulsos
%    3. REORIENTAR -> girar hasta la orientación final del waypoint
% ================================================================
P.r    = 0.008;                 % radio rueda [m]
P.l    = 0.054;                 % distancia entre ruedas [m]
P.dt   = 0.05;                  % paso de tiempo [s]
P.PPR  = 600;                   % pulsos por revolución
P.MMP  = (2*pi*P.r) / P.PPR;    % metros por pulso
P.VLIN = 0.03;                  % m/s
P.VROT = 0.25;                  % rad/s

R_ROBOT_K = 0.0275;
R_RUEDA_K = 0.005;
L_RUEDA_K = 0.010;
CEL_A_M   = 0.02;

wp_metros = [WP_x*CEL_A_M, WP_y*CEL_A_M, TH];

fprintf('\nSimulando Khepera en waypoints...\n');
enc        = [0 0];
reg_pulsos = [];
traj_k     = zeros(0, 3);
fase_k     = false(0, 1);
x_k        = wp_metros(1, :);

for i = 1:NWP-1
    [tramo, x_k, enc, rgp, fs] = ir_wp_k(x_k, wp_metros(i+1, :), i, P, enc);
    traj_k     = [traj_k; tramo];      %#ok<AGROW>
    fase_k     = [fase_k; fs];         %#ok<AGROW>
    reg_pulsos = [reg_pulsos, rgp];    %#ok<AGROW>
end

COL_G = [1.00 0.843 0.000];
COL_A = [0.00 0.902 0.463];
colores = double(fase_k).*COL_G + double(~fase_k).*COL_A;

% --- Reporte de pulsos
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('   REPORTE DE PULSOS — KHEPERA (RRT*)\n');
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
fprintf('  %-28s %10s  %9.1f s\n', 'TIEMPO DE MISION', '', (size(traj_k,1)-1)*P.dt);
fprintf('%s\n', repmat('=', 1, 60));

% --- Escena
fig2 = figure('Name', 'Khepera RRT*', 'Color', BG, 'Position', [120 60 850 900]);
ax2  = axes(fig2);  hold(ax2, 'on');
set(ax2, 'Color', [0.03 0.03 0.03], 'XColor', [0.9 0.9 0.9], 'YColor', [0.9 0.9 0.9], ...
    'GridColor', 'w', 'GridAlpha', 0.15);

capa(ax2, occ,    [0.50 0.50 0.50], 0.35, CEL_A_M);
capa(ax2, M == 1, [0.25 0.45 0.80], 0.70, CEL_A_M);

TXm = TX*CEL_A_M;  TYm = TY*CEL_A_M;
hArb  = plot(ax2, TXm(:), TYm(:), 'Color', [0.35 0.75 1 0.12], 'LineWidth', 0.4);
hRuta = plot(ax2, ruta(:,1)*CEL_A_M, ruta(:,2)*CEL_A_M, 'Color', [1 0.2 0.6 0.5], 'LineWidth', 1);
for k = 1:NWP
    plot(ax2, wp_metros(k,1), wp_metros(k,2), 's', 'Color', 'y', 'MarkerSize', 7);
    text(ax2, wp_metros(k,1)+0.002, wp_metros(k,2)+0.002, num2str(k-1), 'Color', 'y', 'FontSize', 8);
end
hIni  = plot(ax2, wp_metros(1,1),   wp_metros(1,2),   'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
hMeta = plot(ax2, wp_metros(end,1), wp_metros(end,2), 'r*', 'MarkerSize', 14);

axis(ax2, 'equal');
xlim(ax2, [0 ANCHO*CEL_A_M]);  ylim(ax2, [0 ALTO*CEL_A_M]);
xlabel(ax2, 'x [m]');  ylabel(ax2, 'y [m]');  grid(ax2, 'on');
leyenda(ax2, [hArb hRuta hIni hMeta], {'Arbol RRT*', 'Ruta RRT* (sin podar)', 'Inicio', 'Meta'}, 'northeast');
sgtitle(fig2, {'Robot Khepera — RRT* + Cinematica Inversa por Pulsos', ...
               'Amarillo: Orientar   |   Verde: Desplazar   (3 pasos por waypoint)'}, ...
        'Color', 'w', 'FontSize', 11);

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

    set(hBody, 'XData', cx + R_ROBOT_K*cos(ang), 'YData', cy + R_ROBOT_K*sin(ang), 'FaceColor', col);
    set(hDir, 'XData', [cx, cx + 1.7*R_ROBOT_K*cos(th)], ...
              'YData', [cy, cy + 1.7*R_ROBOT_K*sin(th)]);
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
    if fi > 1
        set(hTraj, 'XData', traj_k(1:fi-1,1), 'YData', traj_k(1:fi-1,2), ...
                   'CData', colores(1:fi-1,:));
    end
    if fase_k(fi), txtF = 'Orientando'; else, txtF = 'Desplazando'; end
    hTit.String = sprintf('%s  |  t=%.1fs  |  pos=(%.1f,%.1f)cm  th=%.0f deg  |  EncL=%.0fp  EncR=%.0fp', ...
        txtF, (fi-1)*P.dt, cx*100, cy*100, rad2deg(th), enc(1), enc(2));

    drawnow;
    pause(0.03);
end

%% ================================================================
%  FUNCIONES LOCALES
% ================================================================

% --- Chequeo de colisión de un segmento sobre la grilla de ocupación
function ok = segmentoLibre(p1, p2, occ, res)
    d   = norm(p2 - p1);
    ns  = max(2, ceil(d / res) + 1);
    t   = linspace(0, 1, ns).';
    pts = p1 + t .* (p2 - p1);
    ix  = round(pts(:,1)) + 1;
    iy  = round(pts(:,2)) + 1;
    [W, H] = size(occ);
    if any(ix < 1 | ix > W | iy < 1 | iy > H)
        ok = false;  return;
    end
    ok = ~any(occ(sub2ind([W H], ix, iy)));
end

% --- Poda por línea de vista (shortcut greedy)
function R = podarRuta(ruta, occ, res)
    N = size(ruta, 1);
    R = ruta(1, :);
    i = 1;
    while i < N
        j = N;
        while j > i + 1 && ~segmentoLibre(ruta(i,:), ruta(j,:), occ, res)
            j = j - 1;
        end
        R(end+1, :) = ruta(j, :); %#ok<AGROW>
        i = j;
    end
end

function L = longitud(ruta)
    L = sum(hypot(diff(ruta(:,1)), diff(ruta(:,2))));
end

% --- Cinemática directa (diferencial)
function xn = cinem_dir(x, vL, vR, P)
    v  = P.r * (vR + vL) / 2;
    w  = P.r * (vR - vL) / P.l;
    xn = [x(1) + v*cos(x(3))*P.dt, ...
          x(2) + v*sin(x(3))*P.dt, ...
          x(3) + w*P.dt];
    xn(3) = atan2(sin(xn(3)), cos(xn(3)));
end

function p = m2p(metros, P)
    p = abs(metros) / P.MMP;
end

function p = a2p(rad, P)
    p = abs((P.l/2) * rad) / P.MMP;
end

% --- ORIENTAR (giro en sitio)
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

% --- DESPLAZAR (avance recto)
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