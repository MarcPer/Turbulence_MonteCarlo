% SUB-SCRIPT of Phase_screen_main.m
%   Code for analyzing sampling constraints

close all;
% Set variables
if ~params_flag
    Phase_screen_params;
end
idxg = leng;


rad = Dz*(1 + (k*wn^2/(2*Dz))^2);   % Beam radius of curvature
c = 4;          % Model sensitivity (see pag. 173)
D1p = D1 + c*wvl*Dz/r0sw(idxg);
D2p = D2 + c*wvl*Dz/r0sw(idxg);

d1 = linspace(0, 1.1*wvl*Dz/D2p, 100);
dn = linspace(0, 1.1*wvl*Dz/D1p, 100);
[d1 dn] = meshgrid(d1, dn);

% Constraint 2
hndl = figure;
N2 = log2((wvl * Dz + D1p*dn + D2p*d1) ./ (2 * d1 .* dn));
[C, hh] = contourf(d1, dn, N2);
%clabel(C,hh, 'FontSize', 15, 'Rotation', 0, ...
%    'FontWeight', 'bold');
xlabel('\delta_1 [m]');
ylabel('\delta_n [m]');
colorbar;
hold all;

% Constraint 1
deltan_max = -D2p/D1p*d1 + wvl*Dz/D1p;
plot(d1(1,:), deltan_max(1,:), 'k--', 'Linewidth', 2);
axis([0 d1(end) 0 dn(end)]);
set(gca, 'Color', 'none', 'Layer', 'top');

% Constraint 3
dnmin3 = (1+Dz/rad)*d1 - wvl*Dz/D1;
dnmax3 = (1+Dz/rad)*d1 + wvl*Dz/D1;
plot(d1(1,:), dnmax3(1,:), 'k-.');
set(gca, 'Color', 'none', 'Layer', 'top');
plot(d1(1,:), dnmin3(1,:), 'k-.');
set(gca, 'Color', 'none', 'Layer', 'top');

% Plot values currently in use
plot(delta1, deltan, 'w*', 'MarkerSize', 20)
set(gca, 'Color', 'none', 'Layer', 'top');

hold off;

% Constraint 4
zmax = min([delta1 deltan])^2 * N / wvl;
nmin = ceil(Dz / zmax) + 1;

%% Report
check_c1 = (deltan < -D2p/D1p*delta1 + wvl*Dz/D1p);
Nmin = (wvl * Dz + D1p*deltan + D2p*delta1) ...
    ./ (2 * delta1 .* deltan);
check_c2 = (N > Nmin);
check_c3 = (deltan > (1+Dz/rad)*delta1 - wvl*Dz/D1) & ...
    (deltan < (1+Dz/rad)*delta1 + wvl*Dz/D1);
check_c4 = (max(delta) < zmax);

% Was there any constraint violation?
constr_err = 0;         % Violation flag
if (~check_c1 || ~check_c2 || ~check_c3)
    constr_err = 1;
end

fprintf('Veryfing sampling requirements...\n');
fprintf('Constraint 1: ');
if check_c1
    fprintf('Satisfied\n');
else
    fprintf('Not satisfied\n');
end
fprintf('Constraint 2: ');
if check_c2
    fprintf('Satisfied\n');
else
    fprintf('Not satisfied\n');
end
fprintf('Constraint 3: ');
if check_c3
    fprintf('Satisfied\n');
else
    fprintf('Not satisfied\n');
end
fprintf('Constraint 4: ');
if check_c4
    fprintf('Satisfied\n');
else
    fprintf('Not satisfied\n');
end