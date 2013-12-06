% SUB-SCRIPT of Phase_screen_main.m
%   Setup geometry
%   <LENGTH UNITS: m>

params_flag = 0;

% Define coinc and invers variables if not yet defined
try
    coinc;
catch coin_var_err
    coinc = 0;
    invers = 0;
end
    

%% ------------ [USER-DEFINED PARAMETERS]
% Turbulence parameters
g0 = [0; 4.29; 7.99; 12.1; 16.6; 21.3; 27.9] * 1e-6; % Measured gammas
g = [linspace(1e-6,30e-6,0)'; g0];
%g = [0; 120; 200] * 1e-6;
%g0 = g;
L0 = Inf;               % Turbulence outer scale
l0 = 0;                 % Turbulence inner scale

% Beam and observation plane parameters
wvl = 0.325e-6;         % Wavelength
wn = 60e-6;             % Beam waist (assumed to be at observation plane)

% Slit width
a0 = 50e-6;

% Model parameters
nreals = 500;           % Number of turbulence realizations
z = linspace(0,0.89, ...
    50);                % Propagation planes
zscr_min = 0.44;        % Plane before which there is no turbulence
zscr_max = 0.63;        % Plane after which there is no turbulence
N = 512;                % Grid size in each transverse direction
delta1 = 5.5e-5;        % Grid spacing on the source plane
deltan = 7e-6;          % Grid spacing on the observation plane


%% ------------ [DERIVED QUANTITIES]
k = 2*pi/wvl;           % Wavenumber
g = sort(g);
n = length(z);          % Number of propagation planes
Dz = z(end);            % Total propagation distance
leng = length(g);       % Quantity of turbulence strengths considered
leng0 = length(g0);     % Number of measured turbulence strengths
%alpha = z / Dz;        % Normalized distances
delta = (1-z/Dz) * delta1 + ...
    z/Dz * deltan;      % Grid spacing vector
D2 = 4*wn;              % Region of interest at observation plane
w1 = wn*sqrt(1 + 2*Dz/ ...
    (k*wn^2));          % Beam width at source plane
D1 = 4*w1;              % Region of interest at source plane
qz = -Dz - ...
    1i*k*wn^2/2;        % Beam complex parameter
a = round(a0/deltan);   % Slit width in grid units
nz = length(z);         % Number of propagation planes

% Find indexes of measured gammas (g0) in the enlarged gamma matrix (g)
g0indx = zeros(length(g0),1);
for v = 1 : length(g0)
    idx_vec = find(~(g-g0(v)));
    g0indx(v) = idx_vec(1);
end

% Coherence radius r0i (plane-wave)
zmask = (z > zscr_min & z < zscr_max);
zt = zmask .* z;
zt_idx = find(zt);
nscr = length(zt_idx);
ztmin_idx = zt_idx(1);
ztmax_idx = zt_idx(end);
s = sum((1-zt(ztmin_idx:ztmax_idx)/Dz).^(5/3));

if coinc
    ktemp = k/2;
else
    ktemp = k;
end

r0 = (0.423 * ktemp^2/(7.75 * Dz^(5/3) * s) * g.^2).^(-3/5);
r0 = repmat(r0, 1, n);
r0 = r0 .* repmat(zmask, leng, 1);
r0(isinf(1./r0)) = inf;
r0(isnan(r0)) = inf;

clear s;

% Total coherence radius  r0sw (spherical wave)
aux_mat = repmat(zt/Dz,leng,1);
r0sw = sum(r0.^(-5/3) .* aux_mat.^(5/3), 2).^(-3/5);
clear aux_mat;

%% ------------ [OTHER VARIABLES]
[x1, y1] = meshgrid((-N/2 : N/2-1) * delta1);
sg = exp(-(x1/(0.47*N*delta1)).^16 ... 
    -(y1/(0.47*N*delta1)).^16);     % Super-gaussian filter
Uin = 1/qz*exp(1i*k/(2*qz)*(x1.^2+y1.^2));  % Input beam profile
