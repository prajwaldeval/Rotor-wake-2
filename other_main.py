import matplotlib.pyplot as plt
import numpy as np
from old_defs import *

# applied B.C.
B = 3  # number of blades
U0 = 10  # free stream velocity
TSR = 8  # tip speed ratio
R = 50  # rotor radius
rho = 1.225  # density
Omega = U0 * TSR / R  # angular velocity

a = 0.25  # induction factor for CT = 0.75
U_wake = U0 * (1 - a)  # near wake velocity

## Bound vortex grid
N_cp = 20  # nr. of control points (blade elements)

r_cos,dr_cos = geometry_cosine(N_cp)
r_con,dr_con = geometry_constant(N_cp)

# cosine spacing(0 = disabled; 1 = enabled)
cosspac = 0

if cosspac == 1:
    r = r_cos
else:
    r = r_con

# Iteration parameters
Niter = 1200
error_limit = 0.01
error = 1.0
ConvWeight = 0.3
t = 0

r_bp = r * R
r_cp = (r_bp(2:end) + r_bp(1:end - 1)) / 2
dr_cp = (r_bp(2:end) - r_bp(1:end - 1)) / 2

# input data interpolation
blade_pitch = pitch_WT. * ones(N_cp, 1)
c = interp1(r_cp_WT, c_WT, r_cp)
twist = interp1(r_cp_WT, twist_WT, r_cp)
gamma = interp1(r_cp_WT, gamma_WT, r_cp)

## set up wake geometry and horseshoe vortex grid
# determine total number of filament for each spanwise section
# each section will have two trailing filament and one bound filament

% set
up
wake
length in terms
of
rotor
diameter.
D = 2 * R;
L0 = 2 * D; % reference
length
for wake length
    Lx_wake = 4 * D; % specified
    wake
    length

# compute the required time to discretize the trailing filament.
tend = Lx_wake / U_wake

# set desired discretized number
for trailing filament.
    N_trail = 300 * (Lx_wake / L0)  # number of trailing filament grids

# wake expansion
U_wake = linspace(U0 * (1 - a), U0 * (1 - 2 * a), N_trail - 1)
# U_wake = U0 * (1 - a) * ones(1, N_trail - 1)
r_exp = U_wake(1)/ U_wake

# Finally, we can set up the coordinate(x, y, z) at each discretized
# filament.For(x, y, z), 2 * N_trail + 1 represents two trailing vortices
% and a
bound
filament.
% The
arrangement:
% N_trail + 2, N_trail + 3, N_trail + 4, ..., 2 * N_trail + 1(trailing
votrex)
% N_trail + 1(bound
vortex)
% N_trail, N_trail - 1, N_trail - 2, ..., 1(trailing
votrex)
x_trail = zeros(N_cp, 2 * N_trail + 1, B);
y_trail = zeros(N_cp, 2 * N_trail + 1, B);
z_trail = zeros(N_cp, 2 * N_trail + 1, B);

% t_trail is used
to
compute
the
frozen
wake
coordinate.
% tend is directly
derived
from specified wake

length and wake
velocity.
t_trail = linspace(0, tend, N_trail - 1);

for nn = 1:B
x_trail(:, 1: N_trail - 1, nn)   = repmat(c / 4 + c, [1 N_trail - 1])...
+ flip(t_trail). * flip(U_wake);
x_trail(:, N_trail + 3: end, nn) = flip(x_trail(:, 1: N_trail - 1), 2);
x_trail(:, N_trail + 1, nn)     = c / 4;
x_trail(:, N_trail, nn)       = c / 4;
x_trail(:, N_trail + 2, nn)     = c / 4;

y_trail(:, 1: N_trail - 1, nn)   = flip(r_exp). * repmat(r_cp - dr_cp, [1 N_trail - 1]). * ...
sin(Omega * flip(t_trail) + (360 / B) * (nn - 1) / 180 * pi);
y_trail(:, N_trail + 3: end, nn) = r_exp. * repmat(r_cp + dr_cp, [1 N_trail - 1]). * ...
sin(Omega * t_trail + (360 / B) * (nn - 1) / 180 * pi);
y_trail(:, N_trail + 1, nn)     = r_cp * sin((360 / B) * (nn - 1) / 180 * pi);
y_trail(:, N_trail, nn)       = (r_cp - dr_cp) * sin((360 / B) * (nn - 1) / 180 * pi);
y_trail(:, N_trail + 2, nn)     = (r_cp + dr_cp) * sin((360 / B) * (nn - 1) / 180 * pi);

z_trail(:, 1: N_trail - 1, nn)   = flip(r_exp). * repmat(r_cp - dr_cp, [1 N_trail - 1]). * ...
cos(Omega * flip(t_trail) + (360 / B) * (nn - 1) / 180 * pi);
z_trail(:, N_trail + 3: end, nn) = r_exp. * repmat(r_cp + dr_cp, [1 N_trail - 1]). * ...
cos(Omega * t_trail + (360 / B) * (nn - 1) / 180 * pi);
z_trail(:, N_trail + 1, nn)     = r_cp * cos((360 / B) * (nn - 1) / 180 * pi);
z_trail(:, N_trail, nn)       = (r_cp - dr_cp) * cos((360 / B) * (nn - 1) / 180 * pi);
z_trail(:, N_trail + 2, nn)     = (r_cp + dr_cp) * cos((360 / B) * (nn - 1) / 180 * pi);
end

% % % visualization
of
the
frozen
wake
geometry
figure(1);
for nn = 1:3
for ii = 1:1: N_cp
plot3(x_trail(ii,:, nn), y_trail(ii,:, nn), z_trail(ii,:, nn), 'color', [0, 0.4470, 0.7410]);
hold
on;
grid
on;
plot3([0 0], [r_bp(1) R] * sin((360 / B) * (nn - 1) / 180 * pi), ...
[r_bp(1)
R]*cos((360 / B) * (nn - 1) / 180 * pi), ...
'linewidth', 2, 'color', [0.5 0.5 0.5])
plot3(c, r_cp * sin((360 / B) * (nn - 1) / 180 * pi), ...
r_cp * cos((360 / B) * (nn - 1) / 180 * pi), ...
'linewidth', 2, 'color', [0.5, 0.5, 0.5])
end
end
xlabel('x', 'fontsize', 14);
ylabel('y', 'fontsize', 14);
zlabel('z', 'fontsize', 14);
view(3);
axis([0 50 - 100 100 - 100 100]);

% % set
up
unit
strength
induction
matrix
% initalize and calculate
matrices
% for velocity induced by horseshoe vortex rings

% setting
unit
induction
matrix
using
biot - savart.
unitU_ind = zeros(N_cp, N_cp);
unitV_ind = zeros(N_cp, N_cp);
unitW_ind = zeros(N_cp, N_cp);
for ii = 1:N_cp
for nn = 1:N_cp
for jj = 1:2 * N_trail
for nb = 1:B
[uu, vv, ww] = ...
v3d_biotsavart(1, x_trail(nn, jj, nb), x_trail(nn, jj + 1, nb), ...
y_trail(nn, jj, nb), y_trail(nn, jj + 1, nb), ...
z_trail(nn, jj, nb), z_trail(nn, jj + 1, nb), ...
x_trail(ii, N_trail + 1, 1), ...
y_trail(ii, N_trail + 1, 1), ...
z_trail(ii, N_trail + 1, 1));
uu(isnan(uu)) = 0;
vv(isnan(vv)) = 0;
ww(isnan(ww)) = 0;
unitU_ind(ii, nn) = unitU_ind(ii, nn) + uu;
unitV_ind(ii, nn) = unitV_ind(ii, nn) + vv;
unitW_ind(ii, nn) = unitW_ind(ii, nn) + ww;
end
end
end
end

% % Solution
iteration.
% calculate
solution
through
an
iterative
process
% Specifying
maximum
iteration
number, compute
the
circulation.
for it = 1:3000
radposition = sqrt(x_trail(:, N_trail + 1, 1).^ 2 + ...
y_trail(:, N_trail + 1, 1).^ 2 + ...
z_trail(:, N_trail + 1, 1).^ 2);

% rotational
velocity
tmp = repmat([-Omega 0 0], [N_cp 1]);
V_rot = cross(tmp, [x_trail(:, N_trail + 1, 1) ...
y_trail(:, N_trail + 1, 1) ...
z_trail(:, N_trail + 1, 1)]);

% compute
total
perceived
velocity
u_total = U0 + unitU_ind * gamma + V_rot(:, 1);
v_total = unitV_ind * gamma + V_rot(:, 2);
w_total = unitW_ind * gamma + V_rot(:, 3);

% calculate
azimuthal and axial
velocity
tmp = repmat([-1 0 0], [N_cp 1]). / radposition;
azimdir = cross(tmp, [x_trail(:, N_trail + 1, 1) ...
y_trail(:, N_trail + 1, 1) ...
z_trail(:, N_trail + 1, 1)]);
Vazim = zeros(N_cp, 1);
for i = 1:N_cp
Vazim(i) = dot(azimdir(i,:), [u_total(i) v_total(i) w_total(i)]); % azimuthal
direction
end
Vaxial = u_total; % axial
velocity

% compute
blade
loads, [fnorm, ftan, gamma]
loads = loadBladeElement(Vaxial, Vazim, twist, c, blade_pitch, af, rho);

% new
estimate
of
circulation
for the blade section
gamma_tmp = loads.gamma;

% derive
other
output
variables
a = -(unitU_ind * gamma + V_rot(:, 1)) / U0;
ap = (Vazim. / (radposition. * repmat(Omega, [N_cp 1])) - 1);
Fnorm = loads.fnorm;
Ftan = loads.ftan;

% update
circulation.
gamma = (1 - ConvWeight) * gamma + ConvWeight * gamma_tmp;

% check
convergence
of
solution
% abserr = max(gamma_tmp - gamma);
% if (abserr < error_iterations)
    % break
% end
% scatter(it, max(gamma));
hold
on;
end

% non - dimensional
factor.
Non_dim_gamma = pi * U0 ^ 2 / (B * Omega);
Non_dim_F = 0.5 * rho * U0 ^ 2 * R;

% compute
rotor
performance
CT = sum(2 * dr_cp. * Fnorm * B. / (0.5 * U0 ^ 2 * pi * R * R))
CP = sum(2 * dr_cp. * Ftan. * r_cp * B * Omega. / (0.5 * U0 ^ 3 * pi * R * R))

figure(2)
p0 = plot(r_cp_WT, gamma_WT / Non_dim_gamma);
hold
on;
p1 = plot(r_cp, gamma / Non_dim_gamma);
lg = legend([p0 p1], {'BEM', 'Lifting Line'});

figure(3)
p0 = plot(r_cp, Fnorm / Non_dim_F);
hold
on;
p1 = plot(r_cp, Ftan / Non_dim_F);
hold
on;
p3 = plot(r_cp_WT, Ftan_WT / Non_dim_F, '--');
p2 = plot(r_cp_WT, Fnorm_WT / Non_dim_F, '--');
lg = legend([p0 p1 p2 p3], {'Lifting Line Fn', 'Lifting Line Ft', 'BEM Fn', 'BEM Ft'});