% 1D transient convection--diffusion solver
% Assumption:
% - Uniform cell size
% - fixed fluid properties (density, diffusion coefficient)

clear all;
close all;
clc;

% USER INPUT ===========================

%    geometric
L = 1; % length [m]
n = 50; % no of nodes
A = 1; % cross sectional area

%    material properties
rho = 1.0; % density [kg / m^3]
diff = 0.02; % diffusion coefficient [kg m^-1 s^-1)]



%     convection properties
U = 0.0 % Transport velocity [m / s]


%     convection term discretisation

conv_method = 2; % 1: central; 2: upwind;

%     transient parameters

t_end = 5.0 % [s]
dt = 0.02 % [s]


%    boundary conditions
BC_phi_L = 0; % phi at LHS end
BC_phi_R = 0; % phi at RHS end

%     initial conditions
phi_init = zeros(n,1) % Initial condition
pos_start = 0.4 % normalised start position for initial field value of 1
pos_end = 0.6 % normalised end position for initial field value of 1
phi_init(floor(pos_start*n)+1 : floor(pos_end*n)-1) = 1

% sine wave as initial condition




% PROCESSING BEGIN ======================


% Allocate memory to matrix/arrays
x = zeros(n,1); % node coordinates
phi_old = zeros(n,1) % initial condition
phi = zeros(n,1)
phi_all = zeros(n, size(0:dt:t_end)(1))
coeff = zeros(n,n); % coefficient matrix of governing system of equations
b = zeros(n,1); % RHS matrix of governing system of equations


% Gridding/meshing -------------------------

dx = L / n; %   Node spacing

%   Calculates node coordinates.
x(1) = dx / 2;
for i = 2:n
  x(i) = x(i-1) + dx;
end


phi = phi_init;
phi_all(:,1) = phi;

fig_handle = figure();
set(fig_handle, 'position', [100 100 800 500]);

title({['1-D convection-diffusion']})
xlabel('Spatial co-ordinate, x (m) \rightarrow')
ylabel('phi \rightarrow')
ylim([0 1.1])

% Time increment
for t = dt:dt:t_end
  CFL = U*dt/dx;
  sprintf("t = %4.2f, CFL = %f", t, CFL)


  phi_old = phi;

  % Constructing system of equations [A][u] = [b] ---------------------

  %   coefficient matrix [A]

  for i = 1:n

    D = diff / dx * A;
    F = rho * U * A;
    a_W = max([F, D + F/2, 0])
    a_E = max([-F, D - F/2, 0])
    a_P0 = rho * dx * A / dt;
    s_P = 0;
    s_u = 0;
    a_P = a_W + a_E + a_P0 - s_P;

    if i ~= 1
      coeff(i,i-1) = a_W;
    end

    coeff(i,i) = -a_P;

    if i ~= n
      coeff(i,i+1) = a_E;
    end

    b(i) = -s_u - a_P0*phi_old(i);

    if i == 1
      b(i) = b(i) - a_W*BC_phi_L;
    end
    if i == n
      b(i) = b(i) - a_E*BC_phi_R;
    end
  end

  % Solve the system of equations [u] = [A]^-1 [b]

  phi = coeff\b;

  phi_all(:,round(t/dt + 1)) = phi;

  % Animated plotting
  h=plot(x,phi,'*-');

  ylim([0 1.1])

  title(sprintf("1-D convection-diffusion, t = %4.2f", t));

  pause(dt);


end


% POSTPROCESSING --------------------




%   Plot solution at all times
if 0
  figure;
  h=plot(x,phi_all,'*-');       % plotting the variable profile
  hold on;

  if 0
    %   Compute and plot exact solution for steady convection-diffusion case Eq.(4.39)
    x_exact = linspace(0.0, L, 100);
    phi_exact = BC_phi_L + (BC_phi_R - BC_phi_L) * (exp(rho * U .* x_exact / diff) - 1) ./ (exp(rho * U * L / diff) - 1);
    plot(x_exact, phi_exact, '--');
  end

  %axis([0 2 0 3])
  title({['1-D convection-diffusion']})
  xlabel('Spatial co-ordinate, x (m) \rightarrow')
  ylabel('phi \rightarrow')
end
