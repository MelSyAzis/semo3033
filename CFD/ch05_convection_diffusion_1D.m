clear all;
close all;
clc;

% References:
% - derivation: https://onedrive.live.com/redir?resid=5A238304660466DF%21126&page=Edit&wd=target%28Log202404.one%7Caa9c0a72-11b3-2942-9e05-1e6c8a01e3b7%2F15%20SEMO%203033%20Derivation%201D%20Convection%7Cedd4a673-f8f1-4ab9-8b8a-09f7dd692418%2F%29&wdorigin=703


% USER INPUT ===========================

%    geometric
L = 1; % length [m]
n = 20; % no of nodes
A = 1; % cross sectional area

%    material properties
rho = 1.0; % density [kg / m^3]
diff = 0.1; % diffusion coefficient [kg m^-1 s^-1)]

%    boundary conditions
BC_phi_L = 1; % phi at LHS end
BC_phi_R = 0; % phi at RHS end


%     convection properties
U = -2.5 % Transport velocity [m / s]


%     convection term discretisation

conv_method = 2; % 1: central; 2: upwind;


% PROCESSING BEGIN ======================


% Allocate memory to matrix/arrays
x = zeros(n,1); # node coordinates
coeff = zeros(n,n); # coefficient matrix of governing system of equations
b = zeros(n,1); # RHS matrix of governing system of equations


% Gridding/meshing -------------------------

dx = L / n; %   Node spacing

%   Calculates node coordinates.
x(1) = dx / 2;
for i = 2:n
  x(i) = x(i-1) + dx;
end


% Constructing system of equations [A][u] = [b] ---------------------

%   coefficient matrix [A]

%     define coefficients for internal nodes (row 2 until row n-1).
%     See Eq (5.18) and those in page 138

for i = 2:n-1
  % diffusion
  D_W = diff / dx;
  D_E = diff / dx;
  D = D_W

  % convection
  F_W = rho * U;
  F_E = rho * U;
  F = F_W

  if conv_method == 1 # central (pg 138)
    a_W = D_W + F_W/2;
    a_E = D_E - F_E/2;
  elseif conv_method == 2 # upwind (pg 147, 148)
    a_W = D_W + max(F_W , 0);
    a_E = D_E + max(0, -F_E);
  end
  s_u = 0;
  s_P = 0;
  a_P = a_W + a_E + (F_E - F_W) - s_P;

  coeff(i,i-1) = a_W;
  coeff(i,i) = -a_P;
  coeff(i,i+1) = a_E;

  b(i) = -s_u;
end


%     define coefficients for boundary node LHS (row 1) Eq(4.33)

if conv_method == 1 % central
  a_W = 0;
  a_E = D_E - F_E/2;
  s_P = -(2*D_W + F_W);
  s_u = (2*D_W + F_W) * BC_phi_L;
elseif conv_method == 2 % upwind
  if U >= 0 % pg 148
    a_W = 0;
    a_E = D_E;
    s_P = -(2*D_W + F_W);
    s_u = (2*D_W + F_W) * BC_phi_L
  elseif U < 0
    a_W = 0;
    a_E = -F_E + D_E;
    s_P = -2*D_E;
    s_u = 2*D_E * BC_phi_L
  end
end

a_P = a_W + a_E + (F_E - F_W) - s_P;


coeff(1,1) = -a_P;
coeff(1,2) = a_E;

b(1) = -s_u;


%     define coefficients for boundary node RHS (row n) Eq(4.36)

if conv_method == 1 % central
  a_W = D_W - F_W/2;
  a_E = 0.0;
  s_u = (2*D_W - F_W) * BC_phi_R;
  s_P = -(2*D_W - F_W);
elseif conv_method == 2 % upwind
  if U >= 0 % pg 148
    a_W = F + D;
    a_E = 0.0;
    s_u = 2*D * BC_phi_R;
    s_P = -2*D;
  elseif U < 0
    a_W = D;
    a_E = 0;
    s_u = (-F + 2*D) * BC_phi_R;
    s_P = F - 2*D
  end
end

a_P = a_W + a_E + (F_E - F_W) - s_P;

coeff(n,n) = -a_P;
coeff(n,n-1) = a_W;

b(n) = -s_u;


% Solve the system of equations [u] = [A]^-1 [b]

phi = coeff\b;


% POSTPROCESSING --------------------


%   Compute exact solution Eq.(4.39)
x_exact = linspace(0.0, L, 100)
phi_exact = BC_phi_L + (BC_phi_R - BC_phi_L) * (exp(rho * U .* x_exact / diff) - 1) ./ (exp(rho * U * L / diff) - 1)


%   Plot solution

h=plot(x,phi,'*-');       % plotting the variable profile
hold on;

plot(x_exact, phi_exact, '--');

%axis([0 2 0 3])
title({['1-D convection-diffusion']})
xlabel('Spatial co-ordinate, x (m) \rightarrow')
ylabel('phi \rightarrow')
