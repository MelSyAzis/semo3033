clear all;
close all;
clc;

% USER INPUT ===========================

%    geometric
L = 0.02; % length [m]
n = 4; % no of nodes
A = 1; % cross sectional area

%    material properties
q = 10000e3; % heat generation [ W / m^3 ]
k = 0.5; % thermal conductivity [ W / m.K ]

%    boundary conditions
T_L = 100; % T at LHS end [deg C or K]
T_R = 200; % T at RHS end [deg C or K]


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
%     See Eq(4.30)

for i = 2:n-1
  a_W = k * A / dx;
  a_E = k * A / dx;
  s_u = -q * A * dx;
  s_P = 0;
  a_P = - a_W - a_E + s_P;

  coeff(i,i-1) = a_W;
  coeff(i,i) = a_P;
  coeff(i,i+1) = a_E;

  b(i) = s_u;
end


%     define coefficients for boundary node LHS (row 1) Eq(4.33)

a_W = 0;
a_E = k * A / dx;
s_u = -q * A * dx - 2 * k * A / dx * T_L;
s_P = - 2 * k * A / dx;
a_P = - a_W - a_E + s_P;

coeff(1,1) = a_P;
coeff(1,2) = a_E;

b(1) = s_u;


%     define coefficients for boundary node RHS (row n) Eq(4.36)

a_W = k * A / dx;
a_E = 0.0;
s_u = -q * A * dx - 2 * k * A / dx * T_R;
s_P = - 2 * k * A / dx;
a_P = - a_W - a_E + s_P;

coeff(n,n) = a_P;
coeff(n,n-1) = a_W;

b(n) = s_u;


% Solve the system of equations [u] = [A]^-1 [b]

U = coeff\b;


% POSTPROCESSING --------------------


%   Compute exact solution Eq.(4.39)
x_exact = linspace(0.0, L, 100)
U_exact = ( (T_R - T_L) / L + q/2/k*(L-x_exact) ) .* x_exact + T_L


%   Plot solution

h=plot(x,U,'*-');       %plotting the velocity profile
hold;
plot(x_exact, U_exact, '--');
%axis([0 2 0 3])
title({['1-D heat diffusion (by conduction)']})
xlabel('Spatial co-ordinate, x (m) \rightarrow')
ylabel('Temperature, T (C or K) \rightarrow')
