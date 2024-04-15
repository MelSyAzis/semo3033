clear; clf; clc;

n_elements = 2
n_nodes = 4
E = 30e6
pois = 0.25
thickness = 0.5

post_displacement_scale = 2e3


% nodal coordinates
coord_n_x = [
  3;  % node 1
  3;  % node 2
  0;  % node 3
  0   % node 4
]

coord_n_y = [
  0;  % node 1
  2;  % node 2
  2;  % node 3
  0;  % node 4
]


% element connectivity
conn = [
  1 2 4; % element 1
  3 4 2  % element 2
]


% Force vector

fx = [
  0; % node 1
  0; % node 2
  0; % node 3
  0; % node 4
]

fy = [
  0;        % node 1
  -1000;   % node 2
  0;        % node 3
  0;        % node 4
]


% BOUNDARY CONDITIONS

%   Set whether BC is to be set for each dof.
%     1: ON
%     2: OFF

bc_qx_on = [
  0;  % node 1
  0;  % node 2
  1;  % node 3
  1;  % node 4
]

bc_qy_on = [
  1;  % node 1
  0;  % node 2
  1;  % node 3
  1;  % node 4
]


%   Set BC values only to nodes where BC is turned on
  
bc_qx = [
  0;  % node 1
  0;  % node 2
  0;  % node 3
  0;  % node 4
]

bc_qy = [
  0;  % node 1
  0;  % node 2
  0;  % node 3
  0;  % node 4
]



% =======================================
% PRE-PROCESSING & SOLVER

% No of dof

n_dof = n_nodes * 2

% Allocate memory

J = zeros(2,2, n_elements);
detJ = zeros(n_elements,1);
A = zeros(n_elements,1);
B = zeros(3, 6, n_elements);
D = zeros(3, 3, n_elements);
kloc = zeros(6,6, n_elements);
kglobal = zeros(8,8, n_elements);
q = zeros(n_elements, 1);


% jacobian

for i = 1:n_elements
  J(:,:,i) = [ 
    coord_n_x( conn(i,1) ) - coord_n_x( conn(i,3) ), ...
      coord_n_y( conn(i,1) ) - coord_n_y( conn(i,3) );
    coord_n_x( conn(i,2) ) - coord_n_x( conn(i,3) ), ...
      coord_n_y( conn(i,2) ) - coord_n_y( conn(i,3) )
  ] ;
end

J

for i = 1:n_elements
  detJ(i) = det( J(:,:,i) );
end

detJ


% Area

A = 0.5 * detJ

% B

for i = 1:n_elements
  y23 = coord_n_y ( conn(i, 2) ) - coord_n_y ( conn(i, 3) );
  y31 = coord_n_y ( conn(i, 3) ) - coord_n_y ( conn(i, 1) );
  y12 = coord_n_y ( conn(i, 1) ) - coord_n_y ( conn(i, 2) );
  x23 = coord_n_x ( conn(i, 2) ) - coord_n_x ( conn(i, 3) );
  x13 = coord_n_x ( conn(i, 1) ) - coord_n_x ( conn(i, 3) );
  x21 = coord_n_x ( conn(i, 2) ) - coord_n_x ( conn(i, 1) );
  x32 = coord_n_x ( conn(i, 3) ) - coord_n_x ( conn(i, 2) );
  
  B(:,:,i) = 1/detJ(i) * ...
  [
    y23 0 y31 0 y12 0;
    0 x32 0 x13 0 x21;
    x32 y23 x13 y31 x21 y12
  ];
end

B

% D

for i = 1:n_elements
  D(:,:,i) = E / (1 - pois^2) * ...
  [
    1 pois 0;
    pois 1 0;
    0 0 0.5*(1-pois)
  ];
  
end

D


% k_local

for i = 1:n_elements
  kloc(:, :, i) = thickness * A(i) * B (:, :, i)' * D(:,:,i) * B(:,:,i);
end

kloc


% dof_local

dof_element = zeros(6,n_elements);

for i = 1:n_elements
  for j = 1:3
    dof_element(j*2-1,i) = conn(i,j)*2 - 1;
    dof_element(j*2,i) = conn(i,j)*2;
  end
end


% assemble k_global

k_global = zeros(n_dof, n_dof);

for i = 1:n_elements
  for j = 1:6
    row_global = dof_element(j,i);
    for k = 1:6
      col_global = dof_element(k,i);
      k_global(row_global, col_global) += kloc(j,k,i);
    end
  end
end

k_global


% Assemble BC

bc_q_on = [ bc_qx_on bc_qy_on ]';
bc_q_on = reshape(bc_q_on, [n_dof, 1]);
bc_q = [ bc_qx bc_qy ]';
bc_q = reshape(bc_q, [n_dof, 1]);

bc_q_on
bc_q


% Assemble forces

f = [ fx fy ]';
f = reshape(f, [n_dof, 1]);
f_solve = f(bc_q_on == 0)



% k_global final after applying boundary conditions

k_global_final = k_global(bc_q_on == 0, bc_q_on == 0)

q_solve = k_global_final \ f_solve


% Assemble displacement

q(bc_q_on == 0) = q_solve;
q(bc_q_on == 1) = bc_q(bc_q_on == 1);

q



% POST-PROCESSING VALUES


%   Element values

%     Strain

strain = zeros(3, n_elements);

for i = 1:n_elements
  q_element_id2row = [ conn(i, :)*2-1, conn(i, :)*2 ];
  strain(:, i) = B(:,:,i) * q( q_element_id2row );
endfor

strain


%     Stress

stress = zeros(3, n_elements);

for i = 1:n_elements
  stress(:, i) = D(:,:,i) * strain(:, i);
endfor

stress


%     Principal stresses
stress_pr1 = zeros(n_elements, 1);
stress_pr2 = zeros(n_elements, 1);
stress_shearMax = zeros(n_elements, 1);
stress_normalAvg = zeros(n_elements, 1);

for i = 1:n_elements
  sAvg = 0.5*(stress(1,i) + stress(2,i));
  sDiff = 0.5*(stress(1,i)-stress(2,i));
  sSqrt = sqrt((sDiff).^2 + stress(3,i).^2);
  stress_pr1(i) = sAvg + sSqrt;
  stress_pr2(i) = sAvg - sSqrt;
  stress_shearMax(i) = sSqrt;
  stress_normalAvg(i) = sAvg;
end


%   Plotting

%     Displacement

figure; hold on; axis('equal');

%       Original dimension

for i = 1:n_elements
  for j = 1:3
    
    j2 = j+1;
    if j2 > 3
      j2 = 1;
    end

    node1 = conn(i,j);
    node2 = conn(i,j2);

    x1 = coord_n_x(node1);
    y1 = coord_n_y(node1);

    x2 = coord_n_x(node2);
    y2 = coord_n_y(node2);

    x1_disp = x1 + post_displacement_scale * q(node1*2 -1);
    y1_disp = y1 + post_displacement_scale * q(node1*2);

    x2_disp = x2 + post_displacement_scale * q(node2*2 -1);
    y2_disp = y2 + post_displacement_scale * q(node2*2);
    
    plot([x1 x2],[y1 y2],'k','LineWidth',2);
    plot([x1_disp x2_disp],[y1_disp y2_disp],'r','LineWidth',1);

  end
end
