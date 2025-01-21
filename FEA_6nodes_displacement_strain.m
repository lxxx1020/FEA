% File: FEA_6nodes_displacement_strain.m
% Description: Finite Element Analysis for a 2D thin plate under plane stress conditions
% Author: [Xuan Li]
% 
% This script calculates nodal displacements and element strains for a
% 2D thin plate using 6-node triangular elements. The results are validated
% under different elastic moduli and loading conditions.clc;
clear;

% Material and geometric parameters
E = 40e9;          % Elastic modulus (Pa)
nu = 0.3;          % Poisson's ratio
t = 0.002;         % Thickness (m)

% Node coordinates (unit: m)
nodes_coordinates = [
    0, 0;           % Node 1
    30, 0;          % Node 2
    60, 30;         % Node 3
    30, 30;         % Node 4
    0, 30;          % Node 5
    15, 0;          % Node 6
    0, 15;          % Node 7
    15, 15;         % Node 8
    30, 15;         % Node 9
    45, 15;         % Node 10
    15, 30;         % Node 11
    45, 30          % Node 12
] * 1e-3; % Convert to meters

% Element definitions
elements = [
    1, 4, 5, 8, 11, 7;   % Element 1
    1, 2, 4, 6, 9, 8;    % Element 2
    2, 3, 4, 10, 12, 9;  % Element 3
];

% Gauss integration points and weights
gaussPoints = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
gaussWeights = [1/6; 1/6; 1/6];

% Elastic matrix
D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];

% Initialize global stiffness matrix and external force vector
numNodes = size(nodes_coordinates, 1);
K = zeros(2 * numNodes, 2 * numNodes); % Global stiffness matrix

% Angle array (unit: degrees)
theta = 270;
F_magnitude = 50; % Force magnitude (N)

% Convert angle to radians
theta_rad = deg2rad(theta);

% Calculate force components
Fx = F_magnitude * cos(theta_rad);
Fy = F_magnitude * sin(theta_rad);

fprintf('Current force direction theta = %.1f degrees (Fx = %.2f N, Fy = %.2f N)\n', theta, Fx, Fy);

% Compute force components based on theta
Fx = 50 * cos(theta_rad);  % Force in x-direction (N)
Fy = 50 * sin(theta_rad);  % Force in y-direction (N)

% Initialize the global force vector
F = zeros(2 * numNodes, 1);  % Global external force vector

% Apply forces at Node 3
F(2 * 3 - 1) = Fx;  % Apply Fx at Node 3's x DOF
F(2 * 3) = Fy;      % Apply Fy at Node 3's y DOF

% Initialize the global force vector
F = zeros(2 * numNodes, 1);  % Global external force vector

% Apply forces at Node 3
F(2 * 3 - 1) = Fx;  % Apply Fx at Node 3's x DOF
F(2 * 3) = Fy;      % Apply Fy at Node 3's y DOF

% Assemble the global stiffness matrix
for n = 1:size(elements, 1)
    elemNodes = elements(n, :);
    elemCoords = nodes_coordinates(elemNodes, :);
    Ke = zeros(12); % Element stiffness matrix

    % Gauss integration loop
    for gp = 1:size(gaussPoints, 1)
        L1 = gaussPoints(gp, 1);
        L2 = gaussPoints(gp, 2);
        L3 = 1 - L1 - L2;

        % Shape function derivatives
        dNdL = [
            4 * L1 - 1, 0, -3 + 4 * (L1 + L2), 4 * L2,-4 * L2, 4 - 8 * L1 - 4 * L2;
            0, 4 * L2 - 1, -3 + 4 * (L1 + L2), 4 * L1, 4 - 4 * L1 - 8 * L2, -4 * L1
        ];

        % Jacobian matrix
        J = dNdL * elemCoords;
        detJ = det(J);
        dNdXY = J \ dNdL;

        % Strain-displacement matrix B
        B = zeros(3, 12);
        B(1, 1:2:end) = dNdXY(1, :);
        B(2, 2:2:end) = dNdXY(2, :);
        B(3, 1:2:end) = dNdXY(2, :);
        B(3, 2:2:end) = dNdXY(1, :);

        % Accumulate the element stiffness matrix
        Ke = Ke + (B' * D * B) * detJ * gaussWeights(gp);
    end

    % Assemble into the global stiffness matrix
    globalDofs = reshape([2 * elemNodes - 1; 2 * elemNodes], 1, []);
    K(globalDofs, globalDofs) = K(globalDofs, globalDofs) + Ke;
end

% Solve for nodal displacements
fixedDofs = [1, 2, 9, 10]; % Fixed degrees of freedom (x and y at nodes 1 and 5)
freeDofs = setdiff(1:2 * numNodes, fixedDofs);
Kff = K(freeDofs, freeDofs);
Ff = F(freeDofs);

U = zeros(2 * numNodes, 1);
U(freeDofs) = Kff \ Ff;

% Assemble the global stiffness matrix
for n = 1:size(elements, 1)
    elemNodes = elements(n, :);  % Node indices
    elemCoords = nodes_coordinates(elemNodes, :);  % Node coordinates
    Ke = zeros(12);  % Element stiffness matrix (6 nodes x 2 DOF)

    % Gauss integration loop
    for gp = 1:size(gaussPoints, 1)
        L1 = gaussPoints(gp, 1);
        L2 = gaussPoints(gp, 2);
        L3 = 1 - L1 - L2;

        % Shape function derivatives
        dNdL = [
            4 * L1 - 1,        0,                 -3 + 4 * (L1 + L2),   4 * L2,       -4 * L2,         4 - 8 * L1 - 4 * L2;
            0,                 4 * L2 - 1,        -3 + 4 * (L1 + L2),   4 * L1,        4 - 4 * L1 - 8 * L2, -4 * L1
        ];

        % Jacobian matrix
        J = dNdL * elemCoords;
        detJ = det(J);
        dNdXY = J \ dNdL;  % Transform derivatives from local to global coordinates

        % Strain-displacement matrix B
        B = zeros(3, 12);
        B(1, 1:2:end) = dNdXY(1, :);  % dN/dx
        B(2, 2:2:end) = dNdXY(2, :);  % dN/dy
        B(3, 1:2:end) = dNdXY(2, :);  % dN/dy
        B(3, 2:2:end) = dNdXY(1, :);  % dN/dx

        % Accumulate the element stiffness matrix
        Ke = Ke + (B' * D * B) * detJ * gaussWeights(gp);
    end

    % **Assemble the element stiffness matrix into the global stiffness matrix**
    % Corresponding global DOF
    element_DOF = zeros(1, 12);
    for i = 1:6
        element_DOF(2*i-1) = 2 * elemNodes(i) - 1;  % u_x for node i
        element_DOF(2*i)   = 2 * elemNodes(i);      % u_y for node i
    end

    % Assemble into the global stiffness matrix
    for i = 1:12
        for j = 1:12
            K(element_DOF(i), element_DOF(j)) = ...
                K(element_DOF(i), element_DOF(j)) + Ke(i, j);
        end
    end
end

% Define fixed nodes
fixed_nodes = [1, 5];  % Nodes 1 and 5 are fixed

% Find the global DOF of fixed nodes
fixed_DOF = [];
for i = 1:length(fixed_nodes)
    fixed_DOF = [fixed_DOF, 2 * fixed_nodes(i) - 1, 2 * fixed_nodes(i)];
end

% Free DOF indices (remaining DOF)
free_DOF = setdiff(1:2*numNodes, fixed_DOF);

% Exclude fixed DOF and construct the reduced stiffness matrix and force vector
K_reduced = K(free_DOF, free_DOF);
F_reduced = F(free_DOF);

% Solve for displacements
U_reduced = K_reduced \ F_reduced;

% Construct the complete displacement vector
U = zeros(2 * numNodes, 1);
U(free_DOF) = U_reduced;


% Display displacements
disp('Nodal Displacements (m):');
disp('Node   Ux (m)          Uy (m)');
displacements = reshape(U, 2, [])'; 
for i = 1:size(displacements, 1)
    fprintf('%4d   %12.6e   %12.6e\n', i, displacements(i, 1), displacements(i, 2));
end

% Loop over each element to calculate strain at each vertex
for n = 1:size(elements, 1)
    elemNodes = elements(n, :);  % Nodes of the current element
    elemCoords = nodes_coordinates(elemNodes, :);  % Coordinates of the element nodes
    U_e = zeros(12, 1);  % Displacement vector for this element

    % Extract nodal displacements for this element
    for i = 1:6
        U_e(2*i-1) = U(2*elemNodes(i)-1);  % u_x for node i
        U_e(2*i) = U(2*elemNodes(i));      % u_y for node i
    end

    fprintf('Element %d:\n', n);

    % Loop over each vertex (nodes 1, 2, and 3 of the triangular element)
    for vertex = 1:3
        % Local coordinates of the vertex
        L1 = (vertex == 1);  % L1 = 1 for vertex 1
        L2 = (vertex == 2);  % L2 = 1 for vertex 2
        L3 = (vertex == 3);  % L3 = 1 for vertex 3

        % Compute the shape function derivatives at the vertex
        dNdL = [
            4 * L1 - 1,        0,                 -3 + 4 * (L1 + L2),   4 * L2,       -4 * L2,         4 - 8 * L1 - 4 * L2;
            0,                 4 * L2 - 1,        -3 + 4 * (L1 + L2),   4 * L1,        4 - 4 * L1 - 8 * L2, -4 * L1
        ];

        % Compute Jacobian and derivatives in global coordinates
        J = dNdL * elemCoords;  % Jacobian matrix
        detJ = det(J);
        dNdXY = J \ dNdL;  % Transform derivatives to global coordinates

        % Construct the B-matrix (strain-displacement matrix)
        B = zeros(3, 12);
        B(1, 1:2:end) = dNdXY(1, :);  % dN/dx
        B(2, 2:2:end) = dNdXY(2, :);  % dN/dy
        B(3, 1:2:end) = dNdXY(2, :);  % dN/dy
        B(3, 2:2:end) = dNdXY(1, :);  % dN/dx

        % Compute the strain at this vertex
        strain = B * U_e;

        % Display the strain at the current vertex
        fprintf('  Vertex %d (Node %d):\n', vertex, elemNodes(vertex));
        fprintf('    Strain: [Exx: %.6e, Eyy: %.6e, Exy: %.6e]\n', strain(1), strain(2), strain(3));
    end
end