% File: FEA_3nodes_displacement_strain.m
% Description: Finite Element Analysis for a 2D thin plate under plane stress conditions
% Author: [Xuan Li]
% 
% This script calculates nodal displacements and element strains for a
% 2D thin plate using 3-node triangular elements. The results are validated
% under different elastic moduli and loading conditions.



clc; clear;

% Define elastic modulus array (in Pa)
E_array = [40e9, 10e9, 70e9, 200e9]; % Original value with 3 new values
NU = 0.3; % Poisson's ratio
t = 2e-3; % Thickness (in m)
ID = 1;   % Plane stress problem

% Node coordinates
nodes = [
    0, 0;   % Node1
    30, 0;  % Node2
    60, 30; % Node3
    30, 30; % Node4
    0, 30   % Node5
];

% Element definition (node numbers)
elements = [
    1, 4, 5; % Element 1
    1, 2, 4; % Element 2
    2, 3, 4  % Element 3
];

% Total degrees of freedom
numNodes = size(nodes, 1);
totalDOF = 2 * numNodes; % Each node has 2 degrees of freedom

% Fixed degrees of freedom (Node 1 and Node 5 are fixed)
fixedDOF = [1, 2, 9, 10];

% Initialize results storage
results = struct();

for idx = 1:length(E_array)
    % Get the current elastic modulus
    E = E_array(idx);
    fprintf('Current elastic modulus E = %.1e Pa\n', E);

    % Construct elasticity matrix [D]
    if ID == 1  % Plane stress
        D = (E / (1 - NU^2)) * [1, NU, 0; NU, 1, 0; 0, 0, (1 - NU) / 2];
    elseif ID == 2  % Plane strain
        D = (E / ((1 + NU) * (1 - 2 * NU))) * ...
            [1 - NU, NU, 0; NU, 1 - NU, 0; 0, 0, (1 - 2 * NU) / 2];
    end
    D_inv = inv(D); % Calculate the inverse of [D]

    % Calculate and store element stiffness matrices
    k_all = cell(size(elements, 1), 1); % Store stiffness matrix for each element
    for e = 1:size(elements, 1)
        % Get the node numbers of the current element
        i = elements(e, 1);
        j = elements(e, 2);
        m = elements(e, 3);

        % Get the coordinates of the current element's nodes
        xi = nodes(i, 1); yi = nodes(i, 2);
        xj = nodes(j, 1); yj = nodes(j, 2);
        xm = nodes(m, 1); ym = nodes(m, 2);

        % Call external function to calculate stiffness matrix
        k_all{e} = Triangle2D3Node_Stiffness(E, NU, t, xi, yi, xj, yj, xm, ym, ID);
    end

    % Initialize global stiffness matrix
    KK = zeros(totalDOF, totalDOF);

    % Assemble element stiffness matrices
    for e = 1:size(elements, 1)
        i = elements(e, 1);
        j = elements(e, 2);
        m = elements(e, 3);

        KK = Triangle2D3Node_Assembly(KK, k_all{e}, i, j, m);
    end

    % Apply boundary conditions
    freeDOF = setdiff(1:totalDOF, fixedDOF); % Free degrees of freedom
    K_reduced = KK(freeDOF, freeDOF);

    % Define the load vector (in N)
    P = zeros(totalDOF, 1);
    P(6) = -50; % Force in the -y direction at Node 3
    P_reduced = P(freeDOF);

    % Solve for displacements of free degrees of freedom
    U_reduced = K_reduced \ P_reduced;

    % Construct the full displacement vector
    U_total = zeros(totalDOF, 1);
    U_total(freeDOF) = U_reduced;

    % Calculate stress and strain for each element
    strains = zeros(3, size(elements, 1)); % Initialize strain storage
    for e = 1:size(elements, 1)
        % Get the node numbers of the current element
        i = elements(e, 1);
        j = elements(e, 2);
        m = elements(e, 3);

        % Get the coordinates of the current element's nodes
        xi = nodes(i, 1); yi = nodes(i, 2);
        xj = nodes(j, 1); yj = nodes(j, 2);
        xm = nodes(m, 1); ym = nodes(m, 2);

        % Extract the displacement vector of the current element
        DOF = [2*i-1, 2*i, 2*j-1, 2*j, 2*m-1, 2*m];
        u = U_total(DOF);

        % Call external function to calculate stress
        stresses = Triangle2D3Node_Stress(E, NU, xi, yi, xj, yj, xm, ym, u, ID);

        % Calculate strain from stress using D_inv
        strains(:, e) = D_inv * stresses;
    end

    % Store results
    results(idx).E = E;
    results(idx).U_total = U_total;
    results(idx).strains = strains;
end

% Display calculation results
for idx = 1:length(E_array)
    fprintf('Results for elastic modulus E = %.1e Pa:\n', results(idx).E);
    fprintf('Node displacements (in m):\n');
    disp(reshape(results(idx).U_total, 2, []).'); % Reshape to show [Ux, Uy] per node
    fprintf('Element strains:\n');
    for e = 1:size(elements, 1)
        fprintf('Element %d: Exx = %.2e, Eyy = %.2e, Exy = %.2e\n', ...
                e, results(idx).strains(:, e));
    end
end