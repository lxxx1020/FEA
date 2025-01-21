function z = Triangle2D3Node_Assembly(KK, k, i, j, m)
% Inputs:
%   KK: Current global stiffness matrix
%   k: Stiffness matrix of the current element (6x6)
%   i, j, m: Node numbers of the current element
% Output:
%   z: Updated global stiffness matrix

    % Map the degrees of freedom
    DOF = [
        2*i-1, 2*i, ...
        2*j-1, 2*j, ...
        2*m-1, 2*m
    ];

    % Assemble the stiffness matrix
    for n1 = 1:6
        for n2 = 1:6
            KK(DOF(n1), DOF(n2)) = KK(DOF(n1), DOF(n2)) + k(n1, n2);
        end
    end

    % Return the updated global stiffness matrix
    z = KK;
end
