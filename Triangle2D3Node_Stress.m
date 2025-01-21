function stress = Triangle2D3Node_Stress(E, NU, xi, yi, xj, yj, xm, ym, u, ID)
    % Calculate the area of the element
    A = (xi * (yj - ym) + xj * (ym - yi) + xm * (yi - yj)) / 2;

    % Compute the coefficients of the shape function derivatives
    betai = yj - ym; betaj = ym - yi; betam = yi - yj;
    gammai = xm - xj; gammaj = xi - xm; gammam = xj - xi;

    % Construct the B matrix
    B = [
        betai, 0,     betaj, 0,     betam, 0;
        0,     gammai, 0,     gammaj, 0,     gammam;
        gammai, betai, gammaj, betaj, gammam, betam
    ] / (2 * A);

    % Construct the D matrix
    if ID == 1 % Plane stress
        D = (E / (1 - NU^2)) * [1, NU, 0; NU, 1, 0; 0, 0, (1 - NU) / 2];
    elseif ID == 2 % Plane strain
        D = (E / ((1 + NU) * (1 - 2 * NU))) * ...
            [1 - NU, NU, 0; NU, 1 - NU, 0; 0, 0, (1 - 2 * NU) / 2];
    else
        error('Invalid ID. Use 1 for Plane Stress or 2 for Plane Strain.');
    end

    % Calculate the stress
    stress = D * B * u;
end
