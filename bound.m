function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ  NT  NV NALPHA
% variables
global u v T y_v F_u Alpha solve_fi 

% top and bottom outlet boundary
for I = 1:NPI
    for J = [1 NPJ+1]
        v(I,J) = v(I,J);
    end
end

% right outlet boundary
for J = 1:NPJ+2
    u(NPI+2,J) = u(NPI+1,J);
end

% wall surrounding inlet boundary
for J = 1:NPJ+2
    u(2,J) = 0;
end



end
