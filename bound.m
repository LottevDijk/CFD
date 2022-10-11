function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ  NT  NV NALPHA
% variables
global u v T y_v F_u Alpha solve_fi 

for I = 1:NPI
    for J = [1 NPJ+1]
        u(I,J)=u(I,J+1);
        v(I,J)=v(I,J+1);
        Alpha(I,J)=Alpha(I,J+1);
    end
end
end
