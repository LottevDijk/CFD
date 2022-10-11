function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global XMAX YMAX XMM YMM DMM U_IN NPI NPJ ALPHAIN IinLeft IinRight NT  NV NALPHA
% variables
global u v T y_v F_u Alpha solve_fi 

% Temperature at the walls in Kelvin
T(1:NPI+2,1) = 273.;     % top wall
T(1:NPI+2,NPJ+2) = 273.; % bottom wall
% T(1,1:NPJ+2) = 373.;     %inlet wall


% top and bottom outlet boundary
for I = 1:NPI
    v(I,1) = v(I,2);
    v(I,NPJ+2) = v(I,NPJ+2);
end

% right outlet boundary
for J = 1:NPJ+2
    u(NPI+2,J) = u(NPI+1,J);
end

% wall surrounding inlet boundary
for J = 1:NPJ+2
    u(2,J) = 0;
end

for J = round(IinLeft*NPJ):round(IinRight*NPJ)
    u(2,J) = U_IN;
end

for J = round(IinLeft*NPJ) : round(IinRight*NPJ)
    Alpha(1,J) = ALPHAIN; %inlet gas fraction left
end

% %marshmellow placement
% for I = round(DMM/XMAX*NPI):round((DMM+XMM)/XMAX*NPI)
%     u(I,:) = 0;
% end



end
