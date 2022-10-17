function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN XMAX YMAX XMM YMM DMM JinLeft JinRight Cmu Ti
% variables
global y u v T m_in m_out y_v F_u k eps f time Dt

% Set k and eps at the inlet
k(1,1:NPJ+2)     = 1.5*(U_IN*Ti)^2; % at inlet
eps(1,1:NPJ+2)   = Cmu^0.75 *k(1,1:NPJ+2).^1.5/(0.07*(JinRight-JinLeft)*YMAX*0.5); % at inlet

% Fix temperature at the walls in Kelvin
T(1:NPI+2,1) = 273.; % bottom wall
T(1:NPI+2,NPJ+2) = 273.; % top wall


% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;
    for J = 2:NPJ+1
        m_in = 0;
        m_out = 0;
        j = J;
        AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
        m_in  = m_in  + F_u(2,J)*AREAw;
        m_out = m_out + F_u(NPI+1,J)*AREAw;
    end
% end: globcont()==========================================================

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity

%Outlet top
u(1:NPI+2,NPJ+2) = u(1:NPI+2,NPJ+1);
v(1:NPI+2,NPJ+2) = v(1:NPI+2,NPJ+1);
%Oulet right
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1);
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
% Outlet bottom
u(1:NPI+2,2) = u(1:NPI+2,3);
v(1:NPI+2,2) = v(1:NPI+2,3);

% wall surrounding with inlet boundary
for J = 1:NPJ+2
    u(2,J) = 0;
end
for J = round(JinLeft*NPJ):round(JinRight*NPJ)
    u(2,J) = U_IN;
end

k(NPI+2,2:NPJ+1) = k(NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
T(NPI+2,1:NPJ+2) = T(NPI+1,1:NPJ+2);
end
