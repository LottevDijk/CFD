function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN XMAX YMAX XMM YMM DMM JinLeft JinRight Cmu Ti
% variables
global x y u v T rho m_in m_out m_outr m_outb m_outt x_u y_v F_u F_v k eps f time Dt

% Set k and eps at the inlet
k(1,1:NPJ+2)     = 1.5*(U_IN*Ti)^2; % at inlet
eps(1,1:NPJ+2)   = Cmu^0.75 *k(1,1:NPJ+2).^1.5/(0.07*(JinRight-JinLeft)*YMAX*0.5); % at inlet

% Fix temperature at the walls in Kelvin
T(1:NPI+2,1) = 273.; % bottom wall
T(1:NPI+2,NPJ+2) = 273.; % top wall


% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet

% end: globcont()==========================================================

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity

% wall surrounding with inlet boundary
for J = 1:NPJ+2
    u(2,J) = 0;
end
for J = round(JinLeft*NPJ):round(JinRight*NPJ)
    u(2,J) = U_IN;
end

%% Conservation of mass

convect();
m_in = 0.;
m_out = 0.;
m_outr = 0;
m_outb = 0;
m_outt = 0;

            %Mass in left and right
for J = 1:NPJ+1
    AREAw = y_v(J+1)-y_v(J);
            %Burner inlet
    m_in = m_in + F_u(2,J)*AREAw;
            %Inlet right boundary
    if u(NPI+2,J) < 0
        m_in = m_in - F_u(NPI+2,J) * AREAw;
    end
end
            %Mass in top and bottom
for I = 1:NPI+1
    AREAn = x_u(I+1) - x_u(I);
            %bottom
    if (v(I,2) > 0)
        m_in = m_in + F_v(I,2) * AREAn;
    end
            %top
    if (v(I,NPJ+2) < 0)
        m_in = m_in - F_v(I,NPJ+1) * AREAn;
    end
end

            % mass out rightside
for J = 1:NPJ+1
    if u(NPI+2,J) > 0
        m_outr = m_outr + F_u(NPI+1,J) * AREAw;
    end
end 
            % mass out top and bottom
for I = 1:NPI+1
    %bottom
    if (v(I,2) < 0)
        m_outb = m_outb - F_v(I,2) * AREAn;
    end
    %top
    if (v(I,NPJ+2) > 0)
        m_outt = m_outt + F_v(I,NPJ+1) * AREAn;
    end
end
m_out = m_outr + m_outb + m_outt;

if m_in > 10.0 || m_out > 10.0
    error('Alles gaat helemaal mis')
end

%% Outlet boundaries

%Outlet top
u(3:NPI,NPJ+2) = u(3:NPI,NPJ+1);
for i = 3:NPI
    if v(i,NPJ+2) > 0
        v(i,NPJ+2) = v(i,NPJ+1);
    else
        v(i,NPJ+2) = v(i,NPJ+1) * m_out / m_in;
    end
end
u(NPI+1,NPJ+2) = 0;
v(NPI+1,NPJ+2) = 0;
%Oulet right
for j = 2:NPJ+1
    if u(NPI+2) > 0
        u(NPI+2,j) = u(NPI+1,j);
    else
        u(NPI+2,j) = u(NPI+1,j) * m_out / m_in;
    end
end
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
u(NPI+2,NPJ+1) = 0;
v(NPI+2,NPJ+1) = 0;
% Outlet bottom
u(3:NPI,2) = 0*u(3:NPI,3);
v(3:NPI,2) = 0*v(3:NPI,3) * m_in / m_out;
u(NPI+1,2) = 0;
v(NPI+1,2) = 0;




k(NPI+2,2:NPJ+1) = k(NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
T(NPI+2,1:NPJ+2) = T(NPI+1,1:NPJ+2);
end
