function [] = init()
% Purpose: To initilise all parameters.

% constants
global NPI NPJ LARGE XMAX YMAX Jinmin Jinmax Joutmin Joutmax TOTAL_DUMPS...
       NU NV NALPHA NT NPC NMU NRHO NFIMAX TOTAL_TIME Dt DELTA RHOL RHOG DUMPFREQ
% variables
global x x_u y y_v u v pc p T rho mu Gamma Cp b SP Su d_u d_v omega ...
    SMAX SAVG m_in m_out aP aE aW aN aS F_u F_v u_old v_old pc_old T_old  ...
    dudx dudy dvdx dvdy Alpha rho_old Alpha_old usum vsum psum rhosum ...
    alphasum musum umean vmean pmean rhomean alphamean mumean relax solve_fi print_fi NGAMMA

%% begin: memalloc()
% allocate memory for variables
x   = zeros(1,NPI+2);
x_u = zeros(1,NPI+2);
y   = zeros(1,NPJ+2);
y_v = zeros(1,NPJ+2);

u   = zeros(NPI+2,NPJ+2);
v   = zeros(NPI+2,NPJ+2);
pc  = zeros(NPI+2,NPJ+2);
p   = zeros(NPI+2,NPJ+2);
T   = zeros(NPI+2,NPJ+2);
rho = zeros(NPI+2,NPJ+2);
mu  = zeros(NPI+2,NPJ+2);
Gamma = zeros(NPI+2,NPJ+2);
Cp  = zeros(NPI+2,NPJ+2);
Alpha = zeros(NPI+2,NPJ+2);

u_old  = zeros(NPI+2,NPJ+2);
v_old  = zeros(NPI+2,NPJ+2);
pc_old = zeros(NPI+2,NPJ+2);
T_old  = zeros(NPI+2,NPJ+2);
rho_old  = zeros(NPI+2,NPJ+2);
Alpha_old  = zeros(NPI+2,NPJ+2);

usum = zeros(NPI+2,NPJ+2);
vsum = zeros(NPI+2,NPJ+2);
psum = zeros(NPI+2,NPJ+2);
rhosum = zeros(NPI+2,NPJ+2);
musum = zeros(NPI+2,NPJ+2);
alphasum = zeros(NPI+2,NPJ+2);

umean = zeros(NPI+2,NPJ+2);
vmean = zeros(NPI+2,NPJ+2);
pmean = zeros(NPI+2,NPJ+2);
rhomean = zeros(NPI+2,NPJ+2);
mumean = zeros(NPI+2,NPJ+2);
alphamean = zeros(NPI+2,NPJ+2);

dudx   = zeros(NPI+2,NPJ+2);
dudy   = zeros(NPI+2,NPJ+2);
dvdx   = zeros(NPI+2,NPJ+2);
dvdy   = zeros(NPI+2,NPJ+2);

aP  = zeros(NPI+2,NPJ+2);
aE  = zeros(NPI+2,NPJ+2);
aW  = zeros(NPI+2,NPJ+2);
aN  = zeros(NPI+2,NPJ+2);
aS  = zeros(NPI+2,NPJ+2);
b   = zeros(NPI+2,NPJ+2);
SP  = zeros(NPI+2,NPJ+2);
Su  = zeros(NPI+2,NPJ+2);

F_u = zeros(NPI+2,NPJ+2);
F_v = zeros(NPI+2,NPJ+2);
d_u = zeros(NPI+2,NPJ+2);
d_v = zeros(NPI+2,NPJ+2);

relax = zeros(1,NFIMAX+1);
solve_fi = zeros(1,NFIMAX+1);
print_fi = zeros(1,NFIMAX+1);
% end of memory allocation

%% begin: grid()
% Purpose: Defining the geometrical variables See fig. 6.2-6.4 in ref. 1
% Length of volume element
Dx = XMAX/NPI;
Dy = YMAX/NPJ;
DELTA = sqrt(Dx*Dy);

% Length variable for the scalar points in the x direction
x(1) = 0.;
x(2) = 0.5*Dx;
for I = 3:NPI+1
    x(I) = x(I-1) + Dx;
end
x(NPI+2) = x(NPI+1) + 0.5*Dx;

% Length variable for the scalar points T(i,j) in the y direction
y(1) = 0.;
y(2) = 0.5*Dy;
for J = 3:NPJ+1
    y(J) = y(J-1) + Dy;
end
y(NPJ+2) = y(NPJ+1) + 0.5*Dy;

% Length variable for the velocity components u(i,j) in the x direction
x_u(1) = 0.;
x_u(2) = 0.;
for i = 3:NPI+2
    x_u(i) = x_u(i-1) + Dx;
end

% Length variable for the velocity components v(i,j) in the y direction */
y_v(1) = 0.;
y_v(2) = 0.;
for j = 3:NPJ+2
    y_v(j) = y_v(j-1) + Dy;
end
% end of grid setting

%% begin: gridbound()
% Purpose: Definding the geometrical boundary coordinates
Jinmin   = 2;
Jinmax   = NPJ+1;
Joutmin  = 2;
Joutmax  = NPJ+1;
% end of gridbound setting

%% begin: init()
% Initialising all other variables
omega = 1.0; % Over-relaxation factor for SOR solver

% Initialize convergence parameters at large values
SMAX = LARGE;
SAVG = LARGE;

m_in  = 1.;
m_out = 1.;

% calculate dumpfrequency [1/timesteps], minimum dumpfrequency is equal to 1
DUMPFREQ = max(round( TOTAL_TIME/(Dt*TOTAL_DUMPS)), 1);

u(:,:)     = 0.;       % velocity in x-direction
v(:,:)     = 0.;       % Velocity in y-direction
T(:,:)     = 273.;     % Temperature
Alpha(:,:) = 0.;       % Gas fraction
rho(:,:)   = (1 - Alpha)*RHOL + Alpha*RHOG;  % Density: mixture of liquid and gas
pc(:,:)    = 0.;       % Pressure correction (equivalent to p´ in ref. 1). 
mu(:,:)    = 1.E-3;    % Liquid viscosity
Cp(:,:)    = 1013.;    % J/(K*kg) Heat capacity - assumed constant for this problem
Gamma      = 0.0315./Cp; % Thermal conductivity divided by heat capacity
d_u(:,:)   = 0.;       % Variable d to calculate pc defined in 6.23 
d_v(:,:)   = 0.;       % Variable d to calculate pc defined in 6.23 
b(:,:)     = 0.;	   % The general constant 
SP(:,:)    = 0.;       % Source term 
Su(:,:)    = 0.;	   % Source term 
u_old      = u;         % Velocity in x-direction old timestep
v_old      = v;         % Velocity in y-direction old timestep
pc_old     = pc;        % Pressure correction old timestep
T_old      = T;         % Temperature old timestep
rho_old    = rho;       % density old timestep
Alpha_old(:,:)  = 0.;   % gas fraction old timestep

for I = 1: NPI+2
    for J = 1:NPJ+2
        p(I,J)     = rho(I,J)*9.81*(NPI-I)/NPI*XMAX;      % relative pressure
    end
end

if (solve_fi(NU))
% 		 u(Joutmin:Joutmax,NPJ) = U_IN/1000;
         u(NPI,Joutmin:Joutmax) = 0.;
end
% Important to avoid crash!! Othervise m_out calculated in subroutine globcont 
% would be zero at first iteration=>m_in/m_out =INF
		

% Initialising the logical parameter for which variable fi to solve and to print results for 
solve_fi(1:NFIMAX) = 0; 
print_fi(1:NFIMAX) = 1;
solve_fi(NU)       = 1;
solve_fi(NV)       = 1;
solve_fi(NPC)      = 1;
solve_fi(NALPHA)   = 1;
solve_fi(NRHO)     = 1;
solve_fi(NMU)      = 1;
% solve_fi(NT)       = 1;
% solve_fi(NGAMMA)   = 1;


% Setting the relaxation parameters
relax(NU)   = 0.8;              % See eq. 6.36
relax(NV)   = relax(NU);        % See eq. 6.37
relax(NPC)  = 1.1 - relax(NU);  % See eq. 6.33
relax(NT)   = 1.0;              % Relaxation factor for temperature
relax(NRHO) = 0.1;              % Relaxation factor for density
relax(NALPHA) = 0.25;           % Relaxation factor gas volume fraction

% end of initilization
end

