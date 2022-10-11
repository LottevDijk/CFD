%% Solves: Unsteady, compressible convection-diffusion problems.
% Description:
% This program solves unsteady convection-diffusion problems
% using the transient simple algorithm described in ch. 8.7.1 in "Computational
% Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
% variables follow exactly the notations in this reference, and all
% equations cited are from this reference unless mentioned otherwise.

% Converted from C to Matlab by YTANG
% References: 1. Computational Fluid Dynamics, H.K. Versteeg and W. Malalasekera, Longman Group Ltd, 1995

clear
close all
clc
%% declare all variables and contants
% variables
global  x y u v pc T rho Gamma b SMAX SAVG aP aE aW aN aS Dt Alpha rho_old Cp  current_time NDt d_u d_v
% constants
global NPI NPJ XMAX YMAX XMM YMM DMM LARGE U_IN SMALL Cdrag RHOL RHOG USLIP ALPHAIN CZERO BIT CMUBIT DBUB ...
    NU NV NPC NT NALPHA NP NMU NRHO NGAMMA NFIMAX TOTAL_TIME TOTAL_DUMPS solve_fi IinLeft IinRight


NPI        = 50;       % number of grid cells in x-direction [-]
NPJ        = 50;        % number of grid cells in y-direction [-]
XMAX       = 0.50;       % width of the domain [m]
YMAX       = 0.50;       % height of the domain [m]
XMM        = 0.05;      % width of marshmellow [m]
YMM        = 0.03;      % height of marshmellow [m]
DMM        = 0.20;      % distance of marshmellow to table [m]
MAX_ITER   = 1;        % maximum number of outer iterations [-]
U_ITER     = 5;        % number of Newton iterations for u equation [-]
V_ITER     = 5;         % number of Newton iterations for v equation [-]
PC_ITER    = 30;        % number of Newton iterations for pc equation [-]
T_ITER     = 5;         % number of Newton iterations for T equation [-]
ALPHA_ITER = 5;         % number of Newton iterations for T equation [-]
SMAXneeded = 1E-5;      % maximum accepted error in mass balance [kg/s]
SAVGneeded = 1E-6;      % maximum accepted average error in mass balance [kg/s]
LARGE      = 1E30;      % arbitrary very large value [-]
SMALL      = 1E-30;     % arbitrary very small value [-]
P_ATM      = 101000.;   % athmospheric pressure [Pa]
U_IN       = 0.4;      % in flow velocity [m/s]
Cdrag      = 50E3;

RHOL       = 1000.;      % density of liquid [kg/m3]
RHOG       = 1.2;        % density of gas [kg/m3]
USLIP      = 0.2;        % slip velocity [m/s]
ALPHAIN    = 0.05;        % gas fraction at inlet [-]
IinLeft    = 0.45;       % relative location left side inlet
IinRight   = 0.55;       % relative location right side inlet
CZERO      = 0.1;
BIT        = 1;
CMUBIT     = 0.5;
DBUB       = 0.002;      % bubble diameter
TOTAL_DUMPS = 400;       % total number of dumps

NU         = 1;
NV         = 2;
NPC        = 3;
NP         = 4;
NT         = 5;
NRHO       = 6;
NMU        = 7;
NGAMMA     = 8;
NALPHA     = 9;
NFIMAX     = 9;

NPRINT     = 10;
Dt         = 0.01;      % time step for flowsolver
TOTAL_TIME = 10.;        % total time

%% start main function here
init(); % initialization
bound(); % set boundary values that remain untouched during the iteration process
NDt = 0;
for current_time = Dt:Dt:TOTAL_TIME
    iter = 0;
    NDt = NDt+1;
    
    % outer iteration loop
    while iter < MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded
        
        if solve_fi(NU) || solve_fi(NV)
            derivatives();
        end
        if solve_fi(NU)
            ucoeff();            
            u = solveGS(u, b, aE, aW, aN, aS, aP, U_ITER, 0.25);
        end
        
        if solve_fi(NV)
            vcoeff();            
            v = solveGS(v, b, aE, aW, aN, aS, aP, V_ITER, 0.25);
        end
        
        if solve_fi(NPC)
            bound();
            pccoeff();            
            pc = solve(pc, b, aE, aW, aN, aS, aP, PC_ITER, 0.1);
        end
        
        velcorr(); % Correct pressure and velocity
        
        if solve_fi(NT)
            Tcoeff();
            T = solveGS(T, b, aE, aW, aN, aS, aP, T_ITER, 0.25);
        end
        
        if solve_fi(NALPHA)
            Alphacoeff();            
            Alpha = solveGS(Alpha, b, aE, aW, aN, aS, aP, ALPHA_ITER, 0.1);
        end
        
        % Calculate the density rho(I, J) in the fluid as a function of the ideal gas law.
        % Note: rho at the walls are not needed in this case, and therefore not calculated.
        if solve_fi(NRHO)
            rho_old = rho;
            rho = (1.0-Alpha)*RHOL + Alpha*RHOG;
        end
        
        if solve_fi(NMU)       
            viscosity();
        end
        
        %Calculate the thermal conductivity in the fluid as a function of temperature.
        % Max error in the actual temperature interval is 0.9%
        if solve_fi(NGAMMA)
            for I = 1:NPI+2
                for J=1:NPJ+2
                    Gamma(I,J) = (6.1E-5*T(I,J) + 8.4E-3)/Cp(I,J);
                    if Gamma(I,J) < 0.
                        output();
                        fprintf("Error: Gamma(%d,%d) = %e\n", I, J, Gamma(I,J));
                    end
                end
            end
        end
        
        bound();
        storeresults(); % Store data at current time level in arrays for "old" data
        calcmean();     %Calculate time averaged data fields
        
        iter = iter +1; % increase iteration number
    end  % end of while loop (outer interation)
%       
%             % Store all results in output.txt
%             fp = fopen('output.txt','w');
%             fprintf(fp, "%11.4e\t%12.4e\t%12.4e\t%12.4e\n",current_time,...
%                 u(round(NPI/2),round(NPJ/2)),v(round(NPI/2),round(NPJ/2)), Alpha(round(NPI/2),round(NPJ/2)));
%             fclose(fp);

%Make snaps allong time
if mod(current_time,.5) <= 0.01 
    Frame = figure (1)
    pcolor(x,y,Alpha');
    colorbar;
%     quiver(x,y,u',v');
    fileLocation = "D:\MW courses\1 - 4RM00 - Introduction to computational fluid dynamics\wc5\Images";
    fileName = sprintf('image%d.png',current_time);
    plotName = strcat(fileLocation,fileName); 
    print(Frame,plotName,'-dpng','-r200');  % Save as png file (raster format) '-r200' controls its resolution
end

    % begin: printConv(time,iter)==========================================
    % print convergence to the screen
    if current_time == Dt
        fprintf ('Iter\tTime\t tu(%d,%d)\t tv(%d,%d)\t tAlpha(%d,%d) \t SMAX \t SAVG\n', ...
            round(NPI/2),round(NPJ/2),round(NPI/2),round(NPJ/2),round(NPI/2),round(NPJ/2));
    end
    if rem(NDt, NPRINT) == 0
        du = d_u(round(NPI/2),round(NPJ/2))*(pc(round(NPI/2)-1,round(NPJ/2))-pc(round(NPI/2),round(NPJ/2)));
        dv = d_v(round(NPI/2),round(NPJ/2))*(pc(round(NPI/2),round(NPJ/2)-1)-pc(round(NPI/2),round(NPJ/2)));
        
        fprintf ("%4d\t%11.4e\t%12.4e\t%12.4e\t%12.4e\t%10.2e\t%10.2e\n",iter,current_time,u(round(NPI/2),round(NPJ/2)),...
            v(round(NPI/2),round(NPJ/2)), Alpha(round(NPI/2),round(NPJ/2)), SMAX, SAVG);
    end
    % end: printConv(time, iter)===========================================
    % reset SMAX and SAVG
    SMAX = LARGE;
    SAVG = LARGE;
%     if rem(NDt, 100) == 0
%         f = figure('visible', 'off');
%         pcolor(x,y,Alpha');
%         saveas(f,sprintf('Alpha_%d.png',NDt));
%         close(f);
%     end
end  % end of calculation
%%
figure(1)
quiver(x,y,u',v');
figure(2)
pcolor(x,y,Alpha');
figure(3)
plot(y,u(NPI+2,:))
% figure(3)
% pcolor(x,y,T');

% %% begin: output()
% % Save all results in output.txt
% fp = fopen('output.txt','w');
% for I = 1:NPI+1
%     i = I;
%     for J = 2:NPJ+1
%         j = J;
%         ugrid = 0.5*(u(i,J)+u(i+1,J)); % interpolated horizontal velocity
%         vgrid = 0.5*(v(I,j)+v(I,j+1)); % interpolated vertical velocity
%         fprintf(fp,'%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\t%11.5e\n',...
%             x(I), y(J), ugrid, vgrid, pc(I,J), T(I,J), rho(I,J), mu(I,J), Gamma(I,J), ...
%             k(I,J), eps(I,J), uplus(I,J), yplus(I,J), yplus1(I,J), yplus2(I,J));
%     end
%     fprintf(fp, '\n');
% end
% fclose(fp);
%
% % Save vorticity in vort.txt
% vort = fopen('vort.txt','w');
% for I = 2:NPI+1
%     i = I;
%     for J = 2:NPJ+1
%         j = J;
%         vorticity = (u(i,J) - u(i,J-1)) / (y(J) - y(J-1)) - (v(I,j) - v(I-1,j)) / (x(I) - x(I-1));
%         fprintf(vort, '%11.5e\t%11.5e\t%11.5e\n',x(I), y(J), vorticity);
%     end
%     fprintf(vort,'\n');
% end
% fclose(vort);
%
% % Save streamlines in str.txt
% str = fopen('str.txt', 'w');
% for I = 1:NPI+1
%     i = I;
%     for J = 1:NPJ+1
%         j = J;
%         stream = -0.5*(v(I+1,j)+v(I,j))*(x(I+1)-x(I))+0.5*(u(i,J+1)+u(i,J))*(y(J+1)-y(J));
%         fprintf(str, '%11.5e\t%11.5e\t%11.5e\n',x(I), y(J), stream);
%     end
%     fprintf(str,'\n');
% end
% fclose(str);
%
% % Save horizontal velocity components in velu.txt
% velu = fopen('velu.txt','w');
% for I = 2:NPI+2
%     i = I;
%     for J = 1:NPJ+2
%         fprintf(velu, '%11.5e\t%11.5e\t%11.5e\n',x_u(i), y(J), u(i,J));
%     end
%     fprintf(velu, '\n');
% end
% fclose(velu);
%
% % Save vertical velocity components in velv.txt
% velv = fopen('velv.txt','w');
% for I = 1:NPI+2
%     for J = 2:NPJ+2
%         j = J;
%         fprintf(velv, '%11.5e\t%11.5e\t%11.5e\n',x(I), y_v(j), v(I,j));
%     end
%     fprintf(velv,'\n');
% end
% fclose(velv);
% % end output()


