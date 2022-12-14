function [] = vcoeff()
% Purpose: To calculate the coefficients for the v equation.

% constants
global NPI NPJ  Dt NV LARGE XMAX YMAX XMM YMM DMM
% variables
global x x_u y y_v v p mu SP Su F_u F_v d_v relax v_old Istart Iend ...
    Jstart Jend b aE aW aN aS aP dvdy dudy rho_old

Istart = 2;
Iend = NPI+1;
Jstart = 3;
Jend = NPJ+1;

convect();

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;       
        % Geometrical parameters: Areas of the cell faces
        AREAw = y(J) - y(J-1); % See fig. 6.4
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i);
        AREAn = AREAs;
        
        % eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.
        Fw = 0.5*(F_u(i,J)   + F_u(i,J-1))*AREAw;
        Fe = 0.5*(F_u(i+1,J) + F_u(i+1,J-1))*AREAe;
        Fs = 0.5*(F_v(I,j)   + F_v(I,j-1))*AREAs;
        Fn = 0.5*(F_v(I,j)   + F_v(I,j+1))*AREAn;
        
        % eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b
        % note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition
        Dw = 0.25*(mu(I-1,J-1) + mu(I,J-1) + mu(I-1,J) + mu(I,J))/(x(I) - x(I-1))*AREAw;
        De = 0.25*(mu(I,J-1) + mu(I+1,J-1) + mu(I,J) + mu(I+1,J))/(x(I+1) - x(I))*AREAe;
        Ds = mu(I,J-1)/(y_v(j) - y_v(j-1))*AREAs;
        Dn = mu(I,J)/(y_v(j+1) - y_v(j))*AREAn;
        
        % The source terms
        muw = 0.25*(mu(I-1,J-1) + mu(I,J-1) + mu(I-1,J) + mu(I,J));
        mue = 0.25*(mu(I,J-1) + mu(I+1,J-1) + mu(I,J) + mu(I+1,J));
        SP(I,j) = 0.;
        Su(I,j) = 0.;
        Su(I,j) = (mu(I,J)*dvdy(I,J) - mu(I,J-1)*dvdy(I,J-1)) / (y(J) - y(J-1)) + ...
                    (mue*dudy(i+1,j) - muw*dudy(i,j)) / (x_u(i+1) - x_u(i));
        Su(i,J) =  Su(i,J)*AREAw*AREAs;
        
        %marshmellow placement
        if any(i == round(DMM/XMAX*NPI) : round((DMM+XMM)/XMAX*NPI))
            if any(J == round((YMAX/2-YMM/2)/YMAX*NPJ):round((YMAX/2+YMM/2)/YMAX*NPJ))
                SP(i,J) = -LARGE;
                Su(i,J) = 0;
            end
        end
        
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max([ Fw, Dw + Fw/2, 0.]);
        aE(I,j) = max([-Fe, De - Fe/2, 0.]);
        aS(I,j) = max([ Fs, Ds + Fs/2, 0.]);
        aN(I,j) = max([-Fn, Dn - Fn/2, 0.]);
        aPold   = 0.5*(rho_old(I,J-1) + rho_old(I,J))*AREAe*AREAn/Dt;
        
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,j) = aW(I,j) + aE(I,j) + aS(I,j) + aN(I,j) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % Calculation of d(I,j) = d_v(I,j) defined in eq. 6.23 for use in the
        % equation for pression correction (eq. 6.32) (see subroutine pccoeff).
        d_v(I,j) = AREAs*relax(NV)/aP(I,j);
        
        % Putting the integrated pressure gradient into the source term b(I,j)
        % The reason is to get an equation on the generalised form
        % (eq. 7.7 ) to be solved by the TDMA algorithm.
        % note: In reality b = a0p*fiP + Su = 0.
        b(I,j) = (p(I,J-1) - p(I,J))*AREAs + Su(I,j) + aPold*v_old(I,j);
        
        % Introducing relaxation by eq. 6.37 . and putting also the last
        % term on the right side into the source term b(i,J)
        aP(I,j) = aP(I,j)/relax(NV);
        b(I,j)  = b(I,j) + (1.0 - relax(NV))*aP(I,j)*v(I,j);
        
        % now we have implemented eq. 6.37 in the form of eq. 7.7
        % and the TDMA algorithm can be called to solve it. This is done
        % in the next step of the main program.        
    end
end
end

