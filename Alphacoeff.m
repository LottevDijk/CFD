function [] = Alphacoeff()
% Purpose: To calculate the coefficients for the Alpha equation.

% constants
global NPI NPJ  Dt NALPHA RHOL USLIP IinRight IinLeft ImarshLeft ImasrhRight Jmarshlow Jmarshhigh U_IN LARGE XMAX YMAX XMM YMM DMM LARGE
% variables
global  x_u u v y_v SP Su  relax Alpha_old Alpha Istart Iend  rho ...
    Jstart Jend b aE aW aN aS aP 

Istart = 2;
Iend = NPI+1;
Jstart = 2;
Jend = NPJ+1;

alphae = 0.;
alphaw = 0.;
alphan = 0.;
alphas = 0.;

for I = Istart:Iend
    i = I;
    for J = Jstart:Jend
        j = J;
        % Geometrical parameters
        % Areas of the cell faces
        AREAw = y_v(j+1) - y_v(j); % = A(i,J) See fig. 6.2 or fig. 6.5
        AREAe = AREAw;
        AREAs = x_u(i+1) - x_u(i); % = A(I,j)
        AREAn = AREAs;
        
        % The convective mass flux defined in eq. 5.8a
        % note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition.    
          Fw = (u(i,J)+ USLIP*RHOL*(2. - Alpha(I,J) - Alpha(I-1,J))/(rho(I-1,J) + rho(I,J)))*AREAw;
          Fe = (u(i+1,J)+ USLIP*RHOL*(2. - Alpha(I,J) - Alpha(I+1,J))/(rho(I+1,J) + rho(I,J)))*AREAe;
          Fs = v(I,j)*AREAs;
          Fn = v(I,j+1)*AREAn;
        
%          if I == 2 && J>=round(IinLeft*NPJ) && J<= round(IinRight*NPJ)
%               Fw = U_IN*AREAw;
%          end

        % The source terms
        SP(I,J) = 0.;
        if (i >ImarshLeft&& I<ImarshRight && J > Jmarshlow &&J<Jmarshhigh) 
            SP(i,J) = -LARGE;
        end 
        Su(I,J) = 0.;
        
        %marshmellow placement
%         if any(i == round(DMM/XMAX*NPI) : round((DMM+XMM)/XMAX*NPI))
%             if any(J == round((YMAX/2-YMM/2)/YMAX*NPJ):round((YMAX/2+YMM/2)/YMAX*NPJ))
%                 SP(i,J) = -LARGE;
%                 Su(i,J) = 0;
%             end
%         end
        
        % The coefficients (hybrid differencing scheme)
        aW(I,j) = max(Fw, 0.);
        aE(I,j) = max(-Fe, 0.);
        aS(I,j) = max(Fs, 0.);
        aN(I,j) = max(-Fn, 0.);       
        aPold   = AREAe*AREAn/Dt;
        
   	%deferred correction alphaw (Barton minus Upwind) 
    if I>Istart
        if u(i,J) >= 0.
            alpha1 = 0.5*(Alpha(I-1,J)-Alpha(I-2,J));
            alpha2 = 0.5*(Alpha(I,J)  -Alpha(I-1,J));
            alpha3 = 0.;
            if Alpha(I,J) <= Alpha(I-1,J)
                alphaw = min(alpha3, max(alpha1, alpha2));
            else
                alphaw = max(alpha3, min(alpha1, alpha2));
            end
        else
            alpha1 = 0.5*(Alpha(I,J)  -Alpha(I+1,J));
            alpha2 = 0.5*(Alpha(I-1,J)-Alpha(I,J));
            alpha3 = 0.;
            if Alpha(I,J) <= Alpha(I-1,J)
                alphaw = max(alpha3, min(alpha1, alpha2));
            else
                alphaw = min(alpha3, max(alpha1, alpha2));
            end
        end
    end
    
    %deferred correction alphae (Barton minus Upwind) 
    if I<Iend
        if u(i+1,J) >= 0.
            alpha1 = 0.5*(Alpha(I,J)-Alpha(I-1,J));
            alpha2 = 0.5*(Alpha(I+1,J)-Alpha(I,J));
            alpha3 = 0.;
            if Alpha(I+1,J) <= Alpha(I,J)
                alphae = min(alpha3, max(alpha1, alpha2));
            else
                alphae = max(alpha3, min(alpha1, alpha2));
            end
        else
            alpha1 = 0.5*(Alpha(I+1,J)-Alpha(I+2,J));
            alpha2 = 0.5*(Alpha(I,J)-Alpha(I+1,J));
            alpha3 = 0.;
            if Alpha(I+1,J) <= Alpha(I,J)
                alphae = max(alpha3, min(alpha1, alpha2));
            else
                alphae = min(alpha3, max(alpha1, alpha2));
            end
        end
    end
    
    %deferred correction alphas (Barton minus Upwind) 
    if J>Jstart
        if v(I,j) >= 0.
            alpha1 = 0.5*(Alpha(I,J-1)-Alpha(I,J-2));
            alpha2 = 0.5*(Alpha(I,J)  -Alpha(I,J-1));
            alpha3 = 0.;
            if Alpha(I,J) <= Alpha(I,J-1)
                alphas = min(alpha3, max(alpha1, alpha2));
            else
                alphas = max(alpha3, min(alpha1, alpha2));
            end
        else
            alpha1 = 0.5*(Alpha(I,J)  -Alpha(I,J+1));
            alpha2 = 0.5*(Alpha(I,J-1)-Alpha(I,J));
            alpha3 = 0.;
            if Alpha(I,J) <= Alpha(I,J-1)
                alphas = max(alpha3, min(alpha1, alpha2));
            else
                alphas = min(alpha3, max(alpha1, alpha2));
            end
        end
    end
        
     %deferred correction alphan (Barton minus Upwind) 
    if J<Jend
        if v(I,j+1) >= 0.
            alpha1 = 0.5*(Alpha(I,J)-Alpha(I,J-1));
            alpha2 = 0.5*(Alpha(I,J+1)-Alpha(I,J));
            alpha3 = 0.;
            if Alpha(I,J+1) <= Alpha(I,J)
                alphan = min(alpha3, max(alpha1, alpha2));
            else
                alphan = max(alpha3, min(alpha1, alpha2));
            end
        else
            alpha1 = 0.5*(Alpha(I,J+1)-Alpha(I,J+2));
            alpha2 = 0.5*(Alpha(I,J)-Alpha(I,J+1));
            alpha3 = 0.;
            if Alpha(I,J+1) <= Alpha(I,J)
                alphan = max(alpha3, min(alpha1, alpha2));
            else
                alphan = min(alpha3, max(alpha1, alpha2));
            end
        end
    end
    
        % eq. 8.31 without time dependent terms (see also eq. 5.14):
        aP(I,J) = aW(I,J) + aE(I,J) + aS(I,J) + aN(I,J) + Fe - Fw + Fn - Fs - SP(I,J) + aPold;
        
        % Setting the source term equal to b        
        b(I,J) = Su(I,J) + aPold*Alpha_old(I,J);
        b(I,J) = b(I,J) - (Fe*alphae - Fw*alphaw + Fn*alphan - Fs*alphas);
        
        % Introducing relaxation by eq. 6.36 . and putting also the last
        % term on the right side into the source term b(i,J)
        
        aP(I,J) = aP(I,J)/relax(NALPHA);
        b(I,J)  = b(I,J) + (1.0 - relax(NALPHA))*aP(I,J)*Alpha(I,J);
        
        % now the TDMA algorithm can be called to solve the equation.
        % This is done in the next step of the main program.        
    end
end
end

