function [] = viscosity()
% Purpose: Calculate the viscosity in the fluid as a function of temperature. 
        % Max error in the actual temp. interval is 0.5%
% constants
global NPI NPJ  DELTA CZERO RHOL USLIP DBUB CMUBIT BIT 
% variables
global RHOG  mu  dudy dvdx dudx dvdy Alpha 

for I = 1:NPI+1
    i = I;
    for J=1:NPJ+1
        j = J;
        crossterms = 0.25*(dudy(i,j) + dudy(i+1,j) + dudy(i,j+1) + dudy(i+1,j+1)...
                      + dvdx(i,j) + dvdx(i+1,j) + dvdx(i,j+1) + dvdx(i+1,j+1));
        EijEij  = 2.0*(dudx(I,J)*dudx(I,J) + dvdy(I,J)*dvdy(I,J)) + crossterms*crossterms;
        mu(I,J) = (1.0-Alpha(I,J))*1.e-3 + Alpha(I,J)*1.e-5 + (CZERO*DELTA)^2*sqrt(EijEij)*RHOL;
        if BIT 
             mu(I,J) = mu(I,J)+ CMUBIT * Alpha(I,J) * DBUB * USLIP * RHOG;
        end
     end
end 

end