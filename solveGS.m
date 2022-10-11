function fi = solveGS(fi,b,aE,aW,aN,aS,aP,ITER,ratio)
% Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel
% variables
global Istart Iend Jstart Jend LARGE SMALL

residual = LARGE;
res0 = SMALL;
iter = 0;

% determine initial summed residual 
for J = Jstart:Jend
    for I = Istart:Iend   
        res0 =  res0 + abs(aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) ...
            + aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J) ...
            - aP(I,J)*fi(I,J));
    end
end

% main inner iteration loop (aka Newton loop) 
while iter<ITER && residual/res0 > ratio
	for I = Istart:Iend
        for J = Jstart:Jend
			fi(I,J) = (aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) ...
            + aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J))/aP(I,J);
        end
    end
    
    % reset residual 
	residual =0.;
    
    for J = Jstart:Jend
        for I = Istart:Iend   
        residual =  residual + abs(aE(I,J)*fi(I+1,J) + aW(I,J)*fi(I-1,J) ...
                     + aN(I,J)*fi(I,J+1) + aS(I,J)*fi(I,J-1) + b(I,J) ...
                     - aP(I,J)*fi(I,J));
        end
    end
    iter = iter+1;
end

end
