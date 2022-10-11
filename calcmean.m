function [] = calcmean()
% Purpose: To Creating time averaged result table.

% constants
global NPI NPJ Dt NDt NZERO DUMPFREQ
% variables
global x y u v  Alpha usum vsum rhosum psum musum alphasum p rho mu ...
    umean vmean pmean rhomean mumean alphamean current_time

Nzero = 10./Dt - 1;
 if NDt > NZERO
     for I=2: NPI+1
         i=I;
         for J=2:NPJ+1
           j=J;
           ugrid = 0.5*(u(i,J)+u(i+1,J));
           vgrid = 0.5*(v(I,j)+v(I,j+1));
			usum(I,J)    = usum(I,J) + ugrid;
			vsum(I,J)    = vsum(I,J) + vgrid;
			psum(I,J)    = psum(I,J) + p(I,J);
			rhosum(I,J)  = rhosum(I,J) + rho(I,J);
			musum(I,J)   = musum(I,J) + mu(I,J);
			alphasum(I,J)= alphasum(I,J) + Alpha(I,J);

			umean(I,J)     = usum(I,J)/(NDt-Nzero);
			vmean(I,J)     = vsum(I,J)/(NDt-Nzero);
			pmean(I,J)     = psum(I,J)/(NDt-Nzero);
			rhomean(I,J)   = rhosum(I,J)/(NDt-Nzero);
			mumean(I,J)    = musum(I,J)/(NDt-Nzero);
			alphamean(I,J) = alphasum(I,J)/(NDt-Nzero);
         end
     end
     
     if(mod(NDt,DUMPFREQ) == 0)
        %meandata = path outputdir + meandata.dat */
        fp = fopen('meandata.dat','w');
        fprintf(fp, "#Time  %11.4e\n",current_time);
        fprintf(fp, "#x\t\ty\t\tugrid\t\tvgrid\t\tp\t\trho\t\tmu\t\tAlpha\n");
        for I=2:NPI+1
            for J=2:NPJ+1
                fprintf(fp, "%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\t%11.4e\n",...
				      	x(I), y(J), umean(I,J), vmean(I,J), pmean(I,J), rhomean(I,J), mumean(I,J), alphamean(I,J));
            end
        end
        fclose(fp);
     end
 end
end
	