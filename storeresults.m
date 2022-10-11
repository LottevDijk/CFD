function [] = storeresults()
% Purpose: To store data at current time level in arrays for "old" data

% constants
global solve_fi NU NV NPC NT NALPHA NPI NPJ
% variables
global  u u_old v v_old pc pc_old T T_old Alpha Alpha_old

if solve_fi(NU)
    u_old(3:NPI+1,2:NPJ+1)   = u(3:NPI+1,2:NPJ+1);
end
if solve_fi(NV)
    v_old(2:NPI+1,3:NPJ+1)   = v(2:NPI+1,3:NPJ+1);
end
if solve_fi(NPC)
    pc_old(2:NPI+1,2:NPJ+1)  = pc(2:NPI+1,2:NPJ+1);
end
if solve_fi(NT)
    T_old(2:NPI+1,2:NPJ+1)   = T(2:NPI+1,2:NPJ+1);
end
if solve_fi(NALPHA)
    Alpha_old(2:NPI+1,2:NPJ+1) = Alpha(2:NPI+1,2:NPJ+1);
end
end

