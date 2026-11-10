default(parisizemax, 8000000000);
\p 250
E = ellinit([0,1,1,-2,0]);
ff = factor(389);
n1 = matsize(ff)[1];
for(i=1, n1, p = ff[i,1]; loc = elllocalred(E, p); print("p=", p, "  kodaira=", loc[2], "  c_p=", loc[4]));
Omega1 = E.omega[1];
H = ellheightmatrix(E, [[0,0],[-1,1]]);
Reg = matdet(H);
lf = lfuncreate(E);
L2 = lfun(lf, 1, 2);
print("L2_250=", L2);
print("Omega1_250=", Omega1);
print("Reg_250=", Reg);
print("Sha250 = ", L2/(4*Omega1*Reg));
print("CHECK250-END");
quit
