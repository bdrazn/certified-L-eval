default(parisizemax, 4000000000);
\p 120
E = ellinit([0,0,1,-1,0]);
print("conductor = ", ellglobalred(E)[1]);
print("disc = ", E.disc);
ff = factor(ellglobalred(E)[1]);
n1 = matsize(ff)[1];
Ct = 1;
for(i=1, n1, pp = ff[i,1]; loc = elllocalred(E, pp); print("p=", pp, "  kodaira=", loc[2], "  c_p=", loc[4]); Ct = Ct*loc[4]);
print("TAMAGAWA PRODUCT = ", Ct);
print("torsion = ", elltors(E));
R = ellrank(E);
print("ellrank = ", R);
pts = R[4];
H = ellheightmatrix(E, pts);
Reg = matdet(H);
Omega1 = E.omega[1];
lf = lfuncreate(E);
L0 = lfun(lf, 1);
L1 = lfun(lf, 1, 1);
print("L(1)   = ", L0);
print("L'(1)  = ", L1);
print("Reg = ", Reg, "   omega1 = ", Omega1);
print("Sha_impl (Omega=omega1)   = ", L1/(Omega1*Reg*Ct));
print("Sha_impl (Omega=2*omega1) = ", L1/(2*Omega1*Reg*Ct));
print("CHECK37-END");
quit
