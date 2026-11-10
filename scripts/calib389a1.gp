default(parisizemax, 4000000000);
\p 150
E = ellinit([0,1,1,-2,0]);
red = ellglobalred(E);
print("conductor = ", red[1]);
print("disc = ", E.disc);
ff = factor(red[1]);
Ct = 1;
for(i=1, matsize(ff)[1],
  p = ff[i,1];
  loc = elllocalred(E, p);
  print("bad p=", p, "  f=", loc[1], "  kodaira=", loc[2], "  c_p=", loc[4]);
  Ct = Ct * loc[4];
);
print("TAMAGAWA PRODUCT (elllocalred) = ", Ct);
print("elltamagawa(E) = ", elltamagawa(E));
print("omega member = ", iferr(E.omega, ER, "NO_OMEGA_MEMBER"));
Omega1 = iferr(E.omega[1], ER, 0);
print("Omega1 = ", Omega1);
lf = lfuncreate(E);
L0 = lfun(lf, 1);
L1 = lfun(lf, 1, 1);
L2 = lfun(lf, 1, 2);
print("L(1)   = ", L0);
print("L'(1)  = ", L1);
print("L''(1) = ", L2);
H = ellheightmatrix(E, [[0,0],[-1,1]]);
Reg = matdet(H);
print("REG = ", Reg);
if(Omega1 != 0,
  print("Sha_impl (Omega=omega1)   = ", L2/(2*Omega1*Reg*Ct));
  print("Sha_impl (Omega=omega1/2) = ", L2/(1*Omega1*Reg*Ct));
  print("Sha_impl (Omega=2*omega1) = ", L2/(4*Omega1*Reg*Ct));
, print("no omega member: Omega_est(Sha=1) = L''/(2*Reg*Ct) = ", L2/(2*Reg*Ct)));
\p 250
lf = lfuncreate(E);
L2b = lfun(lf, 1, 2);
H2 = ellheightmatrix(E, [[0,0],[-1,1]]);
Reg2 = matdet(H2);
Omega2 = iferr(E.omega[1], ER, 0);
print("L''(1) at 250 = ", L2b);
print("REG at 250 = ", Reg2);
if(Omega2 != 0,
  print("Sha_impl at 250 (Omega=omega1)   = ", L2b/(2*Omega2*Reg2*Ct));
  print("Sha_impl at 250 (Omega=omega1/2) = ", L2b/(1*Omega2*Reg2*Ct));
  print("Sha_impl at 250 (Omega=2*omega1) = ", L2b/(4*Omega2*Reg2*Ct));
, print("no omega member at 250"));
print("CALIBRATION-END");
quit
