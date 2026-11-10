default(parisizemax, 4000000000);
default(threadsizemax, 1073741824);
\\ Stability check: the ladder's top rung (10,000 digits) re-run WITH threadsizemax set,
\\ after the lfun thread-stack overflow (8MB default) seen in the first pilot.
\p 10000
E = ellinit([0,1,1,-2,0]);
lf = lfuncreate(E);
L2 = lfun(lf, 1, 2);
print("L2_10000 = ", L2);
print("LADDER-DONE");
quit
