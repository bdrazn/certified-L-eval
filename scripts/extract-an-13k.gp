default(threadsizemax, 1073741824);
E = ellinit([0,1,1,-2,0]);
print("389a1 ", ellglobalred(E)[1], " ", ellrootno(E));
for(i = 1, 12600, print(ellak(E, i)));
E = ellinit([0,0,1,-1,0]);
print("37a1 ", ellglobalred(E)[1], " ", ellrootno(E));
for(i = 1, 12600, print(ellak(E, i)));
E = ellinit([1,-1,1,0,0]);
print("53a1 ", ellglobalred(E)[1], " ", ellrootno(E));
for(i = 1, 12600, print(ellak(E, i)));
E = ellinit([0,-1,1,-2,2]);
print("57b ", ellglobalred(E)[1], " ", ellrootno(E));
for(i = 1, 12600, print(ellak(E, i)));
E = ellinit([0,1,0,-1,0]);
print("20a1 ", ellglobalred(E)[1], " ", ellrootno(E));
for(i = 1, 12600, print(ellak(E, i)));
quit
