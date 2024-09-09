(* ::Package:: *)

AHEMultipoleResolve[Sigma_] :=Module[{px,py,pz,C1, C2, C3, C4, fz3, fxz2, fyz2, fzx2y2, fxyz, fxx23y2, fy3x2y2,orbtal, orbtalname, 
SigmaInput, SHfactor,mx, my, mz}, 
px=Sqrt[3/(4Pi)]Sin[\[Theta]] Cos[\[Phi]];
py=Sqrt[3/(4Pi)]Sin[\[Theta]] Sin[\[Phi]];
pz=Sqrt[3/(4Pi)]Cos[\[Theta]];
C1=Sqrt[7/(16Pi)];  C2=Sqrt[21/(32Pi)];
C3=Sqrt[105/(16Pi)];C4=Sqrt[35/(32Pi)];
fz3=C1*(5Cos[\[Theta]]^3-3Cos[\[Theta]]);
fxz2=C2*(5Cos[\[Theta]]^2-1) Sin[\[Theta]] Cos[\[Phi]];
fyz2=C2*(5Cos[\[Theta]]^2-1) Sin[\[Theta]] Sin[\[Phi]];
fzx2y2=C3*Sin[\[Theta]]^2*Cos[\[Theta]]Cos[2 \[Phi]];
fxyz=C3*Sin[\[Theta]]^2*Cos[\[Theta]]Sin[2 \[Phi]];
fxx23y2=C4*Sin[\[Theta]]^3*(Cos[\[Phi]]^2-3 Sin[\[Phi]]^2)Cos[\[Phi]];
fy3x2y2=C4*Sin[\[Theta]]^3*(3 *Cos[\[Phi]]^2- Sin[\[Phi]]^2)Sin[\[Phi]];
orbtal={pz,px,py,fz3,fxz2, fyz2, fzx2y2, fxyz, fxx23y2,fy3x2y2};
orbtalname={"\!\(\*SubscriptBox[\(p\), \(z\)]\)","\!\(\*SubscriptBox[\(p\), \(x\)]\)","\!\(\*SubscriptBox[\(p\), \(y\)]\)",
"\!\(\*SubscriptBox[\(f\), \(z3\)]\)","\!\(\*SubscriptBox[\(f\), SuperscriptBox[\(xz\), \(2\)]]\)", 
"\!\(\*SubscriptBox[\(f\), SuperscriptBox[\(yz\), \(2\)]]\)", "\!\(\*SubscriptBox[\(f\), \(z \((\*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)])\)\)]\)",
 "\!\(\*SubscriptBox[\(f\), \(xyz\)]\)", "\!\(\*SubscriptBox[\(f\), \(x \((\*SuperscriptBox[\(x\), \(2\)] - 3 \*SuperscriptBox[\(y\), \(2\)])\)\)]\)",
 "\!\(\*SubscriptBox[\(f\), \(y \((3 \*SuperscriptBox[\(x\), \(2\)] - \*SuperscriptBox[\(y\), \(2\)])\)\)]\)"};

SigmaInput=Sigma/.{Subscript[n, x]->Sin[\[Theta]]Cos[\[Phi]],  Subscript[n, y] -> Sin[\[Theta]]Sin[\[Phi]], Subscript[n, z]-> Cos[\[Theta]]};

(*Print["The tested function is: ",SigmaInput//MatrixForm]*);  (*For test*)
SHfactor=ConstantArray[0,{3,10}];
For [k=1,k<=3,k++,
  For[m=1,m<=10,m++,
  factor=Integrate[SigmaInput[[k]]*orbtal[[m]]*Sin[\[Theta]], {\[Theta], 0, Pi},{\[Phi], 0, 2 Pi}];
  SHfactor[[k,m]]=factor
  (*Print[SHfactor[[k,m]],"*" orbtalname[[m]]]*)
]];
(*Print[SHfactor//MatrixForm]*) (*For test*);
(*Output the orbits*)
Print["sigma_x= ", Total[SHfactor[[1,;;]]*orbtalname]];
Print["sigma_y= ", Total[SHfactor[[2,;;]]*orbtalname]];
Print["sigma_z= ", Total[SHfactor[[3,;;]]*orbtalname]];
]
