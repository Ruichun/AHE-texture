(* ::Package:: *)

MSGdata[SPNo_] := Module[{MSGDATA, MSGOP, gray,braLatt},

MSGDATA=Import[InputPath<>"\\"<>"MSGData.mx"];
MSGOP=MSGDATA["MSGOP"];  (*The MSGOP key of MSGDATA*)
gray=MSGDATA["gray"]; (*The gray key of MSGDATA*)

braLatt=<|
  "CubicP" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, {{a, 0, 0}, {0, a, 0}, {0, 0, a}}}, 
  "CubicF" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, {{0, a/2, a/2}, {a/2, 0, a/2}, {a/2, a/2, 0}}}, 
  "CubicI" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, a}}, {{-a/2, a/2, a/2}, {a/2, -a/2, a/2}, {a/2, a/2, -a/2}}}, 
  "TetragonalP" -> {{{a, 0, 0}, {0, a, 0}, {0, 0, c}},  {{a, 0, 0}, {0, a, 0}, {0, 0, c}}}, 
  "TetragonalI" ->  {{{a, 0, 0}, {0, a, 0}, {0, 0, c}}, {{-a/2, a/2, c/2}, {a/2, -a/2, c/2}, {a/2, a/2, -c/2}}}, 
  "OrthorhombicP" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {{a, 0, 0}, {0, b, 0}, {0, 0, c}}}, 
  "OrthorhombicF" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {{0, b/2, c/2}, {a/2, 0, c/2}, {a/2, b/2, 0}}}, 
  "OrthorhombicI" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {{-a/2, b/2, c/2}, {a/2, -b/2, c/2}, {a/2, b/2, -c/2}}}, 
  "OrthorhombicC" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}},  {{a/2, b/2, 0}, {-a/2, b/2, 0}, {0, 0, c}}}, 
  "OrthorhombicA" -> {{{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {{a, 0, 0}, {0,-b/2, c/2}, {0, b/2, c/2}}},
  "HexagonalP" -> {{{a, 0, 0}, {-a/2, (Sqrt[3]*a)/2, 0}, {0, 0, c}}, {{a, 0, 0}, {-a/2, (Sqrt[3]*a)/2, 0}, {0, 0, c}}}, 
  "TrigonalR"->{{{Sqrt[3] a,0,0},{-((Sqrt[3] a)/2),(3 a)/2,0},{0,0,3 c}}, {{(Sqrt[3] a)/2,a/2,c},{-((Sqrt[3] a)/2),a/2,c},{0,-a,c}}},
(*  "TrigonalR" -> {{{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}, 
    {{a, 0, 0}, {a*Cos[\[Alpha]], a*Sin[\[Alpha]], 0}, 
     {a*Cos[\[Alpha]], a*(Cos[\[Alpha]] - Cos[\[Alpha]]^2)*Csc[\[Alpha]], 
      a*Sqrt[1 - 3*Cos[\[Alpha]]^2 + 2*Cos[\[Alpha]]^3]*Csc[\[Alpha]]}}}, *)
  "MonoclinicP" -> {{{a,0,0},{0,0,b},{c Cos[\[Beta]],c Sin[\[Beta]],0}}, 
  {{a,0,0},{0,0,b},{c Cos[\[Beta]],c Sin[\[Beta]],0}}},  
  "MonoclinicB" -> {{{a, 0, 0}, {0, 0, b}, {c*Cos[\[Beta]], c*Sin[\[Beta]],0}},
   {{a/2,  0, b/2}, {-a/2, 0, b/2}, {c*Cos[\[Beta]], c*Sin[\[Beta]],0}}}, 
     
  "TriclinicP" -> {{{a, 0, 0}, {b*Cos[\[Gamma]], b*Sin[\[Gamma]], 0}, 
     {c*Cos[\[Beta]], c*(Cos[\[Alpha]] - Cos[\[Beta]]*Cos[\[Gamma]])*
       Csc[\[Gamma]], c*Sqrt[1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 + 
         2*Cos[\[Alpha]]*Cos[\[Beta]]*Cos[\[Gamma]] - Cos[\[Gamma]]^2]*
       Csc[\[Gamma]]}}, {{a, 0, 0}, {b*Cos[\[Gamma]], b*Sin[\[Gamma]], 0}, 
     {c*Cos[\[Beta]], c*(Cos[\[Alpha]] - Cos[\[Beta]]*Cos[\[Gamma]])*
       Csc[\[Gamma]], c*Sqrt[1 - Cos[\[Alpha]]^2 - Cos[\[Beta]]^2 + 
         2*Cos[\[Alpha]]*Cos[\[Beta]]*Cos[\[Gamma]] - Cos[\[Gamma]]^2]*
       Csc[\[Gamma]]}}}|>;
(*Output the data*)      
MPGinput=MSGOP[gray[SPNo]];
Print["Space group(BNS): ",MPGinput["MSG"]];
Print["Lattice: ",MPGinput["BRAV"]];
(*Print["Primitive Lattice Vector: ",braLatt[MPGinput["BRAV"]][[2]]];*)
Print["Conventional Lattice Vector: ",braLatt[MPGinput["BRAV"]][[2]]];

spgop=MPGinput["SymmetryOperation"];
Basis=braLatt[MPGinput["BRAV"]][[2]]/.{a->1, b->1,c->1};
{spgop, Basis} (*The output*)
]; (*The end*)      
       
