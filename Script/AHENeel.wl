(* ::Package:: *)

(* AHE-texture: To determine determine the anomalous Hall vector \[Sigma]H and the corresponding multipoles in the magnetic order space *)
(*See https://arxiv.org/abs/2411.10147*)
(*Rui-Chun Xiao, Anhui University, xiaoruichun@ahu.edu.cn*)
AHENeel[SymmOperation_,Basis_, Magup_] := Module[{OsMatCry, OsVectCry, isSwitch, BasisT,
OsMat, Index, index, Tensor2, Tensor2New, mysolveRi, Tensor4, Tensor4New, Tensor4tmp, 
LengthOfOsMat, m,p,q,l,i,j,k,
\[Sigma]1,\[Sigma]3, Neel},

{Nsym,n}=Dimensions[SymmOperation];

OsMatCry=ConstantArray[0,{Nsym,3,3}];
OsVectCry=ConstantArray[0,{Nsym,3}];
(**To construct the R and t**)
For [i=1,i<=Nsym,i++,
OsMatCry[[i,;;,;;]]=SymmOperation[[i,2]];
OsVectCry[[i,;;]]=SymmOperation[[i,3]];
];

(**To determine whether the aom postion is change or not**)
{Nmagup,m}=Dimensions[Magup];
isSwitch=ConstantArray[-1,Nsym];
For[i=1,i<=Nsym,i++,
 k=1; (*Chech the First atom is OK*)
 newposition=OsMatCry[[i,;;,;;]] . Magup[[k]]+OsVectCry[[i,;;]];
  For[j=1,j<=Nmagup,j++,
     If[ Mod[newposition[[1]]-Magup[[j,1]],1]==0 && Mod[newposition[[2]]-Magup[[j,2]],1]==0&&Mod[newposition[[3]]-Magup[[j,3]],1]==0,
     isSwitch[[i]]=1;
]]];

OsMat=ConstantArray[0,{Nsym,3,3}];
BasisT=Transpose[Basis];
For[i=1,i<=Nsym,i++,
   OsMat[[i,;;,;;]]=BasisT . OsMatCry[[i,;;,;;]] . Inverse[BasisT]
];
LengthOfOsMat = Length[OsMat];  
(*Define the 1st order matrix*)
 Index={"X","Y","Z"}; index={"x","y","z"};
Tensor2=ConstantArray[0,{3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
   Tensor2[[i,j]]=Subscript[T, Index[[i]] <> index[[j]]];
]];
(*Solve the 1st order matrix*)
Tensor2New=Tensor2;
For [r=1,r<=LengthOfOsMat,r++,
 Ri=OsMat[[r]];
  Tensor2rot=isSwitch[[r]]*Ri . Tensor2New . Inverse[Ri];
  mysolveRi=Reduce[Tensor2rot==Tensor2New,Flatten[Tensor2]];
  Tensor2New=Tensor2New/.ToRules[mysolveRi];
];

(*******************************)
(*Define the 3nd order term tensor*)
Tensor4=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
    For[m=1,m<=3,m++,
    AA=Sort[{j,k,m}];
    Tensor4[[i,j,k,m]]=Subscript[T, Index[[i]]<> index[[AA[[1]]]]<> index[[AA[[2]]]]<> index[[AA[[3]]]]];
]]]];
Tensor4New=Tensor4;
For [r=1,r<=LengthOfOsMat,r++,  
Ri=OsMat[[r]];
(* to rotate the matrix*)
Tensor4tmp=ConstantArray[0,{3,3,3,3}];
 For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   For[k=1,k<=3,k++,
    For[l=1,l<=3,l++,
     For[m=1,m<=3,m++,
      For[n=1,n<=3,n++,
       For[p=1,p<=3,p++,
        For[q=1,q<=3,q++,
Tensor4tmp[[i,j,k,l]]=Tensor4tmp[[i,j,k,l]]+isSwitch[[r]]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4New[[ m,n,p,q ]];
]]]]]]]];
mysolveRi=Reduce[Tensor4New==Tensor4tmp,Flatten[Tensor4]];
Tensor4New=Tensor4New/.ToRules[mysolveRi];
];

(***print the \[Sigma] vector***)
Clear [i,j,k,l,m,n];
\[Sigma]1={0, 0, 0}; \[Sigma]3={0, 0, 0}; Neel={Subscript[n, x],Subscript[n, y],Subscript[n, z]};
(***The 1st order term**)
For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   \[Sigma]1[[i]]=\[Sigma]1[[i]]+Tensor2New[[i,j]]*Neel[[j]];
]];
(***The 3rd order term**)
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
   For[m=1,m<=3,m++,
    \[Sigma]3[[i]]=\[Sigma]3[[i]]+Tensor4New[[i,j,k,m]]*Neel[[j]]*Neel[[k]]*Neel[[m]];
]]]];
Print["\!\(\*SubscriptBox[\(\[Sigma]\), \(H\)]\)=", \[Sigma]1//MatrixForm,"+",\[Sigma]3//MatrixForm];
Print["Independent 1st order terms: ",Variables[Tensor2New],", Number:",Length[Variables[Tensor2New]]];
Print["Independent 3rd order terms: ",Variables[Tensor4New],", Number:",Length[Variables[Tensor4New]]];
sigma=\[Sigma]1+\[Sigma]3;
sigma
](*The end*)
