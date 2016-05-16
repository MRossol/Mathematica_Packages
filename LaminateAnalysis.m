(* ::Package:: *)

BeginPackage["LaminateAnalysis`"]
Unprotect@@Names["LaminateAnalysis`*"];
ClearAll@@Names["LaminateAnalysis`*"];


(*Elastic Analysis*)
compositeS::usage = "compositeS[fiberProp_,matrixProp_,Vf_] calculates the compliance matrix for a 0/90 symmetric laminate with the given properties";
outOfPlaneConstants::usage = "outOfPlaneConstants[fiberProp_,matrixProp_,Vf_] calculates EZ, \[Nu]zx and Gyz from the given input properties";
expLaminateElasticPropDiff::usage = "expLaminateElasticPropDiff[{fiberProp_,matrixProp_,Vf_},expProp_] calculates the normalized difference between experimental results and laminate calculations";
laminateElasticProp::usage = "laminateElasticProp[fiberProp_,matrixProp_,Vf_] calculates {E0,\[Nu]0,E45,\[Nu]45} from the laminate analysis";
calculateTowProperties::usage = "calculateTowProperties[{Em_,Gm_,\[Nu]m_},{{Ef11_,Ef22_},{Gf12_,Gf23_},{\[Nu]f12_,\[Nu]f23_},{\[Xi]e22_,\[Xi]g12_,\[Xi]g23_},fTow_}] calculates tow properties from the input parameters, outputs {Axial Tow Properties,Matrix Properties,Transverse Tow Properties}";

(*Thermal Analysis*)
compositeAlpha::usage = "compositeAlpha[fiberProp_,matrixProp_,Vf_] calcualtes the laminate \[Alpha] for a 0/90 symmetric laminate with the given properties";


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Elastic Analysis*)


compositeS[fiberProp_,matrixProp_,Vf_]:=Module[
{Efa,Eft,Gfa,Gft,vfa,vft,Ema,Emt,Gma,Gmt,vma,vmt,fibersubs,Ea2,Et2,Ga2,Gt2,va2,vt2,c2,k2,eta2,matrixsubs,Ea1,Et1,Ga1,Gt1,va1,vt1,c1,k1,m1,eta1,Eac,vac,Gac,kc,mc,Gtr,x,Achr,Bchr,Cchr,sols,Gtc,vtc,Etc,S0,c,s,theta,A,S90,C0,C90,Sall},

If[SameQ[Length[fiberProp],2],
{Efa,Eft,Gfa,Gft,vfa,vft}={fiberProp[[1]],fiberProp[[1]],fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[2]],fiberProp[[2]]};
,
{Efa,Eft,Gfa,Gft,vfa,vft}=fiberProp;
];
If[SameQ[Length[matrixProp],2],
{Ema,Emt,Gma,Gmt,vma,vmt}={matrixProp[[1]],matrixProp[[1]],matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[2]],matrixProp[[2]]};
,
{Ema,Emt,Gma,Gmt,vma,vmt}=matrixProp;
];

(*Fiber Properties*)
fibersubs={Ea2->Efa,Et2->Eft,Ga2->Gfa,Gt2->Gft,va2->vfa,vt2->vft,c2->Vf};
k2=Ea2*Et2/(2*Ea2-4*Et2*va2^2-2*Ea2*vt2) /.fibersubs;
eta2=3-4*1/2*(1-Gt2/k2)/.fibersubs;

(*Matrix Properties*)
matrixsubs={Ea1->Ema,Et1->Emt,Ga1->Gma,Gt1->Gmt,va1->vma,vt1->vmt,c1->1-Vf};
k1=Ea1*Et1/(2*Ea1-4*Et1*va1^2-2*Ea1*vt1) /.matrixsubs;
m1=1+4*k1*va1^2/Ea1/.matrixsubs;
eta1=3-4*1/2*(1-Gt1/k1)/.matrixsubs;

(*Axial Ply Properties Hashin*)
Eac=Ea1*c1+Ea2*c2+4*(va2-va1)^2*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
vac=va1*c1+va2*c2+(va2-va1)*(1/k1-1/k2)*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
Gac=Ga1*(Ga1*c1+Ga2*(1+c2))/(Ga1*(1+c2)+Ga2*c1)/.fibersubs/.matrixsubs;
kc=(k1*(k2+Gt1)*c1+k2*(k1+Gt1)*c2)/((k2+Gt1)*c1+(k1+Gt1)*c2)/.fibersubs/.matrixsubs;
mc=1+4*kc*vac^2/Eac/.fibersubs/.matrixsubs;
Gtr=Gt2/Gt1/.fibersubs/.matrixsubs;

(*Transverse Ply Properties Christianson*)
Achr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(Gtr*eta1+eta2*eta1-(Gtr*eta1-eta2)*c2^3)*(c2*eta1*(Gtr-1)-(Gtr*eta1+1))/.fibersubs/.matrixsubs;
Bchr=-3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+1/2*(eta1*Gtr+(Gtr-1)*c2+1)*((eta1-1)*(Gtr+eta2)-2*(Gtr*eta1-eta2)*c2^3)+c2/2*(eta1+1)*(Gtr-1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
Cchr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(eta1*Gtr+(Gtr-1)*c2+1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
sols=Quiet[Solve[Achr*x^2+2*Bchr*x+Cchr==0,x]];

Gtc=(Gt1*x/.sols[[2]])/.matrixsubs;
vtc=(kc-mc*Gtc)/(kc+mc*Gtc);
Etc=2*(1+vtc)*Gtc;

(*Compliance and Stiffness Matrixes*)
S0={{1/Eac,-vac/Eac,0},{-vac/Eac,1/Etc,0},{0,0,1/(2*Gac)}};
{c,s}={Cos[theta],Sin[theta]};
A={{c^2,s^2,2*s*c},{s^2,c^2,-2*s*c},{-s*c,s*c,c^2-s^2}};
S90=(Inverse[A].S0.A)/.theta->Pi/2;

C0=Inverse[S0];
C90=Inverse[S90];

Sall=Simplify[Inverse[(C0+C90)/2]]
];


outOfPlaneConstants[fiberProp_,matrixProp_,Vf_]:=Module[
{Efa,Eft,Gfa,Gft,vfa,vft,Ema,Emt,Gma,Gmt,vma,vmt,fibersubs,Ea2,Et2,Ga2,Gt2,va2,vt2,c2,k2,eta2,matrixsubs,Ea1,Et1,Ga1,Gt1,va1,vt1,c1,k1,m1,eta1,Eac,vac,Gac,kc,mc,Gtr,x,Achr,Bchr,Cchr,sols,Gtc,vtc,Etc,Ez,nuzx,Gyz},

If[SameQ[Length[fiberProp],2],
{Efa,Eft,Gfa,Gft,vfa,vft}={fiberProp[[1]],fiberProp[[1]],fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[2]],fiberProp[[2]]};
,
{Efa,Eft,Gfa,Gft,vfa,vft}=fiberProp;
];
If[SameQ[Length[matrixProp],2],
{Ema,Emt,Gma,Gmt,vma,vmt}={matrixProp[[1]],matrixProp[[1]],matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[2]],matrixProp[[2]]};
,
{Ema,Emt,Gma,Gmt,vma,vmt}=matrixProp;
];

(*Fiber Properties*)
fibersubs={Ea2->Efa,Et2->Eft,Ga2->Gfa,Gt2->Gft,va2->vfa,vt2->vft,c2->Vf};
k2=Ea2*Et2/(2*Ea2-4*Et2*va2^2-2*Ea2*vt2) /.fibersubs;
eta2=3-4*1/2*(1-Gt2/k2)/.fibersubs;

(*Matrix Properties*)
matrixsubs={Ea1->Ema,Et1->Emt,Ga1->Gma,Gt1->Gmt,va1->vma,vt1->vmt,c1->1-Vf};
k1=Ea1*Et1/(2*Ea1-4*Et1*va1^2-2*Ea1*vt1) /.matrixsubs;
m1=1+4*k1*va1^2/Ea1/.matrixsubs;
eta1=3-4*1/2*(1-Gt1/k1)/.matrixsubs;

(*Axial Ply Properties Hashin*)
Eac=Ea1*c1+Ea2*c2+4*(va2-va1)^2*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
vac=va1*c1+va2*c2+(va2-va1)*(1/k1-1/k2)*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
Gac=Ga1*(Ga1*c1+Ga2*(1+c2))/(Ga1*(1+c2)+Ga2*c1)/.fibersubs/.matrixsubs;
kc=(k1*(k2+Gt1)*c1+k2*(k1+Gt1)*c2)/((k2+Gt1)*c1+(k1+Gt1)*c2)/.fibersubs/.matrixsubs;
mc=1+4*kc*vac^2/Eac/.fibersubs/.matrixsubs;
Gtr=Gt2/Gt1/.fibersubs/.matrixsubs;

(*Transverse Ply Properties Christianson*)
Achr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(Gtr*eta1+eta2*eta1-(Gtr*eta1-eta2)*c2^3)*(c2*eta1*(Gtr-1)-(Gtr*eta1+1))/.fibersubs/.matrixsubs;
Bchr=-3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+1/2*(eta1*Gtr+(Gtr-1)*c2+1)*((eta1-1)*(Gtr+eta2)-2*(Gtr*eta1-eta2)*c2^3)+c2/2*(eta1+1)*(Gtr-1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
Cchr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(eta1*Gtr+(Gtr-1)*c2+1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
sols=Quiet[Solve[Achr*x^2+2*Bchr*x+Cchr==0,x]];

Gtc=(Gt1*x/.sols[[2]])/.matrixsubs;
vtc=(kc-mc*Gtc)/(kc+mc*Gtc);
Etc=2*(1+vtc)*Gtc;

Ez=-((Eac Etc (Eac+Etc+2 Etc vac))/((Etc^2 vac^2+Eac^2 (-1+vtc^2)-Eac Etc (1+2 vac (1+vtc)))));
nuzx=-((Etc (Etc vac^2+Eac (vac+vtc+vac vtc)))/(2 (Etc^2 vac^2+Eac^2 (-1+vtc^2)-Eac Etc (1+2 vac (1+vtc)))));
Gyz=(2 Gac Gtc)/(Gac+Gtc);

(*{{Ez,Eac,Etc},{nuzx,vac,vtc},{Gyz,Gac,Gtc}}*)
{Ez,nuzx,Gyz}
];


expLaminateElasticPropDiff[{fiberProp_,matrixProp_,Vf_},expProp_]:=Module[
{c,s,theta,A,eps45try,eps45,eps45rot,El45,vl45,E45l,v45l,eps0try,eps0,El0,vl0,E0l,v0l,sheartry,epsshear,G0,Gl,E0,E45,G,v0,v45},

{c,s}={Cos[theta],Sin[theta]};
A={{c^2,s^2,2*s*c},{s^2,c^2,-2*s*c},{-s*c,s*c,c^2-s^2}};

(*Simple Loading to calculate Composite Properties*)
eps45try=compositeS[fiberProp,matrixProp,Vf].(1*{1/2,1/2,1/2});
eps45={1/El45,-1*vl45/El45,0};
eps45rot=A.eps45/.theta->-Pi/4;
{E45l,v45l}={El45,vl45}/.First[Quiet[Solve[eps45try==eps45rot,{El45,vl45}]]];

eps0try=compositeS[fiberProp,matrixProp,Vf].(1*{1,0,0});
eps0={1/El0,-vl0/El0,0};
{E0l,v0l}={El0,vl0}/.First[Quiet[Solve[eps0==eps0try,{El0,vl0}]]];

(*sheartry=Simplify[compositeS[{Ef,vf},{Em,vm},Vf].{0,0,1}];
epsshear={0,0,1/(2*G0)};
Gl=First[G0/.Quiet[Solve[sheartry\[Equal]epsshear,G0]]];

{E0,E45,G,v0,v45}=expProp;
((E0-E0l)/E0)^2+((v0-v0l)/v0)^2+((E45-E45l)/E45)^2+((v45-v45l)/v45)^2+((G-Gl)/G)^2
]*)

If[Length[expProp]==2,
{E0,v0}=expProp;
((E0-E0l)/E0)^2+((v0-v0l)/v0)^2
,
{E0,v0,E45,v45}=expProp;
((E0-E0l)/E0)^2+((v0-v0l)/v0)^2+((E45-E45l)/E45)^2+((v45-v45l)/v45)^2
]
];


laminateElasticProp[fiberProp_,matrixProp_,Vf_]:=Module[
{c,s,theta,A,eps45try,eps45,eps45rot,El45,vl45,E45l,v45l,eps0try,eps0,El0,vl0,E0l,v0l,sheartry,epsshear,G0,Gl,E0,E45,G,v0,v45},

{c,s}={Cos[theta],Sin[theta]};
A={{c^2,s^2,2*s*c},{s^2,c^2,-2*s*c},{-s*c,s*c,c^2-s^2}};

(*Simple Loading to calculate Composite Properties*)
eps45try=compositeS[fiberProp,matrixProp,Vf].(1*{1/2,1/2,1/2});
eps45={1/El45,-1*vl45/El45,0};
eps45rot=A.eps45/.theta->-Pi/4;
{E45l,v45l}={El45,vl45}/.First[Quiet[Solve[eps45try==eps45rot,{El45,vl45}]]];

eps0try=compositeS[fiberProp,matrixProp,Vf].(1*{1,0,0});
eps0={1/El0,-vl0/El0,0};
{E0l,v0l}={El0,vl0}/.First[Quiet[Solve[eps0==eps0try,{El0,vl0}]]];

{E0l,v0l,E45l,v45l}
];


calculateTowProperties[{E1_,G1_,\[Nu]1_},{{Ea2_,Et2_},{Ga2_,Gt2_},{va2_,vt2_},c2_}]:=Module[
{Et11,Et22,Gt12,Gt23,\[Nu]t12,\[Nu]t21,\[Nu]t23,AxialTowProperties,TransverseTowProperties,MatrixProperties,c1,k1,k2,Ga1,Gt1,Ea1,Et1,va1,vt1,eta1,eta2,Gtr,Achr,Bchr,Cchr,sols,Eac,vac,vac2,Gac,Gtc,vtc,Etc,kc,mc},

Ga1=G1;
Gt1=G1;
Ea1=E1;
Et1=E1;
va1=\[Nu]1;
vt1=\[Nu]1;
k1=Ea1*Et1/(2*Ea1-4*Et1*va1^2-2*Ea1*vt1) ;
k2=Ea2*Et2/(2*Ea2-4*Et2*va2^2-2*Ea2*vt2) ;
c1= 1 - c2;
eta1=3-4*1/2*(1-Gt1/k1);
eta2=3-4*1/2*(1-Gt2/k2);

Eac=Ea1*c1+Ea2*c2+4*(va2-va1)^2*c1*c2/(c1/k2+c2/k1+1/Gt1);
vac=va1*c1+va2*c2+(va2-va1)*(1/k1-1/k2)*c1*c2/(c1/k2+c2/k1+1/Gt1);
Gac=Ga1*(Ga1*c1+Ga2*(1+c2))/(Ga1*(1+c2)+Ga2*c1);
kc=(k1*(k2+Gt1)*c1+k2*(k1+Gt1)*c2)/((k2+Gt1)*c1+(k1+Gt1)*c2);
mc=1+4*kc*vac^2/Eac;

Gtr=Gt2/Gt1;
Achr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(Gtr*eta1+eta2*eta1-(Gtr*eta1-eta2)*c2^3)*(c2*eta1*(Gtr-1)-(Gtr*eta1+1));
Bchr=-3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+1/2*(eta1*Gtr+(Gtr-1)*c2+1)*((eta1-1)*(Gtr+eta2)-2*(Gtr*eta1-eta2)*c2^3)+c2/2*(eta1+1)*(Gtr-1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3);
Cchr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(eta1*Gtr+(Gtr-1)*c2+1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3);
sols=Solve[Achr*x^2+2*Bchr*x+Cchr==0,x];
Gtc=Gt1*x/.sols[[2]];
vtc=(kc-mc*Gtc)/(kc+mc*Gtc);
Etc=2*(1+vtc)*Gtc;
vac2=(vac/Eac)*Etc;
vtc=Etc/(2*Gtc)-1;

AxialTowProperties={Eac,Etc,Etc,vac,vac,vtc,Gac,Gac,Gtc};
TransverseTowProperties={Etc,Etc,Eac,vtc,vac2,vac2,Gtc,Gtc,Gac};
MatrixProperties={E1,\[Nu]1};

{AxialTowProperties,MatrixProperties,TransverseTowProperties}

];


(* ::Subsection::Closed:: *)
(*Thermal Analysis*)


compositeAlpha[fiberProp_,matrixProp_,Vf_]:=Module[
{Efa,Eft,Gfa,Gft,vfa,vft,afa,aft,Ema,Emt,Gma,Gmt,vma,vmt,ama,amt,fibersubs,Ea2,Et2,Ga2,Gt2,va2,vt2,c2,k2,eta2,matrixsubs,Ea1,Et1,Ga1,Gt1,va1,vt1,c1,k1,m1,eta1,Eac,vac,Gac,kc,mc,Gtr,x,Achr,Bchr,Cchr,sols,Gtc,vtc,Etc,S0,c,s,theta,A,S90,C0,C90,Sall,Sf,Sm,Sbar,P,alphaf,alpham,abara,abart,aa,at,a0,a90,aAll},

If[SameQ[Length[fiberProp],3],
{Efa,Eft,Gfa,Gft,vfa,vft,afa,aft}={fiberProp[[1]],fiberProp[[1]],fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[1]]/(2*(1+fiberProp[[2]])),fiberProp[[2]],fiberProp[[2]],fiberProp[[3]],fiberProp[[3]]};
,
{Efa,Eft,Gfa,Gft,vfa,vft,afa,aft}=fiberProp;
];
If[SameQ[Length[matrixProp],3],
{Ema,Emt,Gma,Gmt,vma,vmt,ama,amt}={matrixProp[[1]],matrixProp[[1]],matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[1]]/(2*(1+matrixProp[[2]])),matrixProp[[2]],matrixProp[[2]],matrixProp[[3]],matrixProp[[3]]};
,
{Ema,Emt,Gma,Gmt,vma,vmt,ama,amt}=matrixProp;
];

(*Fiber Properties*)
fibersubs={Ea2->Efa,Et2->Eft,Ga2->Gfa,Gt2->Gft,va2->vfa,vt2->vft,c2->Vf};
k2=Ea2*Et2/(2*Ea2-4*Et2*va2^2-2*Ea2*vt2) /.fibersubs;
eta2=3-4*1/2*(1-Gt2/k2)/.fibersubs;

(*Matrix Properties*)
matrixsubs={Ea1->Ema,Et1->Emt,Ga1->Gma,Gt1->Gmt,va1->vma,vt1->vmt,c1->1-Vf};
k1=Ea1*Et1/(2*Ea1-4*Et1*va1^2-2*Ea1*vt1) /.matrixsubs;
m1=1+4*k1*va1^2/Ea1/.matrixsubs;
eta1=3-4*1/2*(1-Gt1/k1)/.matrixsubs;

(*Axial Ply Properties Hashin*)
Eac=Ea1*c1+Ea2*c2+4*(va2-va1)^2*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
vac=va1*c1+va2*c2+(va2-va1)*(1/k1-1/k2)*c1*c2/(c1/k2+c2/k1+1/Gt1)/.fibersubs/.matrixsubs;
Gac=Ga1*(Ga1*c1+Ga2*(1+c2))/(Ga1*(1+c2)+Ga2*c1)/.fibersubs/.matrixsubs;
kc=(k1*(k2+Gt1)*c1+k2*(k1+Gt1)*c2)/((k2+Gt1)*c1+(k1+Gt1)*c2)/.fibersubs/.matrixsubs;
mc=1+4*kc*vac^2/Eac/.fibersubs/.matrixsubs;
Gtr=Gt2/Gt1/.fibersubs/.matrixsubs;

(*Transverse Ply Properties Christianson*)
Achr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(Gtr*eta1+eta2*eta1-(Gtr*eta1-eta2)*c2^3)*(c2*eta1*(Gtr-1)-(Gtr*eta1+1))/.fibersubs/.matrixsubs;
Bchr=-3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+1/2*(eta1*Gtr+(Gtr-1)*c2+1)*((eta1-1)*(Gtr+eta2)-2*(Gtr*eta1-eta2)*c2^3)+c2/2*(eta1+1)*(Gtr-1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
Cchr=3*c2*c1^2*(Gtr-1)*(Gtr+eta2)+(eta1*Gtr+(Gtr-1)*c2+1)*(Gtr+eta2+(Gtr*eta1-eta2)*c2^3)/.fibersubs/.matrixsubs;
sols=Quiet[Solve[Achr*x^2+2*Bchr*x+Cchr==0,x]];

Gtc=(Gt1*x/.sols[[2]])/.matrixsubs;
vtc=(kc-mc*Gtc)/(kc+mc*Gtc);
Etc=2*(1+vtc)*Gtc;

(*Compliance and Stiffness Matrixes*)
S0={{1/Eac,-vac/Eac,0},{-vac/Eac,1/Etc,0},{0,0,1/(2*Gac)}};
{c,s}={Cos[theta],Sin[theta]};
A={{c^2,s^2,2*s*c},{s^2,c^2,-2*s*c},{-s*c,s*c,c^2-s^2}};
S90=(Inverse[A].S0.A)/.theta->Pi/2;

C0=Inverse[S0];
C90=Inverse[S90];

Sall=Inverse[(C0+C90)/2];

(*Fiber and Matrix Stiffness*)
Sf={{1/Ea2,-va2/Ea2,0},{-va2/Ea2,1/Et2,0},{0,0,1/(2*Ga2)}}/.fibersubs;
Sm={{1/Ea1,-va1/Ea1,0},{-va1/Ea1,1/Et1,0},{0,0,1/(2*Ga1)}}/.matrixsubs;
Sbar=Sf Vf+Sm (1-Vf);
P=Inverse[(Sf-Sm)];

(*Thermal Properties*)
alphaf={afa,aft,0};
alpham={ama,amt,0};
abara=alphaf[[1]]Vf+alpham[[1]](1-Vf);
abart=alphaf[[1]]Vf+alpham[[2]](1-Vf);

(*Ply CTE*)
aa=abara+(alphaf-alpham).P.(S0[[1]]-Sbar[[1]]);
at=abart+(alphaf-alpham).P.(S0[[2]]-Sbar[[2]]);
a0={aa,at,0};
a90={at,aa,0};

(*laminate CTE*)
aAll=Sall.(C0.a0+C90.a90)/2
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["LaminateAnalysis`*"],
{
}];
EndPackage[ ]
