(* ::Package:: *)

BeginPackage["MarcovChainFunctions`"]
Unprotect@@Names["MarcovChainFunctions`*"];
ClearAll@@Names["MarcovChainFunctions`*"];


(*PreFunctions*)
ptmset::usage = "ptm and svec are 3D arrays with distinct entries for xi_x, xi_y, and xi_z ";
ptmtri::usage = "in notation of paper, b=gamma and c=beta";
ptmtricyclic::usage = "tridiagonal PTM generator that has periodic conditions at extreme values, in notation of paper, b=gamma and c=beta";
varcovar::usage = "find characteristics of PTM";
xigen::usage = "module for generating sequence of random numbers ";
quasicontphi::usage = "modify sequence of angles phi to make quasi-continuous by shifting all following values if change in value exceeds half difference between extrema of mesh 
elements of phi are assumed to be distributed over [0, 2Pi]";
xidiff::usage = "compute average change in xi vs. separation of grid points ";
xidiffabs::usage = "compute average change in Abs[xi] vs. separation of grid points";
dfact::useage = "(* evaluate ratio of double factorials *)";

(*Markov Chain Single Sequence Generator*)
markovRanSeq::usage = "code to generate a sequence of random numbers on a uniform linear grid using the Markov Chain algorithm of Blacklock et al, 2012 input: nn = number of grid points 
sig = RMSD of distribution function for random numbers 
colen = correlation length for random numbers in units of grid interval 
m = 10, 20, or 40 is number of discrete values allowed for random numbers\[IndentingNewLine]maxsig = multiple of RMSD of maximum allowed random number value\[IndentingNewLine]if !cyclic, generated numbers will be normally distributed around zero \[IndentingNewLine]if cyclic, generated numbers can jump from one extreme to the other\[IndentingNewLine]use cyclic=True and maxsig << rmsd to generate numbers that are uniformly distributed and cyclic (e.g., angles on [0,2 Pi]))\[IndentingNewLine]other distributions can be derived by transformation of variables ";
rwStep::usage = "calibrate step size in random walk model using xidiff formed from expt'l data in xi output step size in radians is normalized by grid step\[IndentingNewLine]between input data points and thus has units length^-1 this anticipates possible change in grid step in virtual specimen generator";
phiGen::usage = "module for generating nn angles for virtual specimen via random walk with expectation Abs[change in angle] = measured change in expt'l data\[IndentingNewLine]for any separation of pairs of grid points input step size from rwstep, which has been normalized to expt'l grid spacing input gridstep can be different from grid spacing in expt'l data";
marcovDFT::useage = "marcovDFT[{amp_,ampCoLen_},{phi_,phiCoLen_},useRange_,nn_,OptionsPattern[]] recreates DFT spectram for the given amp and phi inputs in the desired S range using a the marcov chain generator.";


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Prefunctions*)


(* ptm and svec are 3D arrays with distinct entries for xi_x, xi_y, and xi_z *)
ptmset[rmsd_,covar_,m_,a_,alpha_,varc1_,covarc2_,cyclic_]:=Module[
{x,c1,c2,ptm,sol1,eig,svec},
If[!cyclic,
sol1=x/.FindRoot[varc1[x]==rmsd/(a Sqrt[m (m+1)/3]),{x,0.1}][[1]];
c1=((1-alpha)(1-sol1))/2;
];
(* choose c1=gamma for required variance *)
c2=x/.FindRoot[covarc2[x]==covar,{x,8}][[1]]; 
c2=IntegerPart[2^(c2-1)];(* choose c2 = power N to fit correlation  *)

If[!cyclic,
ptm=ptmtri[m,c1,alpha];
,
ptm=ptmtricyclic[m,alpha];
];
ptm=MatrixPower[ptm,c2]; 
svec=varcovar[ptm];
{ptm,svec}
];


(* in notation of paper, b=gamma and c=beta *)
ptmtri[m_,c1_,alpha_]:=Module[
{a,b,c,d,aa},
a=alpha;
b=c1;
c=1-a-b;
d=(b+c)/2;

aa=DiagonalMatrix[ConstantArray[a,2m+1]];
aa=aa+DiagonalMatrix[ConstantArray[b,m]~Join~ConstantArray[c,m],1];
aa=aa+DiagonalMatrix[ConstantArray[c,m]~Join~ConstantArray[b,m],-1];
aa=ReplacePart[aa,{{m+2,m+1}->d,{m,m+1}->d,{2,1}->1-a,{2m,2m+1}->1-a}];
aa
];


(* tridiagonal PTM generator that has periodic conditions at extreme values *)
(* in notation of paper, b=gamma and c=beta *)
ptmtricyclic[m_,alpha_]:=Module[
{a,b,c,d,aa},
a=alpha;(* seek uniform pdf *) 
b=(1-a)/2;
c=b;
d=(b+c)/2;

aa=DiagonalMatrix[ConstantArray[a,2m+1]];
aa=aa+DiagonalMatrix[ConstantArray[b,m]~Join~ConstantArray[c,m],1];
aa=aa+DiagonalMatrix[ConstantArray[c,m]~Join~ConstantArray[b,m],-1];
aa=ReplacePart[aa,{{m+2,m+1}->d,{m,m+1}->d,{2m+1,1}->b,{1,2m+1}->b}];
aa
];


(*find characteristics of PTM*)
varcovar[ptm_]:=Module[
{eg1,stvect},
eg1=Eigensystem[ptm];
stvect=eg1[[2,1]]/(Mean[eg1[[2,1]]]Length[eg1[[2,1]]]);
stvect
];


(* module for generating sequence of random numbers *)
xigen[svec_,a_,ptm_,nn_,m_]:=Module[
{dist,cpd,r,mr},
dist=svec;
Table[
cpd=Table[{-(Length[dist]+1)/2+j,Total[dist[[;;j]]]},{j,Length[dist]}];
r=RandomReal[];
mr=Select[cpd,(r<#[[2]])&][[1,1]];
dist=ReplacePart[ConstantArray[0,2m+1],(mr+m+1)->1];
dist=ptm.dist;
mr a
,{i,nn}]
];


(* modify sequence of angles phi to make quasi-continuous by shifting all following values *)
(* if change in value exceeds half difference between extrema of mesh *)
(* elements of phi are assumed to be distributed over [0, 2Pi] *)
(* brian cox Oct 16 2013 *)
quasicontphi[phi_]:=Module[
{del=0,phiOut={}},
AppendTo[phiOut,phi[[1]]];
Do[If[phi[[i]]-phi[[i-1]]>Pi,del=del-2Pi];
If[phi[[i-1]]-phi[[i]]>Pi,del=del+2Pi];
AppendTo[phiOut,phi[[i]]+del],
{i,2,Length[phi]}];
phiOut
]


(* compute average change in xi vs. separation of grid points *)
xidiff[xi_]:=Module[
{diff},
diff=Table[Sum[xi[[i+k]]-xi[[i]],{i,Length[xi]-k}]/(Length[xi]-k),{k,0,Length[xi]-1}];
diff
]


(* compute average change in Abs[xi] vs. separation of grid points *)
xidiffabs[xi_]:=Module[
{diff},
diff=Table[Sum[Abs[xi[[i+k]]-xi[[i]]],{i,Length[xi]-k}]/(Length[xi]-k),{k,0,Length[xi]-1}];
diff
]


(* evaluate ratio of double factorials *)
dfact[i_]:=Module[
{temp},
If[
i==2IntegerPart[i/2],
temp=(i-1)!!/(i-2)!!,
temp=i!!/(i-1)!!
];
temp
]


(* ::Subsection::Closed:: *)
(*Markov Chain (Single Sequence)*)


Options[markovRanSeq]={
mVal->10,maxsigVal->3,cyclicBoole->False
};


(* code to generate a sequence of random numbers on a uniform linear grid *)
(* using the Markov Chain algorithm of Blacklock et al, 2012 *)
(* input: *)
(* nn = number of grid points *)
(* sig = RMSD of distribution function for random numbers *)
(* colen = correlation length for random numbers in units of grid interval *)
(* m = 10, 20, or 40 is number of discrete values allowed for random numbers *)
(* maxsig = multiple of RMSD of maximum allowed random number value, recommended value is 3 *)
(* if !cyclic, generated numbers will be normally distributed around zero *)
(* if cyclic, generated numbers can jump from one extreme to the other *)
(* use cyclic=True and maxsig << rmsd to generate numbers that are uniformly distributed and cyclic (e.g., angles on [0,2 Pi])) *)
(* other distributions can be derived by transformation of variables *)
(* brian cox August 2013 *)
markovRanSeq[nn_,rmsd_,colen_,OptionsPattern[]]:=Module[{
m,maxsig,cyclic,alpha,covar,a,varc1datam10,varc1datam20,varc1datam40,varc1data,x,varc1,covarc2datam10,covarc2datam20,covarc2datam40,covarc2data,covarc2,temp,ptm,svec,xi,io=False
},
m=OptionValue[mVal];
maxsig=OptionValue[maxsigVal];
cyclic=OptionValue[cyclicBoole];
alpha=0.9;  (* don't change this number *)
(* assigned values for parameters alpha, m, and maxsig determine mesh size, a, for xi *)
covar=1.-1./colen; (* covariance from grid point to next grid point *)
a=maxsig rmsd/m;   (* mesh size for discretized random numbers *)

(* using the following interpolating functions to set the PTM *)
(* is only consistent with choosing m=10, m=20, or m=40 and alpha = 0.9 *) 
(* these functions summarize analysis by Brian Cox of PTM's as described in Blacklock et al 2012 *)
varc1datam10={{0.`,0.9558432735738244`},{0.1000000000000001`,0.718677165671946`},{0.20000000000000007`,0.5183767972589418`},{0.30000000000000004`,0.3795376917139528`},{0.40000000000000013`,0.29071374726361493`},{0.5000000000000001`,0.23343099554923147`},{0.6000000000000002`,0.19461694031507906`},{0.7000000000000002`,0.16682092869070997`},{0.8000000000000002`,0.1459685490160884`},{0.9000000000000002`,0.12974982402495536`},{1.0000000000000002`,0.11677484162422842`}};
varc1datam20={{0.`,0.9765098200076662`},{0.1000000000000001`,0.5275109150343484`},{0.20000000000000007`,0.29700397348957647`},{0.30000000000000004`,0.19916963783861671`},{0.400000000,0.14940310610137422},{0.5000000000000001`,0.11952285716282474`},{0.6000000000000002`,0.09960238409707248`},{0.7000000000000002`,0.08537347209533938`},{0.8000000000000002`,0.07470178808338679`},{0.9000000000000002`,0.06640158940747504`},{1.0000000000000002`,0.05976143046671966`}};
varc1datam40={{0.`,0.9878839173431374`},{0.05000000000000005`,0.5331241514266216`},{0.1000000000000001`,0.3004466843014011`},{0.15000000000000002`,0.20157196126159319`},{0.20000000000000007`,0.151213610087148`},{0.2500000000000001`,0.12097165799548439`},{0.30000000000000004`,0.10080972954031583`},{0.3500000000000001`,0.08640833984078489`},{0.40000000000000013`,0.07560729736294722`},{0.4500000000000001`,0.06720648654468185`},{0.5000000000000001`,0.0604858378916395`},{0.5500000000000002`,0.05498712535480318`},{0.6000000000000002`,0.050404864908168055`},{0.6500000000000001`,0.04652756760818502`},{0.7000000000000002`,0.0432041699224356`},{0.7500000000000002`,0.04032389192709704`},{0.8000000000000002`,0.037803648681445104`},{0.8500000000000002`,0.03557990464153467`},{0.9000000000000002`,0.033603243272746025`},{0.9500000000000002`,0.03183465152160171`},{1.0000000000000002`,0.030242918945456693`}};
varc1=Interpolation[Which[m==10,varc1datam10,m==20,varc1datam20,m==40,varc1datam40]];

covarc2datam10={{1,0.9954892663230123`},{2,0.9910087241796338`},{3,0.9821350141883031`},{4,0.9647184441022717`},{5,0.931096796513396`},{6,0.8681114760971854`},{7,0.7563561701070254`},{8,0.5768454816161662`},{9,0.3379117464303905`},{10,0.11675178085522502`},{11,0.013976956273629975`},{12,0.00020036549475288749`},{13,4.1176016782789425`*^-8},{14,1.7248275481518765`*^-15}};
covarc2datam20={{1,0.9988841200921237`},{2,0.9977704921633825`},{3,0.9955497329541253`},{4,0.991132787340377`},{5,0.982389600250172`},{6,0.9652323762946868`},{7,0.9320973766192114`},{8,0.869951009865906`},{9,0.7594362987511677`},{10,0.5812260652083066`},{11,0.34259734938929115`},{12,0.11973259895518845`},{13,0.014657760649312407`},{14,0.00021971572674392352`}};
covarc2datam40={{1,0.9997210605952394`},{2,0.9994423025492974`},{3,0.9988853020259298`},{4,0.9977732134717437`},{5,0.9955559390346221`},{6,0.991146138470542`},{7,0.9824162827980779`},{8,0.9652832413290935`},{9,0.9321918516195303`},{10,0.8701189793377145`},{11,0.7597071582530956`},{12,0.5815871475736222`},{13,0.34294629225995094`},{14,0.11992985600385801`}};
covarc2=Interpolation[Which[m==10,covarc2datam10,m==20,covarc2datam20,m==40,covarc2datam40]];

(* set up prob trans matrices and steady state distribns for each type of tow *)
{ptm,svec}=ptmset[rmsd,covar,m,a,alpha,varc1,covarc2,cyclic];

(* execute generation of sequence of random numbers *)
xi=xigen[svec,a,ptm,nn,m];
xi

];


Options[rwStep]={
GridStep->1
};


(* calibrate step size in random walk model using xidiff formed from expt'l data in xi *)
(* output step size in radians is normalized by grid step *)
(* between input data points and thus has units length^-1 *)
(* this anticipates possible change in grid step in virtual specimen generator *)
rwStep[phi_,OptionsPattern[]]:=Module[
{phiIn,dphi,w,a},
phiIn=quasicontphi[phi];
(*dphi=xidiff[phi];*)
dphi=xidiffabs[phiIn];
w=Reverse[Range[1,Length[dphi]]];
a=MapIndexed[#1/dfact[First[#2]]&,dphi];
Total[MapThread[#1 #2&,{w,a}]]/(Total[w]/OptionValue[GridStep])
]


Options[phiGen]={
GridStep->1
};


(* module for generating nn angles for virtual specimen via random walk *)
(* with expectation Abs[change in angle] = measured change in expt'l data *)
(* for any separation of pairs of grid points *)
(* input step size from rwstep, which has been normalized to expt'l grid spacing *)
(* input gridstep can be different from grid spacing in expt'l data *)
phiGen[rwstep_,nn_,OptionsPattern[]]:=Module[
{rwstep0,r,phi0,phiOut={}}, 
rwstep0=rwstep OptionValue[GridStep];
r=Random[];
phi0=-Pi+2 Pi r; 
AppendTo[phiOut,phi0];
Do[ 
r=Random[];
If[r<0.5,phi0=phi0-rwstep0,phi0=phi0+rwstep0];
AppendTo[phiOut,phi0];
,{i,2,nn}];
phiOut
]


Options[marcovDFT]={
mVal->10,maxsigVal->3,RandomSeedValue->None,GridStep->1
};


marcovDFT[{amp_,ampCoLen_},{phiRWstep_},useRange_,nn_,OptionsPattern[]]:=Module[
{newDFT,sRange,newAmp,newPhi},

newDFT=ConstantArray[0.+0.I,{nn,Length[amp]}];

sRange=Which[SameQ[useRange,All],All,!ListQ[useRange]&&SameQ[Length[useRange],2],Range[useRange[[1]],useRange[[2]]]~Join~Range[-useRange[[2]],-useRange[[1]]],ListQ[useRange],useRange];

If[NumberQ[OptionValue[RandomSeedValue]],
If[ListQ[ampCoLen],
SeedRandom[OptionValue[RandomSeedValue]];
newAmp=MapThread[(#1+#2)&,{amp[[sRange,1]],MapThread[N[markovRanSeq[nn,#1,#2,mVal->OptionValue[mVal],maxsigVal->OptionValue[maxsigVal]]]&,{amp[[sRange,2]],ampCoLen[[sRange]]}]}];
,
SeedRandom[OptionValue[RandomSeedValue]];
newAmp=MapThread[(#1+#2)&,{amp[[sRange,1]],Map[N[markovRanSeq[nn,#,ampCoLen,mVal->OptionValue[mVal],maxsigVal->OptionValue[maxsigVal]]]&,amp[[sRange,2]]]}];
];
If[ListQ[phiRWstep],
SeedRandom[OptionValue[RandomSeedValue]];
newPhi=Map[phiGen[#,nn,GridStep->OptionValue[GridStep]]&,phiRWstep[[sRange]]];
,
SeedRandom[OptionValue[RandomSeedValue]];
newPhi=ConstantArray[phiGen[phiRWstep,nn],Length[newAmp]];
];
,
SeedRandom[];
If[ListQ[ampCoLen],
newAmp=MapThread[(#1+#2)&,{amp[[sRange,1]],MapThread[N[markovRanSeq[nn,#1,#2,mVal->OptionValue[mVal],maxsigVal->OptionValue[maxsigVal]]]&,{amp[[sRange,2]],ampCoLen[[sRange]]}]}];
,
newAmp=MapThread[(#1+#2)&,{amp[[sRange,1]],Map[N[markovRanSeq[nn,#,ampCoLen,mVal->OptionValue[mVal],maxsigVal->OptionValue[maxsigVal]]]&,amp[[sRange,2]]]}];
];
If[ListQ[phiRWstep],
newPhi=Map[phiGen[#,nn,GridStep->OptionValue[GridStep]]&,phiRWstep[[sRange]]];
,
newPhi=ConstantArray[phiGen[phiRWstep,nn],Length[newAmp]];
];
];

newDFT[[All,sRange]]=Transpose[MapThread[(#1 Cos[#2]+#1 Sin[#2]I)&,{newAmp,newPhi}]];
newDFT
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["MarcovChainFunctions`*"],
{

}];
EndPackage[ ]
