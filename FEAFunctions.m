(* ::Package:: *)

BeginPackage["FEAFunctions`"]
Unprotect@@Names["FEAFunctions`*"];
ClearAll@@Names["FEAFunctions`*"];


(*Modify INP*)
getMaterialOrder::usage = "getMaterialOrder[filesIn_] determines the order of material properties in the given .inp file";
selectMaterialProperties::usage = "selectMaterialProperties[fileIn_,{AxialPropsIn_,TransversePropsIn_,MatrixPropsIn_},fileName_] GUI for creation of new .inp with the input material properties";
changeMaterialProperties::usage = "changeMaterialProperties[fileIn_,PropsIn_,magnitude_,fileName_] changes the material properties in the input .inp and saves a new inp";
\[Eta]::usage = "\[Eta][fiberProp_,matrixProp_,\[Xi]_] function for halpin-tsai";
calculateTowProperties::usage = "calculateTowProperties[{Em_,Gm_,\[Nu]m_},{{Ef11_,Ef22_},{Gf12_,Gf23_},{\[Nu]f12_,\[Nu]f23_},{\[Xi]e22_,\[Xi]g12_,\[Xi]g23_},fTow_}] calculates tow properties from the input parameters, outputs {Axial Tow Properties,Matrix Properties,Transverse Tow Properties}";
calculateTowPropertiesNew::usage = "calculateTowPropertiesNew[{Em_,Gm_,\[Nu]m_},{{Ef11_,Ef22_},{Gf12_,Gf23_},{\[Nu]f12_,\[Nu]f23_},{\[Xi]e22_,\[Xi]g12_,\[Xi]g23_},fTow_}] calculates tow properties from the input parameters, outputs {Axial Tow Properties,Matrix Properties,Transverse Tow Properties}";

(*Import FE Data*)
importNodalValues::usage = "importNodalValues[inpPath_,rptPath_] imports nodal coords, instance labels, and node sets from inp file and obd values from rpt file.  Output is {node value, x, y, z, rpt values, instance, node set}";
importNsets::usage = "importNsets[inpPath_] imports nodal sets with corresponding instance and nodes associated with each";
extractNset::usage = "extractNset[inpValues_,nSetIn_] extracts nSetIn from inpValues";
importRPT::usage = "importRPT[pathIn_] imports .rpt files from abaqus outputing list of experimental value for each part, Output={Legends,Parts}";
importPeriodicLD::usage = "importPeriodicLD[pathIn_,loadDir_] imports U and RF for the loadDir of the Fict node used in periodic models";
extractLoadDisp::usage = "extractLoadDisp[rptIn_,nSetIn_,loadDirection_] extracts and calculates {Mean[U],Total[RF]} from the given RPT and nSet in the loadDirection";
createNodalValues::usage = "createNodalValues[rptIn_,nodesIn_] combines nodal values {node,x,y,z} with rpt values to output {label,part}.";
extractNodeSet::usage = "extractNodeSet[rptIn_,nSetIn_,nodesIn_] extracts nodes in the nSetIn from rptIn and nodesIn. Output is {label,nodeValues}, node values could be segmented by part/section.";

(*Filter Strains*)
gridData::usage = "gridData[dataIn_,stepSize_,] regrids dataIn onto a grid of spacing stepSize.";
forwardDifference::usage = "forwardDiff[gridIn_,vector_] calculates the forward difference for a grid of data in the veritcal{1,0} or horizontal{0,1} direction";
filterStrains::usage = "filterStrains[dataIn_,stepSize_,filterSize_] calculates {\[Epsilon]xx,\[Epsilon]yy,\[Gamma]xy} from the given dispGridIn and filters using a Gaussian Filter of filterSize, outputs";

(*Line Scans*)
LineScans::usage = "oneVarLineScan[dataIn_,varListIn_,{{x1_,y1_},{x2_,y2_}},{var1_}] returns {x,y,var} from all nodes on line from {x1,y1} to {x2,y2}";

(*Get Lines*)
node1::usage = "node with {x1,y1}";
node2::usage = "node with {x2,y2}";
linesOut::usage = "List of coordinates for {node1,node2} in the camera array ({0,0} is upper left corner)";
getLines::usage = "getCoordinates[filesPath_,frame_,varIn_] allows for selection of nodal locations";


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Modify INP*)


getMaterialOrder[fileIn_]:=Module[{textListIn,materialPositions,magnitude,propertiesIn},
textListIn=ReadList[fileIn,String];
materialPositions=#[[1]]&/@Position[Map[StringPosition[#,"*Material,"]&,textListIn],1,Infinity];

Map[textListIn[[#]]&,materialPositions]

];


selectMaterialProperties[fileIn_,{AxialPropsIn_,TransversePropsIn_,MatrixPropsIn_},fileName_]:=Module[{textListIn,materialPositions,magnitude,propertiesIn},
textListIn=ReadList[fileIn,String];
materialPositions=#[[1]]&/@Position[Map[StringPosition[#,"*Material,"]&,textListIn],1,Infinity];
magnitude=ChoiceDialog["Modulus Units",{"Pa"->10.^9,"MPa"->10^3,"GPa"->1}];

Do[
propertiesIn=ChoiceDialog["Properties for " <>textListIn[[materialPositions[[i]]]]<>" "<>textListIn[[materialPositions[[i]]+1]],{"Axial Properties"->AxialPropsIn,"Transverse Properties"->TransversePropsIn,"Matrix Properties"->MatrixPropsIn}];
If[Length[propertiesIn]>2,
textListIn[[materialPositions[[i]]+2]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[1]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[2]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[3]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[4]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[5]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[6]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[7]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[8]]*magnitude],3]];
textListIn[[materialPositions[[i]]+3]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[9]]*magnitude],3]]<>",";
,
textListIn[[materialPositions[[i]]+2]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[1]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[2]]],3]];
];
,{i,3}];

Export[FileNameJoin[{DirectoryName[fileIn],fileName<>".inp"}],textListIn,"Text"]
];


changeMaterialProperties[fileIn_,PropsIn_,magnitude_,fileName_]:=Module[{textListIn,materialPositions,propertiesIn},
textListIn=ReadList[fileIn,String];
materialPositions=#[[1]]&/@Position[Map[StringPosition[#,"*Material,"]&,textListIn],1,Infinity];

Do[
propertiesIn=PropsIn[[i]];
If[Length[propertiesIn]>2,
textListIn[[materialPositions[[i]]+2]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[1]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[2]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[3]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[4]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[5]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[6]]],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[7]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[8]]*magnitude],3]];
textListIn[[materialPositions[[i]]+3]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[9]]*magnitude],3]]<>",";
,
textListIn[[materialPositions[[i]]+2]]=" "<>ToString[SetPrecision[FortranForm[propertiesIn[[1]]*magnitude],3]]<>", "<>ToString[SetPrecision[FortranForm[propertiesIn[[2]]],3]];
];
,{i,Length[PropsIn]}];

Export[FileNameJoin[{DirectoryName[fileIn],fileName<>".inp"}],textListIn,"Text"]
];


\[Eta][fiberProp_,matrixProp_,\[Xi]_]:=N[((fiberProp/matrixProp)-1)/((fiberProp/matrixProp)+\[Xi])];


calculateTowProperties[{Em_,Gm_,\[Nu]m_},{{Ef11_,Ef22_},{Gf12_,Gf23_},{\[Nu]f12_,\[Nu]f23_},{\[Xi]e22_,\[Xi]g12_,\[Xi]g23_},fTow_}]:=Module[
{Et11,Et22,Gt12,Gt23,\[Nu]t12,\[Nu]t21,\[Nu]t23,AxialTowProperties,TransverseTowProperties,MatrixProperties},

Et11=Ef11*fTow+Em*(1-fTow);
Et22=Em*(1+\[Xi]e22*\[Eta][Ef22,Em,\[Xi]e22]*fTow)/(1-\[Eta][Ef22,Em,\[Xi]e22]*fTow);
Gt12=Gm*(1+\[Xi]g12*\[Eta][Gf12,Gm,\[Xi]g12]*fTow)/(1-\[Eta][Gf12,Gm,\[Xi]g12]*fTow);
Gt23=Gm*(1+\[Xi]g23*\[Eta][Gf23,Gm,\[Xi]g23]*fTow)/(1-\[Eta][Gf23,Gm,\[Xi]g23]*fTow);
\[Nu]t12=\[Nu]f12*fTow+\[Nu]m*(1-fTow);
\[Nu]t21=(\[Nu]t12/Et11)*Et22;
\[Nu]t23=Et22/(2*Gt23)-1;

AxialTowProperties={Et11,Et22,Et22,\[Nu]t12,\[Nu]t12,\[Nu]t23,Gt12,Gt12,Gt23};
TransverseTowProperties={Et22,Et22,Et11,\[Nu]t23,\[Nu]t21,\[Nu]t21,Gt23,Gt23,Gt12};
MatrixProperties={Em,\[Nu]m};

{AxialTowProperties,MatrixProperties,TransverseTowProperties}

];


calculateTowPropertiesNew[{E1_,G1_,\[Nu]1_},{{Ea2_,Et2_},{Ga2_,Gt2_},{va2_,vt2_},c2_}]:=Module[
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
(*Import FE Data*)


importNodalValues[inpPath_,rptPath_]:=Module[
{inpIn,rptIn,inpInstances,nodalCoords,rptInstances,nodalRPTvalues,nodalValues,nSets,entryLength,instanceNearest},

inpIn=Import[inpPath,"CSV"];
rptIn=Select[Import[rptPath,"Table"],(Length[#]>1)&];

inpInstances=Map[StringTrim[#,(" name=\""|"\"")]&,Extract[inpIn,Position[inpIn[[All,1]],"*Instance"]][[All,2]]];
nodalCoords=MapThread[inpIn[[First[#1]+1;;First[#2]-1]]&,{Position[inpIn[[All,1]],"*Node"],Position[inpIn[[All,1]],"*Element"]}];

rptInstances=Map[StringSplit[StringReplace[StringTrim[ToString[#],("{"|"}")],","->""],": "][[2]]&,Extract[rptIn,Position[rptIn[[All,1]],"Field"][[2;;]]]];
nodalRPTvalues=Extract[Map[Cases[rptIn[[#[[1]];;#[[2]]-1]],{_Integer,__}][[All,2;;]]&,Partition[Flatten[Position[rptIn[[All,1]],1]],2,1,1,0]],Map[First[Position[StringMatchQ[rptInstances,#,IgnoreCase->True],True]]&,inpInstances]];

nodalValues=MapThread[Join[#1,#2,ConstantArray[{#3},Length[#1]],2]&,{nodalCoords,nodalRPTvalues,inpInstances}];

nodalValues

];


importNsets[inpPath_]:=Module[
{inpIn,nSets,headerBoole},

inpIn=Import[inpPath,"CSV"];

nSets=Select[Map[{inpIn[[First[#]]],#}&,Position[inpIn[[All,1]],"*Nset"]],(Length[#[[1]]]==3&&#[[1,3]]!= " internal")&];
nSets=Map[{StringTrim[#[[1,2]],(" nset=\""|"\"")],StringTrim[#[[1,3]],(" instance=\""|"\"")],#[[2]]}&,nSets];

headerBoole=MapIndexed[{NumberQ[#1[[1]]],#2}&,inpIn];
nSets[[All,3]]=Map[Flatten[inpIn[[#]]]&,Map[#[[1]]+1;;First[Select[headerBoole[[#[[1]]+1;;]],(#[[1]]==False)&][[1,2]]]-1&,nSets[[All,3]]]];

nSets=GatherBy[nSets,(#[[1]])&];
Transpose[{nSets[[All,1,1]],nSets[[All,All,2;;]]}]

];


extractNset[inpValues_,nSetIn_]:=Module[
{nodeValues,instances,nSet},

If[SameQ[Length[Dimensions[inpValues]],2],
nodeValues=GatherBy[inpValues,Last];
instances=nodeValues[[All,1,-1]];
nodeValues=nodeValues[[All,All,1;;-2]];
,
instances=inpValues[[All,1,-1]];
nodeValues=inpValues[[All,All,1;;-2]];
];

If[SameQ[Length[nSetIn],2],
nSet=nSetIn[[2]];
nSet[[All,1]]=Flatten[Map[First[Position[instances,#]]&,nSet[[All,1]]]];
,
nSet=nSetIn;
nSet[[All,1]]=Flatten[Map[First[Position[instances,#]]&,nSet[[All,1]]]];
];

Flatten[Map[Extract[nodeValues[[#[[1]]]],Partition[#[[2]],1]]&,nSet],1]

];


Options[importRPT]={rptLabel->"Node"};


importRPT[pathIn_,OptionsPattern[]]:=Module[
{listIn,labelsPos,labels,tableN,partPos,partLabel,parts,data},

listIn=ReadList[pathIn,String];

Which[StringMatchQ[OptionValue[rptLabel],"N*",IgnoreCase->True]
,
labelsPos=Position[StringPosition[listIn,"Node          RF"],{_,_}][[All,1]];
labels=Map[StringSplit[listIn[[#]]]&,labelsPos];
tableN=Map[Length[#]&,labels];
,
StringMatchQ[OptionValue[rptLabel],"E*",IgnoreCase->True]
,
labelsPos=Position[StringPosition[listIn,"Element             Int"],{_,_}][[All,1]];
labels=Map[StringSplit[listIn[[#]]]&,labelsPos];
tableN=Map[Length[#]&,labels];
];

partPos=Transpose[{Position[StringPosition[listIn,"Field Output reported"],{_,_}][[All,1]],Position[StringPosition[listIn,"Minimum"],{_,_}][[All,1]]-1}];
partLabel=Nearest[partPos[[All,1]]->partPos];
partPos=Map[First[partLabel[#]]&,labelsPos];
parts=Map[StringTrim[StringSplit[listIn[[#]],":"][[-1]]]&,partPos[[All,1]]];
partPos=MapThread[#1+3;;#2&,{labelsPos,partPos[[All,2]]}];

data=Chop[MapThread[ReadList[StringToStream[StringJoin[listIn[[#1]]]],Table[Number,{#2}]]&,{partPos,tableN}]];

If[Length[data]==1,
{{First[parts],First[labels]},First[data]}
,
parts=Map[First[StringSplit[#,"."]]&,parts];
partPos=Map[Position[parts,#]&,Union[parts]];
partPos=Map[If[Length[#]==1,First[#],#]&,partPos];
labels=Map[Flatten[DeleteDuplicates[Extract[labels,#]]]&,partPos];
data=Map[Extract[data,#]&,partPos];

{MapThread[{#1,#2}&,{Union[parts],labels}],data}
]

];


importPeriodicLD[pathIn_,loadDir_]:=Module[
{listIn,fictPos,labelPos,label},

listIn=ReadList[pathIn,String];

fictPos=First[Position[StringPosition[listIn,"FICT-1"],{_,_}]][[1]];
labelPos=First[Position[StringPosition[listIn[[fictPos;;]],"Node"],{_,_}]][[1]]+fictPos-1;
label=StringSplit[listIn[[labelPos]]];

Extract[First[ReadList[StringToStream[StringJoin[listIn[[labelPos+3]]]],Table[Number,{Length[label]}]]],Map[First[Position[label,#]]&,{"U.U"<>ToString[loadDir],"RF.RF"<>ToString[loadDir]}]]
];


extractLoadDisp[dataIn_,nSetIn_]:=Module[
{rptFlat,nSetData},

If[Length[Dimensions[dataIn]]==2,
rptFlat=SortBy[dataIn,First];
,
rptFlat=SortBy[DeleteDuplicates[Flatten[dataIn,1]],First];
];

nSetData=Extract[rptFlat,nSetIn];

{Mean[nSetData[[All,3]]],Total[nSetData[[All,2]]]}

];


createNodalValues[dataIn_,nodesIn_]:=Module[
{},

Map[Join[Extract[nodesIn,Partition[#[[All,1]],1]],#[[All,2;;]],2]&,dataIn]

];


extractNodeSet[dataIn_,nSetIn_,nodesIn_]:=Module[
{nSet,label,rptValues,nodeValues,nodeNearest,partValues},

nSet=Flatten[nSetIn];

rptValues=Map[Join[Extract[nodesIn,Partition[#[[All,1]],1]],#[[All,2;;]],2]&,Select[dataIn,(Length[Intersection[#[[All,1]],nSet]]>0)&]];

If[Length[Dimensions[rptValues]]==3,
nodeValues=First[rptValues];
nodeNearest=Nearest[nodeValues[[All,1]]->nodeValues[[All,2;;]]];
nodeValues=Map[First[nodeNearest[#]]&,nSet];
,
nodeValues=Table[
partValues=rptValues[[i]];
nodeNearest=Nearest[partValues[[All,1]]->partValues[[All,2;;]]];
Map[First[nodeNearest[#]]&,Intersection[partValues[[All,1]],nSet]]
,{i,Length[rptValues]}];
];

nodeValues
];


(* ::Subsection::Closed:: *)
(*Strain Filter*)


gridData[dataIn_,stepSize_]:=Module[
{dataInterp,xyRange,dataOut},

dataInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,dataIn],InterpolationOrder->1];

xyRange=SortBy[dataIn[[All,1;;2]],Total][[{1,-1}]];
dataOut=Table[Quiet[{x,y,dataInterp[x,y]}]
,{y,xyRange[[1,2]],xyRange[[2,2]],stepSize},{x,xyRange[[1,1]],xyRange[[2,1]],stepSize}]

];


forwardDifference[gridIn_,vector_]:=Module[
{dim,paddedRange,paddedGrid,gridOut},
dim=Dimensions[gridIn];

Which[vector=={1,0},
paddedRange={1}~Join~Range[dim[[1]]]~Join~{dim[[1]]};
paddedGrid=gridIn[[paddedRange,All]];

gridOut=((paddedGrid[[3;;-1,All]]-paddedGrid[[2;;-2,All]])-(paddedGrid[[1;;-3,All]]-paddedGrid[[2;;-2,All]]))/2.;
,
vector=={0,1},
paddedRange={1}~Join~Range[dim[[2]]]~Join~{dim[[2]]};
paddedGrid=gridIn[[All,paddedRange]];

gridOut=((paddedGrid[[All,3;;-1]]-paddedGrid[[All,2;;-2]])-(paddedGrid[[All,1;;-3]]-paddedGrid[[All,2;;-2]]))/2.;
];

gridOut
];


filterStrains[dataIn_,stepSize_,filterSize_]:=Module[
{Uinterp,Vinterp,xyRange,gridData,sigmaClusters,holePos,dX,dY,dUdX,dVdY,dUdY,dVdX,strainData,filterR,clipRange,gf,xRange,yRange,xClip,yClip,gfX,gfY,window,weights,out},

Uinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,dataIn],InterpolationOrder->1];
Vinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[4]]}&,dataIn],InterpolationOrder->1];

xyRange=SortBy[dataIn[[All,1;;2]],Total][[{1,-1}]];
gridData=Table[Quiet[{x,y,Uinterp[x,y],Vinterp[x,y]}]
,{y,xyRange[[1,2]],xyRange[[2,2]],stepSize},{x,xyRange[[1,1]],xyRange[[2,1]],stepSize}];

dX=forwardDifference[gridData[[All,All,1]],{0,1}];
dY=forwardDifference[gridData[[All,All,2]],{1,0}];
dUdX=forwardDifference[gridData[[All,All,3]],{0,1}]/dX;
dVdY=forwardDifference[gridData[[All,All,4]],{1,0}]/dY;
dUdY=forwardDifference[gridData[[All,All,3]],{1,0}]/dY;
dVdX=forwardDifference[gridData[[All,All,4]],{0,1}]/dX;

strainData=ConstantArray[{0.,0.,0.,0.,0.},Dimensions[gridData][[1;;2]]];
strainData[[All,All,1;;2]]=gridData[[All,All,1;;2]];
strainData[[All,All,3]]=dUdX;
strainData[[All,All,4]]=dVdY;
strainData[[All,All,5]]=dVdX+dUdY;

filterR=Floor[filterSize/2];

If[filterR==0,
out=strainData;
,
clipRange=Transpose[{{1,1},Dimensions[strainData][[1;;2]]}];
gf=GaussianMatrix[{filterR,filterSize/4.}];

DistributeDefinitions[strainData,filterR,clipRange,gf];
out=ParallelTable[
If[strainData[[yi,xi,-1]]!=-1.,
xRange={xi-filterR,xi+filterR};
yRange={yi-filterR,yi+filterR};
xClip=Clip[xRange,clipRange[[2]]];
yClip=Clip[yRange,clipRange[[1]]];
gfX={1,-1}-(xRange-xClip);
gfY={1,-1}-(yRange-yClip);

window=Flatten[strainData[[yClip[[1]];;yClip[[2]],xClip[[1]];;xClip[[2]],3;;]],1];
weights=Flatten[gf[[gfY[[1]];;gfY[[2]],gfX[[1]];;gfX[[2]]]]];
If[MemberQ[window[[All,-1]],-1.],
weights[[Flatten[Position[window[[All,-1]],-1.]]]]=0;
];

strainData[[yi,xi,1;;2]]~Join~Map[#.weights/Total[weights]&,Transpose[window[[All,1;;3]]]]
,
strainData[[yi,xi]]
]
,{yi,clipRange[[1,1]],clipRange[[1,2]]},{xi,clipRange[[2,1]],clipRange[[2,2]]}];
];

out
];


(* ::Subsection::Closed:: *)
(*Line Scans*)


Options[LineScans]={};


LineScans[dataIn_,linesIn_,OptionsPattern[]]:=Module[{
step,x1,x2,y1,y2,xIn,yIn,lineIn,lineFit,x,dataxyNearest
},

step=EuclideanDistance[#[[1]],#[[2]]]&@@Nearest[dataIn[[All,1;;2]],dataIn[[1,1;;2]],2];

Table[
x1=linesIn[[i,1,1]];
x2=linesIn[[i,2,1]];
y1=linesIn[[i,1,2]];
y2=linesIn[[i,2,2]];

xIn=Sort[{x1,x2}];
yIn=Sort[{y1,y2}];

If[TrueQ[Abs[xIn[[1]]-xIn[[2]]]<Abs[yIn[[1]]-yIn[[2]]]],
lineIn=Reverse[{{x1,y1},{x2,y2}},2];
lineFit[x_]= Fit[lineIn, {1,x},x];
lineIn=Table[{lineFit[i],i},{i,yIn[[1]],yIn[[2]],step/10}];
,
lineIn={{x1,y1},{x2,y2}};
lineFit[x_]= Fit[lineIn, {1,x},x];
lineIn=Table[{i,lineFit[i]},{i,xIn[[1]],xIn[[2]],step/10}];
];

dataxyNearest=Nearest[dataIn[[All,{1,2}]]->Automatic];

Extract[dataIn,DeleteDuplicates[Table[dataxyNearest[lineIn[[i]]],{i,Length[lineIn]}]]]

,{i,1,Length[linesIn]}]
];


(* ::Subsection::Closed:: *)
(*Get Lines*)


node1={0,0};
node2={0,0};
linesOut={};


Options[getLines]={
ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,plotRange->All,color->"Rainbow",size->500,ReduceData->False,NumberofPoints->50000
};


getLines[dataIn_,OptionsPattern[]]:=Module[{
takeNumber,reducedData,badNodes,incr,contourList,colorF,colorFS,contPlt,greyPlt
},

If[OptionValue[ReduceData],
takeNumber=Round[(Length[dataIn])/OptionValue[NumberofPoints]];
reducedData=Take[dataIn,{1,-1,takeNumber}];
,
reducedData=dataIn;
]

If[ListQ[OptionValue[ContourRange]],
incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData,First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

Share[];

Deploy[
Grid[{{
LocatorPane[
Dynamic[{node1,node2}],
Graphics[{
contPlt[[1]],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[{node1,node2}]],
Orange,
Dynamic[Line/@linesOut],
Dynamic[Point/@Flatten[linesOut,1]]
},
PlotRange->OptionValue[plotRange],
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]}]},
{Button["Save coordinates",AppendTo[linesOut,{node1,node2}]]}}]]
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["StitchFunctions`*"],
{
"node1","node2","linesOut"
}];
EndPackage[ ]
