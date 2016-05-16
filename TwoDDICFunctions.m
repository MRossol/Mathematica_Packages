(* ::Package:: *)

BeginPackage["TwoDDICFunctions`"]
Unprotect@@Names["TwoDDICFunctions`*"];
ClearAll@@Names["TwoDDICFunctions*"];


(*Data Import*)
DICdataFiles::usage = "DICdataFiles[pathIn_] returns DIC files and variable list from v6 files folder";
allDataFiles::usage = "allDataFiles[pathIn_] returns DIC files and variable list from v6 files folder and get image files from parent directory";
importDICdata::usage = "importDICdata[pathIn_,OptionsPattern[]] extracts data for all files in pathIn.  Can specify which variables.";

(*Load Displacement Import*)
DICOutput::usage = "DICOutput[pathIn_,fullScale_] returns image number and load from the DIC .csv file";
MTSOutput::usage = "MTSOutput[pathIn_] returns all data from MTS .DAT file with format {Disp,Load,Strain,External Inputs}";

(*Line Scans*)
DIClineScan::usage = "DIClineScan[pathIn_,frames_,lines_,vars_,OptionsPattern[]] returns {x,y,vars} for all nodes on linesIn, with options to return metric or pixel {x,y} values ";
dispScan::usage = "dispScan[pathIn_,frames_,lines_,vars_,OptionsPattern[]] like DIClineSweep with rotation for Begely Peel Test";
DIClineSweep::usage = "DIClineSweep[pathIn_,frames_,lines_,vars_,OptionsPatern[]] returns {x,y,vars} for all lines parallel to linesIn[[1]] sweeping out to the point linesIn[[2]], with options to return metric or pixel {x,y} values";
codTool::usage = "CODdisp[pathIn_,frames_,lines_,vars_,OptionsPattern[]] calculates the change in distance between every point on two parallel lines. ";

(*Contour Plots*)
DICoverLay::usage = "DICoverLay[filesPath_,frame_,varIn_] renders contours of varIn overlayed on the DIC image for the given frame";
DICLegend::usage = "DICLegend[{{sMin_,sMax_},contours_,labels_},labelIn_] renders a legend corresponding to the DIC contour overlay with identical inputs";

(*GUIs*)
nodes::usage = "nodes with {{x1,y1},{x2,y2}}";
coordinates::usage = "List of coordinates for nodes in the camera array ({0,0} is upper left corner)";
coordinatesOut::usage = "List of coordinates for nodes in the real space({0,0} is bottom left corner)";
getCoordinates::usage = "getCoordinates[filesPath_,frame_,varIn_] allows for selection of nodal locations";
guiNodes::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3}}";
dispGUI::usage = "dispGUI[filesPath_,frame_,varIn_] allows for selection of line and end point for disp analysis";
getNodesGUI::usage = "getNodesGUI[filesPath_,frame_,varIn_,OptionsPattern[]]allows for selection of line and end point for disp analysis";
codGUI::usage = "codGUI[filesPath_,frame_,varIn_] allows for selection of line for disp analysis";

(*Extensometers*)
DICextensometer::usage = "DICextensometers[pathIn_,extsIn_,OptionsPattern[]] returns strains or displacemens for the given extsIn, option to return a single frame";

(*COD*)
codAnalysis::usage = "codAnalysis[dataIn_,varListIn_,{{x1_,y1_},{x2_,y2_}},input1_,input2_,method_] gives {{x,y,cod},{x,{COD slope fit},{original y-v data}}},input 1 = number of nodes used to calculate center (~ 2 subsets), inputs 2 = data excluded from discontinuity (~2/3 subset/magnification) ";
initarrays::usage = "initarrays[ydata_,crit1_] for codAnalysis";
findindices::usage = "findindices[ydata_,ycenter_,crit_] for codAnalysis";
slopefind::usage = "slopefind[data_,ycenter_,indices_] finds cod for codAnalysis";
findcent::usage = "findcent[data_,ans_,coeff_,indices_] finds center of crack for codAnalysis";
workhorsesub::usage = "workhorsesub[data_,crit2_,ans_] for codAnalysis";
slopefind2::usage = "slopefind2[data_,fitpoints_,crit_,mypseudo_] finds cod for codAnalysis";


Begin["`Private`"]


Needs["GUIKit`"]


ProgressDialog[]:=GUIRun[ 
Widget["Frame",{ 
WidgetGroup[{
Widget["Label",{"text"->"Percent complete:"},Name->"label"],Widget["ProgressBar",
{"minimum"->0,"maximum"->100,
"preferredSize"->
  Widget["Dimension",{"width"->300,"height"->25}]},Name->"bar"] 
  }, WidgetLayout -> {
"Grouping"-> Column, 
"Border" -> {{15,15},{25,20}}}],
"location"->Widget["Point",{"x"->400,"y"->400}],"title"->"Computation Progress",
"resizable"->False},
Name->"frame"]
];


(* ::Subsection::Closed:: *)
(*Data Import*)


DICdataFiles[pathIn_]:=Module[{
matv6Path,varList,matv6Files
},

matv6Path=FileNames["*v6 files*",pathIn,Infinity];
If[Length[matv6Path]>1,
matv6Path=ChoiceDialog["Multiple data sources found...\n  Pick v6 .mat data location:",(FileBaseName[ParentDirectory[#]]->#)&/@matv6Path];,
matv6Path=matv6Path[[1]]
];

matv6Files=FileNames["*.mat",matv6Path];
varList=StringTrim/@(ReadList[FileNameJoin[{matv6Path,"index_list.txt"}],{Number,String}][[All,2]]);

{matv6Files,varList}
];


allDataFiles[pathIn_]:=Module[{
matv6Path,varList,matv6Files,imgFiles,imgPath
},

matv6Path=FileNames["*v6 files*",pathIn,Infinity];
If[Length[matv6Path]>1,
matv6Path=ChoiceDialog["Multiple data sources found...\n  Pick v6 .mat data location:",(FileBaseName[ParentDirectory[#]]->#)&/@matv6Path];,
matv6Path=matv6Path[[1]]
];

matv6Files=FileNames["*.mat",matv6Path];
imgPath=ParentDirectory[matv6Path];
imgFiles=Flatten[FileNames[FileBaseName[#]<>".tif",imgPath]&/@matv6Files];
varList=StringTrim/@(ReadList[FileNameJoin[{matv6Path,"index_list.txt"}],{Number,String}][[All,2]]);

{matv6Files,varList,imgFiles}
];


 Options[importDICdata]={
GetFiles->All,GetVariables->All
};


importDICdata[pathIn_,OptionsPattern[]]:=Module[{
filesIn,variables},

filesIn=DICdataFiles[pathIn];

If[SameQ[OptionValue[GetVariables],All],
variables=All;
,
variables=Map[Position[filesIn[[2]],#][[1,1]]&,OptionValue[GetVariables]];
];

If[NumberQ[OptionValue[GetFiles]],
MapThread[List,Import[filesIn[[1,OptionValue[GetFiles]]]],2][[All,All,variables]]
,
Map[MapThread[List,Import[#],2][[All,All,variables]]&,filesIn[[1,OptionValue[GetFiles]]]]
]

];


(* ::Subsection::Closed:: *)
(*Load Displacement Import*)


Options[DICOutput]={
Channel->1
};


DICOutput[pathIn_,fullScale_,OptionsPattern[]]:=Module[{
dataIn,i,imageNum,load
},

dataIn=Import[pathIn];
If[NumberQ[dataIn[[1,1]]],i=1,i=2];
imageNum=dataIn[[i;;,1]];
load=dataIn[[i;;,(5+OptionValue[Channel])]]*(fullScale/10);

Transpose[{imageNum,load}]
];


MTSOutput[pathIn_]:=Module[{
dataIn
},
(*
dataIn=Import[pathIn];

Table[DeleteCases[
If[NumberQ[#],#,
If[NumberQ[ToExpression[StringReplace[StringTrim[#,","],"E"->" 10^"]]],ToExpression[StringReplace[StringTrim[#,","],"E"->" 10^"]]]]&/@dataIn[[i]],Null],
{i,Length[dataIn]}]
*)

Import[pathIn,"CSV"]
];


(* ::Subsection::Closed:: *)
(*Line Scans*)


Options[DIClineScan]={
CoordinateSystem->"Pixel",FlipYaxis->False,ForceVertical->False
};


DIClineScan[pathIn_,frames_,lines_,vars_,OptionsPattern[]]:=Module[{
filesIn,dataIn,varsIn,linesIn,varsOut,xyVars,ref,refData,step,nearestX,nearestXPos,nearestY,nearestYPos,xIn,yIn,lineFit,y,x,coords,positions,data,dataOut,yMax
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[vars]==0,
varsIn={vars};,
varsIn=vars;];

If[Length[Dimensions[lines]]==2,
linesIn={lines};,
linesIn=lines;];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x_c","y_c"}~Join~varsIn];,
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}~Join~varsIn];
];

xyVars=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}];

ref=ProgressDialog[];

refData=MapThread[List,Import[dataIn[[1]]],2];

step=refData[[1,2,xyVars[[1]]]]-refData[[1,1,xyVars[[1]]]];

nearestX=Nearest[refData[[1,All,xyVars[[1]]]]];
nearestXPos=Nearest[refData[[1,All,xyVars[[1]]]]->Automatic];
nearestY=Nearest[refData[[All,1,xyVars[[2]]]]];
nearestYPos=Nearest[refData[[All,1,xyVars[[2]]]]->Automatic];

positions=Table[

xIn=Sort[Map[#[[1]]&,linesIn[[j]]]];
yIn=Sort[Map[#[[2]]&,linesIn[[j]]]];

If[TrueQ[First[Differences[xIn]]<First[Differences[yIn]]],
lineFit[y_]= Chop[Fit[Reverse[linesIn[[j]],2], {1,y},y]];
coords=Table[{Round[lineFit[y]],y},
{y,First[nearestY[yIn[[1]]]],First[nearestY[yIn[[2]]]],step}];
,
lineFit[x_]= Chop[Fit[linesIn[[j]], {1,x},x]];
coords=Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xIn[[1]]]],First[nearestX[xIn[[2]]]],step}];
];

Map[{First[nearestYPos[#[[2]]]],First[nearestXPos[#[[1]]]]}&,coords]

,{j,Length[linesIn]}];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

dataOut=Map[Extract[data[[All,All,varsOut]],#]&,positions];

If[MemberQ[vars,"v"]&&OptionValue[FlipYaxis],
dataOut[[All,(Position[vars,"v"][[1,1]]+2)]]=-1.*dataOut[[All,(Position[vars,"v"][[1,1]]+2)]];
,
dataOut=dataOut;
];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
dataOut
,
If[OptionValue[FlipYaxis],
If[OptionValue[ForceVertical],
yMax=2449;,
If[data[[-1,-1,xyVars[[2]]]]>2049,
yMax=2449;,yMax=2049;];];
dataOut[[All,2]]=yMax-dataOut[[All,2]];
dataOut
,
dataOut
]
]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[NumberQ[frames],
If[First[Dimensions[linesIn]]==1,
dataOut[[1,1]]
,
dataOut[[1]]
]
,
If[First[Dimensions[linesIn]]==1,
dataOut[[All,1]]
,
dataOut
]
]

];


Options[dispScan]={
CoordinateSystem->"Pixel",FlipYaxis->False,ForceVertical->False,GetStDev->False
};


dispScan[pathIn_,frames_,lines_,OptionsPattern[]]:=Module[{
filesIn,dataIn,varsIn,linesIn,varsOut,xyVars,ref,dataExport,data,step,nearestX,nearestXPos,nearestY,nearestYPos,line1,angle,midPoint,line2,fitY,rot,fit,y,x,m,b,b1,m1,yRange,xRange,lFit,lCoords,lineFit,coords,positions,dataOut,yMax,coordsOut
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[Dimensions[lines]]==2,
linesIn=lines;,
linesIn=lines[[1]];];

xyVars=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x_c","y_c","u_c","v_c"}];
,
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y","u","v"}];
];

ref=ProgressDialog[];

dataOut=Table[

data=MapThread[List,Import[dataIn[[i]]],2];

step=data[[1,2,xyVars[[1]]]]-data[[1,1,xyVars[[1]]]];

nearestX=Nearest[data[[1,All,xyVars[[1]]]]];
nearestXPos=Nearest[data[[1,All,xyVars[[1]]]]->Automatic];
nearestY=Nearest[data[[All,1,xyVars[[2]]]]];
nearestYPos=Nearest[data[[All,1,xyVars[[2]]]]->Automatic];

line1=Sort[linesIn[[1]]];
angle=ArcTan@@MapThread[#2-#1&,line1]/Degree-90.;
midPoint={Mean[line1[[All,1]]],Mean[line1[[All,2]]]};
line2=Sort[{midPoint,{midPoint[[1]]+10Cos[angle Degree],midPoint[[2]]+10Sin[angle Degree]}}];
fitY=TrueQ[First[Differences[Sort[line2[[All,1]]]]]<First[Differences[Sort[line2[[All,2]]]]]];

If[fitY,
fit=Chop[FindFit[Reverse[line2,2], b+m y,{b,m},y]];
lineFit[y_]=b+m y/.fit;

angle=ArcTan@@MapThread[#2-#1&,line2]/Degree;
If[angle>0,
rot=RotationMatrix[(-90+angle)Degree];
,
rot=RotationMatrix[(90+angle)Degree];
];

yRange=Sort[{midPoint[[2]],linesIn[[2,2]]}];
coords=Table[{Round[lineFit[y]],y},
{y,First[nearestY[yRange[[1]]]],First[nearestY[yRange[[2]]]],step}];
,
fit=Chop[FindFit[line2, b+m x,{b,m},x]];
lineFit[x_]=b+m x/.fit;

angle=ArcTan@@MapThread[#2-#1&,line2]/Degree;
rot=RotationMatrix[angle Degree];

xRange=Sort[{midPoint[[1]],linesIn[[2,1]]}];
coords=Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xRange[[1]]]],First[nearestX[xRange[[2]]]],step}];
];

Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(j/Length[coords])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing line "<>ToString[j]<>" of "<>ToString[Length[coords]]<>" for frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

If[fitY,
lFit={b1->y+m x,m1->-m}/.{x->coords[[j,1]],y->coords[[j,2]],fit[[2]]};
lineFit[x_]=b1+m1 x/.lFit;

xRange=Sort[line1[[All,1]]];
lCoords=Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xRange[[1]]]],First[nearestX[xRange[[2]]]],step}];
,
lFit={b1->x+m y,m1->-m}/.{x->coords[[j,1]],y->coords[[j,2]],fit[[2]]};
lineFit[y_]=b1+m1 y/.lFit;

yRange=Sort[line1[[All,2]]];
lCoords=Table[{Round[lineFit[y]],y},
{y,First[nearestY[yRange[[1]]]],First[nearestY[yRange[[2]]]],step}];
];

positions=Map[{First[nearestYPos[#[[2]]]],First[nearestXPos[#[[1]]]]}&,lCoords];

dataOut=Map[rot.#&,Extract[data[[All,All,varsOut[[3;;4]]]],positions]];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
coordsOut=rot.data[[First[nearestYPos[coords[[j,2]]]],First[nearestXPos[coords[[j,1]]]],varsOut[[1;;2]]]];
If[OptionValue[GetStDev],
dataOut=coordsOut~Join~{Mean[dataOut[[All,1]]],Mean[dataOut[[All,2]]],StandardDeviation[dataOut[[All,1]]],StandardDeviation[dataOut[[All,2]]]}
,
dataOut=coordsOut~Join~{Mean[dataOut[[All,1]]],Mean[dataOut[[All,2]]]}
]
,
coordsOut=rot.coords[[j]];
If[OptionValue[GetStDev],
dataOut=coordsOut~Join~{Mean[dataOut[[All,1]]],Mean[dataOut[[All,2]]],StandardDeviation[dataOut[[All,1]]],StandardDeviation[dataOut[[All,2]]]};
,
dataOut=coordsOut~Join~{Mean[dataOut[[All,1]]],Mean[dataOut[[All,2]]]};
];
If[OptionValue[FlipYaxis],
If[OptionValue[ForceVertical],
yMax=2449;,
If[data[[-1,-1,xyVars[[2]]]]>2049,
yMax=2449;,yMax=2049;];];
dataOut[[2]]=yMax-dataOut[[2]];
dataOut[[4]]=-1.*dataOut[[4]];
dataOut
,
dataOut
]
]

,{j,Length[coords]}]
,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref]; 

If[NumberQ[frames],
dataOut[[1]]
,
dataOut
]

];


Options[DIClineSweep]={
CoordinateSystem->"Pixel",FlipYaxis->False,ForceVertical->False
};


DIClineSweep[pathIn_,frames_,lines_,vars_,OptionsPattern[]]:=Module[{
filesIn,dataIn,varsIn,linesIn,varsOut,xyVars,ref,refData,step,nearestX,nearestXPos,nearestY,nearestYPos,line1,angle,midPoint,line2,fitY,fit,y,x,m,b,b1,m1,yRange,xRange,lFit,lCoords,lineFit,coords,positions,data,dataOut,yMax
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[vars]==0,
varsIn={vars};,
varsIn=vars;];

If[Length[Dimensions[lines]]==2,
linesIn={lines};,
linesIn=lines;];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x_c","y_c"}~Join~varsIn];,
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}~Join~varsIn];
];

xyVars=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}];

ref=ProgressDialog[];


refData=MapThread[List,Import[dataIn[[1]]],2];

step=refData[[1,2,xyVars[[1]]]]-refData[[1,1,xyVars[[1]]]];

nearestX=Nearest[refData[[1,All,xyVars[[1]]]]];
nearestXPos=Nearest[refData[[1,All,xyVars[[1]]]]->Automatic];
nearestY=Nearest[refData[[All,1,xyVars[[2]]]]];
nearestYPos=Nearest[refData[[All,1,xyVars[[2]]]]->Automatic];

positions=Table[

line1=Sort[linesIn[[j,1]]];
angle=ArcTan@@MapThread[#2-#1&,line1]/Degree-90.;
midPoint={Mean[line1[[All,1]]],Mean[line1[[All,2]]]};
line2=Sort[{midPoint,{midPoint[[1]]+10Cos[angle Degree],midPoint[[2]]+10Sin[angle Degree]}}];
fitY=TrueQ[First[Differences[Sort[line2[[All,1]]]]]<First[Differences[Sort[line2[[All,2]]]]]];

If[fitY,
fit=Chop[FindFit[Reverse[line2,2], b+m y,{b,m},y]];
lineFit[y_]=b+m y/.fit;

yRange=Sort[{midPoint[[2]],linesIn[[j,2,2]]}];
coords=Table[{Round[lineFit[y]],y},
{y,First[nearestY[yRange[[1]]]],First[nearestY[yRange[[2]]]],step}];
,
fit=Chop[FindFit[line2, b+m x,{b,m},x]];
lineFit[x_]=b+m x/.fit;

xRange=Sort[{midPoint[[1]],linesIn[[j,2,1]]}];
coords=Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xRange[[1]]]],First[nearestX[xRange[[2]]]],step}];
];

Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(k/Length[coords])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing line "<>ToString[k]<>" of "<>ToString[Length[coords]]<>" in area "<>ToString[j]<>" of "<>ToString[Length[linesIn]]];

If[fitY,
lFit={b1->y+m x,m1->-m}/.{x->coords[[k,1]],y->coords[[k,2]],fit[[2]]};
lineFit[x_]=b1+m1 x/.lFit;

xRange=Sort[line1[[All,1]]];
lCoords=Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xRange[[1]]]],First[nearestX[xRange[[2]]]],step}];
,
lFit={b1->x+m y,m1->-m}/.{x->coords[[k,1]],y->coords[[k,2]],fit[[2]]};
lineFit[y_]=b1+m1 y/.lFit;

yRange=Sort[line1[[All,2]]];
lCoords=Table[{Round[lineFit[y]],y},
{y,First[nearestY[yRange[[1]]]],First[nearestY[yRange[[2]]]],step}];
];

Map[{First[nearestYPos[#[[2]]]],First[nearestXPos[#[[1]]]]}&,lCoords]

,{k,Length[coords]}]
,{j,Length[linesIn]}];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

dataOut=Map[Extract[data[[All,All,varsOut]],#]&,positions,{2}];

If[MemberQ[vars,"v"]&&OptionValue[FlipYaxis],
dataOut[[All,(Position[vars,"v"][[1,1]]+2)]]=-1.*dataOut[[All,(Position[vars,"v"][[1,1]]+2)]];
,
dataOut=dataOut;
];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
dataOut
,
If[OptionValue[FlipYaxis],
If[OptionValue[ForceVertical],
yMax=2449;,
If[data[[-1,-1,xyVars[[2]]]]>2049,
yMax=2449;,yMax=2049;];];
dataOut[[All,2]]=yMax-dataOut[[All,2]];
dataOut
,
dataOut
]
]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[NumberQ[frames],
If[First[Dimensions[linesIn]]==1,
dataOut[[1,1]]
,
dataOut[[1]]
]
,
If[First[Dimensions[linesIn]]==1,
dataOut[[All,1]]
,
dataOut
]
]

];


Options[codTool]={
CoordinateSystem->"Metric"
};


codTool[pathIn_,frames_,lines_,OptionsPattern[]]:=Module[{
filesIn,dataIn,linesIn,varsOut,xyVars,ref,refData,step,nearestX,nearestXPos,nearestY,nearestYPos,xIn,yIn,lineFit,y,x,coords,positions,take,shorten,lo,data,dataOut,yMax
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[Dimensions[lines]]==4,
linesIn=lines[[1]];
,
linesIn=lines;
];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x_c","y_c","u_c","v_c"}];
,
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y","u","v"}];
];

xyVars=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}];

ref=ProgressDialog[];

refData=MapThread[List,Import[dataIn[[1]]],2];

step=refData[[1,2,xyVars[[1]]]]-refData[[1,1,xyVars[[1]]]];

nearestX=Nearest[refData[[1,All,xyVars[[1]]]]];
nearestXPos=Nearest[refData[[1,All,xyVars[[1]]]]->Automatic];
nearestY=Nearest[refData[[All,1,xyVars[[2]]]]];
nearestYPos=Nearest[refData[[All,1,xyVars[[2]]]]->Automatic];

coords=Table[

xIn=Sort[Map[#[[1]]&,linesIn[[j]]]];
yIn=Sort[Map[#[[2]]&,linesIn[[j]]]];

If[TrueQ[First[Differences[xIn]]<First[Differences[yIn]]],
lineFit[y_]= Chop[Fit[Reverse[linesIn[[j]],2], {1,y},y]];
Table[{Round[lineFit[y]],y},
{y,First[nearestY[yIn[[1]]]],First[nearestY[yIn[[2]]]],step}]
,
lineFit[x_]= Chop[Fit[linesIn[[j]], {1,x},x]];
Table[{x,Round[lineFit[x]]},
{x,First[nearestX[xIn[[1]]]],First[nearestX[xIn[[2]]]],step}]
]

,{j,Length[linesIn]}];

If[Length[coords[[1]]]==Length[coords[[2]]],
positions=Map[{First[nearestYPos[#[[2]]]],First[nearestXPos[#[[1]]]]}&,coords,{2}];
,
take=Switch[EuclideanDistance@@coords[[All,1]]<EuclideanDistance@@coords[[All,-1]],True,1;;-2,False,2;;-1];
shorten=Switch[Length[coords[[1]]]<Length[coords[[2]]],True,2,False,1];
coords[[shorten]]=coords[[shorten,take]];
positions=Map[{First[nearestYPos[#[[2]]]],First[nearestXPos[#[[1]]]]}&,coords,{2}];
];

refData=Map[Extract[refData[[All,All,varsOut]],#]&,positions];
refData=Map[{#[[1]]+#[[3]],#[[2]]+#[[4]]}&,refData,{2}];
lo=MapThread[{EuclideanDistance[refData[[1,-1]],#1],EuclideanDistance[#1,#2]}&,refData];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

data=Map[Extract[data[[All,All,varsOut]],#]&,positions];
data=Map[{#[[1]]+#[[3]],#[[2]]+#[[4]]}&,data,{2}];
MapThread[{#3[[1]],(EuclideanDistance[#1,#2]-#3[[2]])}&,data~Join~{lo}]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

dataOut

];


(* ::Subsection::Closed:: *)
(*Contour Plotting Functions*)


Options[DICoverLay]={
ContourOpacity->1,contourShading->Automatic,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,plotRange->All,color->"Rainbow",NumberofPoints->50000,size->500,overlay->True,NormalDisp->False,FixRotation->False,FilterSize->1
};


DICoverLay[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,uNumb,vNumb,UNumb,VNumb,varNumb,sigmaNumb,dataIn,yMax,takeNumber,xCenter,yCenter,uCenterLine,vCenterLine,uLine,vLine,m,b,theta,meanDisp,reducedData,meanSigma,incr,contourList,colorFS,colorF,contPlt,greyPlt
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
uNumb=Position[filesIn[[2]],"u"][[1,1]];
vNumb=Position[filesIn[[2]],"v"][[1,1]];
UNumb=Position[filesIn[[2]],"u_c"][[1,1]];
VNumb=Position[filesIn[[2]],"v_c"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn[[All,All,{xNumb,yNumb}]]=dataIn[[All,All,{xNumb,yNumb}]]+dataIn[[All,All,{uNumb,vNumb}]];

Which[varIn=="u_c",
If[OptionValue[FixRotation],
xCenter=Round[Length[dataIn[[1]]]/2];
yCenter=Round[Length[dataIn]/2];
uCenterLine=Delete[dataIn[[All,xCenter,{yNumb,UNumb}]],Position[dataIn[[All,xCenter,sigmaNumb]],-1.]];
vCenterLine=Delete[dataIn[[yCenter,All,{xNumb,VNumb}]],Position[dataIn[[yCenter,All,sigmaNumb]],-1.]];
uLine=Function[{y},Evaluate[m y +b/.FindFit[uCenterLine,m y+b,{m,b},y]]];
vLine=Function[{x},Evaluate[m x +b/.FindFit[vCenterLine,m x +b,{m,b},x]]];
theta=(Abs[ArcTan[(vLine[vCenterLine[[-1,1]]]-vLine[vCenterLine[[1,1]]])/(vCenterLine[[-1,1]]-vCenterLine[[1,1]])]]+Abs[ArcTan[(uLine[uCenterLine[[-1,1]]]-uLine[uCenterLine[[1,1]]])/(uCenterLine[[-1,1]]-uCenterLine[[1,1]])]])/2;

dataIn[[All,All,UNumb]]=dataIn[[All,All,UNumb]]-theta*dataIn[[All,All,yNumb]];
];

If[OptionValue[CenterDisp],
meanDisp=Mean[Select[Flatten[dataIn,1],(#[[sigmaNumb]]!=-1.)&][[All,UNumb]]];
dataIn[[All,All,UNumb]]=Map[#-meanDisp&,dataIn[[All,All,UNumb]],{2}];
];
,
varIn=="v_c",
If[OptionValue[FixRotation],
xCenter=Round[Length[dataIn[[1]]]/2];
yCenter=Round[Length[dataIn]/2];
uCenterLine=Delete[dataIn[[All,xCenter,{yNumb,UNumb}]],Position[dataIn[[All,xCenter,sigmaNumb]],-1.]];
vCenterLine=Delete[dataIn[[yCenter,All,{xNumb,VNumb}]],Position[dataIn[[yCenter,All,sigmaNumb]],-1.]];
uLine=Function[{y},Evaluate[m y +b/.FindFit[uCenterLine,m y+b,{m,b},y]]];
vLine=Function[{x},Evaluate[m x +b/.FindFit[vCenterLine,m x +b,{m,b},x]]];
theta=(Abs[ArcTan[(vLine[vCenterLine[[-1,1]]]-vLine[vCenterLine[[1,1]]])/(vCenterLine[[-1,1]]-vCenterLine[[1,1]])]]+Abs[ArcTan[(uLine[uCenterLine[[-1,1]]]-uLine[uCenterLine[[1,1]]])/(uCenterLine[[-1,1]]-uCenterLine[[1,1]])]])/2;

dataIn[[All,All,VNumb]]=dataIn[[All,All,VNumb]]+theta*dataIn[[All,All,xNumb]];
];

If[OptionValue[CenterDisp],
meanDisp=Mean[Select[Flatten[dataIn,1],(#[[sigmaNumb]]!=-1.)&][[All,VNumb]]];
dataIn[[All,All,VNumb]]=Map[#-meanDisp&,dataIn[[All,All,VNumb]],{2}];
];
,
varIn=="Z",
If[OptionValue[DeformedZ],
dataIn[[All,All,varNumb]]=dataIn[[All,All,varNumb]]+dataIn[[All,All,WNumb]];
];
];

dataIn=Map[{#[[xNumb]],yMax-#[[yNumb]],#[[varNumb]],#[[sigmaNumb]]}&,dataIn,{2}];

meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
reducedData=Select[Flatten[dataIn,1],(#[[4]]!=-1.)&];

If[OptionValue[NumberofPoints]<Length[reducedData],
takeNumber=Round[(Length[reducedData])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[reducedData,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{Chop[#],OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{Chop[#],OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(First[meanSigma[{#1,#2}]]>0.&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]},
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

If[OptionValue[overlay],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]]
},
PlotRange->OptionValue[plotRange],
PlotRangePadding->None,
AspectRatio->Automatic,
ImageSize->OptionValue[size]]
,
Graphics[{
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]]
},
PlotRange->OptionValue[plotRange],
PlotRangePadding->None,
AspectRatio->Automatic,
ImageSize->OptionValue[size]]
]
];


Options[DICLegend]={
color->"Rainbow",size->{Automatic,500},fontSize->16,textSide->"Right"
};


DICLegend[{{sMin_,sMax_},contours_,lines_,labels_},labelIn_,OptionsPattern[]]:=Module[{
Acolor,
incr=N[(sMax-sMin)/contours],
incr2=N[(sMax-sMin)/lines],
incr3=N[(sMax-sMin)/labels],
AColor,rSet,Amesh,myTxt,myLabel
},

Acolor=ColorData[OptionValue[color]][#]&/@Rescale[Range[sMin+incr,sMax,incr]];
rSet=Rectangle[{0,#},{0.5,#+incr 10/(sMax-sMin)}]&/@Drop[10Rescale[Range[sMin,sMax,incr]],-1];
Amesh=Rescale[Range[sMin+incr2,sMax-incr2,incr2],{sMin,sMax},{0,10}];
If[StringMatchQ[OptionValue[textSide],"r*",IgnoreCase->True],
myLabel=Text[Style[labelIn,FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0,10.5},{-1,0},{1,0}];
myTxt=Text[Style[ToString[Chop[#]],FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0.7,Rescale[#,{sMin,sMax},{0,10}]},{-1,0},{1,0}]&/@Range[sMin,sMax,incr3];
,
myLabel=Text[Style[labelIn,FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0.5,10.5},{1,0},{1,0}];
myTxt=Text[Style[ToString[Chop[#]],FontFamily->"Arial",FontSize->OptionValue[fontSize]],{-0.2,Rescale[#,{sMin,sMax},{0,10}]},{1,0},{1,0}]&/@Range[sMin,sMax,incr3];
];
Acolor=Take[Acolor,Length[rSet]];
rSet=Riffle[Acolor,rSet];

Graphics[{
GrayLevel[0.9],
Rectangle[{0,0},{0.5,10}],
rSet,

Black,AbsoluteThickness[2],
Line[{{0,0},{0.5,0},{0.5,10},{0,10},{0,0}}],
Line[{{0,#},{0.5,#}}]&/@Amesh,

myTxt,
myLabel
},
AspectRatio->Automatic,
ImageSize->OptionValue[size]]
];


(* ::Subsection::Closed:: *)
(*GUIs*)


nodes={{0,0},{0,0}};
coordinates={};
coordinatesOut={};


Options[getCoordinates]={
ContourOpacity->1,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,color->"Rainbow",NumberofPoints->30000,size->600,FilterSize->1
};


getCoordinates[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,varNumb,sigmaNumb,dataIn,yMax,takeNumber,reducedData,meanSigma,incr,contourList,colorF,colorFS,contPlt,greyPlt
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn=Map[{#[[xNumb]],(yMax-#[[yNumb]]),#[[varNumb]],#[[sigmaNumb]]}&,dataIn,{2}];

meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
reducedData=Select[Flatten[dataIn,1],(#[[4]]!=-1.)&];

If[OptionValue[NumberofPoints]<Length[reducedData],
takeNumber=Round[(Length[reducedData])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[reducedData,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,Opacity[0.5,Black]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,None}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(First[meanSigma[{#1,#2}]]>0.&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]},
ColorFunction->colorF,
ColorFunctionScaling->colorFS];

Deploy[
Grid[{{
LocatorPane[
Dynamic[nodes],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[nodes]],
Orange,
Dynamic[Line/@coordinates],
Dynamic[Point/@Flatten[coordinates,1]]
},
PlotRange->{{1,Length[greyPlt[[1]]]},{1,Length[greyPlt]}},
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]}]},
{Grid[{{Button["Clear coordinates",coordinates={};],
Button["Clear coordinatesOut",coordinatesOut={};]}}]},
{Button["Save coordinates",
yMax=Length[greyPlt]+1;
AppendTo[coordinates,nodes];
AppendTo[coordinatesOut,Map[{#[[1]],yMax-#[[2]]}&,nodes]];]}}]]
];


guiNodes={{0,0},{1,1},{2,2}};
coordinates={};
coordinatesOut={};


Options[dispGUI]={
ContourOpacity->1,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,color->"Rainbow",NumberofPoints->20000,size->600,magMuliplier->5.,FilterSize->1
};


dispGUI[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,varNumb,sigmaNumb,dataIn,yMax,step,meanSigma,takeNumber,reducedData,badNodes,incr,contourList,colorF,colorFS,contPlt,greyPlt
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn=Map[{#[[xNumb]],(yMax-#[[yNumb]]),#[[varNumb]],#[[sigmaNumb]]}&,dataIn,{2}];

meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
reducedData=Select[Flatten[dataIn,1],(#[[4]]!=-1.)&];

If[OptionValue[NumberofPoints]<Length[reducedData],
takeNumber=Round[(Length[reducedData])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[reducedData,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,Opacity[0.5,Black]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,None}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(First[meanSigma[{#1,#2}]]>0.&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]},
ColorFunction->colorF,
ColorFunctionScaling->colorFS];

Deploy[
Grid[{{Grid[{{
LocatorPane[
Dynamic[guiNodes],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}]],
Orange,
Dynamic[Line/@coordinates],
Dynamic[Point/@Flatten[coordinates,1]]
},
PlotRange->{{1,Length[greyPlt[[1]]]},{1,Length[greyPlt]}},
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["X",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->24],
Style["X",FontColor->RGBColor[0.4,0.05,0.1],FontSize->24],
Style["X",FontColor->RGBColor[0.05,0.1,0.4],FontSize->24]}]
,
Pane[LocatorPane[
Dynamic[guiNodes],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.0004],
Dynamic[Line[{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}]],
Orange,
Dynamic[Line/@coordinates],
Dynamic[Point/@Flatten[coordinates,1]]
},
PlotRange->{{1,Length[greyPlt[[1]]]},{1,Length[greyPlt]}},
AspectRatio->Automatic,
ImageSize->OptionValue[size]*OptionValue[magMuliplier]],
Appearance->{
Style["X",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->36],
Style["X",FontColor->RGBColor[0.4,0.05,0.1],FontSize->36],
Style["X",FontColor->RGBColor[0.05,0.1,0.4],FontSize->36]}],ImageSize->{OptionValue[size]},
Scrollbars->True,
ScrollPosition->{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]}]}}]},
{Grid[{{Button["Clear coordinates",coordinates={};],
Button["Clear coordinatesOut",coordinatesOut={};]}}]},
{Button["Save coordinates",
yMax=Length[greyPlt]+1;
AppendTo[coordinates,{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}];
AppendTo[coordinatesOut,{frame,{Map[{#[[1]],(yMax-#[[2]])}&,guiNodes[[1;;2]]],{guiNodes[[3,1]],(yMax-guiNodes[[3,2]])}}}];]}}]]
];


Options[getNodesGUI]={
ContourOpacity->1,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,color->"Rainbow",NumberofPoints->20000,size->600,FilterSize->1
};


getNodesGUI[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,varNumb,sigmaNumb,dataIn,yMax,takeNumber,reducedData,meanSigma,incr,contourList,colorF,colorFS,contPlt,greyPlt
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn=Map[{#[[xNumb]],(yMax-#[[yNumb]]),#[[varNumb]],#[[sigmaNumb]]}&,dataIn,{2}];

meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
reducedData=Select[Flatten[dataIn,1],(#[[4]]!=-1.)&];

If[OptionValue[NumberofPoints]<Length[reducedData],
takeNumber=Round[(Length[reducedData])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[reducedData,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,Opacity[0.5,Black]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,None}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(First[meanSigma[{#1,#2}]]>0.&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]},
ColorFunction->colorF,
ColorFunctionScaling->colorFS];

Deploy[
Grid[{{
LocatorPane[
Dynamic[guiNodes],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}]],
Orange,
Dynamic[Line/@coordinates],
Dynamic[Point/@Flatten[coordinates,1]]
},
PlotRange->{{1,Length[greyPlt[[1]]]},{1,Length[greyPlt]}},
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12]
}]},
{Grid[{{Button["Clear coordinates",coordinates={};],
Button["Clear coordinatesOut",coordinatesOut={};]}}]},
{Button["Save coordinates",
yMax=Length[greyPlt]+1;
AppendTo[coordinates,
{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}];
AppendTo[coordinatesOut,{Map[{#[[1]],yMax-#[[2]]}&,guiNodes[[1;;2]]],{guiNodes[[3,1]],yMax-guiNodes[[3,2]]}}];]}}]]
];


Options[codGUI]={
ContourOpacity->1,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,color->"Rainbow",NumberofPoints->20000,size->600,FilterSize->1
};


codGUI[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,varNumb,sigmaNumb,dataIn,yMax,takeNumber,reducedData,meanSigma,incr,contourList,colorF,colorFS,contPlt,greyPlt,v1,v2,v3,vo
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn=Map[{#[[xNumb]],(yMax-#[[yNumb]]),#[[varNumb]],#[[sigmaNumb]]}&,dataIn,{2}];

meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
reducedData=Select[Flatten[dataIn,1],(#[[4]]!=-1.)&];

If[OptionValue[NumberofPoints]<Length[reducedData],
takeNumber=Round[(Length[reducedData])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[reducedData,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,Opacity[0.5,Black]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,None}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(First[meanSigma[{#1,#2}]]>0.&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]},
ColorFunction->colorF,
ColorFunctionScaling->colorFS];

Deploy[
Grid[{{
LocatorPane[
Dynamic[guiNodes],
Graphics[{
Raster[Reverse[greyPlt]],
Opacity[OptionValue[ContourOpacity]],
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.005],
Dynamic[Line[{guiNodes[[1;;2]],
Round[{guiNodes[[1]]+((guiNodes[[3]]-guiNodes[[1]])-((guiNodes[[2]]-guiNodes[[1]]).(guiNodes[[3]]-guiNodes[[1]]))/((guiNodes[[2]]-guiNodes[[1]]).(guiNodes[[2]]-guiNodes[[1]]))(guiNodes[[2]]-guiNodes[[1]])),guiNodes[[2]]+((guiNodes[[3]]-guiNodes[[1]])-((guiNodes[[2]]-guiNodes[[1]]).(guiNodes[[3]]-guiNodes[[1]]))/((guiNodes[[2]]-guiNodes[[1]]).(guiNodes[[2]]-guiNodes[[1]]))(guiNodes[[2]]-guiNodes[[1]]))}]}]],
Orange,
Dynamic[Line/@coordinates],
Dynamic[Point/@Flatten[coordinates,1]]
},
PlotRange->{{1,Length[greyPlt[[1]]]},{1,Length[greyPlt]}},
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12]
}]},
{Grid[{{Button["Clear coordinates",coordinates={};],
Button["Clear coordinatesOut",coordinatesOut={};]}}]},
{Button["Save coordinates",
yMax=Length[greyPlt]+1;
v1=guiNodes[[2]]-guiNodes[[1]];
v2=guiNodes[[3]]-guiNodes[[1]];
v3=(v1.v2)/(v1.v1)v1;
vo=v2-(v1.v2)/(v1.v1)v1;
AppendTo[coordinates,{guiNodes[[1;;2]],Round[{guiNodes[[1]]+vo,guiNodes[[2]]+vo}]}];
AppendTo[coordinatesOut,Map[{#[[1]],yMax-#[[2]]}&,{guiNodes[[1;;2]],Round[{guiNodes[[1]]+vo,guiNodes[[2]]+vo}]},{2}]];
]}}]]
];


(* ::Subsection::Closed:: *)
(*Extensometers*)


Options[DICextensometer]={
Displacements->False,PrintExtLength->False
};


DICextensometer[pathIn_,frames_,exts_,OptionsPattern[]]:=Module[{
filesIn,dataIn,extsIn,xNumb,yNumb,uNumb,vNumb,sigmaNumb,referenceData,refNearestX,refNearestY,points,refXYZ,Lo,ref,Lf,data,dataXYZ
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[Dimensions[extsIn]]==2,
extsIn={exts};,
extsIn=exts;];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
uNumb=Position[filesIn[[2]],"u"][[1,1]];
vNumb=Position[filesIn[[2]],"v"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

referenceData=MapThread[List,Import[dataIn[[1]]],2];
refNearestX=Nearest[referenceData[[1,All,xNumb]]->Automatic];
refNearestY=Nearest[referenceData[[All,1,yNumb]]->Automatic];

points=Map[SortBy[#,Last]&,Map[{First[refNearestX[#[[1]]]],First[refNearestY[#[[2]]]]}&,extsIn,{2}]];

refXYZ=Map[{#[[xNumb]]+#[[uNumb]],#[[yNumb]]+#[[vNumb]]}&,MapThread[{referenceData[[#1[[2]],#1[[1]]]],referenceData[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

Lo=Map[EuclideanDistance@@#&,refXYZ];

If[OptionValue[PrintExtLength],
Print[Lo];
];

ref=ProgressDialog[];

Lf=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing Frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

dataXYZ=Map[{#[[xNumb]]+#[[uNumb]],#[[yNumb]]+#[[vNumb]]}&,MapThread[{data[[#1[[2]],#1[[1]]]],data[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

Map[EuclideanDistance@@#&,dataXYZ]

,{i,1,Length[dataIn]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[OptionValue[Displacements],
Map[MapThread[#1-#2&,{#,Lo}]&,Lf]
,
Map[MapThread[(#1-#2)/#2&,{#,Lo}]&,Lf]
]

];


(* ::Subsection::Closed:: *)
(*COD Functions*)


Options[codAnalysis]={
method->2,ForceVert->False, CrackDirection->"x"
};


(*input1 = # of nodal points used to compute the slope, must be odd, number of pixels should be ~100, input2 = size of the excluded region in mm, typicall 2/3 of subest size*)
codAnalysis[pathIn_,frame_,{{x1_,y1_},{x2_,y2_}},input1_,input2_,OptionsPattern[]]:=Module[
{filesIn,dataIn,varListIn,data,xNumb,yNumb,YNumb,XNumb,dispNumb,xIn,yIn,yS,yE,xS,xE,xyGrid,calFactor,yMax,alldata,ans,codprofileall,dataOut},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frame]];

varListIn=filesIn[[2]];

data=MapThread[List,Import[dataIn],2];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"x_c"][[1,1]];
YNumb=Position[varListIn,"y_c"][[1,1]];

xIn=Sort[{x1,x2}];
yIn=Sort[{y1,y2}];

xS=First[Nearest[Flatten[data[[1,All,{xNumb}]]]->Automatic,xIn[[1]]]];
xE=First[Nearest[Flatten[data[[1,All,{xNumb}]]]->Automatic,xIn[[2]]]];
yS=First[Nearest[Flatten[data[[All,1,{yNumb}]]]->Automatic,yIn[[1]]]];
yE=First[Nearest[Flatten[data[[All,1,{yNumb}]]]->Automatic,yIn[[2]]]];

xyGrid=Take[data,{yS,yE},{xS,xE},All];

calFactor=EuclideanDistance[{xyGrid[[1,1,{XNumb,YNumb}]]},{xyGrid[[-1,-1,{XNumb,YNumb}]]}]/EuclideanDistance[{xyGrid[[1,1,{xNumb,yNumb}]]},{xyGrid[[-1,-1,{xNumb,yNumb}]]}];

If[OptionValue[ForceVert],yMax=2449,If[data[[-1,-1,yNumb]]>2048,yMax=2449,yMax=2049]];

If[TrueQ[OptionValue[CrackDirection]=="x"],

dispNumb=Position[varListIn,"v_c"][[1,1]];

alldata=Transpose[Map[{#[[1]]*calFactor,(yMax-#[[2]])*calFactor,#[[3]] }&,xyGrid[[All,All,{xNumb,yNumb,dispNumb}]],{2}],{2,1,3}];
ans = Switch[OptionValue[method],1,initarrays[alldata[[1,All,2]],input1],2,PseudoInverse[Table[{i,1},{i,1,input1}]]];
codprofileall=Switch[OptionValue[method],1,Map[{#[[1,1]],workhorsesub[#,input2,ans]}&,alldata],2,Map[{#[[1,1]],slopefind2[#,input1,input2,ans]}&,alldata]];

dataOut={Map[{#[[1]],#[[2,3]],#[[2,2,4]]}&,codprofileall],MapIndexed[{#1[[1]],#1[[2,1]],alldata[[#2,All,2;;3]][[1]]}&,codprofileall]};,

dispNumb=Position[varListIn,"u_c"][[1,1]];

alldata=Map[{(yMax-#[[1]])*calFactor,#[[2]]*calFactor,#[[3]] }&,xyGrid[[All,All,{yNumb,xNumb,dispNumb}]],{2}];
ans = Switch[OptionValue[method],1,initarrays[alldata[[1,All,2]],input1],2,PseudoInverse[Table[{i,1},{i,1,input1}]]];
codprofileall=Switch[OptionValue[method],1,Map[{#[[1,1]],workhorsesub[#,input2,ans]}&,alldata],2,Map[{#[[1,1]],slopefind2[#,input1,input2,ans]}&,alldata]];

dataOut={Map[{#[[2,3]],#[[1]],#[[2,2,4]]}&,codprofileall],MapIndexed[{#1[[1]],#1[[2,1]],alldata[[#2,All,2;;3]][[1]]}&,codprofileall]};
];

dataOut
];
(*Output:{{x,y,cod},{center,{COD slope fit},{original position-disp data}}}*)


initarrays[ydata_,crit1_]:=Module[
{ymin,ymax,ycenter,yhalfspan,critnew,myA2,myA3},

ymin=Min[ydata];
ymax=Max[ydata];
ycenter = (ymax+ymin)/2;
yhalfspan=(ymax-ymin)/2;
critnew=yhalfspan-crit1;
myA2=Map[{#,1}&,ydata];
myA3=Map[{#,1,1}&,ydata];

{myA2,myA3,ycenter,critnew}
];


findindices[ydata_,ycenter_,crit_]:=Module[
{},

Flatten[Position[ydata,_?(Abs[#-ycenter]<crit&)]]
];


slopefind[data_,ycenter_,indices_]:=Module[
{datanew,myA,coeff},

datanew= Join[data[[1;;(indices[[1]]-1)]],data[[(indices[[-1]]+1);;-1]]];
myA =Map[If[#[[2]]>ycenter,{0,#[[2]],1,1},{#[[2]],0,1,0}]&,datanew];
coeff = PseudoInverse[myA].datanew[[All,3]];

{Transpose[{datanew[[All,2]],myA.coeff}],coeff}
];


findcent[data_,ans_,coeff_,indices_]:=Module[
{extrap1,extrap2,ydata,vdata,dev1,dev2,devmin},

extrap1=(ans[[1]].coeff[[{1,3}]])[[indices]];
extrap2=(ans[[2]].coeff[[2;;4]])[[indices]];
ydata=(data[[All,2]])[[indices]];
vdata=(data[[All,3]])[[indices]];
dev1=Abs[vdata-extrap1];
dev2=Abs[vdata-extrap2];
devmin=MapThread[Min,{dev1,dev2}];

ydata[[Ordering[devmin,-1]]][[1]]
];


workhorsesub[data_,crit2_,ans_]:=Module[
{indexexcludecent,slopeans,ycenternew},

indexexcludecent=findindices[data[[All,2]],ans[[3]],ans[[4]]];
slopeans=slopefind[data,ans[[3]],indexexcludecent];
Do[
ycenternew=findcent[data,ans,slopeans[[2]],indexexcludecent];
indexexcludecent=findindices[data[[All,2]],ycenternew,crit2];
slopeans=slopefind[data,ycenternew,indexexcludecent];,
{n,3}];

{slopeans[[1]],(slopeans[[2,2]]-slopeans[[2,1]])*ycenternew+slopeans[[2,4]]}
];


slopefind2[data_,fitpoints_,crit_,mypseudo_]:=Module[
{mypart,ansy,ans,center,datanew,myA,coeff},

mypart = Transpose[Partition[data[[All,{2,3}]],fitpoints,1]];
ansy=mypart[[(fitpoints+1)/2,All,1]];
ans=(mypseudo.mypart[[All,All,2]])[[1]];
center=ansy[[Ordering[Abs[ans],-1]]][[1]];
datanew= Select[data,Abs[#[[2]]-center]>crit&];
myA =Map[If[#[[2]]>center,{0,#[[2]]-center,1,1},{#[[2]]-center,0,1,0}]&,datanew];
coeff = PseudoInverse[myA].datanew[[All,3]];

{Transpose[{datanew[[All,2]],myA.coeff}],coeff,center}
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["TwoDDICFunctions`*"],
{
"nodes","coordinates","coordinatesOut","guiNodes"
}];
EndPackage[ ]
