(* ::Package:: *)

BeginPackage["ThreeDDICFunctions`"]
Unprotect@@Names["ThreeDDICFunctions`*"];
ClearAll@@Names["ThreeDDICFunctions`*"];


(*Data Import*)
DICdataFiles::usage = "DICdataFiles[pathIn_] returns DIC files and variable list from v6 files folder";
allDataFiles::usage = "allDataFiles[pathIn_] returns DIC files and variable list from v6 files folder and get image files from parent directory";
importDICdata::usage = "importDICdata[pathIn_,OptionsPattern[]] extracts data for all files in pathIn.  Can specify which variables.";
replaceBadData::usage = "replaceBadData[pathIn_,frame_] replaces all bad data using interpolation for file specified.";
cropData::usage = "cropData[dataIn_,pointsIn_] crop DIC data to the smallest array defined by pointsIn.";

(*Load Displacement Import*)
DICOutput::usage = "DICOutput[pathIn_,fullScale_] returns image number and load from the DIC .csv file";
MTSOutput::usage = "MTSOutput[pathIn_] returns all data from MTS .DAT file with format {Disp,Load,Strain,External Inputs}";

(*DIC Analysis*)
getMagnification::usage = "getMagnification[filesPath_,frames_] returns {mm/pixels, pixels/mm} for frame or an average of all frames";

(*Error Analysis*)
dispError::usage = "dispError[pathIn_,frame_] returns {{u stDev,v st Dev},{U stDev, V stDev, W st Dev}} for input file";
strainError::usage = "strainError[pathIn_,frame_] returns {{exx stDev,eyy stDev, exy stDeV},{exx mean,eyy mean, exy mean}} for input file";

(*Line Scans*)
DIClineScan::usage = "DIClineScan[pathIn_,frames_,lines_,vars_,OptionsPattern[]] returns {x,y,vars} for all nodes on linesIn, with options to return metric or pixel {x,y} values ";
DIClineSweep::usage = "DIClineSweep[pathIn_,frames_,lines_,vars_,OptionsPattern[]] returns {x,y,vars} for all lines parallel to linesIn[[1]] sweeping out to the point linesIn[[2]], with options to return metric or pixel {x,y} values";
codTool::usage = "CODdisp[pathIn_,frames_,lines_,vars_,OptionsPattern[]] calculates the change in distance between every point on two parallel lines. ";
getDICgrid::usage = "getDICgrid[pathIn_,frames_,grids_,vars_] extracts rectangular grids of data bounded by gridsIn for frames with the given variables vars.";

(*Extensometers*)
DICextensometer::usage = "DICextensometers[pathIn_,extsIn_,OptionsPattern[]] returns strains or displacemens for the given extsIn, option to return a single frame";
extAngle::usage = "extAngle[pathIn_,frames_,exts_,OptionsPattern[]] returns the angle between extIn and the x-axis, option to return a single frame";
tangentModulus::usage = "tangentModulus[dataIn_,windowSize_] returns {strain, stress, t-modulus} for given windowSize";
aveStrains::usage = "aveStrain[pathIn_,frames_] calculates the average {exx,eyy,exy} for all good nodes in frames";
aveSlopes::usage = "aveSlopes[pathIn_,frames_] calculates the average {du/dx,dv/dy,(du/dy+dv/dx)/2} for a planar fit to all the u and v data";

(*Disp and Strains*)
dataGrid::usage = "dataGrid[dataIn_] uses interpolation to create an orthogonal grid of the input data. Options to locate holes/notches. Input is{X,Y,data,sigma}";
dispGrid::usage = "dispGrid[dataIn_] uses interpolation to create an orthogonal grid of U and V values centered in the grid. Options to locate holes/notches. Input is{X,Y,U,V,sigma}.";
forwardDiff::usage = "forwardDiff[gridIn_,vector_] calculates the forward difference for a grid of data in the veritcal{1,0} or horizontal{0,1} direction";
strainGrid::usage = "strainGrid[dispGridIn_,filterSize_] calculates {\[Epsilon]xx,\[Epsilon]yy,\[Gamma]xy} from the given dispGridIn and filters using a Gaussian Filter of filterSize, outputs";
calculateStrains::usage = "calculateStrains[pathIn_,frame_,filterSize_] calculates strains for the given frame using the give filter size. Bad nodes are excluded during filtering.";

(*Varun COD Analysis*)
codAnalysis::usage = "codAnalysis[pathIn_,frame_,{{x1_,y1_},{x2_,y2_}},centerCalSize_,disconSize_] gives {{x,y,cod},{x,{COD slope fit},{original y-v data}}},centerCalSize = number of nodes used to calculate center (~ 2 subsets), disconSize = data excluded from discontinuity (~2/3 subset/magnification), Output:{{x,y,cod},{center,{COD slope fit},{original position-disp data}}}";
initarrays::usage = "initarrays[ydata_,crit1_] for codAnalysis";
findindices::usage = "findindices[ydata_,ycenter_,crit_] for codAnalysis";
slopefind::usage = "slopefind[data_,ycenter_,indices_] finds cod for codAnalysis";
findcent::usage = "findcent[data_,ans_,coeff_,indices_] finds center of crack for codAnalysis";
workhorsesub::usage = "workhorsesub[data_,crit2_,ans_] for codAnalysis";
slopefind2::usage = "slopefind2[data_,fitpoints_,crit_,mypseudo_] finds cod for codAnalysis";
averageFilter::usage = "averageFilter[dataIn_,filterLength_] averages groups of data equal to filter length";

(*COD Slope Analysis*)
linearRegression::usage = "linearRegression[dataIn_] calculates the linear regression for dataIn, output is {b,m}";
slopeFit::usage = "slopeFit[dataIn_,w_] fits the strain peak associated with a crack to a parabola and outputs (center,COD)";
measureCOD::usage = "measureCOD[lineIn_,w_,OptionsPattern[]] analyzes any cracks on lineIn using the slope of the line and outputs ({x,y,cod},{line Slopes,line Displacements}";
codSlopeAnalysis::usage = "codSlopeAnalysis[pathIn_,frames_,grids_,hSub_,OptionsPattern[]] calculates the CODs and crack centers (x,y,COD) for any cracks on every line scan in every gridIn. Output is {{Crack Locations},{linePos,{slopes},{displacements}}";

(*Get Coordinates*)
nodes::usage = "nodes with {{x1,y1},{x2,y2}}";
coordinates::usage = "List of coordinates for {node1,node2} in real space({0,0} is bottom left corner)";
coordinatesOut::usage = "List of coordinates for {node1,node2} in the camera array ({0,0} is upper left corner)";
getCoordinates::usage = "getCoordinates[filesPath_,frame_,varIn_] allows for selection of nodal locations";
guiNodes::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3}}";
getNodesGUI::usage = "getNodesGUI[filesPath_,frame_,varIn_,OptionsPattern[]]allows for selection of line and end point for disp analysis";
codGUI::usage = "codGUI[filesPath_,frame_,varIn_] allows for selection of line for disp analysis";

(*Contour Overlays*)
getDeformedData::"getDeformedData[pathIn_,frame_,varIn_] extracts an array of {X,Y,varIn,sigma} in the deformed or reference frame with options for centering and rotating U&V or plotting deformed or reference Z";
DICoverLay::usage = "DICoverLay[filesPath_,frame_,varIn_] renders contours of varIn overlayed on the DIC image for the given frame";

(*Legends*)
DICLegend::usage = "DICLegend[{{sMin_,sMax_},contours_,labels_},labelIn_] renders a legend corresponding to the DIC contour overlay with identical inputs";


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


importDICdata[pathIn_,files_,variablesIn_]:=Module[{
filesIn,variables},

filesIn=DICdataFiles[pathIn];


If[SameQ[variablesIn,All],
variables=All;
,
If[Length[variablesIn]==0,
variables=Position[filesIn[[2]],variablesIn][[1,1]];
,
variables=Map[Position[filesIn[[2]],#][[1,1]]&,variablesIn];
];
];

If[NumberQ[files],
MapThread[List,Import[filesIn[[1,files]]],2][[All,All,variables]]
,
Map[MapThread[List,Import[#],2][[All,All,variables]]&,filesIn[[1,files]]]
]

];


replaceBadData[pathIn_,frame_]:=Module[
{dataIn,dataFlat,xRange,yRange,interpX,interpY,interpZ,badNodesPos,badNodesxy},

dataIn=importDICdata[pathIn,frame,{"x","y","X","Y","Z","sigma"}];
dataFlat=Select[Flatten[dataIn,1],(#[[6]]!=-1.)&][[All,1;;5]];

xRange={Min[dataFlat[[All,1]]],Max[dataFlat[[All,1]]]};
yRange={Min[dataFlat[[All,2]]],Max[dataFlat[[All,2]]]};

interpX=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,dataFlat],InterpolationOrder->1];
interpY=Interpolation[Map[{{#[[1]],#[[2]]},#[[4]]}&,dataFlat],InterpolationOrder->1];
interpZ=Interpolation[Map[{{#[[1]],#[[2]]},#[[5]]}&,dataFlat],InterpolationOrder->1];

xRange=Sort[Map[First[Nearest[dataIn[[1,All,1]]->Automatic,#]]&,xRange]];
yRange=Sort[Map[First[Nearest[dataIn[[All,1,2]]->Automatic,#]]&,yRange]];

dataIn=dataIn[[yRange[[1]];;yRange[[2]],xRange[[1]];;xRange[[2]]]];

badNodesPos=Position[dataIn[[All,All,6]],-1.];
badNodesxy=Map[Extract[dataIn,#][[1;;2]]&,badNodesPos];

dataIn=dataIn[[All,All,3;;5]];0

Do[
dataIn[[badNodesPos[[i,1]],badNodesPos[[i,2]]]]={interpX@@badNodesxy[[i]],interpY@@badNodesxy[[i]],interpZ@@badNodesxy[[i]]};
,{i,Length[badNodesPos]}];

dataIn

];


cropData[dataIn_,pointsIn_]:=Module[
{xRange,yRange,dim,xyRange,xPos,yPos},

xRange=Sort[pointsIn[[All,1]]];
yRange=Sort[pointsIn[[All,2]]];

dim=Round[Dimensions[dataIn][[1;;2]]/2];
xyRange=Table[
xPos=First[Nearest[dataIn[[dim[[1]],All,1]]->Automatic,xRange[[i]]]];
yPos=First[Nearest[dataIn[[All,dim[[2]],2]]->Automatic,yRange[[i]]]];
xPos=First[Nearest[dataIn[[yPos,All,1]]->Automatic,xRange[[i]]]];
yPos=First[Nearest[dataIn[[All,xPos,2]]->Automatic,yRange[[i]]]];
{xPos,yPos}, {i,Length[xRange]}];
xRange=Sort[xyRange[[All,1]]];
yRange=Sort[xyRange[[All,2]]];

dataIn[[yRange[[1]];;yRange[[2]],xRange[[1]];;xRange[[2]]]]

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
(*DIC Analysis*)


getMagnification[filesPath_,frames_]:=Module[{
filesIn,dataIn,varListIn,xNumb,yNumb,XNumb,YNumb,UNumb,VNumb,WNumb,sigmaNumb,ref,calFactor,data,goodData
},

filesIn=DICdataFiles[filesPath];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

varListIn=filesIn[[2]];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

ref=ProgressDialog[];

calFactor=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/(Length[filesIn[[1]]]))]];
ref @ SetPropertyValue[{"label", "text"}, "Analyzing File "<>ToString[(i)]<>" of "<>ToString[(Length[filesIn[[1]]])]];

data=MapThread[List,Import[dataIn[[i]]],2];

goodData=SortBy[Delete[Flatten[data,1],Position[Flatten[data[[All,All,sigmaNumb]],1],-1.]],(#[[xNumb]]+#[[yNumb]])&];

EuclideanDistance[goodData[[-1,{XNumb,YNumb}]],goodData[[1,{XNumb,YNumb}]]]/EuclideanDistance[goodData[[-1,{xNumb,yNumb}]],goodData[[1,{xNumb,yNumb}]]]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

Share[];

{1/Mean[calFactor],Mean[calFactor]}
];


(* ::Subsection::Closed:: *)
(*Error Analysis*)


dispError[filesPath_,frame_]:=Module[{
filesIn,dataIn,varListIn,uNumb,vNumb,UNumb,VNumb,WNumb,sigmaNumb,data,goodData,mag
},

filesIn=DICdataFiles[filesPath];
dataIn=filesIn[[1,frame]];

varListIn=filesIn[[2]];

uNumb=Position[varListIn,"u"][[1,1]];
vNumb=Position[varListIn,"v"][[1,1]];
UNumb=Position[varListIn,"U"][[1,1]];
VNumb=Position[varListIn,"V"][[1,1]];
WNumb=Position[varListIn,"W"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

data=Flatten[MapThread[List,Import[dataIn],2],1];

goodData=Delete[data,Position[data[[All,sigmaNumb]],-1.]];

(*Output is {{u stDev,v st Dev},{U stDev, V stDev, W st Dev}}*)
{{StandardDeviation[goodData[[All,uNumb]]],StandardDeviation[goodData[[All,vNumb]]]},{StandardDeviation[goodData[[All,UNumb]]],StandardDeviation[goodData[[All,VNumb]]],StandardDeviation[goodData[[All,WNumb]]]}}
];


strainError[filesPath_,frame_]:=Module[{
filesIn,dataIn,varListIn,eyyNumb,exxNumb,exyNumb,sigmaNumb,data,goodData
},

filesIn=DICdataFiles[filesPath];
dataIn=filesIn[[1,frame]];

varListIn=filesIn[[2]];

exxNumb=Position[varListIn,"exx"][[1,1]];
eyyNumb=Position[varListIn,"eyy"][[1,1]];
exyNumb=Position[varListIn,"exy"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

data=Flatten[MapThread[List,Import[dataIn],2][[All,All,{exxNumb,eyyNumb,exyNumb,sigmaNumb}]],1];

goodData=Delete[data,Position[data[[All,4]],-1.]];

(*Output is {{exx stDev,eyy stDev, exy stDeV},{exx mean,eyy mean, exy mean}}*)
{{StandardDeviation[goodData[[All,1]]],StandardDeviation[goodData[[All,2]]],StandardDeviation[goodData[[All,3]]]},{Mean[goodData[[All,1]]],Mean[goodData[[All,2]]],Mean[goodData[[All,3]]]}}
];


(* ::Subsection::Closed:: *)
(*Line Scans*)


Options[DIClineScan]={
CoordinateSystem->"Metric",ForceVertical->False
};


DIClineScan[pathIn_,frames_,lines_,vars_,OptionsPattern[]]:=Module[{
filesIn,dataIn,varsIn,linesIn,varsOut,xyVars,ref,refData,step,nearestX,nearestXPos,nearestY,nearestYPos,xIn,yIn,lineFit,y,x,m,b,coords,positions,data,dataOut,yMax
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
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"X","Y"}~Join~varsIn];
,
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

If[i==1,
data=refData;
,
data=MapThread[List,Import[dataIn[[i]]],2];
];

dataOut=Map[Extract[data[[All,All,varsOut]],#]&,positions];

If[MemberQ[vars,"v"],
dataOut[[All,(Position[vars,"v"][[1,1]]+2)]]=-1.*dataOut[[All,(Position[vars,"v"][[1,1]]+2)]];
,
dataOut=dataOut;
];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
dataOut
,
If[OptionValue[ForceVertical],
yMax=2449;
,
If[data[[-1,-1,xyVars[[2]]]]>2049,
yMax=2449;
,
yMax=2049;
];
];
dataOut[[All,2]]=yMax-dataOut[[All,2]];
dataOut
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


Options[DIClineSweep]={
CoordinateSystem->"Pixel",ForceVertical->False
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
varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"X","Y"}~Join~varsIn];,
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

If[i==1,
data=refData;
,
data=MapThread[List,Import[dataIn[[i]]],2];
];

dataOut=Map[Extract[data[[All,All,varsOut]],#]&,positions,{2}];

If[MemberQ[vars,"v"],
dataOut[[All,(Position[vars,"v"][[1,1]]+2)]]=-1.*dataOut[[All,(Position[vars,"v"][[1,1]]+2)]];
,
dataOut=dataOut;
];

If[StringMatchQ[OptionValue[CoordinateSystem],"m*",IgnoreCase->True],
dataOut
,
If[OptionValue[ForceVertical],
yMax=2449;
,
If[data[[-1,-1,xyVars[[2]]]]>2049,
yMax=2449;
,
yMax=2049;
];
];
dataOut[[All,2]]=yMax-dataOut[[All,2]];
dataOut
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


codTool[pathIn_,frames_,lines_,OptionsPattern[]]:=Module[{
filesIn,dataIn,linesIn,varsOut,xyVars,ref,refData,step,nearestX,nearestXPos,nearestY,nearestYPos,xIn,yIn,lineFit,y,x,coords,positions,lo,data,dataOut,yMax
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

varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,{"X","Y","Z","U","V","W"}];

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

refData=Map[Extract[refData[[All,All,varsOut]],#]&,positions];
refData=Map[{#[[1]]+#[[4]],#[[2]]+#[[5]],#[[3]]+#[[6]]}&,refData,{2}];
lo=MapThread[{EuclideanDistance[refData[[1,-1]],#1],EuclideanDistance[#1,#2]}&,refData];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

If[i==1,
data=refData;
,
data=MapThread[List,Import[dataIn[[i]]],2];
];

data=Map[Extract[data[[All,All,varsOut]],#]&,positions];
data=Map[{#[[1]]+#[[4]],#[[2]]+#[[5]],#[[3]]+#[[6]]}&,data,{2}];
MapThread[{#3[[1]],(EuclideanDistance[#1,#2]-#3[[2]])}&,data~Join~{lo}]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

dataOut

];


getDICgrid[pathIn_,frames_,grids_,vars_]:=Module[{
filesIn,dataIn,varsIn,gridsIn,varsOut,xyVars,ref,refData,xNearest,yNearest,gridPos,x1,x2,y1,y2,xIn,yIn,dataOut,data,gridsPos
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[vars]==0,
varsIn={vars};,
varsIn=vars;];

If[Length[Dimensions[grids]]==2,
gridsIn={grids};,
gridsIn=grids;];

varsOut=Map[Position[filesIn[[2]],#][[1,1]]&,varsIn];
xyVars=Map[Position[filesIn[[2]],#][[1,1]]&,{"x","y"}];

ref=ProgressDialog[];

ref @ SetPropertyValue[{"label", "text"}, "Calculating grid positions"];

refData=MapThread[List,Import[dataIn[[1]]],2];
xNearest=Nearest[refData[[1,All,xyVars[[1]]]]->Automatic];
yNearest=Nearest[refData[[All,1,xyVars[[2]]]]->Automatic];

gridsPos=Table[
{{x1,y1},{x2,y2}}=gridsIn[[i]];

xIn=Sort[{x1,x2}];
yIn=Sort[{y1,y2}];

{First[xNearest[xIn[[1]]]];;First[xNearest[xIn[[2]]]],First[yNearest[yIn[[1]]]];;First[yNearest[yIn[[2]]]]}
,{i,Length[gridsIn]}];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

If[i==1,
data=refData;
,
data=MapThread[List,Import[dataIn[[i]]],2];
];

Map[data[[#[[2]],#[[1]],varsOut]]&,gridsPos]

,{i,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[NumberQ[frames],
If[First[Dimensions[gridsIn]]==1,
dataOut[[1,1]]
,
dataOut[[1]]
]
,
If[First[Dimensions[gridsIn]]==1,
dataOut[[All,1]]
,
dataOut
]
]

];


(* ::Subsection::Closed:: *)
(*Extensometers*)


Options[DICextensometer]={
Displacements->False,PrintExtLength->False
};


DICextensometer[pathIn_,frames_,exts_,OptionsPattern[]]:=Module[{
filesIn,dataIn,extsIn,xNumb,yNumb,XNumb,YNumb,ZNumb,UNumb,VNumb,WNumb,sigmaNumb,referenceData,refNearestX,refNearestY,points,refXYZ,Lo,ref,Lf,data,dataXYZ
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[Dimensions[extsIn]]==2,
extsIn={exts};,
extsIn=exts;];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
ZNumb=Position[filesIn[[2]],"Z"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
WNumb=Position[filesIn[[2]],"W"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

referenceData=MapThread[List,Import[dataIn[[1]]],2];
refNearestX=Nearest[referenceData[[1,All,xNumb]]->Automatic];
refNearestY=Nearest[referenceData[[All,1,yNumb]]->Automatic];

points=Map[SortBy[#,Last]&,Map[{First[refNearestX[#[[1]]]],First[refNearestY[#[[2]]]]}&,extsIn,{2}]];

refXYZ=Map[{#[[XNumb]]+#[[UNumb]],#[[YNumb]]+#[[VNumb]],#[[ZNumb]]+#[[WNumb]]}&,MapThread[{referenceData[[#1[[2]],#1[[1]]]],referenceData[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

Lo=Map[EuclideanDistance@@#&,refXYZ];

If[OptionValue[PrintExtLength],
Print[Lo];
];

ref=ProgressDialog[];

Lf=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing Frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

dataXYZ=Map[{#[[XNumb]]+#[[UNumb]],#[[YNumb]]+#[[VNumb]],#[[ZNumb]]+#[[WNumb]]}&,MapThread[{data[[#1[[2]],#1[[1]]]],data[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

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


Options[extAngle]={
getRotation->False
};


extAngle[pathIn_,frames_,exts_,OptionsPattern[]]:=Module[{
filesIn,dataIn,extsIn,xNumb,yNumb,XNumb,YNumb,ZNumb,UNumb,VNumb,WNumb,sigmaNumb,referenceData,refNearestX,refNearestY,points,refXYZ,phi,ref,phiOut,data,dataXYZ
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[Dimensions[extsIn]]==2,
extsIn={exts};,
extsIn=exts;];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
ZNumb=Position[filesIn[[2]],"Z"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
WNumb=Position[filesIn[[2]],"W"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

referenceData=MapThread[List,Import[dataIn[[1]]],2];
refNearestX=Nearest[referenceData[[1,All,xNumb]]->Automatic];
refNearestY=Nearest[referenceData[[All,1,yNumb]]->Automatic];

points=Map[SortBy[#,Last]&,Map[{First[refNearestX[#[[1]]]],First[refNearestY[#[[2]]]]}&,extsIn,{2}]];

refXYZ=Map[{#[[XNumb]]+#[[UNumb]],#[[YNumb]]+#[[VNumb]],#[[ZNumb]]+#[[WNumb]]}&,MapThread[{referenceData[[#1[[2]],#1[[1]]]],referenceData[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

phi=Map[N[ArcTan[(#[[1,2]]-#[[2,2]])/(#[[1,1]]-#[[2,1]])]/Degree]&,refXYZ];

ref=ProgressDialog[];

phiOut=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing Frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=MapThread[List,Import[dataIn[[i]]],2];

dataXYZ=Map[{#[[XNumb]]+#[[UNumb]],#[[YNumb]]+#[[VNumb]],#[[ZNumb]]+#[[WNumb]]}&,MapThread[{data[[#1[[2]],#1[[1]]]],data[[#2[[2]],#2[[1]]]]}&,{points[[All,1]],points[[All,2]]}],{2}];

Map[N[ArcTan[(#[[1,2]]-#[[2,2]])/(#[[1,1]]-#[[2,1]])]/Degree]&,dataXYZ]

,{i,2,Length[dataIn]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[OptionValue[getRotation],
Map[MapThread[#2-#1&,{#,phi}]&,{phi}~Join~phiOut]
,
{phi}~Join~phiOut
]
];


tangentModulus[dataIn_,windowSize_]:=Module[{rangeS,rangeE,shift,m,x,b},
rangeS[s_]:=Switch[s<1,True,1,False,s];
rangeE[e_]:=Switch[e>Length[dataIn],True,Length[dataIn],False,e];
shift=Round[windowSize/2];

Table[dataIn[[i]]~Join~{m}/.FindFit[dataIn[[rangeS[i-shift];;rangeE[i+shift]]],m*x+b,{m,b},x]
,{i,1,Length[dataIn]}]
];


aveStrains[pathIn_,frames_]:=Module[{
filesIn,dataIn,exxNumb,eyyNumb,exyNumb,sigmaNumb,ref,aveStrain,data
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

exxNumb=Position[filesIn[[2]],"exx"][[1,1]];
eyyNumb=Position[filesIn[[2]],"eyy"][[1,1]];
exyNumb=Position[filesIn[[2]],"exy"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

ref=ProgressDialog[];

aveStrain=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing Frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=Flatten[MapThread[List,Import[dataIn[[i]]],2],1][[All,{exxNumb,eyyNumb,exyNumb,sigmaNumb}]];
data=Select[data,(#[[4]]!= -1.)&];

Mean[data[[All,1;;3]]]

,{i,1,Length[dataIn]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

aveStrain

];


aveSlopes[pathIn_,frames_]:=Module[{
filesIn,dataIn,XNumb,YNumb,UNumb,VNumb,sigmaNumb,ref,aveStrain,data,du,ax,mx,x,bx,y,ay,by,my,dv
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

ref=ProgressDialog[];

aveStrain=Table[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing Frame "<>ToString[i]<>" of "<>ToString[Length[dataIn]]];

data=Flatten[MapThread[List,Import[dataIn[[i]]],2],1][[All,{XNumb,YNumb,UNumb,VNumb,sigmaNumb}]];
data=Select[data,(#[[5]]!=-1.)&];

du={mx,bx}/.FindFit[data[[All,{1,2,3}]],ax +mx x+bx y,{ax,mx,bx},{x,y}];
dv={by,my}/.FindFit[data[[All,{1,2,4}]],ay +by x+my y,{ay,by,my},{x,y}];

{du[[1]],dv[[2]],(du[[2]]+dv[[1]])/2.}

,{i,1,Length[dataIn]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

aveStrain

];


(* ::Subsection::Closed:: *)
(*Grid, Disp, and Strains*)


Options[dataGrid]={
BadRegions->0,StepSize->Automatic,XYRange->Automatic
};


dataGrid[dataIn_,OptionsPattern[]]:=Module[
{dim,xStep,yStep,step,data,cleanData,dataInterp,sigmaInterp,xyRange,gridData,sigmaClusters,holePos},

If[SameQ[OptionValue[StepSize],Automatic],
dim=Round[Dimensions[dataIn][[1;;2]]/2];
xStep=Abs[Mean[Differences[Select[dataIn[[dim[[1]],All,{1,4}]],(#[[2]]!=-1.)&][[All,1]]]]];
yStep=Abs[Mean[Differences[Select[dataIn[[All,dim[[2]],{2,4}]],(#[[2]]!=-1.)&][[All,1]]]]];
step=Mean[{xStep,yStep}];
,
step=OptionValue[StepSize];
];

If[SameQ[OptionValue[BadRegions],0],
data=dataIn;
cleanData=Select[Flatten[data,1],(#[[-1]]!=-1.)&];

dataInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,cleanData],InterpolationOrder->1];

If[SameQ[OptionValue[XYRange],Automatic],
xyRange=SortBy[cleanData[[All,1;;2]],Total][[{1,-1}]];
,
xyRange=OptionValue[XYRange];
];

gridData=Table[Quiet[{x,y,dataInterp[x,y]}]
,{y,xyRange[[1,2]],xyRange[[2,2]],step},{x,xyRange[[1,1]],xyRange[[2,1]],step}];
gridData
,
data=dataIn[[All,All,{1,2,3,4,4}]];
data[[All,All,5]]=MeanFilter[data[[All,All,5]],1];
cleanData=Select[Flatten[data,1],(#[[4]]!=-1.)&][[All,{1,2,3,5}]];

dataInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,cleanData],InterpolationOrder->1];
sigmaInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[4]]}&,cleanData],InterpolationOrder->1];

If[SameQ[OptionValue[XYRange],Automatic],
xyRange=SortBy[cleanData[[All,1;;2]],Total][[{1,-1}]];
,
xyRange=OptionValue[XYRange];
];
gridData=Table[Quiet[{x,y,dataInterp[x,y],sigmaInterp[x,y],0.1}]
,{y,xyRange[[1,2]],xyRange[[2,2]],step},{x,xyRange[[1,1]],xyRange[[2,1]],step}];

sigmaClusters=MorphologicalComponents[Map[Boole[#[[4]]<0]&,gridData,{2}]];
holePos=Flatten[Map[Position[sigmaClusters,#]&,SortBy[ComponentMeasurements[sigmaClusters,"Count"],Last][[-IntegerPart[OptionValue[BadRegions]];;,1]]],1];

gridData[[All,All,5]]=ReplacePart[gridData[[All,All,5]],holePos->-1.];
gridData[[All,All,{1,2,3,5}]]
]
];


Options[dispGrid]={
BadRegions->0,StepSize->Automatic,XYRange->Automatic,DeformedXY->False
};


dispGrid[dataIn_,OptionsPattern[]]:=Module[
{dim,xStep,yStep,step,data,cleanData,meanV,meanU,deformedData,Uinterp,Vinterp,sigmaInterp,xyRange,gridData,sigmaClusters,holePos},

If[SameQ[OptionValue[StepSize],Automatic],
dim=Round[Dimensions[dataIn][[1;;2]]/2];
xStep=Abs[Mean[Differences[Select[dataIn[[dim[[1]],All,{1,5}]],(#[[2]]!=-1.)&][[All,1]]]]];
yStep=Abs[Mean[Differences[Select[dataIn[[All,dim[[2]],{2,5}]],(#[[2]]!=-1.)&][[All,1]]]]];
step=Mean[{xStep,yStep}];
,
step=OptionValue[StepSize];
];

If[SameQ[OptionValue[BadRegions],0],
data=dataIn;
cleanData=Select[Flatten[data,1],(#[[-1]]!=-1.)&];
meanU=Mean[cleanData[[All,3]]];
meanV=Mean[cleanData[[All,4]]];

If[OptionValue[DeformedXY],
deformedData=Map[{#[[1]]+#[[3]]-meanU,#[[2]]+#[[4]]-meanV,#[[3]]-meanU,#[[4]]-meanV}&,cleanData];
,
deformedData=Map[{#[[1]],#[[2]],#[[3]]-meanU,#[[4]]-meanV}&,cleanData];
];

Uinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,deformedData],InterpolationOrder->1];
Vinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[4]]}&,deformedData],InterpolationOrder->1];

If[SameQ[OptionValue[XYRange],Automatic],
xyRange=SortBy[deformedData[[All,1;;2]],Total][[{1,-1}]];
,
xyRange=OptionValue[XYRange];
];
gridData=Table[Quiet[{x,y,Uinterp[x,y],Vinterp[x,y]}]
,{y,xyRange[[1,2]],xyRange[[2,2]],step},{x,xyRange[[1,1]],xyRange[[2,1]],step}];
gridData
,
data=dataIn[[All,All,{1,2,3,4,5,5}]];
data[[All,All,6]]=MeanFilter[data[[All,All,6]],1];
cleanData=Select[Flatten[data,1],(#[[5]]!=-1.)&];
meanU=Mean[cleanData[[All,3]]];
meanV=Mean[cleanData[[All,4]]];

If[OptionValue[DeformedXY],
deformedData=Map[{#[[1]]+#[[3]]-meanU,#[[2]]+#[[4]]-meanV,#[[3]]-meanU,#[[4]]-meanV,#[[6]]}&,cleanData];
,
deformedData=Map[{#[[1]],#[[2]],#[[3]]-meanU,#[[4]]-meanV,#[[6]]}&,cleanData];
];

Uinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,deformedData],InterpolationOrder->1];
Vinterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[4]]}&,deformedData],InterpolationOrder->1];
sigmaInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[5]]}&,deformedData],InterpolationOrder->1];

If[SameQ[OptionValue[XYRange],Automatic],
xyRange=SortBy[deformedData[[All,1;;2]],Total][[{1,-1}]];
,
xyRange=OptionValue[XYRange];
];
gridData=Table[Quiet[{x,y,Uinterp[x,y],Vinterp[x,y],sigmaInterp[x,y],0.1}]
,{y,xyRange[[1,2]],xyRange[[2,2]],step},{x,xyRange[[1,1]],xyRange[[2,1]],step}];

sigmaClusters=MorphologicalComponents[Map[Boole[#[[5]]<0]&,gridData,{2}]];
holePos=Flatten[Map[Position[sigmaClusters,#]&,SortBy[ComponentMeasurements[sigmaClusters,"Count"],Last][[-IntegerPart[OptionValue[BadRegions]];;,1]]],1];

gridData[[All,All,6]]=ReplacePart[gridData[[All,All,6]],holePos->-1.];
gridData[[All,All,{1,2,3,4,6}]]
]
];


forwardDiff[gridIn_,vector_]:=Module[
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


strainGrid[gridData_,filterSize_]:=Module[
{dX,dY,dUdX,dVdY,dUdY,dVdX,sigmaClusters,holePos,strainData,filterR,clipRange,gf,xRange,yRange,xClip,yClip,gfX,gfY,window,weights,out},

dX=forwardDiff[gridData[[All,All,1]],{0,1}];
dY=forwardDiff[gridData[[All,All,2]],{1,0}];
dUdX=forwardDiff[gridData[[All,All,3]],{0,1}]/dX;
dVdY=forwardDiff[gridData[[All,All,4]],{1,0}]/dY;
dUdY=forwardDiff[gridData[[All,All,3]],{1,0}]/dY;
dVdX=forwardDiff[gridData[[All,All,4]],{0,1}]/dX;

filterR=Floor[filterSize/2];

If[Dimensions[gridData][[3]]==4,
strainData=ConstantArray[{0.,0.,0.,0.,0.},Dimensions[gridData][[1;;2]]];
strainData[[All,All,1;;2]]=gridData[[All,All,1;;2]];
strainData[[All,All,3]]=dUdX;
strainData[[All,All,4]]=dVdY;
strainData[[All,All,5]]=dVdX+dUdY;

If[filterR==0,
out=strainData;
,
clipRange=Transpose[{{1,1},Dimensions[strainData][[1;;2]]}];
gf=GaussianMatrix[{filterR,filterSize/4.}];

DistributeDefinitions[strainData,filterR,clipRange,gf];
out=ParallelTable[

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

,{yi,clipRange[[1,1]],clipRange[[1,2]]},{xi,clipRange[[2,1]],clipRange[[2,2]]}];
];
,
strainData=ConstantArray[{0.,0.,0.,0.,0.,0.},Dimensions[gridData][[1;;2]]];
strainData[[All,All,1;;2]]=gridData[[All,All,1;;2]];
strainData[[All,All,3]]=ReplacePart[dUdX,holePos->"NaN"];
strainData[[All,All,4]]=ReplacePart[dVdY,holePos->"NaN"];
strainData[[All,All,5]]=ReplacePart[dVdX+dUdY,holePos->"NaN"];
strainData[[All,All,6]]=gridData[[All,All,5]];

If[filterR==0,
out=strainData;
,
clipRange=Transpose[{{1,1},Dimensions[strainData][[1;;2]]}];
gf=GaussianMatrix[{filterR,filterSize/4.}];

DistributeDefinitions[strainData,filterR,clipRange,gf];
out=ParallelTable[
If[strainData[[yi,xi,6]]!=-1.,
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

strainData[[yi,xi,1;;2]]~Join~Map[#.weights/Total[weights]&,Transpose[window[[All,1;;3]]]]~Join~strainData[[yi,xi,{6}]]
,
strainData[[yi,xi]]
]
,{yi,clipRange[[1,1]],clipRange[[1,2]]},{xi,clipRange[[2,1]],clipRange[[2,2]]}];
];
];

out
];


calculateStrains[pathIn_,frame_,filterSize_]:=Module[
{filesIn,XNumb,YNumb,UNumb,VNumb,sigmaNumb,dataIn,sigma,sigmaPos,dX,dY,dUdX,dVdY,dUdY,dVdX,strainData,filterR,clipRange,gf,xRange,yRange,xClip,yClip,gfX,gfY,window,weights,out},

filesIn=DICdataFiles[pathIn];

XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
dataIn=dataIn[[All,All,{XNumb,YNumb,UNumb,VNumb,sigmaNumb}]];
dataIn[[All,All,1;;2]]=dataIn[[All,All,1;;2]]+dataIn[[All,All,3;;4]];
sigma=Map[If[#<0,-1.,#]&,MeanFilter[dataIn[[All,All,5]],1],{2}];

sigmaPos=Position[sigma,-1.];
dX=ReplacePart[forwardDiff[dataIn[[All,All,1]],{0,1}],sigmaPos->Indeterminate];
dY=ReplacePart[forwardDiff[dataIn[[All,All,2]],{1,0}],sigmaPos->Indeterminate];
dUdX=forwardDiff[dataIn[[All,All,3]],{0,1}]/dX;
dVdY=forwardDiff[dataIn[[All,All,4]],{1,0}]/dY;
dUdY=forwardDiff[dataIn[[All,All,3]],{1,0}]/dY;
dVdX=forwardDiff[dataIn[[All,All,4]],{0,1}]/dX;

strainData=ConstantArray[{0.,0.,0.,0.,0.,0.},Dimensions[dataIn][[1;;2]]];
strainData[[All,All,1;;2]]=dataIn[[All,All,1;;2]];
strainData[[All,All,3]]=ReplacePart[dUdX,sigmaPos->"NaN"];
strainData[[All,All,4]]=ReplacePart[dVdY,sigmaPos->"NaN"];
strainData[[All,All,5]]=ReplacePart[dVdX+dUdY,sigmaPos->"NaN"];
strainData[[All,All,6]]=sigma;

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

window=Flatten[strainData[[yClip[[1]];;yClip[[2]],xClip[[1]];;xClip[[2]],3;;6]],1];
weights=Flatten[gf[[gfY[[1]];;gfY[[2]],gfX[[1]];;gfX[[2]]]]];
If[MemberQ[window[[All,-1]],-1.],
weights[[Flatten[Position[window[[All,-1]],-1.]]]]=0;
];

strainData[[yi,xi,1;;2]]~Join~Map[#.weights/Total[weights]&,Transpose[window[[All,1;;3]]]]~Join~strainData[[yi,xi,{6}]]
,
strainData[[yi,xi]]
]
,{yi,clipRange[[1,1]],clipRange[[1,2]]},{xi,clipRange[[2,1]],clipRange[[2,2]]}];
];

out
];


(* ::Subsection::Closed:: *)
(*Varun COD Analysis*)


Options[codAnalysis]={
method->2,ForceVert->False,CrackDirection->"x"
};


(*centerCalSize = number of nodes used to calculate center (1.5-2 subsets), disconSize = data excluded from discontinuity (~2/3 subset/magnification)*)
codAnalysis[pathIn_,frame_,{{x1_,y1_},{x2_,y2_}},centerCalSize_,disconSize_,OptionsPattern[]]:=Module[
{filesIn,dataIn,varListIn,data,xNumb,yNumb,YNumb,XNumb,dispNumb,xIn,yIn,xNearest,yNearest,yS,yE,xS,xE,xyGrid,calFactor,yMax,alldata,ans,codprofileall,dataOut},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frame]];

varListIn=filesIn[[2]];

data=MapThread[List,Import[dataIn],2];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];

xIn=Sort[{x1,x2}];
yIn=Sort[{y1,y2}];

xNearest=Nearest[data[[1,All,xNumb]]->Automatic];
yNearest=Nearest[data[[All,1,yNumb]]->Automatic];

xS=First[xNearest[xIn[[1]]]];
xE=First[xNearest[xIn[[2]]]];
yS=First[yNearest[yIn[[1]]]];
yE=First[yNearest[yIn[[2]]]];

xyGrid=data[[yS;;yE,xS;;xE]];

calFactor=EuclideanDistance[{xyGrid[[1,1,{XNumb,YNumb}]]},{xyGrid[[-1,-1,{XNumb,YNumb}]]}]/EuclideanDistance[{xyGrid[[1,1,{xNumb,yNumb}]]},{xyGrid[[-1,-1,{xNumb,yNumb}]]}];

If[OptionValue[ForceVert],
yMax=2449
,
If[data[[-1,-1,yNumb]]>2049,
yMax=2449
,
yMax=2049]];

If[StringMatchQ[OptionValue[CrackDirection],"x",IgnoreCase->True],
dispNumb=Position[varListIn,"V"][[1,1]];

alldata=Transpose[Map[{#[[1]]*calFactor,(yMax-#[[2]])*calFactor,#[[3]] }&,xyGrid[[All,All,{xNumb,yNumb,dispNumb}]],{2}],{2,1,3}];
ans = Switch[OptionValue[method],1,initarrays[alldata[[1,All,2]],centerCalSize],2,PseudoInverse[Table[{i,1},{i,1,centerCalSize}]]];
codprofileall=Switch[OptionValue[method],1,Map[{#[[1,1]],workhorsesub[#,disconSize,ans]}&,alldata],2,Map[{#[[1,1]],slopefind2[#,centerCalSize,disconSize,ans]}&,alldata]];

dataOut={Map[{#[[1]],#[[2,3]],#[[2,2,4]]}&,codprofileall],MapIndexed[{#1[[1]],#1[[2,1]],alldata[[#2,All,2;;3]][[1]]}&,codprofileall]};
,
dispNumb=Position[varListIn,"U"][[1,1]];

alldata=Map[{(yMax-#[[1]])*calFactor,#[[2]]*calFactor,#[[3]] }&,xyGrid[[All,All,{yNumb,xNumb,dispNumb}]],{2}];
ans = Switch[OptionValue[method],1,initarrays[alldata[[1,All,2]],centerCalSize],2,PseudoInverse[Table[{i,1},{i,1,centerCalSize}]]];
codprofileall=Switch[OptionValue[method],1,Map[{#[[1,1]],workhorsesub[#,disconSize,ans]}&,alldata],2,Map[{#[[1,1]],slopefind2[#,centerCalSize,disconSize,ans]}&,alldata]];

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


averageFilter[dataIn_,filterLength_]:=Module[{
filterStep,dataOut={}
},

Do[
filterStep=Floor[Length[dataIn[[i]]]/(Abs[dataIn[[i,-1,1]]-dataIn[[i,1,1]]]/(filterLength))];
AppendTo[dataOut,
Table[{Mean[dataIn[[i,j;;j+filterStep-1,1]]],Mean[dataIn[[i,j;;j+filterStep-1,2]]]},
{j,1,Floor[Length[dataIn[[i]]]/filterStep]*filterStep,filterStep}]],
{i,1,Length[dataIn]}];

dataOut
];


(* ::Subsection::Closed:: *)
(*COD Slope Analysis*)


linearRegression[dataIn_]:=Module[
{Y,X,B},

Y=dataIn[[All,2]];
X=Transpose[{ConstantArray[1,Length[dataIn]],dataIn[[All,1]]}];
B=Inverse[Transpose[X].X].Transpose[X].Y

];


slopeFit[dataIn_,w_]:=Module[
{fit,max,x},
fit=Fit[dataIn,{1,x,x^2},x];
max=FindMaximum[{fit,dataIn[[1,1]]<=x<=dataIn[[-1,1]]},{x,Mean[dataIn[[All,1]]]}];

{x,max[[1]]*w}/.max[[2]]
];


Options[measureCOD]={
CrackDirection->"x",CrackThreshold->1.75,GlobalStrain->Automatic
};


measureCOD[lineIn_,w_,OptionsPattern[]]:=Module[
{lineData,linePos,lineInterp,globalSlope,slope,slopeNearest,globalStrain,slopeThreshold,slopePos,crackSlopes,COD},

If[StringMatchQ[OptionValue[CrackDirection],"x*",IgnoreCase->True],
lineData=SortBy[Select[lineIn,(#[[-1]]!=-1)&][[All,1;;4]],(#[[2]])&];

linePos=Mean[lineData[[All,1]]];
lineData=lineData[[All,{2,4}]];
lineInterp[y_]:=Interpolation[lineData][y];
slope=Map[{#[[1]],lineInterp'[#[[1]]]}&,lineData];
slopeNearest=Nearest[slope[[All,2]]->slope];

If[NumberQ[OptionValue[GlobalStrain]],
globalStrain=OptionValue[GlobalStrain];
,
globalStrain=linearRegression[lineData][[2]];
];

slopeThreshold=OptionValue[CrackThreshold]*globalStrain;
slopePos=Flatten[MapIndexed[#1*First[#2]&,Split[ Map[Boole[#>slopeThreshold]&,slope[[All,2]]]],2]];
slopePos=Select[Map[Position[slopePos,#]&,Union[slopePos]],(Length[#]>2)&];

If[Length[slopePos]>1,
crackSlopes=Map[Extract[slope,#]&,slopePos[[2;;]]];
COD=Reverse[SortBy[Map[{linePos}~Join~slopeFit[#,w]&,crackSlopes],Last]];
,
COD={{linePos,Null,Null}};
];
,
lineData=SortBy[Select[lineIn,(#[[-1]]!=-1)&][[All,1;;4]],(#[[1]])&];

linePos=Mean[lineData[[All,2]]];
lineData=lineData[[All,{1,3}]];
lineInterp[x_]:=Interpolation[lineData][x];
slope=Map[{#[[1]],lineInterp'[#[[1]]]}&,lineData];
slopeNearest=Nearest[slope[[All,2]]->slope];

If[NumberQ[OptionValue[GlobalStrain]],
globalStrain=OptionValue[GlobalStrain];
,
globalStrain=linearRegression[lineData][[2]];
];

slopeThreshold=OptionValue[CrackThreshold]*globalStrain;
slopePos=Flatten[MapIndexed[#1*First[#2]&,Split[ Map[Boole[#>slopeThreshold]&,slope[[All,2]]]],2]];
slopePos=Select[Map[Position[slopePos,#]&,Union[slopePos]],(Length[#]>2)&];

If[Length[slopePos]>1,
crackSlopes=Map[Extract[slope,#]&,slopePos[[2;;]]];
COD=Reverse[SortBy[Map[{linePos}~Join~slopeFit[#,w]&,crackSlopes][[All,{2,1,3}]],Last]];
,
crackSlopes={Null};
COD={{Null,linePos,Null}};
];
];

{COD,{linePos,slope,crackSlopes,lineData}}
];


Options[codSlopeAnalysis]={
CrackDirection->"x",CrackThreshold->1.5,ReferenceStrain->"Global"
};


codSlopeAnalysis[pathIn_,frames_,grids_,hSub_,OptionsPattern[]]:=Module[{
filesIn,dataIn,gridsIn,xNumb,yNumb,XNumb,YNumb,UNumb,VNumb,sigmaNumb,dataOut,ref,refData,xNearest,yNearest,gridsPos,x1,y1,x2,y2,xIn,yIn,goodData,w,data,globalStrain,x,y,a,b,m,gridData,codProfiles
},

filesIn=DICdataFiles[pathIn];
dataIn=filesIn[[1,frames]];

If[Length[dataIn]==0,
dataIn={dataIn};,
dataIn=dataIn;];

If[Length[Dimensions[grids]]==2,
gridsIn={grids};,
gridsIn=grids;];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

ref=ProgressDialog[];

ref @ SetPropertyValue[{"label", "text"}, "Calculating grid positions"];

refData=MapThread[List,Import[dataIn[[1]]],2];
xNearest=Nearest[refData[[1,All,xNumb]]->Automatic];
yNearest=Nearest[refData[[All,1,yNumb]]->Automatic];

gridsPos=Table[
{{x1,y1},{x2,y2}}=gridsIn[[i]];

xIn=Sort[{x1,x2}];
yIn=Sort[{y1,y2}];

{First[xNearest[xIn[[1]]]];;First[xNearest[xIn[[2]]]],First[yNearest[yIn[[1]]]];;First[yNearest[yIn[[2]]]]}
,{i,Length[gridsIn]}];

goodData=SortBy[Select[Flatten[refData,1],(#[[sigmaNumb]]!=-1.)&],(#[[xNumb]]+#[[yNumb]])&];
w=0.682*hSub*EuclideanDistance[goodData[[-1,{XNumb,YNumb}]],goodData[[1,{XNumb,YNumb}]]]/EuclideanDistance[goodData[[-1,{xNumb,yNumb}]],goodData[[1,{xNumb,yNumb}]]];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(frame/Length[dataIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing frame "<>ToString[frame]<>" of "<>ToString[Length[dataIn]]];

If[frame==1,
data=refData;
goodData=goodData;
,
data=MapThread[List,Import[dataIn[[frame]]],2];
goodData=Select[Flatten[data,1],(#[[sigmaNumb]]!=-1.)&];
];

Table[
gridData=data[[gridsPos[[grid,2]],gridsPos[[grid,1]],{XNumb,YNumb,UNumb,VNumb,sigmaNumb}]];

If[StringMatchQ[OptionValue[CrackDirection],"x*",IgnoreCase->True],
Which[NumberQ[OptionValue[ReferenceStrain]],
globalStrain=OptionValue[ReferenceStrain];
,
StringMatchQ[OptionValue[ReferenceStrain],"g*",IgnoreCase->True],
globalStrain=m/.FindFit[goodData[[All,{XNumb,YNumb,VNumb}]],a+b*x+m*y,{a,b,m},{x,y}];
,
StringMatchQ[OptionValue[ReferenceStrain],"l*",IgnoreCase->True],
globalStrain="Local";
];

codProfiles=Map[measureCOD[#,w,CrackDirection->"x",CrackThreshold->OptionValue[CrackThreshold],GlobalStrain->globalStrain]&,Transpose[gridData,{2,1,3}]];
,
Which[NumberQ[OptionValue[ReferenceStrain]],
globalStrain=OptionValue[ReferenceStrain];
,
StringMatchQ[OptionValue[ReferenceStrain],"g*",IgnoreCase->True],
globalStrain=m/.FindFit[goodData[[All,{XNumb,YNumb,UNumb}]],a+m*x+b*y,{a,m,b},{x,y}];
,
StringMatchQ[OptionValue[ReferenceStrain],"l*",IgnoreCase->True],
globalStrain="Local";
];

codProfiles=Map[measureCOD[#,w,CrackDirection->"y",CrackThreshold->OptionValue[CrackThreshold],GlobalStrain->globalStrain]&,gridData];
];

{codProfiles[[All,1]],codProfiles[[All,2]]}
,{grid,Length[gridsPos]}]

,{frame,Length[dataIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[NumberQ[frames],
If[First[Dimensions[gridsIn]]==1,
dataOut[[1,1]]
,
dataOut[[1]]
]
,
If[First[Dimensions[gridsIn]]==1,
dataOut[[All,1]]
,
dataOut
]
]

];


(* ::Subsection::Closed:: *)
(*Get Coordinates*)


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
yMax=Length[greyPlt];
AppendTo[coordinates,nodes];
AppendTo[coordinatesOut,Map[{#[[1]],yMax-#[[2]]}&,nodes]];]}}]]
];


guiNodes={{0,0},{1,1},{2,2}};
coordinates={};
coordinatesOut={};


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
Style["\[FilledCircle]",FontColor->RGBColor[0.5,0.1,0.4],FontSize->12]
}]},
{Grid[{{Button["Clear coordinates",coordinates={};],
Button["Clear coordinatesOut",coordinatesOut={};]}}]},
{Button["Save coordinates",
yMax=Length[greyPlt];
AppendTo[coordinates,{guiNodes[[1;;2]],{{Mean[guiNodes[[1;;2,1]]],Mean[guiNodes[[1;;2,2]]]},guiNodes[[3]]}}];
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
yMax=Length[greyPlt];
v1=guiNodes[[2]]-guiNodes[[1]];
v2=guiNodes[[3]]-guiNodes[[1]];
v3=(v1.v2)/(v1.v1)v1;
vo=v2-(v1.v2)/(v1.v1)v1;
AppendTo[coordinates,{guiNodes[[1;;2]],Round[{guiNodes[[1]]+vo,guiNodes[[2]]+vo}]}];
AppendTo[coordinatesOut,Map[{#[[1]],yMax-#[[2]]}&,{guiNodes[[1;;2]],Round[{guiNodes[[1]]+vo,guiNodes[[2]]+vo}]},{2}]];
]}}]]
];


(* ::Subsection::Closed:: *)
(*Contour Overlays*)


Options[getDeformedData]={
DeformedXY->True,CenterDisp->True,FixRotation->False,DeformedZ->False
};


getDeformedData[pathIn_,frame_,varIn_,OptionsPattern[]]:=Module[
{filesIn,XNumb,ZNumb,YNumb,UNumb,VNumb,WNumb,varNumb,sigmaNumb,dataIn,xCenter,yCenter,uCenterLine,vCenterLine,uLine,vLine,m,b,theta,meanDisp},

filesIn=DICdataFiles[pathIn];

XNumb=Position[filesIn[[2]],"X"][[1,1]];
YNumb=Position[filesIn[[2]],"Y"][[1,1]];
ZNumb=Position[filesIn[[2]],"Z"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
WNumb=Position[filesIn[[2]],"W"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];

If[OptionValue[DeformedXY],
dataIn[[All,All,{XNumb,YNumb}]]=dataIn[[All,All,{XNumb,YNumb}]]+dataIn[[All,All,{UNumb,VNumb}]];
];

Which[varIn=="U",
If[OptionValue[FixRotation],
xCenter=Round[Length[dataIn[[1]]]/2];
yCenter=Round[Length[dataIn]/2];
uCenterLine=Delete[dataIn[[All,xCenter,{YNumb,UNumb}]],Position[dataIn[[All,xCenter,sigmaNumb]],-1.]];
vCenterLine=Delete[dataIn[[yCenter,All,{XNumb,VNumb}]],Position[dataIn[[yCenter,All,sigmaNumb]],-1.]];
uLine=Function[{y},Evaluate[m y +b/.FindFit[uCenterLine,m y+b,{m,b},y]]];
vLine=Function[{x},Evaluate[m x +b/.FindFit[vCenterLine,m x +b,{m,b},x]]];
theta=(Abs[ArcTan[(vLine[vCenterLine[[-1,1]]]-vLine[vCenterLine[[1,1]]])/(vCenterLine[[-1,1]]-vCenterLine[[1,1]])]]+Abs[ArcTan[(uLine[uCenterLine[[-1,1]]]-uLine[uCenterLine[[1,1]]])/(uCenterLine[[-1,1]]-uCenterLine[[1,1]])]])/2;

dataIn[[All,All,UNumb]]=dataIn[[All,All,UNumb]]-theta*dataIn[[All,All,YNumb]];
];

If[OptionValue[CenterDisp],
meanDisp=Mean[Select[Flatten[dataIn,1],(#[[sigmaNumb]]!=-1.)&][[All,UNumb]]];
dataIn[[All,All,UNumb]]=Map[#-meanDisp&,dataIn[[All,All,UNumb]],{2}];
];
,
varIn=="V",
If[OptionValue[FixRotation],
xCenter=Round[Length[dataIn[[1]]]/2];
yCenter=Round[Length[dataIn]/2];
uCenterLine=Delete[dataIn[[All,xCenter,{YNumb,UNumb}]],Position[dataIn[[All,xCenter,sigmaNumb]],-1.]];
vCenterLine=Delete[dataIn[[yCenter,All,{XNumb,VNumb}]],Position[dataIn[[yCenter,All,sigmaNumb]],-1.]];
uLine=Function[{y},Evaluate[m y +b/.FindFit[uCenterLine,m y+b,{m,b},y]]];
vLine=Function[{x},Evaluate[m x +b/.FindFit[vCenterLine,m x +b,{m,b},x]]];
theta=(Abs[ArcTan[(vLine[vCenterLine[[-1,1]]]-vLine[vCenterLine[[1,1]]])/(vCenterLine[[-1,1]]-vCenterLine[[1,1]])]]+Abs[ArcTan[(uLine[uCenterLine[[-1,1]]]-uLine[uCenterLine[[1,1]]])/(uCenterLine[[-1,1]]-uCenterLine[[1,1]])]])/2;

dataIn[[All,All,VNumb]]=dataIn[[All,All,VNumb]]+theta*dataIn[[All,All,XNumb]];
];

If[OptionValue[CenterDisp],
meanDisp=Mean[Select[Flatten[dataIn,1],(#[[sigmaNumb]]!=-1.)&][[All,VNumb]]];
dataIn[[All,All,VNumb]]=Map[#-meanDisp&,dataIn[[All,All,VNumb]],{2}];
];
,
varIn=="Z",
If[OptionValue[DeformedZ],
dataIn[[All,All,ZNumb]]=dataIn[[All,All,ZNumb]]+dataIn[[All,All,WNumb]];
];
];

dataIn[[All,All,{XNumb,YNumb,varNumb,sigmaNumb}]]

];


Options[DICoverLay]={
ContourOpacity->1,contourShading->Automatic,ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,plotRange->All,color->"Rainbow",NumberofPoints->50000,size->500,overlay->True,CenterDisp->True,FixRotation->False,DeformedZ->False,FilterSize->1
};


DICoverLay[filesPath_,frame_,varIn_,OptionsPattern[]]:=Module[{
filesIn,imagePathIn,xNumb,yNumb,uNumb,vNumb,UNumb,VNumb,WNumb,varNumb,sigmaNumb,dataIn,yMax,boole,takeNumber,xCenter,yCenter,uCenterLine,vCenterLine,uLine,vLine,m,b,theta,meandDisp,reducedData,meanSigma,incr,contourList,colorFS,colorF,contPlt,greyPlt
},

filesIn=allDataFiles[filesPath];
imagePathIn=filesIn[[3,frame]];

greyPlt=Import[imagePathIn,"GrayLevels"];

xNumb=Position[filesIn[[2]],"x"][[1,1]];
yNumb=Position[filesIn[[2]],"y"][[1,1]];
uNumb=Position[filesIn[[2]],"u"][[1,1]];
vNumb=Position[filesIn[[2]],"v"][[1,1]];
UNumb=Position[filesIn[[2]],"U"][[1,1]];
VNumb=Position[filesIn[[2]],"V"][[1,1]];
WNumb=Position[filesIn[[2]],"W"][[1,1]];
varNumb=Position[filesIn[[2]],varIn][[1,1]];
sigmaNumb=Position[filesIn[[2]],"sigma"][[1,1]];

dataIn=MapThread[List,Import[filesIn[[1,frame]]],2];
yMax=Length[greyPlt]+1;

dataIn[[All,All,{xNumb,yNumb}]]=dataIn[[All,All,{xNumb,yNumb}]]+dataIn[[All,All,{uNumb,vNumb}]];

Which[varIn=="U",
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
varIn=="V",
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


(* ::Subsection::Closed:: *)
(*Legends*)


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


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["ThreeDDICFunctions`*"],
{
"nodes","coordinates","coordinatesOut","guiNodes"
}];
EndPackage[ ]
