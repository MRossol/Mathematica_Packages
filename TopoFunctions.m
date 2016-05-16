(* ::Package:: *)

BeginPackage["TopoFunctions`"]
Unprotect@@Names["TopoFunctions`*"];
ClearAll@@Names["TopoFunctions`*"];


(*Data Import*)
getDICfiles::usage = "getDICfiles[pathIn_] returns DIC files and variable list from v6 files folder and get image files from parent directory";
getStitchData::usage = "getStitchData[filesPath_,frames_] returns a list of DICdata of {X,Y,Z,sigma} for the given frames";

(*DIC Analysis*)
getMagnificationandStep::usage = "getMagnification[filesPath_] returns {{pixels/mm, mm/pixels},{step(pixels),step (mm)}} for dataIn";

(*Stitching*)
stitchValues::usage = "stitchValues[dataIn_,vars_] calculates the sum of squares distance between ref points and data points";
stitchDICdata::usage = "stitchDICdata[dataIn_,disps_, stitches DIC dataIn using given individual disps in horizontal (right or left) or vertical(up or down) direction with interpolation";
stitchVertical::usage = "stitchVertical1[dataIn_,disps_] stitches vertical set of DIC dataIn together using exact disps given through interpolatino functions in metric space to remove perspective";
stitchHorizontal::usage = "stitchHorizontal1[dataIn_,disps_] stitches horizontal set of DIC dataIn together using exact disps given through interpolatino functions in metric space to remove perspective";
autoPlaneFit::usage = "autoPlaneFit[dataIn] auto plane fits the input data";

(*Stitch Post Processing*)
zBarLocal::usage = "zBarLocal[dataIn_,unitCell_,step_] calculates the average Z value of the surrounding unit cell at each point";
getThickness::usage = "getThickness[{frontIn_,backIn_},aveThickness] calculates the local thickenss at every point";
getCurvatures::usage = "getCurvatures[dataIn_] calculates the principle curvatures and mean curvature for the surfaceIn outputs {{k1,k2},H}";

(*Thresholding, Clustering, and Segmentation*)
arrayClusters::usage = "arrayClustering[dataIn_,{thresholdMin_,thresholdMax_},OptionsPattern[]] uses MorphologicalComponents to find clusters in arrays of data.";
arrayClusterWindow::usage = "arrayClusterWindow[dataIn_,{thresholdMin_,thresholdMax_},{center_,sizes_},OptionsPattern[]] applies arrayClusters to a window (subset) of dataIn with given center and sizes";
imageClusters::usage = "imageClusters[dataIn_] uses MorphologicalComponents to find clusters in a binarized image or binarized image data.";
thresholdWindow::usage = "thresholdWindow[dataIn_,{thresholdMin_,thresholdMax_},{center_,sizes_},OptionsPattern[]] thresholds the given window";
getClusters::usage = "getClusters[dataIn, cutOff] finds clusters of data using the cut of distance cutOff, can iterate over the adcajency matrix if needed to improve connectivity";
gatherClusters::usage = "gatherClusters[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},OptionsPattern[]] gathers all clusters in dataIn (array)";
gatherListClusters::usage = "gatherListClusters[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},OptionsPattern[]] gather all clusters for dataIn as a list";
imageListClusters::usage = "imageListClusters[dataIn_,cutOff_,OptionsPattern[]] gather all clusters for image dataIn as a list";
clusterWindow::usage = "clusterWindow[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},{center_,size_},OptionsPattern[]] extracts clusters in the window of given size with given center";
segregateTows::usage = "segregates tows based on ratio of steDev[x]/stDev[y]";
combineClusters::usage = "combineClusters[dataIn1_,dataIn2_] combines clusters in dataIn2 with existing clusters in dataIn1";

(*Tow Analysis*)
x::usage = "equation variable";
y::usage = "equation variable";
th::usage = "rotation variable";
rSquared::usage = "rSquared[dataIn_,fitIn_] calculates the R^2 for the input data and fit";
findMax::usage = "findMax[dataIn_,model_,parameters_,variables_] fits dataIn to the given model with given independant variables and unknown parameter and used the fit to find the local Maximum";
findMin::usage = "findMin[dataIn_,model_,parameters_,variables_] fits dataIn to the given model with given independant variables and unknown parameter and used the fit to find the local Minimum";
findCrest::usage = "findCrest[dataIn_] finds the crest of the tow, the raw data and fit data can be used";
findTrough::usage = "findTrough[dataIn_] finds the Trough of the tow, the raw data and fit data can be used";

(*Rotation and Translation*)
flipData::usage = "Flips data about specified axis";
rotateData::usage = "rotateData[dataIn_,rot_] applies the given rotation matrix to dataIn";
shiftData::usage = "shiftData[dataIn_,shift_] shifts data by shift={\[CapitalDelta]x,\[CapitalDelta]y,\[CapitalDelta]z}";
rotateShiftData::usage = "rotateShiftData[dataIn_,rot_,shift_] rotates data by the given rotation matrix then shifts data by the given {\[CapitalDelta]x,\[CapitalDelta]y,\[CapitalDelta]z}";
deltaH::usage = "deltaH[dataIn_,vars_] is minimized with respect to {z,thy,thx} to determine rotation and shift of back needed to align out of plane with front.";

(*Autocorrelation and Variances*)
spatialAC::usage = "spatialAC[dataIn_,bin_,OptionsPattern[]] performs spatial autocorrelation on dataIn using bin size bin";
posVariances::usage = "posVariances[dataIn_,bin_,OptionsPattern[]] extracts positional variances from dataIn output is {x,y,var_x,var_xy,var_y,n}";
gatherACclusters::usage = "gatherACclusters[dataIn_,cutOff_,OptionsPattern[]] gather all clusters for AC or posVariance data in (gatherListClusters with out the thresholding)";

(*Angular Analysis*)
angleMatrix::usage = "angleMatrix[dataIn_,OptionsPattern[]] calculate rotation matrix for dataIn, or matrix of angle between every point in dataIn";
angularHist::usage = "angularHist[dataIn_,bin_,OptionsPattern[]] converts angleMatrix data to histogram values with given bin size";
peakLocation::usage = "peakLocation[dataIn_,bin_,range_,theta_] calculates the angle difference positive and negative peaks or between peak and PeakValue";
peakFitting::usage = "peakFitting[dataIn_,bin_,range_,theta_] calculates the angle difference but outputs {{{anglePos,angleNeg},angleDiff},{posInterp,negInterp},{posData,negData}}";
peakInfo::usage = "peakInfo[dataIn_,peaksIn_] finds the maximum in the vicinity of peaks in and returns (xMax,yMax,width) for peaksIn";
getShearAngle::usage = "getShearAngle[anglesIn_,unitCellAngles] calculates the pure shear angle from anglesIn = {ideal\[Theta], d\[Theta]} and the ideal unit cell angles = {ideal\[Theta],{ux,uy}}";

(*Displacement and Strain Analysis*)
dispFromIdeal::usage = "dispFromIdeal[dataIn_,idealIn_] calculates displacements {dx,dy} or {u,v} between dataIn and idealIn";
strainFilter::usage = "strainFilter[dataIn_,cutOffR_] calculates strains from dispIn using neighbors within cutOffR = r or (rx,ry)";

(*Neighbors*)
getNN::usgae = "getNN[dispIn_,cutOffR_] gets neighbors within cutOffR = r or (rx,ry)";
getNeighbors::usage = "getNeighbors[dataIn_,angleRange_,distanceRange_,OptionsPattern[]]gets neighbors within angles and distances";
neighborDistance::usage = "neighborDistance[neighborsIn_,vars_] extracts the sum of square of the distances between neighborsIn";
alignmentDistance::usage = "alignmentDistance[neighborsIn_,vars_] extracts the sum of square of the distances between points to be aligned";

(*Line Scans*)
lineScan::usage = "lineScan[dataIn_,lines_] extracts lines in";

(*Get Reposition Displacements*)
nodalCoords1::usage = "nodes with {{x1,y1},{x2,y2}}";
nodalCoords2::usage = "nodes with {{x1,y1},{x2,y2}}";
RepositionDispsOut::usage = "List of reposition displacements in mm taken from the mean of end points on lines on TIFF image 1 and image 2";
getRepositionDisp::usage = "getRepositionDisp[filesPath_,{Frame1_,Frame2_}] allows for selection of corresponding points on two sequential TIFF images";
StitchCoordsOut::usage = "List of reposition Coords for corresponding points on data sets to be stitched ";
getStitchCoords::usage = "getRepositionDisp[filesPath_,{Frame1_,Frame2_}] allows for selection of corresponding points on two sequential TIFF images";


(*Contour Overlays*)
DICcontourPlot::usage = "DICcontourPlot[dataIn_ renders a surface topography plot of the input stitch data with options to AutoPlaneFit dataIn";
DICcontourPlot2::usage = "DICcontourPlot2[dataIn_ renders a surface topography plot of the input stitch data with options to AutoPlaneFit dataIn, and to plot with or with out bad data (region function)";

(*Selection GUIs*)
points::usage = "dynamic points {{x1,y1},{x2,y2}} for GUI getPoints";
pointsOut::usage = "saved Points from GUI getPoints";
getPoints::usage = "getPoints[pltIn_] select points on pltIn using GUI";


Begin["`Private`"]


Needs["GUIKit`"]
Needs["HierarchicalClustering`"]


ProgressDialog[]:=GUIRun[ 
Widget["Frame",{ 
WidgetGroup[{
Widget["Label",{"text"->"Percent complete:"},
Name->"label"],
Widget["ProgressBar",
{"minimum"->0,"maximum"->100,
"preferredSize"->
  Widget["Dimension",{"width"->300,"height"->25}]},
Name->"bar"] 
  }, 
WidgetLayout -> {
"Grouping"-> Column, 
"Border" -> {{15,15},{25,20}}}],
"location"->Widget["Point",{"x"->400,"y"->400}],"title"->"Computation Progress",
"resizable"->False},
Name->"frame"]
];


(* ::Subsection::Closed:: *)
(*Data Import*)


getDICfiles[pathIn_]:=Module[{
matv6Path,varList,matv6Files,imgFiles,imgPath
},

matv6Path=FileNames["*v6 files*",pathIn,Infinity];
If[Length[matv6Path]>1,
matv6Path=ChoiceDialog["Multiple data sources found...\n  Pick v6 .mat data location:",(FileBaseName[ParentDirectory[#]]->#)&/@matv6Path];,
matv6Path=matv6Path[[1]]
];

matv6Files=FileNames["*.mat",matv6Path];
imgPath=ParentDirectory[matv6Path];
imgFiles=Flatten[FileNames[FileBaseName[#]<>".tif",imgPath]&/@matv6Files];varList=StringTrim/@(ReadList[FileNameJoin[{matv6Path,"index_list.txt"}],{Number,String}][[All,2]]);

{matv6Files,varList,imgFiles}
];


getStitchData[filesPath_,frames_]:=Module[
{filesIn,varListIn,xNumb,yNumb,XNumb,YNumb,ZNumb,sigmaNumb
},

filesIn=getDICfiles[filesPath];

varListIn=filesIn[[2]];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];
ZNumb=Position[varListIn,"Z"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

Map[MapThread[List,Import[filesIn[[1,#]]],2][[All,All,{XNumb,YNumb,ZNumb,sigmaNumb}]]&,frames]
];


(* ::Subsection::Closed:: *)
(*DIC Analysis*)


getMagnificationandStep[filesPath_]:=Module[{
filesIn,varListIn,xNumb,yNumb,XNumb,YNumb,UNumb,VNumb,WNumb,sigmaNumb,ref,calFactor={},step={},data,goodData
},

$HistoryLength=1;

filesIn=getDICfiles[filesPath];

varListIn=filesIn[[2]];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

ref=ProgressDialog[];

Do[

ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/(Length[filesIn[[1]]]))]];
ref @ SetPropertyValue[{"label", "text"}, "Analyzing File "<>ToString[(i)]<>" of "<>ToString[(Length[filesIn[[1]]])]];

data=MapThread[List,Import[filesIn[[1,i]]],2];

goodData=SortBy[Delete[Flatten[data,1],Position[Flatten[data[[All,All,sigmaNumb]],1],-1.]],(#[[xNumb]]+#[[yNumb]])&];

AppendTo[calFactor,EuclideanDistance[goodData[[-1,{XNumb,YNumb}]],goodData[[1,{XNumb,YNumb}]]]/EuclideanDistance[goodData[[-1,{xNumb,yNumb}]],goodData[[1,{xNumb,yNumb}]]]];
AppendTo[step,Abs[data[[1,1,xNumb]]-data[[1,2,xNumb]]]];

,{i,Length[filesIn[[1]]]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

Share[];

{{1/Mean[calFactor],Mean[calFactor]},{Mean[step],Mean[step]*Mean[calFactor]}}
];


(* ::Subsection::Closed:: *)
(*Stitching*)


stitchValues[dataIn_,vars_]:=Module[
{ref,data,rot},

ref=dataIn[[1]];
data=dataIn[[2]];

rot=RotationMatrix[vars[[3]]Degree];

data=Transpose[rot.Transpose[data]];
data[[All,1]]=data[[All,1]]+vars[[1]];
data[[All,2]]=data[[All,2]]+vars[[2]];

Sum[(EuclideanDistance[ref[[i]],data[[i]]]^2),{i,Length[ref]}]
];


Options[stitchDICdata]={
StitchDirection->"Down",AutoPlaneFit->False,UseMeanSlope->False,RetainEdges->False,SortInterp->True
};


stitchDICdata[dataIn_,disps_,stepSize_,OptionsPattern[]]:=Module[{
dataOut,xo,yo,planeFit,a,b,c,x,y,th,xdisp=0,ydisp=0,theta=0,dispsOut,ref,data,rot
},

$HistoryLength=1;
If[Length[Dimensions[disps]]==2,
If[StringMatchQ[OptionValue[StitchDirection],"d*",IgnoreCase->True]||StringMatchQ[OptionValue[StitchDirection],"r*",IgnoreCase->True],
dispsOut=Table[
xdisp=xdisp+disps[[i,1]];
ydisp=ydisp+disps[[i,2]];
{xdisp,ydisp,0},
{i,1,Length[disps]}];
,
dispsOut=Table[
xdisp=xdisp-disps[[i,1]];
ydisp=ydisp-disps[[i,2]];
{xdisp,ydisp,0},
{i,1,Length[disps]}];
];
,
If[StringMatchQ[OptionValue[StitchDirection],"d*",IgnoreCase->True]||StringMatchQ[OptionValue[StitchDirection],"r*",IgnoreCase->True],
dispsOut=Table[
ref=disps[[i,1]];
data=disps[[i,2]];
rot=RotationMatrix[theta Degree];
ref=Transpose[rot.Transpose[ref]];
ref[[All,1]]=ref[[All,1]]+xdisp;
ref[[All,2]]=ref[[All,2]]+ydisp;

{xdisp,ydisp,theta}={x,y,th}/.Minimize[stitchValues[{ref,data},{x,y,th}],{x,y,th}][[2]];
{xdisp,ydisp,theta},
{i,1,Length[disps]}];
,
dispsOut=Table[
ref=disps[[i,2]];
data=disps[[i,1]];
rot=RotationMatrix[theta Degree];
ref=Transpose[rot.Transpose[ref]];
ref[[All,1]]=ref[[All,1]]+xdisp;
ref[[All,2]]=ref[[All,2]]+ydisp;

{xdisp,ydisp,theta}={x,y,th}/.Minimize[stitchValues[{ref,data},{x,y,th}],{x,y,th}][[2]];
{xdisp,ydisp,theta},
{i,1,Length[disps]}];
];
];

If[StringMatchQ[OptionValue[StitchDirection],"u*",IgnoreCase->True]||StringMatchQ[OptionValue[StitchDirection],"d*",IgnoreCase->True],
dataOut=stitchVertical[dataIn,dispsOut,stepSize,UseMeanSlope->OptionValue[UseMeanSlope],RetainEdges->OptionValue[RetainEdges],SortInterp->OptionValue[SortInterp]];
,
dataOut=stitchHorizontal[dataIn,dispsOut,stepSize,UseMeanSlope->OptionValue[UseMeanSlope],RetainEdges->OptionValue[RetainEdges],SortInterp->OptionValue[SortInterp]];
];

Share[dataOut];

If[OptionValue[AutoPlaneFit],
planeFit[x_,y_]=a+b*x+c*y/.FindFit[Flatten[dataOut[[1]],1],a+b*x+c*y,{a,b,c},{x,y}];
xo=dataOut[[1,-1,1,1]];
yo=dataOut[[1,-1,1,2]];

{Map[{(#[[1]]-xo),
(#[[2]]-yo),(#[[3]]-planeFit[#[[1]],#[[2]]])}&,dataOut[[1]],{2}],dataOut[[2]]}
,
dataOut
]

];


Options[stitchVertical]={
UseMeanSlope->False,RetainEdges->False,SortInterp->True
};


stitchVertical[dataIn_,disps_,step_,OptionsPattern[]]:=Module[{
xMins={},xMaxs={},yStarts={},yCenters={},zShift=0,interpsOut={},zShiftFits={},stitchDisp,ref,dataIn1,dataIn2,rot,yShiftCenter,yCenter,centerLine,data1F,data2F,xMin,xMax,interp1,interp2,edge,centerLine1,centerLine2,slope,m1,m2,x,b,b1,b2,line1,line2,fitLine1,fitLine2,xCoords,centerNearest,yCoords,dataOut
},

$HistoryLength=1;

stitchDisp={{0,0,0}}~Join~disps;

ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/(Length[dataIn]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing data pair "<>ToString[(i-1)]<>" of "<>ToString[(Length[dataIn]-1)]];

dataIn1=dataIn[[i-1]];
dataIn2=dataIn[[i]];

rot=RotationMatrix[stitchDisp[[i-1,3]] Degree];
dataIn1[[All,All,1;;2]]=Map[rot.#&,dataIn1[[All,All,1;;2]],{2}];
dataIn1[[All,All,1]]=dataIn1[[All,All,1]]+stitchDisp[[i-1,1]];
dataIn1[[All,All,2]]=dataIn1[[All,All,2]]+stitchDisp[[i-1,2]];
dataIn1[[All,All,3]]=dataIn1[[All,All,3]]+zShift;

rot=RotationMatrix[stitchDisp[[i,3]] Degree];
dataIn2[[All,All,1;;2]]=Map[rot.#&,dataIn2[[All,All,1;;2]],{2}];
dataIn2[[All,All,1]]=dataIn2[[All,All,1]]+stitchDisp[[i,1]];
dataIn2[[All,All,2]]=dataIn2[[All,All,2]]+stitchDisp[[i,2]];

If[Length[dataIn1[[1,1]]]<4.,
data1F=Flatten[dataIn1,1];
,
data1F=Delete[Flatten[dataIn1,1],Position[Flatten[dataIn1,1][[All,4]],-1.]];
];

If[Length[dataIn2[[1,1]]]<4.,
data2F=Flatten[dataIn2,1];
,
data2F=Delete[Flatten[dataIn2,1],Position[Flatten[dataIn2,1][[All,4]],-1.]];
];

If[Positive[Mean[disps][[2]]],
AppendTo[yStarts,Min[data1F[[All,2]]]];
yShiftCenter=(Min[data2F[[All,2]]]+(Max[data1F[[All,2]]]-Min[data2F[[All,2]]])/2);
,
AppendTo[yStarts,Max[data1F[[All,2]]]];
yShiftCenter=(Min[data1F[[All,2]]]+(Max[data2F[[All,2]]]-Min[data1F[[All,2]]])/2);
];

AppendTo[yCenters,yShiftCenter];

If[OptionValue[RetainEdges],
xMin=Min[{Min[data1F[[All,1]]],Min[data2F[[All,1]]]}];
xMax=Max[{Max[data1F[[All,1]]],Max[data2F[[All,1]]]}];
,
xMin=Max[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{Min[data1F[[All,1]]],yShiftCenter}]][[1]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{Min[data2F[[All,1]]],yShiftCenter}]][[1]]}];
xMax=Min[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{Max[data1F[[All,1]]],yShiftCenter}]][[1]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{Max[data2F[[All,1]]],yShiftCenter}]][[1]]}];
];

AppendTo[xMins,xMin];
AppendTo[xMaxs,xMax];

xMin=Max[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{Min[data1F[[All,1]]],yShiftCenter}]][[1]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{Min[data2F[[All,1]]],yShiftCenter}]][[1]]}];
xMax=Min[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{Max[data1F[[All,1]]],yShiftCenter}]][[1]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{Max[data2F[[All,1]]],yShiftCenter}]][[1]]}];

If[OptionValue[SortInterp],
interp1=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data1F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data2F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
,
interp1=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data1F],InterpolationOrder->1];
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data2F],InterpolationOrder->1];
];

edge=Round[Length[Range[xMin,xMax,step]]/100.];

centerLine1=Table[{i,interp1[i,yShiftCenter]},{i,xMin+step*edge,xMax-step*edge,step}];
centerLine2=Table[{i,interp2[i,yShiftCenter]},{i,xMin+step*edge,xMax-step*edge,step}];

If[OptionValue[UseMeanSlope],
slope=Mean[{m1,m2}]/.FindFit[centerLine1,m1*x+b,{m1,b},x]/.FindFit[centerLine2,m2*x+b,{m2,b},x];
zShift=(b1-b2)/.FindFit[centerLine1,slope*x+b1,{b1},x]/.FindFit[centerLine2,slope*x+b2,{b2},x];

line1[x_]=slope*x+b/.FindFit[centerLine1,slope*x+b,{b},x];
line2[x_]=slope*x+b/.FindFit[centerLine2,slope*x+b,{b},x];

fitLine1=Table[{i,line1[i]},{i,xMin+step*edge,xMax-step*edge,step}];
fitLine2=Table[{i,line2[i]},{i,xMin+step*edge,xMax-step*edge,step}];

AppendTo[zShiftFits,{centerLine1,centerLine2,fitLine1,fitLine2}];
,
zShift=(Mean[centerLine1[[All,2]]]-Mean[centerLine2[[All,2]]]);

AppendTo[zShiftFits,{zShift,centerLine1,centerLine2}];
];

AppendTo[interpsOut,interp1];

,{i,2,Length[dataIn]}];

data2F[[All,3]]=data2F[[All,3]]+zShift;

If[OptionValue[SortInterp],
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data2F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
,
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data2F],InterpolationOrder->1];
];

AppendTo[interpsOut,interp2];

If[OptionValue[RetainEdges],
xCoords={Min[xMins],Max[xMaxs]};
,
xCoords={Max[xMins],Min[xMaxs]};
];

If[Positive[Mean[disps][[2]]],
centerNearest=Nearest[Table[i,{i,yStarts[[1]],Max[data2F[[All,2]]],step}]];
yCoords=Reverse[Partition[{yStarts[[1]]}~Join~Flatten[Map[{First[centerNearest[#-step]],First[centerNearest[#]]}&,yCenters]]~Join~{Max[data2F[[All,2]]]},2],2];

dataOut=Table[{i,j,interpsOut[[1]][i,j]},
{i,yCoords[[1,1]], yCoords[[1,2]],-1.*step},
{j,xCoords[[1]],xCoords[[2]],step}];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(k/(Length[interpsOut]))]];
ref @ SetPropertyValue[{"label", "text"}, "Compiling stitch "<>ToString[(k)]<>" of "<>ToString[(Length[interpsOut])]];

dataOut=Join[Table[{i,j,interpsOut[[k]][i,j]},
{j,yCoords[[k,1]], yCoords[[k,2]],-1.*step},
{i,xCoords[[1]],xCoords[[2]],step}],dataOut];

,{k,2,Length[interpsOut]}];
,
centerNearest=Nearest[Table[i,{i,Min[data2F[[All,2]]],yStarts[[1]],step}]];
yCoords=Partition[{yStarts[[1]]}~Join~Flatten[Map[{First[centerNearest[#+step]],First[centerNearest[#]]}&,yCenters]]~Join~{Min[data2F[[All,2]]]},2];

dataOut=Table[{i,j,interpsOut[[1]][i,j]},
{j,yCoords[[1,1]], yCoords[[1,2]],-1.*step},
{i,xCoords[[1]],xCoords[[2]],step}];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(k/(Length[interpsOut]))]];
ref @ SetPropertyValue[{"label", "text"}, "Compiling stitch "<>ToString[(k)]<>" of "<>ToString[(Length[interpsOut])]];

dataOut=Join[dataOut,Table[{i,j,interpsOut[[k]][i,j]},
{j,yCoords[[k,1]], yCoords[[k,2]],-1.*step},
{i,xCoords[[1]],xCoords[[2]],step}]];

,{k,2,Length[interpsOut]}];
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

Share[];

{dataOut,zShiftFits}

];


Options[stitchHorizontal]={
UseMeanSlope->False,RetainEdges->False,SortInterp->True
};


stitchHorizontal[dataIn_,disps_,step_,OptionsPattern[]]:=Module[{
xStarts={},xCenters={},yMins={},yMaxs={},zShift=0,interpsOut={},zShiftFits={},stitchDisp,ref,dataIn1,dataIn2,rot,data1F,data2F,xShiftCenter,xCenter,centerLine,yMin,yMax,interp1,interp2,edge,centerLine1,centerLine2,slope,m1,m2,x,b,b1,b2,line1,line2,fitLine1,fitLine2,yCoords,centerNearest,xCoords,dataOut
},

$HistoryLength=1;

stitchDisp={{0,0,0}}~Join~disps;

ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/(Length[dataIn]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing data pair "<>ToString[(i-1)]<>" of "<>ToString[(Length[dataIn]-1)]];

dataIn1=dataIn[[i-1]];
dataIn2=dataIn[[i]];

rot=RotationMatrix[stitchDisp[[i-1,3]] Degree];
dataIn1[[All,All,1;;2]]=Map[rot.#&,dataIn1[[All,All,1;;2]],{2}];
dataIn1[[All,All,1]]=dataIn1[[All,All,1]]+stitchDisp[[i-1,1]];
dataIn1[[All,All,2]]=dataIn1[[All,All,2]]+stitchDisp[[i-1,2]];
dataIn1[[All,All,3]]=dataIn1[[All,All,3]]+zShift;

rot=RotationMatrix[stitchDisp[[i,3]] Degree];
dataIn2[[All,All,1;;2]]=Map[rot.#&,dataIn2[[All,All,1;;2]],{2}];
dataIn2[[All,All,1]]=dataIn2[[All,All,1]]+stitchDisp[[i,1]];
dataIn2[[All,All,2]]=dataIn2[[All,All,2]]+stitchDisp[[i,2]];

If[Length[dataIn1[[1,1]]]<4.,
data1F=Flatten[dataIn1,1];
,
data1F=Delete[Flatten[dataIn1,1],Position[Flatten[dataIn1,1][[All,4]],-1.]];
];
If[Length[dataIn2[[1,1]]]<4.,
data2F=Flatten[dataIn2,1];
,
data2F=Delete[Flatten[dataIn2,1],Position[Flatten[dataIn2,1][[All,4]],-1.]];
];

If[Positive[Mean[disps][[1]]],
AppendTo[xStarts,Min[data1F[[All,1]]]];
xShiftCenter=(Min[data2F[[All,1]]]+(Max[data1F[[All,1]]]-Min[data2F[[All,1]]])/2);
,
AppendTo[xStarts,Max[data1F[[All,1]]]];
xShiftCenter=(Min[data1F[[All,1]]]+(Max[data2F[[All,1]]]-Min[data1F[[All,1]]])/2);
];
AppendTo[xCenters,xShiftCenter];

If[OptionValue[RetainEdges],
yMin=Min[{Min[data1F[[All,2]]],Min[data2F[[All,2]]]}];
yMax=Max[{Max[data1F[[All,2]]],Max[data2F[[All,2]]]}];
,
yMin=Max[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{xShiftCenter,Min[data1F[[All,2]]]}]][[2]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{xShiftCenter,Min[data2F[[All,2]]]}]][[2]]}];
yMax=Min[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{xShiftCenter,Max[data1F[[All,2]]]}]][[2]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{xShiftCenter,Max[data2F[[All,2]]]}]][[2]]}];
];

AppendTo[yMins,yMin];
AppendTo[yMaxs,yMax];

yMin=Max[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{xShiftCenter,Min[data1F[[All,2]]]}]][[2]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{xShiftCenter,Min[data2F[[All,2]]]}]][[2]]}];
yMax=Min[{First[Nearest[Flatten[dataIn1[[All,All,{1,2}]],1],{xShiftCenter,Max[data1F[[All,2]]]}]][[2]],First[Nearest[Flatten[dataIn2[[All,All,{1,2}]],1],{xShiftCenter,Max[data2F[[All,2]]]}]][[2]]}];

If[OptionValue[SortInterp],
interp1=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data1F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data2F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
,
interp1=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data1F],InterpolationOrder->1];
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data2F],InterpolationOrder->1];
];

edge=Round[Length[Range[yMin,yMax,step]]/100.];

centerLine1=Table[{i,interp1[xShiftCenter,i]},{i,yMin+step*edge,yMax-step*edge,step}];
centerLine2=Table[{i,interp2[xShiftCenter,i]},{i,yMin+step*edge,yMax-step*edge,step}];

If[OptionValue[UseMeanSlope],
slope=Mean[{m1,m2}]/.FindFit[centerLine1,m1*x+b,{m1,b},x]/.FindFit[centerLine2,m2*x+b,{m2,b},x];
zShift=(b1-b2)/.FindFit[centerLine1,slope*x+b1,{b1},x]/.FindFit[centerLine2,slope*x+b2,{b2},x];

line1[x_]=slope*x+b/.FindFit[centerLine1,slope*x+b,{b},x];
line2[x_]=slope*x+b/.FindFit[centerLine2,slope*x+b,{b},x];

fitLine1=Table[{i,line1[i]},{i,yMin+step*edge,yMax-step*edge,step}];
fitLine2=Table[{i,line2[i]},{i,yMin+step*edge,yMax-step*edge,step}];

AppendTo[zShiftFits,{centerLine1,centerLine2,fitLine1,fitLine2}];
,
zShift=(Mean[centerLine1[[All,2]]]-Mean[centerLine2[[All,2]]]);

AppendTo[zShiftFits,{zShift,centerLine1,centerLine2}];
];

AppendTo[interpsOut,interp1];

,{i,2,Length[dataIn]}];

data2F[[All,3]]=data2F[[All,3]]+zShift;

If[OptionValue[SortInterp],
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,SortBy[data2F,Total[{#[[1]],#[[2]]}]&]],InterpolationOrder->1];
,
interp2=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data2F],InterpolationOrder->1];
];

AppendTo[interpsOut,interp2];

If[OptionValue[RetainEdges],
yCoords={Min[yMins],Max[yMaxs]};
,
yCoords={Max[yMins],Min[yMaxs]};
];

If[Positive[Mean[disps][[1]]],
centerNearest=Nearest[Table[i,{i,xStarts[[1]],Max[data2F[[All,1]]],step}]];
xCoords=Partition[{xStarts[[1]]}~Join~Flatten[Map[{First[centerNearest[#-step]],First[centerNearest[#]]}&,xCenters]]~Join~{Max[data2F[[All,1]]]},2];

dataOut=Table[{i,j,interpsOut[[1]][i,j]},
{j,yCoords[[2]], yCoords[[1]],-1.*step},
{i,xCoords[[1,1]],xCoords[[1,2]],step}];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(k/(Length[dataIn]))]];
ref @ SetPropertyValue[{"label", "text"}, "Compiling stitch "<>ToString[(k)]<>" of "<>ToString[(Length[dataIn])]];

dataOut=Join[dataOut,Table[{i,j,interpsOut[[k]][i,j]},
{j,yCoords[[2]], yCoords[[1]],-1.*step},
{i,xCoords[[k,1]],xCoords[[k,2]],step}],2];

,{k,2,Length[interpsOut]}];
,
centerNearest=Nearest[Table[i,{i,Min[data2F[[All,1]]],xStarts[[1]],step}]];
xCoords=Reverse[Partition[{xStarts[[1]]}~Join~Flatten[Map[{First[centerNearest[#+step]],First[centerNearest[#]]}&,xCenters]]~Join~{Min[data2F[[All,1]]]},2],2];

dataOut=Table[{i,j,interpsOut[[1]][i,j]},
{j,yCoords[[2]], yCoords[[1]],-1.*step},
{i,xCoords[[1,1]],xCoords[[1,2]],step}];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(k/(Length[dataIn]))]];
ref @ SetPropertyValue[{"label", "text"}, "Compiling stitch "<>ToString[(k)]<>" of "<>ToString[(Length[dataIn])]];

dataOut=Join[Table[{i,j,interpsOut[[k]][i,j]},
{j,yCoords[[2]], yCoords[[1]],-1.*step},
{i,xCoords[[k,1]],xCoords[[k,2]],step}],dataOut,2];

,{k,2,Length[interpsOut]}];
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

Share[];
{dataOut,zShiftFits}

];


Options[autoPlaneFit]={
GetPlane->False,PlaneIn->None,ShiftOrigin->True
};


autoPlaneFit[dataIn_,OptionsPattern[]]:=Module[
{data,planeFit,xo,yo,a,b,c,x,y,dataOut,vars,shift},

If[Length[Dimensions[dataIn]]==3,
data=Flatten[dataIn[[All,All,1;;3]],1];
,
data=dataIn[[All,1;;3]];
];

If[OptionValue[GetPlane],
{a,b,c}/.FindFit[data,a*x+b*y+c,{a,b,c},{x,y}]
,
If[ListQ[OptionValue[PlaneIn]],
If[Length[OptionValue[PlaneIn]]<3,
planeFit={OptionValue[PlaneIn][[1]],OptionValue[PlaneIn][[2]],c}/.FindFit[data,OptionValue[PlaneIn][[1]]*x+OptionValue[PlaneIn][[2]]*y+c,{c},{x,y}];
,
planeFit={OptionValue[PlaneIn][[1]],OptionValue[PlaneIn][[2]],OptionValue[PlaneIn][[3]]};
];
,
planeFit={a,b,c}/.FindFit[data,a*x+b*y+c,{a,b,c},{x,y}];
];

dataOut=dataIn;
If[Length[Dimensions[dataOut]]==3,
vars=dataOut[[All,All,1;;3]];
vars[[All,All,3]]=1;
dataOut[[All,All,3]]=dataOut[[All,All,3]]-vars.planeFit;
,
vars=dataOut[[All,1;;3]];
vars[[All,3]]=1;
dataOut[[All,3]]=dataOut[[All,3]]-vars.planeFit;
];

If[OptionValue[ShiftOrigin],
If[Length[Dimensions[dataIn]]==3,
xo=dataIn[[-1,1,1]];
yo=dataIn[[-1,1,2]];

shift=ConstantArray[{xo,yo},Dimensions[dataOut][[1;;2]]];
dataOut[[All,All,1;;2]]=dataOut[[All,All,1;;2]]-shift;
,
xo=SortBy[dataIn,Total][[1,1]];
yo=SortBy[dataIn,Total][[1,2]];

shift=ConstantArray[{xo,yo},Length[dataOut]];
dataOut[[All,1;;2]]=dataOut[[All,1;;2]]-shift;
];
dataOut
,
dataOut
]
]
];


(* ::Subsection::Closed:: *)
(*Stitch Post Processing*)


Options[zBarLocal]={
SquareFilter->True
};


zBarLocal[dataIn_,unitCell_,step_,OptionsPattern[]]:=Module[{
filterRadius,dataOut
},

If[OptionValue[SquareFilter],
filterRadius=Round[Mean[Map[#/2/step&,unitCell]]];
,
filterRadius=Round[Map[#/2/step&,unitCell]];
];

dataOut=dataIn;
dataOut[[All,All,3]]=MeanFilter[dataOut[[All,All,3]],filterRadius];

dataOut
];


getThickness[frontIn_,backIn_]:=Module[
{backFlat,backNearestXY},

backFlat=Flatten[backIn,1];

backNearestXY=Nearest[backFlat[[All,1;;2]]->Automatic];

Map[{#[[1]],#[[2]],#[[3]]-backFlat[[First[backNearestXY[{#[[1]],#[[2]]}]],3]]}&,frontIn,{2}]
];


getCurvatures[surfaceIn_]:=Module[
{surface,firstD,F1,n,secondD,F2,k},

surface={x,y,Normal[surfaceIn]};

firstD=Transpose[D[surface,{{x,y}}]];
F1=Simplify[{{#[[1]].#[[1]],#[[1]].#[[2]]},{#[[2]].#[[1]],#[[2]].#[[2]]}}&@firstD];

n=Simplify[Cross@@firstD/Norm[Cross@@firstD]];
secondD=Map[Transpose[#]&,D[firstD,{{x,y}}]];
F2=Simplify[{{#[[1,1]].n,#[[1,2]].n},{#[[2,1]].n,#[[2,2]].n}}&@secondD];

k=Eigenvalues[Inverse[F1].F2];

{k,Mean[k]}
];


(* ::Subsection::Closed:: *)
(*Thresholding*)


Options[arrayClusters]={
ThresholdVariable->4,VariablesOut->1;;3
};


arrayClusters[dataIn_,{thresholdMin_,thresholdMax_},OptionsPattern[]]:=Module[
{threshold,peakL,markers,peakPos,out},

threshold=Map[Boole[thresholdMin<#<thresholdMax]&,dataIn[[All,All,OptionValue[ThresholdVariable]]],{2}];

peakL=MorphologicalComponents[threshold];
markers=Union[Flatten[peakL]][[2;;-1]];
DistributeDefinitions[peakL,markers];
peakPos=ParallelMap[Position[peakL,#]&,markers];

DistributeDefinitions[peakPos,dataIn];
out=ParallelMap[Extract[dataIn,#]&,peakPos];

out[[All,All,OptionValue[VariablesOut]]]

];


Options[arrayClusterWindow]={
ThresholdVariable->4,VariablesOut->1;;3
};


arrayClusterWindow[dataIn_,{thresholdMin_,thresholdMax_},{center_,sizes_},OptionsPattern[]]:=Module[
{sizeX,sizeY,nearestX,nearestY,xRange,yRange},

If[ListQ[sizes],
sizeX=sizes[[1]];
sizeY=sizes[[2]];
,
sizeX=sizes;
sizeY=sizes;
];

nearestX=Nearest[dataIn[[Round[Dimensions[dataIn][[1]]/2.],All,1]]->Automatic];
nearestY=Nearest[dataIn[[All,Round[Dimensions[dataIn][[2]]/2.],2]]->Automatic];

xRange=Sort[Map[First[nearestX[#]]&,{center[[1]]-sizeX/2.,center[[1]]+sizeX/2.}]];
yRange=Sort[Map[First[nearestY[#]]&,{center[[2]]-sizeY/2.,center[[2]]+sizeY/2.}]];

arrayClusters[dataIn[[yRange[[1]];;yRange[[2]],xRange[[1]];;xRange[[2]]]],{thresholdMin,thresholdMax},ThresholdVariable->OptionValue[ThresholdVariable],VariablesOut->OptionValue[VariablesOut]]

];


imageClusters[dataIn_]:=Module[
{threshold,peakL,markers},

If[ImageQ[dataIn],
threshold=ImageData[dataIn];
,
threshold=dataIn;
];

peakL=MorphologicalComponents[threshold];
markers=Union[Flatten[peakL]][[2;;-1]];

DistributeDefinitions[peakL,markers];
Reverse[ParallelMap[Position[peakL,#]&,markers],{3}]


];


Options[thresholdWindow]={
ThresholdVariable->4
};


thresholdWindow[dataIn_,{thresholdMin_,thresholdMax_},{center_,sizes_},OptionsPattern[]]:=Module[
{data,sizeX,sizeY},

If[Length[Dimensions[dataIn]]==3,
data=Flatten[dataIn,1];
,
data=dataIn;
];

If[ListQ[sizes],
sizeX=sizes[[1]];
sizeY=sizes[[2]];
,
sizeX=sizes;
sizeY=sizes;
];

Select[data,(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax&&(center[[1]]-sizeX/2.)<#[[1]]<(center[[1]]+sizeX/2.)&&(center[[2]]-sizeY/2.)<#[[2]]<(center[[2]]+sizeY/2.))&][[All,{1,2,3}]]

];


Options[getClusters]={
Iterations->5
};


getClusters[dataIn_,cutOff_,OptionsPattern[]]:=Module[
{distanceMatrix,adjacencyMatrix},

distanceMatrix=DistanceMatrix[dataIn,DistanceFunction->EuclideanDistance];

If[OptionValue[Iterations]>1,
adjacencyMatrix=SparseArray[Map[Boole[#<=cutOff]&,distanceMatrix,{2}]];

Do[
adjacencyMatrix=adjacencyMatrix.adjacencyMatrix;
adjacencyMatrix=Map[Boole[#!=0]&,adjacencyMatrix,{2}];
Share[adjacencyMatrix];
,
{i,OptionValue[Iterations]}];
,
adjacencyMatrix=Map[Boole[#<=cutOff]&,distanceMatrix,{2}];
];

adjacencyMatrix=Normal[adjacencyMatrix];
SetSharedVariable[dataIn];
SetSharedVariable[adjacencyMatrix];

Share[];
ParallelMap[dataIn[[#]]&,ParallelMap[Flatten[Position[adjacencyMatrix,#]]&,Gather[adjacencyMatrix][[All,1]]]]

];


Options[gatherClusters]={
Iterations->5,OverlapFraction->4,xWindows->25,yWindows->25,ThresholdVariable->4
};


gatherClusters[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},OptionsPattern[]]:=Module[
{ref,frame=1,yMax,xMax,yRange,xRange,yStep,xStep,yPairs,xPairs,clusters,k,add,clusters2,test},

yMax=Dimensions[dataIn][[1]];
xMax=Dimensions[dataIn][[2]];

yRange=Round[Range[1,yMax,yMax/OptionValue[yWindows]]];
xRange=Round[Range[1,xMax,xMax/OptionValue[xWindows]]];

yStep=Round[(yRange[[2]]-yRange[[1]])/OptionValue[OverlapFraction]];
xStep=Round[(xRange[[2]]-xRange[[1]])/OptionValue[OverlapFraction]];

yPairs=Table[{yRange[[i]];;yRange[[i+1]]+yStep},{i,Length[yRange]-1}]~Join~{{yRange[[-1]];;yMax}};
xPairs=Table[{xRange[[i]];;xRange[[i+1]]+xStep},{i,Length[xRange]-1}]~Join~{{xRange[[-1]];;xMax}};

clusters=getClusters[Select[Flatten[dataIn[[First[yPairs[[1]]],First[xPairs[[1]]]]],1],(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax)&][[All,{1,2,3}]],cutOff,Iterations->OptionValue[Iterations]];
ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(frame/(Length[yPairs]*Length[xPairs]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing window "<>ToString[frame]<>" of "<>ToString[(Length[yPairs]*Length[xPairs])]];
frame++;

clusters2=getClusters[Select[Flatten[dataIn[[First[yPairs[[k]]],First[xPairs[[j]]]]],1],(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax)&][[All,{1,2,3}]],cutOff,Iterations->OptionValue[Iterations]];

add={};

Do[
test=Map[Boole[Length[Intersection[#,clusters2[[i]]]]!=0]&,clusters];
If[MemberQ[test,1],
If[Length[Position[test,1]]>1,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[Flatten[clusters[[Flatten[Position[test,1]]]],1]~Join~clusters2[[i]]];
clusters=Delete[clusters,Position[test,1][[2;;-1]]];
,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[clusters[[First[Flatten[Position[test,1]]]]]~Join~clusters2[[i]]];
];
,
AppendTo[add,clusters2[[i]]];
];
Share[];
,{i,1,Length[clusters2]}
];

clusters=clusters~Join~add;

Share[];
,{k,1,Length[yPairs]},{j,1,Length[xPairs]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];
Share[];

clusters
];


Options[gatherListClusters]={
Iterations->5,OverlapFraction->4,xWindows->25,yWindows->25,ThresholdVariable->4
};


gatherListClusters[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},OptionsPattern[]]:=Module[
{ref,frame=1,yRange,xRange,yStep,xStep,yPairs,xPairs,clusters,k,add,clusters2,test},

xRange={Min[dataIn[[All,1]]],Max[dataIn[[All,1]]]};
yRange={Min[dataIn[[All,2]]],Max[dataIn[[All,2]]]};

xRange=Range[xRange[[1]],xRange[[2]],First[Differences[xRange]]/OptionValue[xWindows]];
yRange=Range[yRange[[1]],yRange[[2]],First[Differences[yRange]]/OptionValue[yWindows]];

xStep=(xRange[[2]]-xRange[[1]])/OptionValue[OverlapFraction];
yStep=(yRange[[2]]-yRange[[1]])/OptionValue[OverlapFraction];

xPairs=Table[{xRange[[i]],xRange[[i+1]]+xStep},{i,Length[xRange]-1}];
yPairs=Table[{yRange[[i]],yRange[[i+1]]+yStep},{i,Length[yRange]-1}];

clusters=getClusters[Select[dataIn,(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax&&xPairs[[1,1]]<#[[1]]<xPairs[[1,2]]&&yPairs[[1,1]]<#[[2]]<yPairs[[1,2]])&][[All,{1,2,3}]],cutOff,Iterations->OptionValue[Iterations]];
ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(frame/(Length[yPairs]*Length[xPairs]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing window "<>ToString[frame]<>" of "<>ToString[(Length[yPairs]*Length[xPairs])]];
frame++;

clusters2=getClusters[Select[dataIn,(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax&&xPairs[[j,1]]<#[[1]]<xPairs[[j,2]]&&yPairs[[k,1]]<#[[2]]<yPairs[[k,2]])&][[All,{1,2,3}]],cutOff,Iterations->OptionValue[Iterations]];

add={};

Do[
test=Map[Boole[Length[Intersection[#,clusters2[[i]]]]!=0]&,clusters];
If[MemberQ[test,1],
If[Length[Position[test,1]]>1,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[Flatten[clusters[[Flatten[Position[test,1]]]],1]~Join~clusters2[[i]]];
clusters=Delete[clusters,Position[test,1][[2;;-1]]];
,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[clusters[[First[Flatten[Position[test,1]]]]]~Join~clusters2[[i]]];
];
,
AppendTo[add,clusters2[[i]]];
];
Share[];
,{i,1,Length[clusters2]}
];

clusters=clusters~Join~add;

Share[];

,{k,1,Length[yPairs]},{j,1,Length[xPairs]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];
Share[];
clusters
];


Options[imageListClusters]={
Iterations->5,OverlapFraction->4,xWindows->25,yWindows->25
};


imageListClusters[dataIn_,cutOff_,OptionsPattern[]]:=Module[
{ref,frame=1,yRange,xRange,yStep,xStep,yPairs,xPairs,clusters,k,add,clusters2,test},

xRange={Min[dataIn[[All,1]]],Max[dataIn[[All,1]]]};
yRange={Min[dataIn[[All,2]]],Max[dataIn[[All,2]]]};

xRange=Range[xRange[[1]],xRange[[2]],First[Differences[xRange]]/OptionValue[xWindows]];
yRange=Range[yRange[[1]],yRange[[2]],First[Differences[yRange]]/OptionValue[yWindows]];

xStep=(xRange[[2]]-xRange[[1]])/OptionValue[OverlapFraction];
yStep=(yRange[[2]]-yRange[[1]])/OptionValue[OverlapFraction];

xPairs=Table[{xRange[[i]],xRange[[i+1]]+xStep},{i,Length[xRange]-1}];
yPairs=Table[{yRange[[i]],yRange[[i+1]]+yStep},{i,Length[yRange]-1}];

clusters=getClusters[Select[dataIn,(xPairs[[1,1]]<#[[1]]<xPairs[[1,2]]&&yPairs[[1,1]]<#[[2]]<yPairs[[1,2]])&][[All,{1,2}]],cutOff,Iterations->OptionValue[Iterations]];
ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(frame/(Length[yPairs]*Length[xPairs]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing window "<>ToString[frame]<>" of "<>ToString[(Length[yPairs]*Length[xPairs])]];
frame++;

clusters2=getClusters[Select[dataIn,(xPairs[[j,1]]<#[[1]]<xPairs[[j,2]]&&yPairs[[k,1]]<#[[2]]<yPairs[[k,2]])&][[All,{1,2}]],cutOff,Iterations->OptionValue[Iterations]];

add={};

Do[
test=Map[Boole[Length[Intersection[#,clusters2[[i]]]]!=0]&,clusters];
If[MemberQ[test,1],
If[Length[Position[test,1]]>1,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[Flatten[clusters[[Flatten[Position[test,1]]]],1]~Join~clusters2[[i]]];
clusters=Delete[clusters,Position[test,1][[2;;-1]]];
,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[clusters[[First[Flatten[Position[test,1]]]]]~Join~clusters2[[i]]];
];
,
AppendTo[add,clusters2[[i]]];
];
Share[];
,{i,1,Length[clusters2]}
];

clusters=clusters~Join~add;

Share[];

,{k,1,Length[yPairs]},{j,1,Length[xPairs]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];
Share[];
clusters
];


Options[clusterWindow]={
Iterations->5,ThresholdVariable->4
};


clusterWindow[dataIn_,cutOff_,{thresholdMin_,thresholdMax_},{center_,size_},OptionsPattern[]]:=Module[
{data},

If[Length[Dimensions[dataIn]]==3,
data=Flatten[dataIn,1];
,
data=dataIn;
];

getClusters[Select[data,(thresholdMin<#[[OptionValue[ThresholdVariable]]]<thresholdMax&&(center[[1]]-size/2.)<#[[1]]<(center[[1]]+size/2.)&&(center[[2]]-size/2.)<#[[2]]<(center[[2]]+size/2.))&][[All,{1,2,3}]],cutOff,Iterations->OptionValue[Iterations]]

];


segregateTows[dataIn_,{stDevMin_,stDevMax_}]:=Module[
{},

Select[Map[{#,Abs[StandardDeviation[#[[All,1]]]/StandardDeviation[#[[All,2]]]]}&,dataIn],(stDevMin<#[[2]]<stDevMax)&][[All,1]]

];


combineClusters[dataIn1_,dataIn2_]:=Module[
{clusters,add,test},

clusters=dataIn1;
add={};

SetSharedVariable[clusters];
SetSharedVariable[add];
DistributeDefinitions[dataIn2];

ParallelDo[
test=Map[Boole[Length[Intersection[#,dataIn2[[i]]]]!=0]&,clusters];
If[MemberQ[test,1],
If[Length[Position[test,1]]>1,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[Flatten[clusters[[Flatten[Position[test,1]]]],1]~Join~dataIn2[[i]]];
clusters=Delete[clusters,Position[test,1][[2;;-1]]];
,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[clusters[[First[Flatten[Position[test,1]]]]]~Join~dataIn2[[i]]];
];
,
AppendTo[add,dataIn2[[i]]];
];
,{i,1,Length[dataIn2]}
];

clusters~Join~add

];


(* ::Subsection::Closed:: *)
(*Tow Analysis*)


rSquared[fitIn_]:=Module[
{ssERR,ssTOT},

ssERR=Total[Map[#^2&,fitIn["FitResiduals"]]];
ssTOT=Total[Map[(#[[3]]-Mean[fitIn["Data"][[All,3]]])^2&,fitIn["Data"]]];

1-(ssERR/ssTOT)
];


Options[findMax]={
GetRSquared->False
};


findMax[dataIn_,model_,parameters_,variables_,OptionsPattern[]]:=Module[
{ref,add,fit,max},

DistributeDefinitions[dataIn,model,parameters,variables,OptionValue[GetRSquared]];

ParallelTable[
fit=NonlinearModelFit[Sort[dataIn[[i]]],model,parameters,variables];
max=Quiet[FindMaximum[{fit["BestFit"],Max[dataIn[[i,All,1]]]>x>Min[dataIn[[i,All,1]]]&&Max[dataIn[[i,All,2]]]>y>Min[dataIn[[i,All,2]]]},{{x,Mean[dataIn[[i,All,1]]]},{y,Mean[dataIn[[i,All,2]]]}}]];
If[OptionValue[GetRSquared],
{{x,y,max[[1]]},{fit["RSquared"],rSquared[fit]}}/.max[[2]]
,
{x,y,max[[1]]}/.max[[2]]
]
,{i,Length[dataIn]}]
];


Options[findMin]={
GetRSquared->False
};


findMin[dataIn_,model_,parameters_,variables_,OptionsPattern[]]:=Module[
{ref,fit,min},

DistributeDefinitions[dataIn,model,parameters,variables,OptionValue[GetRSquared]];

ParallelTable[
fit=NonlinearModelFit[Sort[dataIn[[i]]],model,parameters,variables];
min=Quiet[FindMinimum[{fit["BestFit"],Max[dataIn[[i,All,1]]]>x>Min[dataIn[[i,All,1]]]&&Max[dataIn[[i,All,2]]]>y>Min[dataIn[[i,All,2]]]},{{x,Mean[dataIn[[i,All,1]]]},{y,Mean[dataIn[[i,All,2]]]}}]];
If[OptionValue[GetRSquared],
{{x,y,min[[1]]},{fit["RSquared"],rSquared[fit]}}/.min[[2]]
,
{x,y,min[[1]]}/.min[[2]]
]
,{i,Length[dataIn]}]
];


Options[findCrest]={
GetRSquared->False,CrestDirection->"y",Interpolate->False
};


findCrest[dataIn_,OptionsPattern[]]:=Module[
{data,model,parameters,variables,crests,fit,yValues,max,xValues,step,xRange,yRange,maxes,dataInterp},

If[StringMatchQ[OptionValue[CrestDirection],"y*",IgnoreCase->True],

If[Length[dataIn]>1,
data=dataIn[[1]];
model=dataIn[[2]];
parameters=dataIn[[3]];
variables=dataIn[[4]];

crests=Table[
fit=NonlinearModelFit[data[[i]],model,parameters,variables];
yValues=DeleteDuplicates[Sort[data[[i,All,2]]]];

If[OptionValue[GetRSquared],
DistributeDefinitions[fit,yValues];
{ParallelTable[
max=Quiet[FindMaximum[{fit["BestFit"]/.y->yValues[[i]],x>0},{x}]];

{x,yValues[[i]],max[[1]]}/.max[[2]]

,{i,Length[yValues]}],{fit["RSquared"],rSquared[fit]}}
,
DistributeDefinitions[fit,yValues];
ParallelTable[
max=Quiet[FindMaximum[{fit["BestFit"]/.y->yValues[[i]],x>0},{x}]];

{x,yValues[[i]],max[[1]]}/.max[[2]]

,{i,Length[yValues]}]
]
,{i,Length[data]}];

,
data=dataIn[[1]];
DistributeDefinitions[data];

If[OptionValue[Interpolate],

crests=ParallelTable[

step=Round[EuclideanDistance[data[[i,1]],data[[i,2]]],.001];
yValues=Range[Min[data[[i,All,2]]],Max[data[[i,All,2]]],step];
xRange={Min[data[[i,All,1]]],Max[data[[i,All,1]]]};

dataInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data[[i]]],InterpolationOrder->1];

maxes=Map[Quiet[Maximize[{dataInterp[x,#],xRange[[1]]<=x<=xRange[[2]]},x]]&,yValues];
MapThread[{x,#2,#1[[1]]}/.#1[[2]]&,{maxes,yValues}]

,{i,Length[data]}];

,

crests=ParallelTable[
Map[Extract[#,First[Position[#[[All,3]],Max[#[[All,3]]]]]]&,GatherBy[data[[i]],#[[2]]&]]
,{i,Length[data]}];
];
];
,

If[Length[dataIn]>1,
data=dataIn[[1]];
model=dataIn[[2]];
parameters=dataIn[[3]];
variables=dataIn[[4]];

crests=Table[
fit=NonlinearModelFit[data[[i]],model,parameters,variables];
xValues=DeleteDuplicates[Sort[data[[i,All,1]]]];

If[OptionValue[GetRSquared],
DistributeDefinitions[fit,xValues];
{ParallelTable[
max=Quiet[FindMaximum[{fit["BestFit"]/.x->xValues[[i]],y>0},{y}]];

{xValues[[i]],y,max[[1]]}/.max[[2]]

,{i,Length[xValues]}],{fit["RSquared"],rSquared[fit]}}
,
DistributeDefinitions[fit,xValues];
ParallelTable[
max=Quiet[FindMaximum[{fit["BestFit"]/.x->xValues[[i]],x>0},{x}]];

{xValues[[i]],y,max[[1]]}/.max[[2]]

,{i,Length[xValues]}]
]
,{i,Length[data]}];
,

data=dataIn[[1]];
DistributeDefinitions[data];

If[OptionValue[Interpolate],

crests=ParallelTable[

step=Round[EuclideanDistance[data[[i,1]],data[[i,2]]],.001];
xValues=Range[Min[data[[i,All,1]]],Max[data[[i,All,1]]],step];
yRange={Min[data[[i,All,2]]],Max[data[[i,All,2]]]};

dataInterp=Interpolation[Map[{{#[[1]],#[[2]]},#[[3]]}&,data[[i]]],InterpolationOrder->1];

maxes=Map[Quiet[Maximize[{dataInterp[#,y],yRange[[1]]<=y<=yRange[[2]]},y]]&,xValues];
MapThread[{#2,y,#1[[1]]}/.#1[[2]]&,{maxes,xValues}]

,{i,Length[data]}];

,

crests=ParallelTable[
Map[Extract[#,First[Position[#[[All,3]],Max[#[[All,3]]]]]]&,GatherBy[data[[i]],#[[1]]&]]
,{i,Length[data]}];
];
];
];

crests

];


Options[findTrough]={
GetRSquared->False,CrestDirection->"y"
};


findTrough[dataIn_,OptionsPattern[]]:=Module[
{data,model,parameters,variables,crests,fit,yValues,min,xValues},

If[StringMatchQ[OptionValue[CrestDirection],"y*",IgnoreCase->True],

If[Length[dataIn]>1,
data=dataIn[[1]];
model=dataIn[[2]];
parameters=dataIn[[3]];
variables=dataIn[[4]];

crests=Table[
fit=NonlinearModelFit[data[[i]],model,parameters,variables];
yValues=DeleteDuplicates[Sort[data[[i,All,2]]]];

If[OptionValue[GetRSquared],
DistributeDefinitions[fit,yValues];
{ParallelTable[
min=Quiet[FindMinimum[{fit["BestFit"]/.y->yValues[[i]],x>0},{x}]];

{x,yValues[[i]],min[[1]]}/.min[[2]]

,{i,Length[yValues]}],{fit["RSquared"],rSquared[fit]}}
,
DistributeDefinitions[fit,yValues];
ParallelTable[
min=Quiet[FindMinimum[{fit["BestFit"]/.y->yValues[[i]],x>0},{x}]];

{x,yValues[[i]],min[[1]]}/.min[[2]]

,{i,Length[yValues]}]
]
,{i,Length[data]}];
,
data=dataIn[[1]];
DistributeDefinitions[data];

crests=ParallelTable[
Map[Extract[#,First[Position[#[[All,3]],Min[#[[All,3]]]]]]&,GatherBy[data[[i]],#[[2]]&]]
,{i,Length[data]}];
];
,

If[Length[dataIn]>1,
data=dataIn[[1]];
model=dataIn[[2]];
parameters=dataIn[[3]];
variables=dataIn[[4]];

crests=Table[
fit=NonlinearModelFit[data[[i]],model,parameters,variables];
xValues=DeleteDuplicates[Sort[data[[i,All,1]]]];

If[OptionValue[GetRSquared],
DistributeDefinitions[fit,xValues];
{ParallelTable[
min=Quiet[FindMinimum[{fit["BestFit"]/.x->xValues[[i]],y>0},{y}]];

{xValues[[i]],y,min[[1]]}/.min[[2]]

,{i,Length[xValues]}],{fit["RSquared"],rSquared[fit]}}
,
DistributeDefinitions[fit,xValues];
ParallelTable[
min=Quiet[FindMinimum[{fit["BestFit"]/.x->xValues[[i]],x>0},{x}]];

{xValues[[i]],y,min[[1]]}/.min[[2]]

,{i,Length[xValues]}]
]
,{i,Length[data]}];
,
data=dataIn[[1]];
DistributeDefinitions[data];

crests=ParallelTable[
Map[Extract[#,First[Position[#[[All,3]],Min[#[[All,3]]]]]]&,GatherBy[data[[i]],#[[1]]&]]
,{i,Length[data]}];
];
];

crests

];


(* ::Subsection::Closed:: *)
(*Rotation and Translation*)


Options[flipData]={
FlipAbout->"y",Invert->False
};


flipData[dataIn_,OptionsPattern[]]:=Module[{
dataInRev,Start
},

If[StringMatchQ[OptionValue[FlipAbout],"y*",IgnoreCase->True],
dataInRev=Reverse[dataIn,2];
Start=dataInRev[[1,1,1]];
dataInRev[[All,All,1]]=Start-dataInRev[[All,All,1]];
,
dataInRev=Reverse[dataIn];
Start=dataInRev[[-1,-1,2]];
dataInRev[[All,All,2]]=Start-dataInRev[[All,All,2]];
];

If[OptionValue[Invert],
dataInRev[[All,All,3]]=(-1.)*dataInRev[[All,All,3]];
dataInRev
,
dataInRev]
];


Options[rotateData]={
RotationCenter->None
};


rotateData[dataIn_,rot_,OptionsPattern[]]:=Module[{
dataOut
},
dataOut=dataIn;

If[ListQ[OptionValue[RotationCenter]],
If[SameQ[Length[Dimensions[dataOut]],2],
dataOut[[All,1]]=dataOut[[All,1]]-OptionValue[RotationCenter][[1]];
dataOut[[All,2]]=dataOut[[All,2]]-OptionValue[RotationCenter][[2]];
dataOut[[All,3]]=dataOut[[All,3]]-OptionValue[RotationCenter][[3]];

dataOut[[All,1;;3]]=Transpose[rot.Transpose[dataOut[[All,1;;3]]]];

dataOut[[All,1]]=dataOut[[All,1]]+OptionValue[RotationCenter][[1]];
dataOut[[All,2]]=dataOut[[All,2]]+OptionValue[RotationCenter][[2]];
dataOut[[All,3]]=dataOut[[All,3]]+OptionValue[RotationCenter][[3]];
,
dataOut[[All,All,1]]=dataOut[[All,All,1]]-OptionValue[RotationCenter][[1]];
dataOut[[All,All,2]]=dataOut[[All,All,2]]-OptionValue[RotationCenter][[2]];
dataOut[[All,All,3]]=dataOut[[All,All,3]]-OptionValue[RotationCenter][[3]];

dataOut[[All,All,1;;3]]=Partition[Transpose[rot.Transpose[Flatten[dataOut[[All,All,1;;3]],1]]],Dimensions[dataOut][[2]]];

dataOut[[All,All,1]]=dataOut[[All,All,1]]+OptionValue[RotationCenter][[1]];
dataOut[[All,All,2]]=dataOut[[All,All,2]]+OptionValue[RotationCenter][[2]];
dataOut[[All,All,3]]=dataOut[[All,All,3]]+OptionValue[RotationCenter][[3]];
];

,

If[SameQ[Length[Dimensions[dataOut]],2],
dataOut[[All,1;;3]]=Transpose[rot.Transpose[dataOut[[All,1;;3]]]];
,
dataOut[[All,All,1;;3]]=Partition[Transpose[rot.Transpose[Flatten[dataOut[[All,All,1;;3]],1]]],Dimensions[dataOut][[2]]];
];
];

dataOut

];


shiftData[dataIn_,shift_]:=Module[
{dataOut},

dataOut=dataIn;

If[Length[Dimensions[dataOut]]==2,
dataOut[[All,1]]=dataOut[[All,1]]+shift[[1]];
dataOut[[All,2]]=dataOut[[All,2]]+shift[[2]];
dataOut[[All,3]]=dataOut[[All,3]]+shift[[3]];
,
dataOut[[All,All,1]]=dataOut[[All,All,1]]+shift[[1]];
dataOut[[All,All,2]]=dataOut[[All,All,2]]+shift[[2]];
dataOut[[All,All,3]]=dataOut[[All,All,3]]+shift[[3]];
];

dataOut

];


Options[rotateShiftData]={
RotationCenter->None
};


rotateShiftData[dataIn_,rot_,shift_]:=Module[
{dataOut},

shiftData[rotateData[dataIn,rot,RotationCenter->OptionValue[RotationCenter]],shift]

];


deltaH[dataIn_,vars_]:=Module[
{front, back, hmeasured,rot,hdata},
front=dataIn[[1,All,3]];
back=dataIn[[2]];
hmeasured=dataIn[[3]];

rot=RotationMatrix[vars[[2]]Degree,{0,1,0}].RotationMatrix[vars[[3]]Degree,{1,0,0}];
back=Transpose[rot.Transpose[back]];
back=back[[All,3]]+vars[[1]];
hdata=MapThread[{#1-#2}&,{front,back}];

Sum[(hmeasured[[i]]-hdata[[i]])^2,{i,Length[hdata]}]

];


(* ::Subsection::Closed:: *)
(*Autocorrelation and Variances*)


Options[spatialAC]={
range->Automatic,BothSides->False,GetCounts->False,CutOffR->None,Iterations->1,OverlapFraction->10,Windows->1
};


spatialAC[dataIn_,bin_,OptionsPattern[]]:=Module[
{deltaXY,rangeX,rangeY,countXY,peakL,peakPos,peakN,centersX,centersY,cutOff,countPos},

deltaXY=ConstantArray[dataIn[[All,1;;2]],Length[dataIn[[All,1;;2]]]]-Transpose[ConstantArray[dataIn[[All,1;;2]],Length[dataIn[[All,1;;2]]]]];

If[OptionValue[BothSides],
deltaXY=ReplacePart[deltaXY,Position[Map[Boole[#==0]&,ConstantArray[dataIn[[All,4]],Length[dataIn[[All,4]]]]-Transpose[ConstantArray[dataIn[[All,4]],Length[dataIn[[All,4]]]]],{2}],0]->{Infinity,Infinity}];
,
deltaXY=deltaXY;
];

If[SameQ[OptionValue[range],Automatic],
rangeX={Round[-Median[dataIn[[All,1]]],bin]-bin/2.,Round[Median[dataIn[[All,1]]],bin]+bin/2.,bin};
rangeY={Round[-Median[dataIn[[All,2]]],bin]-bin/2.,Round[Median[dataIn[[All,2]]],bin]+bin/2.,bin};
,
rangeX={Round[OptionValue[range][[1,1]],bin]-bin/2.,Round[OptionValue[range][[1,2]],bin]+bin/2.,bin};
rangeY={Round[OptionValue[range][[2,1]],bin]-bin/2.,Round[OptionValue[range][[2,2]],bin]+bin/2.,bin};
];

centersX=Chop[Range[rangeX[[1]]+bin/2.,rangeX[[2]]-bin/2.,bin]];
centersY=Chop[Range[rangeY[[1]]+bin/2.,rangeY[[2]]-bin/2.,bin]];

countXY=BinCounts[Flatten[deltaXY,1],rangeX,rangeY];

If[OptionValue[GetCounts],
MapThread[{#1,#2,#3}&,{Transpose[ConstantArray[centersX,Length[centersY]]],ConstantArray[centersY,Length[centersX]],countXY},2]
,
If[NumberQ[OptionValue[CutOffR]],
cutOff=OptionValue[CutOffR]/bin;

countPos=Select[Flatten[MapIndexed[#2~Join~{#1}&,countXY,{2}],1],(#[[3]]!= 0.)&][[All,1;;2]];
If[SameQ[OptionValue[Windows],1],
peakL=ReplacePart[ConstantArray[0,Dimensions[countXY]],MapIndexed[#1->First[#2]&,getClusters[countPos,cutOff,Iterations->OptionValue[Iterations]]]];
,
peakL=ReplacePart[ConstantArray[0,Dimensions[countXY]],MapIndexed[#1->First[#2]&,gatherACclusters[countPos,cutOff,Iterations->OptionValue[Iterations],OverlapFraction->OptionValue[OverlapFraction],Windows->OptionValue[Windows]]]];
];
,
peakL=MorphologicalComponents[countXY];
];

DistributeDefinitions[peakL];
peakPos=ParallelMap[Position[peakL,#]&,Union[Flatten[peakL]]];
peakN=Map[#->Total[Extract[countXY,#]]&,peakPos];
peakN[[1,2]]=1;
peakN=ReplacePart[peakL,peakN];
countXY=N[countXY/peakN];

MapThread[{#1,#2,#3}&,{Transpose[ConstantArray[centersX,Length[centersY]]],ConstantArray[centersY,Length[centersX]],countXY},2]
]
];


Options[posVariances]={
range->Automatic,BothSides->False,CutOffR->None,Iterations->1,OverlapFraction->10,Windows->1,PrincipalVariances->False
};


posVariances[dataIn_,bin_,OptionsPattern[]]:=Module[
{deltaXY,rangeX,rangeY,listXY,countXY,cutOff,countPos,distanceMatrix,adjacencyMatrix,peakL,peakPos},

deltaXY=ConstantArray[dataIn[[All,1;;2]],Length[dataIn[[All,1;;2]]]]-Transpose[ConstantArray[dataIn[[All,1;;2]],Length[dataIn[[All,1;;2]]]]];

If[OptionValue[BothSides],
deltaXY=ReplacePart[deltaXY,Position[Map[Boole[#==0]&,ConstantArray[dataIn[[All,4]],Length[dataIn[[All,4]]]]-Transpose[ConstantArray[dataIn[[All,4]],Length[dataIn[[All,4]]]]],{2}],0]->{Infinity,Infinity}];
,
deltaXY=deltaXY;
];

If[SameQ[OptionValue[range],Automatic],
rangeX={Round[-Median[dataIn[[All,1]]],bin]-bin/2.,Round[Median[dataIn[[All,1]]],bin]+bin/2.,bin};
rangeY={Round[-Median[dataIn[[All,2]]],bin]-bin/2.,Round[Median[dataIn[[All,2]]],bin]+bin/2.,bin};
,
rangeX={Round[OptionValue[range][[1,1]],bin]-bin/2.,Round[OptionValue[range][[1,2]],bin]+bin/2.,bin};
rangeY={Round[OptionValue[range][[2,1]],bin]-bin/2.,Round[OptionValue[range][[2,2]],bin]+bin/2.,bin};
];

listXY=BinLists[Flatten[deltaXY,1],rangeX,rangeY];
countXY=Map[Length[#]&,listXY,{2}];

If[NumberQ[OptionValue[CutOffR]],
cutOff=OptionValue[CutOffR]/bin;

countPos=Select[Flatten[MapIndexed[#2~Join~{#1}&,countXY,{2}],1],(#[[3]]!= 0.)&][[All,1;;2]];
If[SameQ[OptionValue[Windows],1],
peakL=ReplacePart[ConstantArray[0,Dimensions[countXY]],MapIndexed[#1->First[#2]&,getClusters[countPos,cutOff,Iterations->OptionValue[Iterations]]]];
,
peakL=ReplacePart[ConstantArray[0,Dimensions[countXY]],MapIndexed[#1->First[#2]&,gatherACclusters[countPos,cutOff,Iterations->OptionValue[Iterations],OverlapFraction->OptionValue[OverlapFraction],Windows->OptionValue[Windows]]]];
];
,
peakL=MorphologicalComponents[countXY];
];

DistributeDefinitions[peakL];
peakPos=ParallelMap[Position[peakL,#]&,Union[Flatten[peakL]][[2;;-1]]];

If[OptionValue[PrincipalVariances],
Map[Join[Mean[Flatten[#,1]],Eigenvalues[Covariance[Flatten[#,1]]],{ArcTan[(Eigenvectors[Covariance[Flatten[#,1]]][[1,2]]/Eigenvectors[Covariance[Flatten[#,1]]][[1,1]])]/Degree},{Length[Flatten[#,1]]}]&,Select[Map[Extract[listXY,#]&,peakPos],(Dimensions[Flatten[#,1]][[1]]>1)&]]
,
Map[Join[Mean[Flatten[#,1]],Flatten[Covariance[Flatten[#,1]]][[{1,2,4}]],{Length[Flatten[#,1]]}]&,Select[Map[Extract[listXY,#]&,peakPos],(Dimensions[Flatten[#,1]][[1]]>1)&]]
]
];


Options[gatherACclusters]={
Iterations->1,OverlapFraction->10,Windows->16
};


gatherACclusters[dataIn_,cutOff_,OptionsPattern[]]:=Module[
{ref,frame=1,yRange,xRange,yStep,xStep,yPairs,xPairs,clusters,add,clusters2,test},

xRange={Min[dataIn[[All,1]]],Max[dataIn[[All,1]]]};
yRange={Min[dataIn[[All,2]]],Max[dataIn[[All,2]]]};

xRange=Range[xRange[[1]],xRange[[2]],First[Differences[xRange]]/OptionValue[Windows]];
yRange=Range[yRange[[1]],yRange[[2]],First[Differences[yRange]]/OptionValue[Windows]];

xStep=(xRange[[2]]-xRange[[1]])/OptionValue[OverlapFraction];
yStep=(yRange[[2]]-yRange[[1]])/OptionValue[OverlapFraction];

xPairs=Table[{xRange[[i]],xRange[[i+1]]+xStep},{i,Length[xRange]-1}];
yPairs=Table[{yRange[[i]],yRange[[i+1]]+yStep},{i,Length[yRange]-1}];

clusters=getClusters[Select[dataIn,(xPairs[[1,1]]<#[[1]]<xPairs[[1,2]]&&yPairs[[1,1]]<#[[2]]<yPairs[[1,2]])&],cutOff,Iterations->OptionValue[Iterations]];
ref=ProgressDialog[];

Do[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(frame/(Length[yPairs]*Length[xPairs]))]];
ref @ SetPropertyValue[{"label", "text"}, "Processing window "<>ToString[frame]<>" of "<>ToString[(Length[yPairs]*Length[xPairs])]];
frame++;

clusters2=getClusters[Select[dataIn,(xPairs[[j,1]]<#[[1]]<xPairs[[j,2]]&&yPairs[[k,1]]<#[[2]]<yPairs[[k,2]])&],cutOff,Iterations->OptionValue[Iterations]];

add={};

Do[
test=Map[Boole[Length[Intersection[#,clusters2[[i]]]]!=0]&,clusters];
If[MemberQ[test,1],
If[Length[Position[test,1]]>1,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[Flatten[clusters[[Flatten[Position[test,1]]]],1]~Join~clusters2[[i]]];
clusters=Delete[clusters,Position[test,1][[2;;-1]]];
,
clusters[[First[Flatten[Position[test,1]]]]]=DeleteDuplicates[clusters[[First[Flatten[Position[test,1]]]]]~Join~clusters2[[i]]];
];
,
AppendTo[add,clusters2[[i]]];
];
Share[];
,{i,1,Length[clusters2]}
];

clusters=clusters~Join~add;

Share[];

,{k,1,Length[yPairs]},{j,1,Length[xPairs]}
];

CloseGUIObject[ref];
ReleaseGUIObject[ref];
Share[];
clusters
];


(* ::Subsection::Closed:: *)
(*Angular Analysis*)


Options[angleMatrix]={
AngleFromXaxis->True,IgnoreQuadrant->False,CutOffDistance->None
};


angleMatrix[dataIn_,OptionsPattern[]]:=Module[
{deltaX,deltaY,dataOut,pos},

deltaX=ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]-Transpose[ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]];
deltaY=ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]-Transpose[ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]];

If[OptionValue[IgnoreQuadrant],
If[OptionValue[AngleFromXaxis],
dataOut=Quiet[N[ArcTan[deltaY/deltaX]/Degree]];
,
dataOut=Quiet[N[ArcTan[deltaX/deltaY]/Degree]];
];
,
If[OptionValue[AngleFromXaxis],
dataOut=Quiet[N[ArcTan[deltaX,deltaY]/Degree]];
,
dataOut=Quiet[N[ArcTan[deltaY,deltaX]/Degree]];
];
];

If[SameQ[First[DeleteDuplicates[Diagonal[dataOut]]],Indeterminate],
dataOut
,
dataOut=ReplacePart[dataOut,{i_,i_}->Indeterminate];
];

If[SameQ[OptionValue[CutOffDistance],None],
dataOut
,
If[NumberQ[OptionValue[CutOffDistance]],
pos=Position[Map[Boole[#>=OptionValue[CutOffDistance]]&,Sqrt[deltaX^2+deltaY^2],{2}],1];
ReplacePart[dataOut,pos->Indeterminate]
,
pos=Position[MapThread[Boole[(#1+#2)!=0]&,{Map[Boole[Abs[#]>=OptionValue[CutOffDistance][[1]]]&,deltaX,{2}],Map[Boole[Abs[#]>=OptionValue[CutOffDistance][[2]]]&,deltaY,{2}]},2],1];
ReplacePart[dataOut,pos->Indeterminate]
]
]

];


Options[angularHist]={
AngleFromXaxis->True,IgnoreQuadrant->False,CutOffDistance->None,UnitArea->False
};


angularHist[dataIn_,bin_,OptionsPattern[]]:=Module[
{dataOut,total},

dataOut=Tally[Sort[Round[DeleteCases[Flatten[angleMatrix[dataIn,AngleFromXaxis->OptionValue[AngleFromXaxis],IgnoreQuadrant->OptionValue[IgnoreQuadrant],CutOffDistance->OptionValue[CutOffDistance]]],Indeterminate],bin]]];

If[OptionValue[UnitArea],
total=Total[dataOut[[All,2]]];
dataOut[[All,2]]=dataOut[[All,2]]/(total*bin);
dataOut
,
dataOut
]

];


Options[peakLocation]={
AngleFromXaxis->True,IgnoreQuadrant->False,CutOffDistance->None,UnitArea->False,PeakValue->None,DistanceOnly->False
};


peakLocation[dataIn_,bin_,range_,theta_,OptionsPattern[]]:=Module[{
rot,dataRot,tallyTheta,model,a,\[Mu],\[Sigma],peak,normalFit,posPeak,negPeak,normalFitP,normalFitN},

rot=RotationMatrix[theta Degree];
dataRot=Map[rot.#&,dataIn[[All,1;;2]]];

tallyTheta=angularHist[dataRot,bin,AngleFromXaxis->OptionValue[AngleFromXaxis],IgnoreQuadrant->OptionValue[IgnoreQuadrant],CutOffDistance->OptionValue[CutOffDistance],UnitArea->OptionValue[UnitArea]];
model=(a*PDF[NormalDistribution[\[Mu],\[Sigma]],x]);

If[NumberQ[OptionValue[PeakValue]],

peak=Select[tallyTheta,(range[[1]]<#[[1]]<range[[2]])&];

normalFit=FindFit[peak,{model,{range[[1]]<\[Mu]<range[[2]]}},{a,\[Mu],\[Sigma]},x,Method->NMinimize];

If[OptionValue[DistanceOnly],
Abs[OptionValue[PeakValue]-\[Mu]/.normalFit]
,
{{\[Mu],OptionValue[PeakValue]-\[Mu]},model,peak}/.normalFit
]
,
posPeak=Select[tallyTheta,(range[[1]]<#[[1]]<range[[2]])&];
negPeak=Select[tallyTheta,(-range[[2]]<#[[1]]<-range[[1]])&];

normalFitP=FindFit[posPeak,{model,{range[[1]]<\[Mu]<range[[2]]}},{a,\[Mu],\[Sigma]},x,Method->NMinimize];
normalFitN=FindFit[negPeak,{model,{-range[[2]]<\[Mu]<-range[[1]]}},{a,\[Mu],\[Sigma]},x,Method->NMinimize];

If[OptionValue[DistanceOnly],
Abs[(\[Mu]/.normalFitP)-(\[Mu]/.normalFitN)]
,
{{{\[Mu]/.normalFitP,\[Mu]/.normalFitN},(\[Mu]/.normalFitP)-(\[Mu]/.normalFitN)},{model/.normalFitP,model/.normalFitN},{posPeak,negPeak}}
]
]

];


Options[peakFitting]={
AngleFromXaxis->True,PeakValue->None
};


peakFitting[dataIn_,bin_,range_,theta_,OptionsPattern[]]:=Module[{
rot,dataRot,tallyTheta,peak,interp,posInterp,negInterp,xp,xn,max,posPeak,negPeak},

rot=RotationMatrix[theta Degree];
dataRot=Map[rot.#&,dataIn[[All,1;;2]]];

tallyTheta=Tally[Sort[Round[Flatten[angleMatrix[dataRot,AngleFromXaxis->OptionValue[AngleFromXaxis]]],bin]]];

If[NumberQ[OptionValue[PeakValue]],

peak=Select[tallyTheta,(range[[1]]<#[[1]]<range[[2]])&];

interp=Interpolation[peak,InterpolationOrder->2, Method->"Spline"];

max=xp/.Flatten[Last[Maximize[{interp[xp],range[[1]]<xp<range[[2]]},xp]]];

{{max,Abs[OptionValue[PeakValue]-max]},interp,peak}
,
posPeak=Select[tallyTheta,(range[[1]]<#[[1]]<range[[2]])&];
negPeak=Select[tallyTheta,(-range[[2]]<#[[1]]<-range[[1]])&];

posInterp=Interpolation[posPeak,InterpolationOrder->2, Method->"Spline"];
negInterp=Interpolation[negPeak,InterpolationOrder->2, Method->"Spline"];

max={xp,xn}/.Flatten[{Last[Maximize[{posInterp[xp],range[[1]]<xp<range[[2]]},xp]],Last[Maximize[{negInterp[xn],-range[[2]]<xn<-range[[1]]},xn]]}];

{{max,max[[1]]+max[[2]]},{posInterp,negInterp},{posPeak,negPeak}}]

];


Options[peakInfo]={
Window->10
};


peakInfo[dataIn_,peaksIn_,OptionsPattern[]]:=Module[
{peaks,peakWindows,data,normalFit,model,a,\[Mu],\[Sigma]},

If[Length[peaksIn]==0,
peaks={peaksIn};
,
peaks=peaksIn;
];

If[SameQ[Length[OptionValue[Window]],Length[peaks]],
peakWindows=MapThread[{#1-#2/2.,#1+#2/2.}&,{peaks,OptionValue[Window]}];
,
peakWindows=Map[{#-OptionValue[Window]/2.,#+OptionValue[Window]/2.}&,peaks];
];
model=(a*PDF[NormalDistribution[\[Mu],\[Sigma]],x]);

Table[
data=Select[dataIn,(peakWindows[[i,1]]<#[[1]]<peakWindows[[i,2]])&];
normalFit=FindFit[data,{model,{peakWindows[[i,1]]<\[Mu]<peakWindows[[i,2]],0<a}},{a,\[Mu],\[Sigma]},x,Method->NMinimize];

{a,\[Mu],\[Sigma]}/.normalFit
,{i,Length[peakWindows]}]

];


getShearAngle[anglesIn_,unitCellAngles_]:=Module[
{nearestCrown,inputs},

nearestCrown=Nearest[unitCellAngles[[All,1]]->unitCellAngles[[All,2]]];
inputs=Map[{#[[1]],(#[[1]]+#[[2]]),First[nearestCrown[#[[1]]]]}&,anglesIn];

Cases[Quiet[Map[{#[[1]],N[ArcTan[(#[[3,1]]*Tan[#[[2]]Degree]-#[[3,2]])/#[[3,1]]]/Degree]}&,inputs]],{_Real,_Real}]

];


(* ::Subsection::Closed:: *)
(*Displacement and Strain Analysis*)


Options[dispFromIdeal]={
IdealPositions->False
};


dispFromIdeal[dataIn_,idealIn_,OptionsPattern[]]:=Module[
{nearestIdeal},

If[OptionValue[IdealPositions],
If[Length[Dimensions[dataIn]]==2,
nearestIdeal=Nearest[idealIn];
Map[First[nearestIdeal[#]]~Join~(#-First[nearestIdeal[#]])&,dataIn]
,
Flatten[MapThread[MapThread[#2~Join~(#1-#2)&,{#1,#2}]&,{dataIn,idealIn}],1]
]
,
If[Length[Dimensions[dataIn]]==2,
nearestIdeal=Nearest[idealIn];
Map[#~Join~(#-First[nearestIdeal[#]])&,dataIn]
,
Flatten[MapThread[MapThread[#1~Join~(#1-#2)&,{#1,#2}]&,{dataIn,idealIn}],1]
]
]

];


strainFilter[dispIn_,cutOffR_]:=Module[
{deltaX,deltaY,adj,adjx,adjy,dataF,dUdx,dUdy,dVdx,dVdy,b,mx,my,x,y},
deltaX=ConstantArray[dispIn[[All,1]],Length[dispIn[[All,1]]]]-Transpose[ConstantArray[dispIn[[All,1]],Length[dispIn[[All,1]]]]];
deltaY=ConstantArray[dispIn[[All,2]],Length[dispIn[[All,2]]]]-Transpose[ConstantArray[dispIn[[All,2]],Length[dispIn[[All,2]]]]];

If[NumberQ[cutOffR],
adj=Map[Boole[#<=cutOffR]&,Sqrt[deltaX^2+deltaY^2],{2}];
,
adjx=Map[Boole[Abs[#]<=cutOffR[[1]]]&,deltaX,{2}];
adjy=Map[Boole[Abs[#]<=cutOffR[[2]]]&,deltaY,{2}];
adj=MapThread[Boole[(#1+#2)==2]&,{adjx,adjy},2];
];

dataF=MapIndexed[{dispIn[[First[#2],1;;2]],dispIn[[#1]]}&,Map[Flatten[Position[#,1]]&,adj]];

Table[
{dUdx,dUdy}={mx,my}/.FindFit[dataF[[i,2,All,{1,2,3}]],b+mx*x+my*y,{b,mx,my},{x,y}];
{dVdx,dVdy}={mx,my}/.FindFit[dataF[[i,2,All,{1,2,4}]],b+mx*x+my*y,{b,mx,my},{x,y}];

dataF[[i,1]]~Join~{dUdx,dVdy,1/2*(dUdy+dVdx),dUdy,dVdx}
,{i,Length[dataF]}
]

];


(* ::Subsection::Closed:: *)
(*Neighbors*)


getNN[dataIn_,cutOffR_]:=Module[
{deltaX,deltaY,adj,adjx,adjy,dataF},

deltaX=ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]-Transpose[ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]];
deltaY=ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]-Transpose[ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]];

If[NumberQ[cutOffR],
adj=Map[Boole[#<=cutOffR]&,Sqrt[deltaX^2+deltaY^2],{2}];
,
adjx=Map[Boole[Abs[#]<=cutOffR[[1]]]&,deltaX,{2}];
adjy=Map[Boole[Abs[#]<=cutOffR[[2]]]&,deltaY,{2}];
adj=MapThread[Boole[(#1+#2)==2]&,{adjx,adjy},2];
];

MapIndexed[{dataIn[[First[#2]]],dataIn[[#1]]}&,Map[Flatten[Position[#,1]]&,adj]]

];


Options[getNeighbors]={
AngleFromXaxis->True,AllNeighbors->False
};


getNeighbors[dataIn_,angleRange_,distanceRange_,OptionsPattern[]]:=Module[
{deltaX,deltaY,angles,distances,neighbors,test},

If[Length[Dimensions[dataIn]]==2,
deltaX=ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]-Transpose[ConstantArray[dataIn[[All,1]],Length[dataIn[[All,1]]]]];
deltaY=ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]-Transpose[ConstantArray[dataIn[[All,2]],Length[dataIn[[All,2]]]]];
,
deltaX=ConstantArray[dataIn[[2,All,1]],Length[dataIn[[1,All,1]]]]-Transpose[ConstantArray[dataIn[[1,All,1]],Length[dataIn[[2,All,1]]]]];
deltaY=ConstantArray[dataIn[[2,All,2]],Length[dataIn[[1,All,2]]]]-Transpose[ConstantArray[dataIn[[1,All,2]],Length[dataIn[[2,All,2]]]]];];

If[OptionValue[AngleFromXaxis],
angles=Quiet[N[ArcTan[deltaY/deltaX]/Degree]];
,
angles=Quiet[N[ArcTan[deltaX/deltaY]/Degree]];
];

distances=Sqrt[deltaX^2+deltaY^2];

distances=ReplacePart[distances,Position[Map[Boole[angleRange[[1]]<#<angleRange[[2]]]&,angles,{2}],0]->Infinity];

If[OptionValue[AllNeighbors],
neighbors=Select[MapIndexed[{#2}~Join~Position[#1,1]&,Map[Boole[distanceRange[[1]]<#<distanceRange[[2]]]&,distances,{2}]],Length[#]==3&];
test=Length[neighbors];
While[Length[Gather[neighbors,Length[Intersection[#1,#2]]>0&]]<test,
neighbors=Map[DeleteDuplicates[Flatten[#,1]]&,Gather[neighbors,Length[Intersection[#1,#2]]>0&]];test=Length[neighbors]];
neighbors
,
Select[MapIndexed[{#2}~Join~Position[#1,1]&,Map[Boole[distanceRange[[1]]<#<distanceRange[[2]]]&,distances,{2}]],Length[#]==3&]
]
];


Options[neighborDistance]={
Distance->None
};


neighborDistance[neighborsIn_,vars_,OptionsPattern[]]:=Module[
{front,back,rot},

front=neighborsIn[[All,1]];
back=neighborsIn[[All,2]];

rot=RotationMatrix[vars[[3]]Degree];

back=Transpose[rot.Transpose[back]];

back[[All,1]]=back[[All,1]]+vars[[1]];
back[[All,2]]=back[[All,2]]+vars[[2]];

If[NumberQ[OptionValue[Distance]],
Sum[(OptionValue[Distance]-EuclideanDistance[front[[i,1]],back[[i]]])^2+(OptionValue[Distance]-EuclideanDistance[front[[i,2]],back[[i]]])^2,{i,Length[front]}]
,
Sum[(EuclideanDistance[front[[i,1]],back[[i]]]^2.+EuclideanDistance[front[[i,2]],back[[i]]]^2.),{i,Length[front]}]
]
];


alignmentDistance[neighborsIn_,vars_,OptionsPattern[]]:=Module[
{ref,data,rot},

ref=neighborsIn[[1]];
data=neighborsIn[[2]];

If[Length[vars]==3,
rot=RotationMatrix[vars[[3]]Degree];

data=Transpose[rot.Transpose[data]];

data[[All,1]]=data[[All,1]]+vars[[1]];
data[[All,2]]=data[[All,2]]+vars[[2]];
,
rot=RotationMatrix[vars[[4]]Degree,{1,0,0}].RotationMatrix[vars[[5]]Degree,{0,1,0}].RotationMatrix[vars[[6]]Degree,{0,0,1}];

data=Transpose[rot.Transpose[data]];

data[[All,1]]=data[[All,1]]+vars[[1]];
data[[All,2]]=data[[All,2]]+vars[[2]];
data[[All,3]]=data[[All,3]]+vars[[3]];
];

Sum[(EuclideanDistance[ref[[i]],data[[i]]]^2),{i,Length[ref]}]

];


(* ::Subsection::Closed:: *)
(*Line Scans*)


Options[lineScan]={
CalculateDistance->False
};


lineScan[dataIn_,lines_,OptionsPattern[]]:=Module[{
data,linesIn,ref,nearestXYPos,step,xIn,yIn,lineFit,coords,positions,dataOut,out
},

If[Length[Dimensions[dataIn]]==3,
data=Flatten[dataIn,1];,
data=dataIn;];

If[Length[Dimensions[lines]]==2,
linesIn={lines};,
linesIn=lines;];

ref=ProgressDialog[];

nearestXYPos=Nearest[data[[All,1;;2]]->Automatic];

dataOut=Table[
ref @ SetPropertyValue[{"bar", "value"}, Round[100*(i/Length[linesIn])]];
ref @ SetPropertyValue[{"label", "text"}, "Procesing line "<>ToString[i]<>" of "<>ToString[Length[linesIn]]];

xIn=Sort[Map[#[[1]]&,linesIn[[i]]]];
yIn=Sort[Map[#[[2]]&,linesIn[[i]]]];

If[TrueQ[First[Differences[xIn]]<First[Differences[yIn]]],
lineFit[y_]= Chop[Fit[Reverse[linesIn[[i]],2], {1,y},y]];
step=EuclideanDistance@@DeleteDuplicates[SortBy[data,#[[2]]&][[All,2]]][[1;;2]];
coords=Table[{lineFit[y],y},
{y,yIn[[1]],yIn[[2]],step/2.}];
,
lineFit[x_]= Chop[Fit[linesIn[[i]], {1,x},x]];
step=EuclideanDistance@@DeleteDuplicates[SortBy[data,#[[1]]&][[All,1]]][[1;;2]];
coords=Table[{x,lineFit[x]},
{x,xIn[[1]],xIn[[2]],step/2.}];
];

positions=DeleteDuplicates[Map[nearestXYPos[#]&,coords]];

out=Extract[data,positions];

If[OptionValue[CalculateDistance],
Map[{#[[1]],#[[2]],EuclideanDistance[out[[1,1;;2]],#[[1;;2]]],#[[3]]}&,out]
,
out
]

,{i,Length[linesIn]}];

CloseGUIObject[ref];
ReleaseGUIObject[ref];

If[First[Dimensions[linesIn]]==1,
dataOut[[1]]
,
dataOut
]

];


(* ::Subsection::Closed:: *)
(*Get Stitch Displacments or Coords*)


nodalCoords1={{0,0},{0,0},{0,0}};
nodalCoords2={{0,0},{0,0},{0,0}};
RepositionDispsOut={};
StitchCoordsOut={};


Options[getRepositionDisp]={
RepositionDirection->"Down",size->500,BoxSize->75
};


getRepositionDisp[filesPath_,{Frame1_,Frame2_},OptionsPattern[]]:=Module[{
filesIn,imagePathIn,varListIn,xNumb,yNumb,XNumb,YNumb,ZNumb,sigmaNumb,DICdata1,DICdata2,nearestxy1,nearestxy2,subsetSize,plotRange1,plotRange2,greyPlt1,greyPlt2,image1,image2,coordsOut,imageRanges,matches
},

$HistoryLength=5;

filesIn=getDICfiles[filesPath];
imagePathIn=filesIn[[3,{Frame1,Frame2}]];

varListIn=filesIn[[2]];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];
ZNumb=Position[varListIn,"Z"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

DICdata1=MapThread[List,Import[filesIn[[1,Frame1]]],2];
DICdata2=MapThread[List,Import[filesIn[[1,Frame2]]],2];

DICdata1=Delete[Flatten[DICdata1,1],Position[Flatten[DICdata1,1][[All,sigmaNumb]],-1.]];
DICdata2=Delete[Flatten[DICdata2,1],Position[Flatten[DICdata2,1][[All,sigmaNumb]],-1.]];

nearestxy1=Nearest[DICdata1[[All,{xNumb,yNumb}]]->DICdata1[[All,{XNumb,YNumb}]]];
nearestxy2=Nearest[DICdata2[[All,{xNumb,yNumb}]]->DICdata2[[All,{XNumb,YNumb}]]];

greyPlt1=Import[imagePathIn[[1]],"GrayLevels"];
greyPlt2=Import[imagePathIn[[2]],"GrayLevels"];

image1=Import[imagePathIn[[1]]];
image2=Import[imagePathIn[[2]]];

plotRange1={Sort[DICdata1[[All,xNumb]]][[{1,-1}]],(Length[greyPlt1]+1-Sort[DICdata1[[All,yNumb]]][[{-1,1}]])};
plotRange2={Sort[DICdata2[[All,xNumb]]][[{1,-1}]],(Length[greyPlt2]+1-Sort[DICdata2[[All,yNumb]]][[{-1,1}]])};

If[ListQ[OptionValue[BoxSize]],
subsetSize=OptionValue[BoxSize];
,
subsetSize={OptionValue[BoxSize],OptionValue[BoxSize]};
];

Share[];

If[StringMatchQ[OptionValue[RepositionDirection],"u*",IgnoreCase->True]||StringMatchQ[OptionValue[RepositionDirection],"d*",IgnoreCase->True],
Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[nodalCoords1],
Graphics[{
Raster[Reverse[greyPlt1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords1]]]
},
PlotRange->plotRange1,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]},
{LocatorPane[
Dynamic[nodalCoords2],
Graphics[{
Raster[Reverse[greyPlt2]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords2]]]
},
PlotRange->plotRange2,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear Reposition Displacement",RepositionDispsOut={};]},
{Button["Calculate and Save Reposition Displacement",
AppendTo[RepositionDispsOut,
coordsOut={Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,nodalCoords1],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,nodalCoords2]};
imageRanges=IntegerPart[Map[{{#[[1]],#[[1]]+subsetSize[[1]]},{#[[2]]-subsetSize[[2]],#[[2]]}}&,coordsOut,{2}]];
matches=MapThread[ImageCorrespondingPoints[ImageTake[image1,#1[[2]],#1[[1]]],ImageTake[image2,#2[[2]],#2[[1]]],"Transformation"->"Rigid"]&,imageRanges];
matches=Table[
{Map[{coordsOut[[1,i,1]]+#[[1]],coordsOut[[1,i,2]]-#[[2]]}&,matches[[i,1]]],Map[{coordsOut[[2,i,1]]+#[[1]],coordsOut[[2,i,2]]-#[[2]]}&,matches[[i,2]]]}
,{i,Length[coordsOut[[1]]]}];
Mean[Map[Mean[MapThread[First[nearestxy1[#1]]-First[nearestxy2[#2]]&,#]]&,matches]]
];]
}}]]
,
Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[nodalCoords1],
Graphics[{
Raster[Reverse[greyPlt1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords1]]]
},
PlotRange->plotRange1,
AspectRatio->Automatic,
ImageSize->{Automatic,OptionValue[size]}
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]
,
LocatorPane[
Dynamic[nodalCoords2],
Graphics[{
Raster[Reverse[greyPlt2]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords2]]]
},
PlotRange->plotRange2,
AspectRatio->Automatic,
ImageSize->{Automatic,OptionValue[size]}
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear Reposition Displacement",RepositionDispsOut={};]},
{Button["Calculate and Save Reposition Displacement",
AppendTo[RepositionDispsOut,
coordsOut={Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,nodalCoords1],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,nodalCoords2]};
imageRanges=IntegerPart[Map[{{#[[1]],#[[1]]+subsetSize[[1]]},{#[[2]]-subsetSize[[2]],#[[2]]}}&,coordsOut,{2}]];
matches=MapThread[ImageCorrespondingPoints[ImageTake[image1,#1[[2]],#1[[1]]],ImageTake[image2,#2[[2]],#2[[1]]],"Transformation"->"Rigid"]&,imageRanges];
matches=Table[
{Map[{coordsOut[[1,i,1]]+#[[1]],coordsOut[[1,i,2]]-#[[2]]}&,matches[[i,1]]],Map[{coordsOut[[2,i,1]]+#[[1]],coordsOut[[2,i,2]]-#[[2]]}&,matches[[i,2]]]}
,{i,Length[coordsOut[[1]]]}];
Mean[Map[Mean[MapThread[First[nearestxy1[#1]]-First[nearestxy2[#2]]&,#]]&,matches]]
];]
}}]]
]

];


Options[getStitchCoords]={
StitchDirection->"Down",size->500,BoxSize->75
};


getStitchCoords[filesPath_,{Frame1_,Frame2_},OptionsPattern[]]:=Module[{
filesIn,imagePathIn,varListIn,xNumb,yNumb,XNumb,YNumb,ZNumb,sigmaNumb,DICdata1,DICdata2,nearestxy1,nearestxy2,subsetSize,plotRange1,plotRange2,greyPlt1,greyPlt2,image1,image2,coordsOut,imageRanges,matches
},

$HistoryLength=5;

filesIn=getDICfiles[filesPath];
imagePathIn=filesIn[[3,{Frame1,Frame2}]];

varListIn=filesIn[[2]];

xNumb=Position[varListIn,"x"][[1,1]];
yNumb=Position[varListIn,"y"][[1,1]];
XNumb=Position[varListIn,"X"][[1,1]];
YNumb=Position[varListIn,"Y"][[1,1]];
ZNumb=Position[varListIn,"Z"][[1,1]];
sigmaNumb=Position[varListIn,"sigma"][[1,1]];

DICdata1=MapThread[List,Import[filesIn[[1,Frame1]]],2];
DICdata2=MapThread[List,Import[filesIn[[1,Frame2]]],2];

DICdata1=Delete[Flatten[DICdata1,1],Position[Flatten[DICdata1,1][[All,sigmaNumb]],-1.]];
DICdata2=Delete[Flatten[DICdata2,1],Position[Flatten[DICdata2,1][[All,sigmaNumb]],-1.]];

nearestxy1=Nearest[DICdata1[[All,{xNumb,yNumb}]]->DICdata1[[All,{XNumb,YNumb}]]];
nearestxy2=Nearest[DICdata2[[All,{xNumb,yNumb}]]->DICdata2[[All,{XNumb,YNumb}]]];

greyPlt1=Import[imagePathIn[[1]],"GrayLevels"];
greyPlt2=Import[imagePathIn[[2]],"GrayLevels"];

image1=Import[imagePathIn[[1]]];
image2=Import[imagePathIn[[2]]];

plotRange1={Sort[DICdata1[[All,xNumb]]][[{1,-1}]],(Length[greyPlt1]+1-Sort[DICdata1[[All,yNumb]]][[{-1,1}]])};
plotRange2={Sort[DICdata2[[All,xNumb]]][[{1,-1}]],(Length[greyPlt2]+1-Sort[DICdata2[[All,yNumb]]][[{-1,1}]])};

If[ListQ[OptionValue[BoxSize]],
subsetSize=OptionValue[BoxSize];
,
subsetSize={OptionValue[BoxSize],OptionValue[BoxSize]};
];

Share[];

If[StringMatchQ[OptionValue[StitchDirection],"u*",IgnoreCase->True]||StringMatchQ[OptionValue[StitchDirection],"d*",IgnoreCase->True],
Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[nodalCoords1],
Graphics[{
Raster[Reverse[greyPlt1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords1]]]
},
PlotRange->plotRange1,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]},
{LocatorPane[
Dynamic[nodalCoords2],
Graphics[{
Raster[Reverse[greyPlt2]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords2]]]
},
PlotRange->plotRange2,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear Stitch Coords",StitchCoordsOut={};]},
{Button["Calculate and Save Stitch Coords",
AppendTo[StitchCoordsOut,
coordsOut={Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,nodalCoords1],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,nodalCoords2]};
imageRanges=IntegerPart[Map[{{#[[1]],#[[1]]+subsetSize[[1]]},{#[[2]]-subsetSize[[2]],#[[2]]}}&,coordsOut,{2}]];
matches=MapThread[ImageCorrespondingPoints[ImageTake[image1,#1[[2]],#1[[1]]],ImageTake[image2,#2[[2]],#2[[1]]],"Transformation"->"Rigid"]&,imageRanges];
matches=Table[
{Map[{coordsOut[[1,i,1]]+#[[1]],coordsOut[[1,i,2]]-#[[2]]}&,matches[[i,1]]],Map[{coordsOut[[2,i,1]]+#[[1]],coordsOut[[2,i,2]]-#[[2]]}&,matches[[i,2]]]}
,{i,Length[coordsOut[[1]]]}];
{Map[Mean[Map[First[nearestxy1[#]]&,#]]&,matches[[All,1]]],Map[Mean[Map[First[nearestxy2[#]]&,#]]&,matches[[All,2]]]}
];]
}}]]
,
Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[nodalCoords1],
Graphics[{
Raster[Reverse[greyPlt1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords1]]]
},
PlotRange->plotRange1,
AspectRatio->Automatic,
ImageSize->{Automatic,OptionValue[size]}
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]
,
LocatorPane[
Dynamic[nodalCoords2],
Graphics[{
Raster[Reverse[greyPlt2]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,nodalCoords2]]]
},
PlotRange->plotRange2,
AspectRatio->Automatic,
ImageSize->{Automatic,OptionValue[size]}
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear Stitch Coords",StitchCoordsOut={};]},
{Button["Calculate and Save Stitch Coords",
AppendTo[StitchCoordsOut,
coordsOut={Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,nodalCoords1],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,nodalCoords2]};
imageRanges=IntegerPart[Map[{{#[[1]],#[[1]]+subsetSize[[1]]},{#[[2]]-subsetSize[[2]],#[[2]]}}&,coordsOut,{2}]];
matches=MapThread[ImageCorrespondingPoints[ImageTake[image1,#1[[2]],#1[[1]]],ImageTake[image2,#2[[2]],#2[[1]]],"Transformation"->"Rigid"]&,imageRanges];
matches=Table[
{Map[{coordsOut[[1,i,1]]+#[[1]],coordsOut[[1,i,2]]-#[[2]]}&,matches[[i,1]]],Map[{coordsOut[[2,i,1]]+#[[1]],coordsOut[[2,i,2]]-#[[2]]}&,matches[[i,2]]]}
,{i,Length[coordsOut[[1]]]}];
{Map[Mean[Map[First[nearestxy1[#]]&,#]]&,matches[[All,1]]],Map[Mean[Map[First[nearestxy2[#]]&,#]]&,matches[[All,2]]]}
];]
}}]]
]

];


(* ::Subsection::Closed:: *)
(*Contour Overlays*)


Options[DICcontourPlot]={
ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,interpolationOrder->None,plotRange->All,color->"Rainbow",size->500,NumberofPoints->50000,AutoPlaneFit->False
};


DICcontourPlot[dataIn_,OptionsPattern[]]:=Module[{
planeFit,step,takeNumber,reducedData,incr,contourList,colorFS,colorF,contPlt,greyPlt
},

takeNumber=Round[Sqrt[(Length[dataIn]*Length[dataIn[[1]]])/OptionValue[NumberofPoints]]];
reducedData=Flatten[Take[dataIn,{1,-1,takeNumber},{1,-1,takeNumber}],1];
step=Abs[reducedData[[2,1]]-reducedData[[1,1]]];

If[OptionValue[AutoPlaneFit],
planeFit[x_,y_]=a+b*x+c*y/.FindFit[reducedData,a+b*x+c*y,{a,b,c},{x,y}];
reducedData=Map[{(#[[1]]),
(#[[2]]),(#[[3]]-planeFit[#[[1]],#[[2]]])}&,reducedData];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],
incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{Chop[#],OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{Chop[#],OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
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
InterpolationOrder->OptionValue[interpolationOrder],
ContourShading->OptionValue[contourShading]];

Share[];

Graphics[{
contPlt[[1]]
},
PlotRange->OptionValue[plotRange],
PlotRangePadding->None,
AspectRatio->Automatic,
ImageSize->OptionValue[size]]
];


Options[DICcontourPlot2]={
ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,plotRange->All,color->"Rainbow",size->500,NumberofPoints->50000,UseRegionFunction->True,AutoPlaneFit->False,RemoveBadNodes->False
};


DICcontourPlot2[dataIn_,OptionsPattern[]]:=Module[{
planeFit,step,takeNumber,reducedData,badNodes,incr,contourList,colorFS,colorF,contPlt,greyPlt
},

takeNumber=Round[Sqrt[(Length[dataIn]*Length[dataIn[[1]]])/OptionValue[NumberofPoints]]];
reducedData=Flatten[Take[dataIn,{1,-1,takeNumber},{1,-1,takeNumber}],1];
step=Abs[reducedData[[2,1]]-reducedData[[1,1]]];
badNodes=Nearest[Map[reducedData[[#,{1,2}]]&,Flatten[Position[reducedData[[All,4]],-1.],1]]];

If[OptionValue[AutoPlaneFit],
planeFit[x_,y_]=a+b*x+c*y/.FindFit[Delete[Flatten[reducedData,1],Position[Flatten[reducedData,1][[All,4]],-1.]][[All,{1,2,3}]],a+b*x+c*y,{a,b,c},{x,y}];
reducedData=Map[{(#[[1]]+1-dataIn[[1,1,1]]),
(dataIn[[-1,-1,2]]+1-#[[2]]),(#[[3]]-planeFit[#[[1]],#[[2]]]),#[[4]]}&,reducedData,{2}];
,
reducedData=reducedData;
];

If[OptionValue[RemoveBadNodes],
reducedData=Delete[reducedData,Position[reducedData[[All,4]],-1.]];
,
reducedData=reducedData;
];

If[ListQ[OptionValue[ContourRange]],
incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{Chop[#],OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{Chop[#],OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;];


If[OptionValue[UseRegionFunction],
contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
RegionFunction->(EuclideanDistance[{#1,#2},First[badNodes[{#1,#2}]]]>step&),
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],ClippingStyle->None,ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];
,
contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],ClippingStyle->None,ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];
];

Share[];

Graphics[{
contPlt[[1]]
},
PlotRange->OptionValue[plotRange],
PlotRangePadding->None,
AspectRatio->Automatic,
ImageSize->OptionValue[size]]
];


(* ::Subsection::Closed:: *)
(*Selection GUIs*)


points={{0,0},{0,0}};
pointsOut={};


Options[getPoints]={
plotRange->All,size->750
};


getPoints[pltIn_,OptionsPattern[]]:=Module[{
},

Deploy[
Grid[{{
LocatorPane[
Dynamic[points],
Graphics[{
pltIn[[1]],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[points]],
Orange,
Dynamic[Line/@pointsOut],
Dynamic[Point/@Flatten[pointsOut,1]]
}
,
PlotRange->OptionValue[plotRange],
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]}]},
{Button["Clear coordinates",pointsOut={};]},
{Button["Save coordinates",AppendTo[pointsOut,points];]}
}]]
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["TopoFunctions`*"],
{
"nodalCoords1","nodalCoords2","RepositionDispsOut","StitchCoordsOut","points","pointsOut"
}];
EndPackage[ ]
