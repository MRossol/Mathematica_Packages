(* ::Package:: *)

BeginPackage["VEECOfunctions`"]
Unprotect@@Names["VEECOfunctions`*"];
ClearAll@@Names["VEECOfunctions`*"];


(*Data Import*)
loadVEECOdata::usage = "loadVEECOdata[pathIn_] loads VEECO data and processes as text";
convertVEECOdata::usage = "convertVEECOdata[pathIn_,wavelength_] converts VEECO data to {x,y,z} in um and removes bad points";

(*Data Manipulation*)
autoPlaneFit::usage = "autoPlaneFit[dataIn] auto plane fits the input data";
getPixelsize::usage = "getPixelsize[pathIn_] determines pixel size from .txt file";
createBWImg::usage = "createBWImg[dataIn_,step_] creates a BW image of VEECO data with white regions designating bad data";
arrayCrop::usage = "arrayCrop[dataIn_,{xRange_,yRange_}] selects data in the x and y ranges given while retaining an array structure";
arraySelect::usage = "arraySelect[dataIn_,{xRange_,yRange_}] selects data in the x and y ranges given from an non-orthogonal array";
arrayThreshold::usage = "arrayThreshold[dataIn_,{thresholdMin_,thresholdMax_},OptionsPattern[]] uses MorphologicalComponents to find clusters in arrays of data.";

(*Delta H*)
localDeltaH::usage = "localDeltaH[dataIn_,OptionsPatter[]] calculates the mean height change for a local region dataIn";
deltaH::usage = "deltaH[dataIn_,OptionsPattern[]] calculates the local height change for dataIn";
subRegionDeltaH::usage = "subRegionDeltaH[dataIn_,aoi_] calculates the local heigh change for a subregion of data using all matrix values in the region";
meanDeltaHscan::usage = "meanDeltaHscan[dataIn_,scanWidth_] scans along x calculating {Mean[x],Mean[y],Mean[\[CapitalDelta]H],StDev[\[CapitalDelta]H],# of points} for  strips of width = scanWidth";
dHandVf::usage = "dHandVf[dataIn_,rangeIn_,OptionsPattern[]] calculates {Mean X, Mean Y, Mean dH, StDev dH, Vf} for the region of dataIn specified";
surfaceDisp::usage = "surfaceDisp[refData_,dispData_] calculates z displacement between ref and disp surfaces";

(*Line Scans*)
lineScan::usage = "lineScan[dataIn_,lines_] extracts lines in";

(*Alignment*)
nodalCoords1::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3},{x4,y4}}";
nodalCoords2::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3},{x4,y4}}";
AlignmentCoordsOut::usage = "List of coordinates corresponding to identical regions in image 1 and 2";
getAlignmentCoords::usage = "getAlignmentCoords[img1,img2,step] allows for selection of corresponding points on two BW images";
BWcoords1::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3}}";
BWcoords2::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3}}";
BWcoordsOut::usage = "List of coordinates corresponding to identical regions in image 1 and 2";
getBWcoords::usage = "getBWcoords[img1,img2,step] allows for selection of corresponding points on two BW images";
refNode::usaage = "node with {x1,y1} for reference surface";
expNode::usage = "node with {x1,y1} for exposed surface";
WindowAlignmentOut::usage = "List of {refNode,expNode} in 3D to align surfaces";
getWindowAlignment::usage = "getWindowAlignment[{refDataIn_,expDataIn_},windowC_,windowS_,OptionsPattern[]] allows selection of identical on subsections of ref and exposed surface";
refNodes::usaage = "nodes with {{x1,y1},{x2,y2},{x3,y3}} for reference surface";
expNodes::usage = "nodes with {{x1,y1},{x2,y2},{x3,y3}} for exposed surface";
FullAlignmentOut::usage = "List of {refNode,expNode} in 3D to align surfaces";
getFullAlignment::usage = "getFullAlignment[{refDataIn_,expDataIn_},windowC_,windowS_,OptionsPattern[]] allows selection of identical on subsections of ref and exposed surface";
edgeAlign::usage = "edgeAlign[dataIn_,coordsIn_,OptionsPattern[]] calculates rotation and translation then repositions dataIn to align edge with y-axis";
alignmentParameters::usage = "alignmentParameters[dataIn_,vars_] calculates {xDisp,yDips,rot} from AlignmentCoordsOut";
sampleAlign::usage = "sampleAlign[dataIn_,coordsIn_] calculates rotation and translation then repositions dataIn to align dataIn with pristine data";

(*Selection GUI*)
points::usage = "dynamic points {{x1,y1},{x2,y2}} for GUI getPoints";
pointsOut::usage = "saved Points from GUI getPoints";
getPoints::usage = "getPoints[pltIn_] select points on pltIn using GUI";
nodes::usage = "nodes with {{x1,y1},{x2,y2}}";
coordinates::usage = "List of coordinates for {node1,node2} in real space({0,0} is bottom left corner)";
coordinatesOut::usage = "List of coordinates for {node1,node2} in the camera array ({0,0} is upper left corner)";
getCoordinates::usage = "getCoordinates[filesPath_,frame_,varIn_] allows for selection of nodal locations";


Begin["`Private`"]


Needs["GUIKit`"]
Needs["HierarchicalClustering`"]


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


Options[loadVEECOdata]={
ExportData->False,RemoveBadData->True,ConvertToArray->False
};


loadVEECOdata[pathIn_,OptionsPattern[]]:=Module[
{fileIn,waveLength,end,data,badPos,pixelSize,pixelData},

fileIn=Import[pathIn,"CSV"];

waveLength=fileIn[[First[Flatten[Position[fileIn[[All,1]],"Wavelength"]]],-1]]/1000.;

end=Flatten[Position[fileIn[[All,1]],"Intensity"]];

If[Dimensions[end]=={1},
data=Cases[fileIn[[1;;end[[1]]]],{_Real,__}];
,
data=Cases[fileIn,{_Real,__}];
];

If[OptionValue[RemoveBadData],
data=Cases[data,{_Real,_Real,_Real}];
,
If[SameQ[Length[Dimensions[data]],2],
badPos=Cases[Position[data[[All,3]],Except[_Real]][[2;;]],{_Integer}];
data[[All,3]]=ReplacePart[data[[All,3]],badPos->0.];
data=Join[data,ReplacePart[ConstantArray[0.1,{Length[data],1}],badPos->{-1.}],2];
,
data=Join[PadRight[data],ReplacePart[ConstantArray[0.1,{Length[data],1}],Position[Map[Length[#]&,data],2]->{-1.}],2];
];
];

data[[All,3]]=data[[All,3]]*waveLength;

If[OptionValue[ConvertToArray],
pixelSize=EuclideanDistance@@data[[1;;2,1;;2]];
pixelData=Round[data[[All,{2,1}]]/pixelSize]+1;
data=ReplacePart[ConstantArray[0,{Max[pixelData[[All,1]]],Max[pixelData[[All,2]]]}],MapThread[#1->#2&,{pixelData,data}]];
,
data=data;
];

If[OptionValue[ExportData],
Export[StringTrim[pathIn,".txt"]<>".MAT",data]
,
data]

];


Options[convertVEECOdata]={
ExportData->False
};


convertVEECOdata[pathIn_,waveLength_,OptionsPattern[]]:=Module[
{data},

data=Cases[Import[pathIn,"CSV"],{_Real,_Real,_Real}];
data[[All,3]]=data[[All,3]]*waveLength;

If[OptionValue[ExportData],
Export[StringTrim[pathIn,".txt"]<>".MAT",data]
,
data]

];


(* ::Subsection::Closed:: *)
(*Data Manipulation*)


Options[autoPlaneFit]={
GetPlane->False,PlaneIn->None,ShiftOrigin->True,RemoveBadData->True
};


autoPlaneFit[dataIn_,OptionsPattern[]]:=Module[
{data,planeFit,xo,yo,a,b,c,x,y,dataOut,ones,shift},

If[SameQ[Length[Dimensions[dataIn]],3],
data=Flatten[dataIn,1];
,
data=dataIn;
];

If[OptionValue[RemoveBadData],
data=Select[data,(#[[4]]!=-1.)&][[All,1;;3]];
,
data=data[[All,1;;3]];
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
ones=dataOut[[All,All,1;;3]];
ones[[All,All,3]]=1;
dataOut[[All,All,3]]=dataOut[[All,All,3]]-ones.planeFit;
,
ones=dataOut[[All,1;;3]];
ones[[All,3]]=1;
dataOut[[All,3]]=dataOut[[All,3]]-ones.planeFit;
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


getPixelsize[pathIn_]:=Module[
{textListIn,textPos,waveLength,data,Bad},

textListIn=ReadList[pathIn,String];

textPos=DeleteDuplicates[Position[Map[StringPosition[#,"Pixel_size"]&,textListIn],_Integer,Infinity][[All,1]]];

ToExpression[StringSplit[textListIn[[textPos[[1]]]],","][[-1]]]

];


createBWImg[dataIn_,step_]:=Module[
{pixelData},

pixelData=dataIn;
pixelData[[All,1;;2]]=Round[pixelData[[All,{2,1}]]/step]+1;

Reverse[ReplacePart[ConstantArray[0,{Max[pixelData[[All,1]]],Max[pixelData[[All,2]]]}],Select[pixelData,(#[[4]]==-1.)&][[All,1;;2]]->1]]
]


arrayCrop[dataIn_,{xRange_,yRange_}]:=Module[
{xPos,yPos},
xPos=Sort[Map[First[Nearest[dataIn[[1,All,1]]->Automatic,#]]&,xRange]];
yPos=Sort[Map[First[Nearest[dataIn[[All,1,2]]->Automatic,#]]&,yRange]];

dataIn[[yPos[[1]];;yPos[[2]],xPos[[1]];;xPos[[2]]]]
]


Options[arraySelect]={
RemoveBadData->True
};


arraySelect[dataIn_,{xRange_,yRange_},OptionsPattern[]]:=Module[
{xyRange,xWindow,yWindow,xPos,yPos},
xyRange={Sort[xRange],Sort[yRange]};

xWindow=First[Differences[xyRange[[1]]]]/4.;
xWindow={-xWindow,xWindow};
xWindow=xyRange[[1]]+xWindow;

yWindow=First[Differences[xyRange[[2]]]]/4.;
yWindow={-yWindow,yWindow};
yWindow=xyRange[[2]]+yWindow;

xPos=Sort[Map[First[Nearest[dataIn[[1,All,1]]->Automatic,#]]&,xWindow]];
yPos=Sort[MapThread[First[Nearest[dataIn[[All,#2,2]]->Automatic,#1]]&,{yWindow,xPos}]];
xPos=Sort[MapThread[First[Nearest[dataIn[[#2,All,1]]->Automatic,#1]]&,{xWindow,yPos}]];
yPos=Sort[MapThread[First[Nearest[dataIn[[All,#2,2]]->Automatic,#1]]&,{yWindow,xPos}]];
xPos=Sort[MapThread[First[Nearest[dataIn[[#2,All,1]]->Automatic,#1]]&,{xWindow,yPos}]];

If[OptionValue[RemoveBadData],
Select[Flatten[dataIn[[yPos[[1]];;yPos[[2]],xPos[[1]];;xPos[[2]]]],1],(xyRange[[1,1]]<#[[1]]<xyRange[[1,2]]&&xyRange[[2,1]]<#[[2]]<xyRange[[2,2]]&&#[[4]]!= -1.)&]
,
Select[Flatten[dataIn[[yPos[[1]];;yPos[[2]],xPos[[1]];;xPos[[2]]]],1],(xyRange[[1,1]]<#[[1]]<xyRange[[1,2]]&&xyRange[[2,1]]<#[[2]]<xyRange[[2,2]])&]
]
]


Options[arrayThreshold]={
ThresholdVariable->3,VariablesOut->1;;4
};


arrayThreshold[dataIn_,{thresholdMin_,thresholdMax_},OptionsPattern[]]:=Module[
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


(* ::Subsection::Closed:: *)
(*Delta H*)


Options[localDeltaH]={
AutoPlaneFit->True,RemoveBadData->True
};


localDeltaH[dataIn_,OptionsPattern[]]:=Module[
{data,vars,fit,mp},

If[SameQ[Length[Dimensions[dataIn]],2],
data=dataIn;
,
data=Flatten[dataIn,1];
];

If[OptionValue[AutoPlaneFit],
If[OptionValue[RemoveBadData],
data=Select[data,(#[[4]]!=-1.)&][[All,1;;3]];
,
data=data[[All,1;;3]];
];
vars=data;
vars[[All,3]]=1;
fit=PseudoInverse[vars].data[[All,3]];
data=data[[All,3]]-vars.fit;
,
data=data[[All,3]]
];

mp=Mean[{Mean[Select[data,(0.<=#)&]],Mean[Select[data,(0.>=#)&]]}];
Mean[Select[data,(mp<#)&]]-Mean[Select[data,(#<mp)&]]
];


Options[deltaH]={
StepSize->1,ZeroMatrix->True
};


deltaH[dataIn_,windowS_,OptionsPattern[]]:=Module[
{windowX,windowY,yLength,xLength,start,xEnd,yEnd,hStep,zMatrixBoole,point,window,pointPos,vars,fit,mp,dh},

If[ListQ[windowS],
{windowX,windowY}=Round[windowS/2];
,
windowX=Round[windowS/2];
windowY=Round[windowS/2];
];

{yLength,xLength}=Dimensions[dataIn][[1;;2]];
start[in_]:=Switch[in<1,True,1,False,in];
xEnd[xIn_]:=Switch[xIn>xLength,True,xLength,False,xIn];
yEnd[yIn_]:=Switch[yIn>yLength,True,yLength,False,yIn];

hStep=OptionValue[StepSize];
zMatrixBoole=OptionValue[ZeroMatrix];

DistributeDefinitions[start,xEnd,yEnd,yLength,xLength,windowX,windowY,hStep,zMatrixBoole];
ParallelTable[
point=dataIn[[y,x]];

If[SameQ[point[[4]],-1.],
point[[3]]=-1.;
,
window=Select[Flatten[dataIn[[start[y-windowY];;yEnd[y+windowY],start[x-windowX];;xEnd[x+windowX]]],1],(#[[4]]!= -1.)&][[All,3]];
If[SameQ[Length[window],0],
point[[3]]=-1.;
,
mp=Mean[window];
If[zMatrixBoole&&point[[3]]>mp,
point[[3]]=0.;
,
dh=Mean[Select[window,(mp<#)&]]-point[[3]];
If[NumberQ[dh],
point[[3]]=dh;
,
point[[3]]=0.;
];
];
];
];
point
,{y,1,yLength,hStep},{x,1,xLength,hStep}]

];


subRegionDeltaH[dataIn_,aoi_]:=Module[
{dataSubset,badPos,goodData,mp,meanM},

dataSubset=autoPlaneFit[arraySelect[dataIn,aoi,RemoveBadData->False]];
badPos=Position[dataSubset[[All,4]],-1.];

goodData=Delete[dataSubset,badPos][[All,3]];
mp=Mean[goodData];
meanM=Mean[Select[goodData,(mp<#)&]];

dataSubset[[All,3]]=ReplacePart[Map[If[#>mp,0,meanM-#]&,dataSubset[[All,3]]],badPos->-1.];

dataSubset

];


meanDeltaHscan[dataIn_,scanWidth_]:=Module[
{scanRange,strip,dataOut={}},

scanRange=Floor[Sort[dataIn[[All,1]]][[{1,-1}]]/scanWidth];

Do[
strip=Select[dataIn,(scanWidth*i<#[[1]]<scanWidth*(i+1)&&#[[3]]>0)&];
If[Length[strip]>1,
AppendTo[dataOut,Mean[strip]~Join~{StandardDeviation[strip[[All,3]]],Length[strip]}]]
,{i,scanRange[[1]],scanRange[[2]]}];
dataOut
];


Options[surfaceDisp]={
AutoPlaneFit->True
};


surfaceDisp[refDataIn_,dispDataIn_,OptionsPattern[]]:=Module[
{refData,dispData,refNearest,disp},

If[SameQ[Length[Dimensions[refDataIn]],3],
refData=Flatten[refDataIn,1];
,
refData=refDataIn;
];
dispData=dispDataIn;

refNearest=Nearest[refData[[All,1;;2]]->refData[[All,3]]];
disp[point_]:=Switch[SameQ[point[[4]],-1.],True,0,False,(point[[3]]-First[refNearest[point[[1;;2]]]])];

If[SameQ[Length[Dimensions[dispData]],3],
dispData[[All,All,3]]=Map[disp[#]&,dispData,{2}];
,
dispData[[All,3]]=Map[disp[#]&,dispData];
];

If[OptionValue[AutoPlaneFit],
autoPlaneFit[dispData,GetPlane->False,PlaneIn->None,ShiftOrigin->False,RemoveBadData->True]
,
dispData
]

];


Options[dHandVf]={
MeanZThreshold->True
};


dHandVf[dataIn_,rangeIn_,OptionsPattern[]]:=Module[
{regionData,total,goodData,meanZ},

regionData=arraySelect[dataIn,rangeIn,RemoveBadData->False];
total=Length[regionData];
goodData=Select[regionData,(#[[4]]!=-1.)&][[All,1;;3]];

If[OptionValue[MeanZThreshold],
meanZ=Mean[goodData[[All,3]]];
goodData=Select[goodData,(#[[3]]>meanZ)&];
Mean[goodData]~Join~{StandardDeviation[goodData[[All,3]]],N[Length[goodData]/total]}
,
goodData=Select[goodData,(#[[3]]>0)&];
Mean[goodData]~Join~{StandardDeviation[goodData[[All,3]]],N[Length[goodData]/total]}
]
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
(*Alignment*)


nodalCoords1={{0,0},{0,0},{0,0},{0,0}};
nodalCoords2={{0,0},{0,0},{0,0},{0,0}};
AlignmentCoordsOut={};


Options[getAlignmentCoords]={
NumberofPoints->100000,RemoveBadData->True,size->500
};


getAlignmentCoords[dataIn1_,dataIn2_,OptionsPattern[]]:=Module[{
data1,data2,takeNumber,reducedData,plt1,plt2},

If[Length[Dimensions[dataIn1]]==2,
data1=dataIn1;
,
data1=Flatten[dataIn1,1];
];
If[Length[Dimensions[dataIn2]]==2,
data2=dataIn2;
,
data2=Flatten[dataIn2,1];
];

If[OptionValue[RemoveBadData],
data1=Select[data1,(#[[4]]!= -1.)&];
data2=Select[data2,(#[[4]]!= -1.)&];
,
data1=data1;
data2=data2;
];

If[OptionValue[NumberofPoints]<Length[data1],
takeNumber=Round[(Length[data1])/OptionValue[NumberofPoints]];
data1=Take[SortBy[data1,(#[[1]]+#[[2]])&],{1,-1,takeNumber}][[All,1;;2]];
,
data1=data1[[All,1;;2]];
];

If[OptionValue[NumberofPoints]<Length[data2],
takeNumber=Round[(Length[data2])/OptionValue[NumberofPoints]];
data2=Take[SortBy[data2,(#[[1]]+#[[2]])&],{1,-1,takeNumber}][[All,1;;2]];
,
data2=data2[[All,1;;2]];
];

plt1=ListPlot[data1,
AspectRatio->Automatic,
Axes->False,
PlotStyle->Black
];

plt2=ListPlot[data2,
AspectRatio->Automatic,
Axes->False,
PlotStyle->Black
];

Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[nodalCoords1],
Graphics[{
plt1[[1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[nodalCoords1~Join~{nodalCoords1[[1]]}]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0,1,0],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0,0,1],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[1,0,0],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[1,1,0],FontSize->12]
}]},
{LocatorPane[
Dynamic[nodalCoords2],
Graphics[{
plt2[[1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[nodalCoords2~Join~{nodalCoords2[[1]]}]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0,1,0],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0,0,1],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[1,0,0],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[1,1,0],FontSize->12]
}]}}]},
{Button["Clear Alignnment Coords",AlignmentCoordsOut={};]},
{Button["Save Alignment Coords",
AppendTo[AlignmentCoordsOut,{nodalCoords1,nodalCoords2}];]
}}]]
];


BWcoords1={{0,0},{0,0},{0,0}};
BWcoords2={{0,0},{0,0},{0,0}};
BWcoordsOut={};


Options[getBWcoords]={
size->500,BoxSize->200
};


getBWcoords[img1_,img2_,step_,OptionsPattern[]]:=Module[{
greyPlt1,greyPlt2,image1,image2,subsetSize,coordsOut,imageRanges,matches
},

greyPlt1=img1;
greyPlt2=img2;

image1=Image[img1];
image2=Image[img2];

If[ListQ[OptionValue[BoxSize]],
subsetSize=OptionValue[BoxSize];
,
subsetSize={OptionValue[BoxSize],OptionValue[BoxSize]};
];

Share[];

Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[BWcoords1],
Graphics[{
Raster[Reverse[greyPlt1]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,BWcoords1]]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]},
{LocatorPane[
Dynamic[BWcoords2],
Graphics[{
Raster[Reverse[greyPlt2]],
GrayLevel[0],Thickness[0.003],Red,
Dynamic[Line[Map[{#,#+{subsetSize[[1]],0},#+subsetSize,#+{0,subsetSize[[2]]},#}&,BWcoords2]]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear BW Coords",BWcoordsOut={};]},
{Button["Save BW Coords",
AppendTo[BWcoordsOut,
coordsOut={Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,BWcoords1],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,BWcoords2]};
imageRanges=IntegerPart[Map[{{#[[1]],#[[1]]+subsetSize[[1]]},{#[[2]]-subsetSize[[2]],#[[2]]}}&,coordsOut,{2}]];
matches=MapThread[ImageCorrespondingPoints[ImageTake[image1,#1[[2]],#1[[1]]],ImageTake[image2,#2[[2]],#2[[1]]],"Transformation"->"Rigid"]&,imageRanges];
matches=Table[
{Map[{coordsOut[[1,i,1]]+#[[1]],coordsOut[[1,i,2]]-#[[2]]}&,matches[[i,1]]],Map[{coordsOut[[2,i,1]]+#[[1]],coordsOut[[2,i,2]]-#[[2]]}&,matches[[i,2]]]}
,{i,Length[coordsOut[[1]]]}];
{Map[{#[[1]],Length[greyPlt1]+1-#[[2]]}&,Map[Mean[#]&,matches[[All,1]]]],Map[{#[[1]],Length[greyPlt2]+1-#[[2]]}&,Map[Mean[#]&,matches[[All,2]]]]}*step
];]
}}]]
];


refNode={0,0};
expNode={0,0};
WindowAlignmentOut={};


Options[getWindowAlignment]={
ContourRange->{-5.,5.},MinorContours->1000,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->Opacity[0.25,Black],contourShading->Automatic,plotRange->All,color->"Rainbow",size->500,NumberofPoints->50000
};


getWindowAlignment[{refDataIn_,expDataIn_},windowC_,windowS_,OptionsPattern[]]:=Module[{
refC,expC,xyRange,xPos,yPos,refData,refNearest,expData,expNearest,meanZ,incr,contourList,colorFS,colorF,refPlt,expPlt
},

If[SameQ[Length[Dimensions[windowC]],2],
{refC,expC}=windowC;
,
refC=windowC;
expC=windowC;
];

If[SameQ[Length[Dimensions[refDataIn]],3],
xyRange=Transpose[{refC-{windowS/1.5,windowS/1.5},refC+{windowS/1.5,windowS/1.5}}];
xPos=Sort[Map[First[Nearest[refDataIn[[1,All,1]]->Automatic,#]]&,xyRange[[1]]]];
yPos=Sort[MapThread[First[Nearest[refDataIn[[All,#2,2]]->Automatic,#1]]&,{xyRange[[2]],xPos}]];
xPos=Sort[MapThread[First[Nearest[refDataIn[[#2,All,1]]->Automatic,#1]]&,{xyRange[[1]],yPos}]];
yPos=Sort[MapThread[First[Nearest[refDataIn[[All,#2,2]]->Automatic,#1]]&,{xyRange[[2]],xPos}]];
xPos=Sort[MapThread[First[Nearest[refDataIn[[#2,All,1]]->Automatic,#1]]&,{xyRange[[1]],yPos}]];
refData=Select[Flatten[refDataIn[[yPos[[1]];;yPos[[2]],xPos[[1]];;xPos[[2]]]],1],(refC[[1]]-windowS/2.<#[[1]]<refC[[1]]+windowS/2.&&refC[[2]]-windowS/2.<#[[2]]<refC[[2]]+windowS/2.&&#[[4]]!=-1.)&];
,
refData=Select[refDataIn,(refC[[1]]-windowS/2.<#[[1]]<refC[[1]]+windowS/2.&&refC[[2]]-windowS/2.<#[[2]]<refC[[2]]+windowS/2.&&#[[4]]!=-1.)&];
];

refNearest=Nearest[refData[[All,1;;2]]->refData[[All,1;;3]]];
If[OptionValue[NumberofPoints]<Length[refData],
refData=Take[refData,{1,-1,Round[(Length[refData])/OptionValue[NumberofPoints]]}];
,
refData=refData;
];

If[SameQ[Length[Dimensions[expDataIn]],3],
xyRange=Transpose[{expC-{windowS/1.5,windowS/1.5},expC+{windowS/1.5,windowS/1.5}}];
xPos=Sort[Map[First[Nearest[expDataIn[[1,All,1]]->Automatic,#]]&,xyRange[[1]]]];
yPos=Sort[MapThread[First[Nearest[expDataIn[[All,#2,2]]->Automatic,#1]]&,{xyRange[[2]],xPos}]];
xPos=Sort[MapThread[First[Nearest[expDataIn[[#2,All,1]]->Automatic,#1]]&,{xyRange[[1]],yPos}]];
yPos=Sort[MapThread[First[Nearest[expDataIn[[All,#2,2]]->Automatic,#1]]&,{xyRange[[2]],xPos}]];
xPos=Sort[MapThread[First[Nearest[expDataIn[[#2,All,1]]->Automatic,#1]]&,{xyRange[[1]],yPos}]];
expData=Select[Flatten[expDataIn[[yPos[[1]];;yPos[[2]],xPos[[1]];;xPos[[2]]]],1],(expC[[1]]-windowS/2.<#[[1]]<expC[[1]]+windowS/2.&&expC[[2]]-windowS/2.<#[[2]]<expC[[2]]+windowS/2.&&#[[4]]!=-1.)&];
,
expData=Select[expDataIn,(expC[[1]]-windowS/2.<#[[1]]<expC[[1]]+windowS/2.&&expC[[2]]-windowS/2.<#[[2]]<expC[[2]]+windowS/2.&&#[[4]]!=-1.)&];
];
expNearest=Nearest[expData[[All,1;;2]]->expData[[All,1;;3]]];
If[OptionValue[NumberofPoints]<Length[expData],
expData=Take[expData,{1,-1,Round[(Length[expData])/OptionValue[NumberofPoints]]}];
,
expData=expData;
];

meanZ=Mean[refData[[All,3]]];
refData[[All,3]]=refData[[All,3]]-meanZ;
meanZ=Mean[expData[[All,3]]];
expData[[All,3]]=expData[[All,3]]-meanZ;


If[ListQ[OptionValue[ContourRange]],
incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;
,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;
];

refPlt=ListContourPlot[SortBy[refData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

expPlt=ListContourPlot[SortBy[expData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[refNode],
Graphics[{
refPlt[[1]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0,0,1],FontSize->12]
}]},
{LocatorPane[
Dynamic[expNode],
Graphics[{
expPlt[[1]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0,0,1],FontSize->12]
}]}}]},
{Button["Clear 3D Alignment",WindowAlignmentOut={};]},
{Button["Save 3D Alignment",
AppendTo[WindowAlignmentOut,{First[refNearest[refNode]],First[expNearest[expNode]]}];]
}}]]
];


refNodes={{0,0},{0,0},{0,0}};
expNodes={{0,0},{0,0},{0,0}};
FullAlignmentOut={};


Options[getFullAlignment]={
ContourRange->{-5.,5.},MinorContours->100,MajorContours->20,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,plotRange->All,color->"Rainbow",size->500,NumberofPoints->50000
};


getFullAlignment[{refDataIn_,expDataIn_},OptionsPattern[]]:=Module[{
refData,refNearest,expData,expNearest,incr,contourList,colorFS,colorF,refPlt,expPlt
},

If[SameQ[Length[Dimensions[refDataIn]],3],
refData=Select[Flatten[refDataIn,1],(#[[4]]!=-1.)&];
,
refData=Select[refDataIn,(#[[4]]!=-1.)&];
];

refNearest=Nearest[refData[[All,1;;2]]->refData[[All,1;;3]]];
If[OptionValue[NumberofPoints]<Length[refData],
refData=Take[refData,{1,-1,Round[(Length[refData])/OptionValue[NumberofPoints]]}];
,
refData=refData;
];

If[SameQ[Length[Dimensions[expDataIn]],3],
expData=Select[Flatten[expDataIn,1],(#[[4]]!=-1.)&];
,
expData=Select[expDataIn,(#[[4]]!=-1.)&];
];
expNearest=Nearest[expData[[All,1;;2]]->expData[[All,1;;3]]];
If[OptionValue[NumberofPoints]<Length[expData],
expData=Take[expData,{1,-1,Round[(Length[expData])/OptionValue[NumberofPoints]]}];
,
expData=expData;
];

If[ListQ[OptionValue[ContourRange]],
incr=(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours];
contourList=Sort[DeleteDuplicates[Map[{#,OptionValue[MajorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MajorContours]]]~Join~Map[{#,OptionValue[MinorContourColor]}&,Range[OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]],(OptionValue[ContourRange][[2]]-OptionValue[ContourRange][[1]])/OptionValue[MinorContours]]],#1[[1]]==#2[[1]]&]];
colorF=(ColorData[{OptionValue[color],{OptionValue[ContourRange][[1]],OptionValue[ContourRange][[2]]}}][#1]&);
colorFS=False;
,
contourList=OptionValue[MinorContours];
colorF=OptionValue[color];
colorFS=True;
];

refPlt=ListContourPlot[SortBy[refData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

expPlt=ListContourPlot[SortBy[expData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

Deploy[
Grid[{{
Grid[{{
LocatorPane[
Dynamic[refNodes],
Graphics[{
refPlt[[1]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]},
{LocatorPane[
Dynamic[expNodes],
Graphics[{
expPlt[[1]]
},
PlotRange->All,
AspectRatio->Automatic,
ImageSize->OptionValue[size]
],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.05,0.1,0.4],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]
}]}}]},
{Button["Clear 3D Alignment",FullAlignmentOut={};]},
{Button["Save 3D Alignment",
AppendTo[FullAlignmentOut,{Map[First[refNearest[#]]&,refNodes],Map[First[expNearest[#]]&,expNodes]}];]
}}]]
];


Options[edgeAlign]={
XYmutiplier->1000.
};


edgeAlign[dataIn_,coordsIn_,OptionsPattern[]]:=Module[
{start,edgeV,rot,x,dataOut},

start={x,Fit[coordsIn[[2]], {1,x},x]}/.First[Solve[Fit[coordsIn[[1]], {1,x},x]==Fit[coordsIn[[2]], {1,x},x],x]];
edgeV=First[Differences[SortBy[coordsIn[[2]],Total]]];
rot=RotationMatrix[ArcTan[edgeV[[1]]/edgeV[[2]]]];

dataOut=dataIn;
If[SameQ[Length[Dimensions[dataOut]],2],
dataOut[[All,1;;2]]=(Transpose[rot.Transpose[dataOut[[All,1;;2]]]]-ConstantArray[start,Length[dataOut]])*OptionValue[XYmutiplier];
,
dataOut[[All,All,1;;2]]=(Partition[Transpose[rot.Transpose[Flatten[dataOut[[All,All,1;;2]],1]]],Dimensions[dataOut][[2]]]-ConstantArray[start,Dimensions[dataOut][[1;;2]]])*OptionValue[XYmutiplier];
];

dataOut
];


alignmentParameters[dataIn_,vars_]:=Module[
{ref,data,rot},

If[SameQ[Length[dataIn],1],
ref=dataIn[[1,1]];
data=dataIn[[1,2]];
,
ref=dataIn[[1]];
data=dataIn[[2]];
];

If[SameQ[Length[vars],6],
rot=RotationMatrix[vars[[4]]Degree,{1,0,0}].RotationMatrix[vars[[5]]Degree,{0,1,0}].RotationMatrix[vars[[6]]Degree,{0,0,1}];

data=Transpose[rot.Transpose[data]];
data[[All,1]]=data[[All,1]]-vars[[1]];
data[[All,2]]=data[[All,2]]-vars[[2]];
data[[All,3]]=data[[All,3]]-vars[[3]];
,
rot=RotationMatrix[vars[[3]]Degree];

data=Transpose[rot.Transpose[data]];
data[[All,1]]=data[[All,1]]-vars[[1]];
data[[All,2]]=data[[All,2]]-vars[[2]];
];

Sum[(EuclideanDistance[ref[[i]],data[[i]]]^2),{i,Length[ref]}]
];


sampleAlign[dataIn_,coordsIn_]:=Module[
{xStart,yStart,zStart,theta,thetaX,thetaY,thetaZ,x,y,z,th,thX,thY,thZ,rot,start,dataOut},

If[SameQ[Dimensions[coordsIn][[-1]],3],
{xStart,yStart,zStart,thetaX,thetaY,thetaZ}={x,y,z,thX,thY,thZ}/.Minimize[alignmentParameters[coordsIn,{x,y,z,thX,thY,thZ}],{x,y,z,thX,thY,thZ}][[2]];
rot=RotationMatrix[thetaX Degree,{1,0,0}].RotationMatrix[thetaY Degree,{0,1,0}].RotationMatrix[thetaZ Degree,{0,0,1}];
start={xStart,yStart,zStart};

dataOut=dataIn;
If[SameQ[Length[Dimensions[dataOut]],2],
dataOut[[All,1;;3]]=(Transpose[rot.Transpose[dataOut[[All,1;;3]]]]-ConstantArray[start,Length[dataOut]]);
,
dataOut[[All,All,1;;3]]=(Partition[Transpose[rot.Transpose[Flatten[dataOut[[All,All,1;;3]],1]]],Dimensions[dataOut][[2]]]-ConstantArray[start,Dimensions[dataOut][[1;;2]]]);
];
,
{xStart,yStart,theta}={x,y,th}/.Minimize[alignmentParameters[coordsIn,{x,y,th}],{x,y,th}][[2]];
rot=RotationMatrix[theta Degree];
start={xStart,yStart};

dataOut=dataIn;
If[SameQ[Length[Dimensions[dataOut]],2],
dataOut[[All,1;;2]]=(Transpose[rot.Transpose[dataOut[[All,1;;2]]]]-ConstantArray[start,Length[dataOut]]);
,
dataOut[[All,All,1;;2]]=(Partition[Transpose[rot.Transpose[Flatten[dataOut[[All,All,1;;2]],1]]],Dimensions[dataOut][[2]]]-ConstantArray[start,Dimensions[dataOut][[1;;2]]]);
];
];

dataOut
];


(* ::Subsection::Closed:: *)
(*Selection GUI*)


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


nodes={{0,0},{0,0}};
coordinatesOut={};


Options[getCoordinates]={
ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,plotRange->All,color->"Rainbow",size->750,NumberofPoints->50000
};


getCoordinates[dataIn_,OptionsPattern[]]:=Module[{
takeNumber,reducedData,incr,contourList,colorFS,colorF,contPlt
},

If[OptionValue[NumberofPoints]<Length[dataIn],
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

contPlt=ListContourPlot[SortBy[reducedData[[All,1;;3]],First],
PlotRange->{All,All,OptionValue[ContourRange]},
Contours->contourList,
If[TrueQ[OptionValue[contourShading]==None],
ClippingStyle->None,
ClippingStyle->{ColorData[OptionValue[color]][0],ColorData[OptionValue[color]][1]}
],
ColorFunction->colorF,
ColorFunctionScaling->colorFS,
ContourShading->OptionValue[contourShading]];

Deploy[
Grid[{{
LocatorPane[
Dynamic[nodes],
Graphics[{
contPlt[[1]],
Opacity[1],
GrayLevel[0],Thickness[0.002],
Dynamic[Line[nodes]],
Orange,
Dynamic[Line/@coordinatesOut],
Dynamic[Point/@Flatten[coordinatesOut,1]]
},
PlotRange->OptionValue[plotRange],
AspectRatio->Automatic,
ImageSize->OptionValue[size]],
Appearance->{
Style["\[FilledCircle]",FontColor-> RGBColor[0.1,0.4,0.05],FontSize->12],
Style["\[FilledCircle]",FontColor->RGBColor[0.4,0.05,0.1],FontSize->12]}]},
{Button["Clear coordinatesOut",coordinatesOut={};]},
{Button["Save coordinates",
AppendTo[coordinatesOut,nodes];]}}]]
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["VEECOfunctions`*"],
{
"nodalCoords1","nodalCoords2","AlignmentCoordsOut","BWcoords1","BWcoords2","BWcoordsOut","refNode","expNode","FullAlignmentOut","refNodes","expNodes","WindowAlignmentOut","points","pointsOut","nodes","coordinatesOut"
}];
EndPackage[ ]
