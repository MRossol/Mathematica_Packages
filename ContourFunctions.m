(* ::Package:: *)

BeginPackage["ContourFunctions`"]
Unprotect@@Names["ContourFunctions`*"];
ClearAll@@Names["ContourFunctions`*"];


(*Plot Functions*)
abbrevColors::usage = "abbrevColors2[{colorSytle_,{zMin_,zMax_},frac_}], take fraction frac of colorStyle over the range zMin to zMax.";
colorF::usage = "colorF[{zMin_,zMax_},colors_,OptionsPattern[]] creates color function for z range {zMin,zMax} with # of colors given";
contourList::usage = "contourList[ContourRange_,MajorContours_,MinorContours_,OptionsPattern[]] creates a contour list for the given contourRange, Major and Minor contours";
contourPlot::usage = "contourPlot[dataIn_,OptionsPattern[]] renders a contour plot of dataIn (List) with control over contours";
reduceData::usage = "reduceData[dataIn_,OptionsPattern[]] reduces data if there are more than the given number of points";
sigmaRegionFunction::usage = "sigmaRegionFunction[dataIn_] creates region function from MeanFilter of Sigma values";
frameTicks::usage = "frameTicks[ticksIn_] outputs custom frame ticks with options for tick size";

(*Legends*)
contourLegend::usage = "contourLegend[{{sMin_,sMax_},contours_,lines_,labels_},labelIn_] renders a legend corresponding to the DIC contour plot with identical inputs";



Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Plot Functions*)


abbrevColors[{colorStyle_,{zMin_,zMax_},frac_}]:=Module[
{colorData,colors=100,colorScheme,colorFun},

colorData=ColorData[{colorStyle,{zMin,zMax}}];
colorScheme=Map[colorData[#]&,Range[zMin,zMax,(zMax-zMin)/colors]];

If[Negative[frac],
colorScheme=colorScheme[[;;-IntegerPart[colors-colors*Abs[frac]+1]]];
,
colorScheme=colorScheme[[IntegerPart[colors-colors*Abs[frac]+1];;]];
];

colorFun=Nearest[Range[zMin,zMax,(zMax-zMin)/(colors*Abs[frac])]->colorScheme];
colorFun[#]&
];


Options[colorF]={
color->"Rainbow"
};


colorF[{zMin_,zMax_},colors_,OptionsPattern[]]:=Module[
{colorScheme},
colorScheme=Map[ColorData[{OptionValue[color],{zMin,zMax}}][#]&,Range[zMin,zMax,(zMax-zMin)/colors]];
colorScheme[[First[Nearest[Range[zMin,zMax,(zMax-zMin)/colors]->Automatic,#]]]]&
];


Options[contourList]={
MajorContourColor->Opacity[0.5,Black],MinorContourColor->None
};


contourList[ContourRange_,MinorContours_,MajorContours_,OptionsPattern[]]:=Module[
{},
If[ListQ[ContourRange],
Sort[DeleteDuplicates[Map[{Chop[#],OptionValue[MajorContourColor]}&,Range[ContourRange[[1]],ContourRange[[2]],(ContourRange[[2]]-ContourRange[[1]])/MajorContours]]~Join~Map[{Chop[#],OptionValue[MinorContourColor]}&,Range[ContourRange[[1]],ContourRange[[2]],(ContourRange[[2]]-ContourRange[[1]])/MinorContours]],#1[[1]]==#2[[1]]&]]
,
Automatic
]
];


Options[reduceData]={
NumberofPoints->50000,RemoveBadData->False,SortTake->False
};


reduceData[dataIn_,OptionsPattern[]]:=Module[{
data,takeNumber,reducedData
},

If[Length[Dimensions[dataIn]]==2,
data=dataIn;
,
data=Flatten[dataIn,1];
];

If[OptionValue[RemoveBadData],
data=Select[data,(#[[4]]!= -1.)&];
,
data=data;
];

If[OptionValue[NumberofPoints]<Length[data],
takeNumber=Round[(Length[data])/OptionValue[NumberofPoints]];
If[OptionValue[SortTake],
reducedData=Take[SortBy[data,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=Take[data,{1,-1,takeNumber}];
];
,
reducedData=data;
];

SortBy[reducedData[[All,1;;3]],First]

];


Options[contourPlot]={
ContourRange->Automatic,MinorContours->Automatic,MajorContours->10,MajorContourColor->Opacity[0.5,Black],MinorContourColor->None,contourShading->Automatic,interpolationOrder->None,plotRange->All,color->"Rainbow",size->500,NumberofPoints->None
};


contourPlot[dataIn_,OptionsPattern[]]:=Module[{
takeNumber,reducedData,incr,contourList,colorFS,colorF,contPlt
},

If[NumberQ[OptionValue[NumberofPoints]],
takeNumber=Round[(Length[dataIn])/OptionValue[NumberofPoints]];
reducedData=Take[SortBy[dataIn,(#[[1]]+#[[2]])&],{1,-1,takeNumber}];
,
reducedData=dataIn;
]

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


Options[sigmaRegionFunction]={
FilterSize->1
};


sigmaRegionFunction[dataIn_,OptionsPattern[]]:=Module[{
meanSigma},
 
If[Length[Dimensions[dataIn]]==2,
meanSigma=Nearest[dataIn[[All,1;;2]]->MeanFilter[dataIn[[All,4]],OptionValue[FilterSize]]];
,
meanSigma=Nearest[Flatten[dataIn[[All,All,1;;2]],1]->Flatten[MeanFilter[dataIn[[All,All,4]],OptionValue[FilterSize]],1]];
];

(First[meanSigma[{#1,#2}]]>0.&)

];


Options[frameTicks]={
TickSize->{0.015,0}
};


frameTicks[xTicks_,yTicks_,OptionsPattern[]]:=Module[
{ticks},

If[SameQ[OptionValue[TickSize],Automatic],
ticks[in_]:=Switch[Length[Dimensions[in]],1,{Map[{#,#}&,in],Map[{#,Null}&,in]},2,Map[Map[{#,#}&,#]&,in]];
,
ticks[in_]:=Switch[Length[Dimensions[in]],1,{Map[{#,#,OptionValue[TickSize]}&,in],Map[{#,Null,OptionValue[TickSize]}&,in]},2,Map[Map[{#,#,OptionValue[TickSize]}&,#]&,in]];
];

Map[ticks[#]&,{yTicks,xTicks}]

];


(* ::Subsection::Closed:: *)
(*Legends*)


Options[contourLegend]={
color->"Rainbow",size->{Automatic,500},fontSize->16,textSide->"Right",reverse->False
};


contourLegend[{{sMin_,sMax_},contours_,lines_,labels_},labelIn_,OptionsPattern[]]:=Module[{
Acolor,
incr=N[(sMax-sMin)/contours],
incr2=N[(sMax-sMin)/lines],
incr3=N[(sMax-sMin)/labels],
AColor,rSet,Amesh,myTxt,myLabel
},

If[OptionValue[reverse],
If[StringQ[OptionValue[color]],
Acolor=ColorData[OptionValue[color]][1-#]&/@Rescale[Range[sMin+incr,sMax,incr]];
,
Acolor=First[OptionValue[color][1-#]]&/@Range[sMin+incr,sMax,incr];
];
,
If[StringQ[OptionValue[color]],
Acolor=ColorData[OptionValue[color]][#]&/@Rescale[Range[sMin+incr,sMax,incr]];
,
Acolor=First[OptionValue[color][#]]&/@Range[sMin+incr,sMax,incr];
];
];

rSet=Rectangle[{0,#},{0.5,#+incr 10/(sMax-sMin)}]&/@Drop[10Rescale[Range[sMin,sMax,incr]],-1];
Amesh=Rescale[Range[sMin+incr2,sMax-incr2,incr2],{sMin,sMax},{0,10}];

If[StringMatchQ[OptionValue[textSide],"r*",IgnoreCase->True],
myLabel=Text[Style[labelIn,FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0,10.5},{-1,-1},{1,0}];
myTxt=Text[Style[ToString[Chop[#]],FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0.7,Rescale[#,{sMin,sMax},{0,10}]},{-1,0},{1,0}]&/@Range[sMin,sMax,incr3];
,
myLabel=Text[Style[labelIn,FontFamily->"Arial",FontSize->OptionValue[fontSize]],{0.5,10.5},{1,-1},{1,0}];
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
Protect @@ Complement[Names["ContourFunctions`*"],
{
}];
EndPackage[ ]
