(* ::Package:: *)

BeginPackage["Voxelization`"]
Unprotect@@Names["Voxelization`*"];
ClearAll@@Names["Voxelization`*"];


(*Import 3D GG*)
convertSurfacePolygons::usage = "convertSurfacePolygons[pathIn_] Imports surface polygons from the 3D GG, caps off the tows, outputs {{surfaceNodesX,surfaceElementsX,surfacePolygonsX},{surfaceNodesY,surfaceElementsY,surfacePolygonsY}}";
modelSize::usage = "modelSize[xTowsIn_,yTowsIn_] determines {min,max} for {x,y,z} from the 3D GG surface polygons";

(*Voxelize*)
voxelSize::usage = "voxelSize[{{xMin_,xMax_,xVoxels_},{yMin_,yMax_,yVoxels_},{zMin_,zMax_,zVoxels_}}] calculates total number of elements and the element size for inputs to voxelize";
voxelNumbers::usage = "voxelNumbers[{{xMin_,xMax_,dx_},{yMin_,yMax_,dy_},{zMin_,zMax_,dz_}}] calculates the model size and number of voxels for a specified element size.";
vertexCheck::usage = "vertexCheck[facetsIn_vertexIn_,rayDir_] checks to see if rays passing through vertexes pass through the facet";
facetCheck::usage = "facetCheck[facetsIn_,ray_,rayDir_] checks to see if rays pass through possible facets, check similar to vertex check is performed if an edge is intersected";
facetIntersect::usage = "facetIntersect[facet_,ray_,rayDir_] determines the interesection point of the ray and the facets they intersect";
voxelize::usage = "voxelize[towsIn_,{{xMin_,xMax_,xVoxels_},{yMin_,yMax_,yVoxels_},{zMin_,zMax_,zVoxels_}}] creates a voxel mesh of the given dimensions and voxels in each direction with the imbeded towsIn, output is {3D Mesh Array, Mesh Vectors, Element Center Ranges}";
meshData::usage = "meshData[{gridIn_,nodeVectors_}] extracts {nodes, elementList, elementSets, faceSets}";
voxelCenters::usage = "voxelCenters[{gridIn_,voxelVectors_}] extracts the voxelCenters and groups them by element Set";

(*Calculate Mechanical Properties*)
towProperties::usage = "calculateTowProperties[{Em_,Gm_,\[Nu]m_},{{Ef11_,Ef22_},{Gf12_,Gf23_},{\[Nu]f12_,\[Nu]f23_},{\[Xi]e22_,\[Xi]g12_,\[Xi]g23_},fTow_}] calculates tow properties from the input parameters, outputs {Axial Tow Properties,Matrix Properties}";
rebarProps::usage = "rebarProps[elemDims_,Vf_] calculate A and S for rebar elements oriented in the XZ and YZ plane.";

(*Create Abaqus Files*)
localOrient::usage = "localOrient[lociPath_,nodeList_,elementList_,towElemSets_] determines the local unit vectors {x',y'} for each element in towElemSets.";
rebarOrient::usage = "rebarOrient[lociPath_,nodeListIn_,elementListIn_,towElemSets_,OptionsPattern[]] calculates the local rebar orientation for each element in towElemSets.";
rebarLayers::usage = "rebarData[lociPath_,nodeList_,elementList_,elemSets_,OptionsPattern[]] creates surface elements inside the tow elements and creates new node and element lists as well as new element sets and the local orientation of the surface elements.";
twoPlyRebarLayer::usage = "twoPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,OptionsPattern[]] creates nodeLists, elementLists, Sets, and Orientation for a two ply symmetric RebarLayer model of the single ply given by gridIn, nodeVector, and lociPath.";
fourPlyRebarLayer::usage = "fourPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,OptionsPattern[]] creates nodeLists, elementLists, Sets, and Orientation for a four ply symmetric RebarLayer model with layers 1 and 2 dissimilar to 3 and 4 of the single ply given by gridIn, nodeVector, and lociPath.";
nPlyRebarLayer::usage = "nPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,plies_,OptionsPattern[]] creates nodeLists, elementLists, Sets, and Orientation for an n=plies ply symmetric RebarLayer model of the single ply given by gridIn, nodeVector, and lociPath.";
periodicFaceSets::usage = "periodicFaceSets[faceSetsIn_,nodeListIn_,OptionsPattern[]] removes conflict nodes corresponding to Set1Dir and Set2Dir. Conflict nodes from Set1Plus are removed from Set2Minus and Set2Plus.";
exportAbaqusFiles::usage = "exportAbaqusFiles[pathOut_,{nodeList_,elementList_,elementSets_,faceSets_,localOrient_}] exports the requisite files for the abaqus model: Node List, Element List, Element Sets, Local Orientation Table.";

(*Rebar INP*)
rebarInp::usage = "rebarInp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_] write the abaqus inp file for a continuum mesh reinforced with *rebar.";

(*RebarLayer INP*)
rebarLayerInp::usage = "rebarLayerInp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_] write the abaqus inp file for a continuum mesh reinforced with *rebarLayer.";
rebarLayerP2Inp::usage = "rebarLayerP2Inp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_] write the abaqus inp file for a continuum mesh reinforced with *rebarLayer and periodic boundary conditions at +/- X and Y.";
rebarLayerP3Inp::usage = "rebarLayerP3Inp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_] write the abaqus inp file for a continuum mesh reinforced with *rebarLayer and periodic boundary conditions in all three directions.";

(*MultiPart INP*)
elasticInp::usage = "elasticInp[pathOut_,modelName_,tows_,{towProp_,matrixProp_},modelDims_,disp_] write the abaqus inp file for a linear elastic model with links to the Node List, Element List, Element Sets, Local Orientation files created by exportAbaqusFiles.";
elasticUMATInp::usage = "elasticUMATInp[pathOut_,modelName_,tows_,{towProp_,matrixProp_},disp_] write the abaqus inp file with UMAT for the matrix with links to the Node List, Element List, Element Sets, Local Orientation files created by exportAbaqusFiles.";



Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Import 3D GG*)


convertSurfacePolygons[pathIn_]:=Module[
{surfacePolygonsX,surfacePolygonsY,surfaceNodesX,surfaceNodesY,towEndsX,towEndsY,endCenterX,endCenterY,cap1,cap2,endCapX,endCapY,surfaceElementsX,surfaceElementsY,elements,elementPos},

surfacePolygonsX=Map[Map[Partition[#,3]&,Import[#]]&,FileNames[pathIn<>"xTOWpolygon_*.csv"]];
surfacePolygonsY=Map[Map[Partition[#,3]&,Import[#]]&,FileNames[pathIn<>"yTOWpolygon_*.csv"]];

surfaceNodesX=Map[DeleteDuplicates[Flatten[#,1]]&,surfacePolygonsX];
surfaceNodesY=Map[DeleteDuplicates[Flatten[#,1]]&,surfacePolygonsY];

towEndsX=Map[GatherBy[#,(#[[1]])&][[{1,-1}]]&,surfaceNodesX];
towEndsY=Map[GatherBy[#,(#[[2]])&][[{1,-1}]]&,surfaceNodesY];

endCenterX=Map[Mean[#]&,towEndsX,{2}];
endCenterY=Map[Mean[#]&,towEndsY,{2}];

surfaceNodesX=MapThread[#2[[{1}]]~Join~#1~Join~#2[[{2}]]&,{surfaceNodesX,endCenterX}];
surfaceNodesY=MapThread[#2[[{1}]]~Join~#1~Join~#2[[{2}]]&,{surfaceNodesY,endCenterY}];

endCapX=Table[
cap1=Partition[SortBy[towEndsX[[i,1]],(ArcTan@@(#[[{2,3}]]-endCenterX[[i,1,{2,3}]]))&],2,1,1];
cap1=Join[cap1,ConstantArray[{endCenterX[[i,1]]},Length[cap1]],2];

cap2=Partition[SortBy[towEndsX[[i,2]],(ArcTan@@(#[[{2,3}]]-endCenterX[[i,2,{2,3}]]))&],2,1,1];
cap2=Join[cap2,ConstantArray[{endCenterX[[i,2]]},Length[cap2]],2];

cap1~Join~cap2
,{i,Length[towEndsX]}];
endCapY=Table[
cap1=Partition[SortBy[towEndsY[[i,1]],(ArcTan@@(#[[{1,3}]]-endCenterY[[i,1,{1,3}]]))&],2,1,1];
cap1=Join[cap1,ConstantArray[{endCenterY[[i,1]]},Length[cap1]],2];

cap2=Partition[SortBy[towEndsY[[i,2]],(ArcTan@@(#[[{1,3}]]-endCenterY[[i,2,{1,3}]]))&],2,1,1];
cap2=Join[cap2,ConstantArray[{endCenterY[[i,2]]},Length[cap2]],2];

cap1~Join~cap2
,{i,Length[towEndsY]}];

surfacePolygonsX=Join[surfacePolygonsX,endCapX,2];
surfacePolygonsY=Join[surfacePolygonsY,endCapY,2];

DistributeDefinitions[surfacePolygonsX,surfaceNodesX];
surfaceElementsX=ParallelTable[
elements=surfacePolygonsX[[i]];
ReplaceAll[elements,MapThread[#1->#2&,{surfaceNodesX[[i]],Range[1,Length[surfaceNodesX[[i]]]]}]]
,{i,Length[surfacePolygonsX]}];
DistributeDefinitions[surfacePolygonsY,surfaceNodesY];
surfaceElementsY=ParallelTable[
elements=surfacePolygonsY[[i]];
ReplaceAll[elements,MapThread[#1->#2&,{surfaceNodesY[[i]],Range[1,Length[surfaceNodesY[[i]]]]}]]
,{i,Length[surfacePolygonsY]}];

{{surfaceNodesX,surfaceElementsX,surfacePolygonsX},{surfaceNodesY,surfaceElementsY,surfacePolygonsY}}
];



modelSize[xTowsIn_,yTowsIn_]:=Module[
{allNodes},

allNodes=Flatten[xTowsIn~Join~yTowsIn,2];

Map[{Min[#],Max[#]}&,Transpose[allNodes]]
];


(* ::Subsection::Closed:: *)
(*Voxelization*)


voxelSize[{{xMin_,xMax_,xVoxels_},{yMin_,yMax_,yVoxels_},{zMin_,zMax_,zVoxels_}}]:=Module[{quadMesh,meshSpacing},
quadMesh={{xMin,xMax,xVoxels},{yMin,yMax,yVoxels},{zMin,zMax,zVoxels}};
meshSpacing=Flatten[Map[Differences[#]&,quadMesh[[All,1;;2]]]/quadMesh[[All,3]]];

Print["Number of Elements = "<>ToString[xVoxels*yVoxels*zVoxels]];
Print["Element Size {dx,dy,dz}(mm) = "<>ToString[meshSpacing]];

meshSpacing
];


voxelNumbers[{{xMin_,xMax_,dx_},{yMin_,yMax_,dy_},{zMin_,zMax_,dz_}}]:=Module[{quadMesh,numbofVoxels,quadMeshRange,quadMeshOut},
quadMesh={{xMin,xMax,dx},{yMin,yMax,dy},{zMin,zMax,dz}};
numbofVoxels=Ceiling[Flatten[Map[Differences[#]&,quadMesh[[All,1;;2]]]/quadMesh[[All,3]]]];
quadMeshRange=Map[{Mean[#],Mean[#]}&,quadMesh[[All,1;;2]]]+MapThread[{-1.,1.}*#1*#2/2.&,{quadMesh[[All,3]],numbofVoxels}];
quadMeshOut=MapThread[#1~Join~{#2}&,{quadMeshRange,numbofVoxels}];

Print["Number of Elements = "<>ToString[Total[numbofVoxels]]];
Print["Number of Voxels {x,y,z} = "<>ToString[numbofVoxels]];
Print["Model Range {{xMin, xMax},{yMin,yMax},{zMin,zMax}} = "<>ToString[quadMeshRange]];

quadMeshOut
];


vertexCheck[facetsIn_,vertexIn_,rayDir_]:=Module[{posEdgePos,negEdgePos,ordering,vertexPos,facetSign,facet,edgePos1,neighborPos1,edgePos2,neighborPos2,facetNormals,normDir},

posEdgePos[pos_]:=Switch[pos,1,2,2,3,3,1];
negEdgePos[pos_]:=Switch[pos,1,3,2,1,3,2];
ordering[edgeDiff_]:=Switch[edgeDiff==1||edgeDiff==-2,True,1,False,-1];

vertexPos=Position[facetsIn,vertexIn];
facetSign=ConstantArray[1,Length[vertexPos]];

Do[
facet=facetsIn[[i]];

edgePos1=posEdgePos[vertexPos[[i,2]]];
neighborPos1=First[DeleteCases[Position[facetsIn,facet[[edgePos1]]],{i,edgePos1}]];

facetSign[[neighborPos1[[1]]]]=-1*ordering[(neighborPos1[[2]]-vertexPos[[neighborPos1[[1]],2]])]*facetSign[[i]];

edgePos2=negEdgePos[vertexPos[[i,2]]];
neighborPos2=First[DeleteCases[Position[facetsIn,facet[[edgePos2]]],{i,edgePos2}]];

facetSign[[neighborPos2[[1]]]]=ordering[(neighborPos2[[2]]-vertexPos[[neighborPos2[[1]],2]])]*facetSign[[i]];

,{i,Length[facetSign]}];

Which[rayDir=="x",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn]][[All,1]];
,
rayDir=="y",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn]][[All,2]];
,
rayDir=="z",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn]][[All,3]];
];

normDir=facetNormals*facetSign;

Max[normDir]<0||Min[normDir]>0
];


facetCheck[facetsIn_,ray_,rayDir_]:=Module[{facets,ordering,facet,remainingFacets,edgesOrder,edgeTest,check,edgeOrder1,adjacentFacet,edgeOrder2,facetSign,facetNormals,normDir},

Which[rayDir=="x",
facets=facetsIn[[All,All,{2,3}]];
,
rayDir=="y",
facets=facetsIn[[All,All,{1,3}]];
,
rayDir=="z",
facets=facetsIn[[All,All,{1,2}]];
];

ordering[edgeDiff_]:=Switch[edgeDiff==1||edgeDiff==-2,True,1,False,-1];

Table[
facet=facets[[i]];
remainingFacets=DeleteCases[Range[Length[facetsIn]],i];

edgesOrder={{1,2,3},{2,3,1},{3,1,2}};
edgeTest=Sign[Map[((facet[[#[[3]]]]-facet[[#[[1]]]])-Projection[facet[[#[[3]]]]-facet[[#[[1]]]],facet[[#[[2]]]]-facet[[#[[1]]]]]).((ray-facet[[#[[1]]]])-Projection[ray-facet[[#[[1]]]],facet[[#[[2]]]]-facet[[#[[1]]]]])&,edgesOrder]];

If[Total[edgeTest]==3,
check=True;
,
If[MemberQ[edgeTest,0],

edgeOrder1=edgesOrder[[First[Flatten[Position[edgeTest,0]]],1;;2]];
adjacentFacet=Flatten[Map[Position[facetsIn[[remainingFacets]],#]&,facetsIn[[i,edgeOrder1]]],1];
edgeOrder2=SortBy[GatherBy[adjacentFacet,First],(Length[#])&][[-1,All,2]];

If[Length[edgeOrder2]==2,
facetSign={-1,ordering[First[Differences[edgeOrder2]]]};

Which[rayDir=="x",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn[[{i,remainingFacets[[adjacentFacet[[1,1]]]]}]]]][[All,1]];
,
rayDir=="y",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn[[{i,remainingFacets[[adjacentFacet[[1,1]]]]}]]]][[All,2]];
,
rayDir=="z",
facetNormals=Map[#/Norm[#]&,Map[Cross[#[[2]]-#[[1]],#[[3]]-#[[1]]]&,facetsIn[[{i,remainingFacets[[adjacentFacet[[1,1]]]]}]]]][[All,3]];
];

normDir=facetNormals*facetSign;

check=Max[normDir]<0||Min[normDir]>0;
,
check=False;
];
,
check=False;
];
];
check
,{i,Length[facetsIn]}]
];


facetIntersect[facet_,ray_,rayDir_]:=Module[{A,B,C,D,out},
A=Chop[facet[[1,2]](facet[[2,3]]-facet[[3,3]])+facet[[2,2]](facet[[3,3]]-facet[[1,3]])+facet[[3,2]](facet[[1,3]]-facet[[2,3]])];
B=Chop[facet[[1,3]](facet[[2,1]]-facet[[3,1]])+facet[[2,3]](facet[[3,1]]-facet[[1,1]])+facet[[3,3]](facet[[1,1]]-facet[[2,1]])];
C=Chop[facet[[1,1]](facet[[2,2]]-facet[[3,2]])+facet[[2,1]](facet[[3,2]]-facet[[1,2]])+facet[[3,1]](facet[[1,2]]-facet[[2,2]])];
D=Chop[-facet[[1,1]](facet[[2,2]]facet[[3,3]]-facet[[3,2]]facet[[2,3]])-facet[[2,1]](facet[[3,2]]facet[[1,3]]-facet[[1,2]]facet[[3,3]])-facet[[3,1]](facet[[1,2]]facet[[2,3]]-facet[[2,2]]facet[[1,3]])];

Which[rayDir=="x",
If[A==0,
out=0;
,
out=(-D -B ray[[1]]-C ray[[2]])/A;
];
,
rayDir=="y",
If[B==0,
out=0;
,
out=(-D-A ray[[1]]-C ray[[2]])/B;
];
,
rayDir=="z",
If[C==0,
out=0;
,
out=(-D-A ray[[1]]-B ray[[2]])/C;
];
];

out
];


voxelize[towsIn_,{{xMin_,xMax_,xVoxels_},{yMin_,yMax_,yVoxels_},{zMin_,zMax_,zVoxels_}}]:=Module[
{tows,quadMesh,meshSpacing,meshRanges,elementCenters,logicArray,tow,towMinMax,facetMinMax,towArrayPos,towArray,towCenters,positivePos,ray,possibleFacets,positiveFacetsPos,vertexFacets,vertexFacetPos,vertexes,positiveFacets,zRange,zPos,yRange,yPos,xRange,xPos,voxelsOut,arrayOut},

quadMesh={{xMin,xMax,xVoxels},{yMin,yMax,yVoxels},{zMin,zMax,zVoxels}};

meshSpacing=Map[Differences[#]&,quadMesh[[All,1;;2]]]/quadMesh[[All,3]];
meshRanges=MapThread[Range@@Join[#1,#2]&,{quadMesh[[All,1;;2]],meshSpacing}];
elementCenters=Map[Mean[#]&,Map[Partition[#,2,1]&,meshRanges],{2}];

logicArray=ConstantArray[0,quadMesh[[{3,1,2},3]]];

Which[Length[Dimensions[towsIn]]==1,
tows=towsIn;
,
Length[Dimensions[towsIn]]==2,
tows=towsIn[[1]]~Join~towsIn[[2]];
,
Length[Dimensions[towsIn]]==3,
tows={towsIn};
,
Length[Dimensions[towsIn]]==4,
tows=towsIn;
];

Do[
tow=tows[[t]];
towMinMax=Map[{Min[#],Max[#]}&,Transpose[Flatten[tow,1]]];
facetMinMax=Map[Map[{Min[#],Max[#]}&,Transpose[#]]&,tow];

towArrayPos=MapThread[Clip[First[Position[Sign[#2-#1[[1]]],1]][[1]]-1,{1,#3[[3]]}];;Clip[Last[Position[Sign[#2-#1[[2]]],-1]][[1]]+1,{1,#3[[3]]}]&,{towMinMax,elementCenters,quadMesh}];

towArray=ConstantArray[0,Dimensions[logicArray[[towArrayPos[[3]],towArrayPos[[1]],towArrayPos[[2]]]]]];
towCenters=MapThread[#1[[#2]]&,{elementCenters,towArrayPos}];

DistributeDefinitions[tow,towCenters,facetMinMax,vertexCheck,facetCheck,facetIntersect];

(*xy plane*)
positivePos=Flatten[ParallelTable[
ray={towCenters[[1,x]],towCenters[[2,y]]};

possibleFacets=Position[Map[Boole[(#[[1,1]]<=ray[[1]]<=#[[1,2]]&&#[[2,1]]<=ray[[2]]<=#[[2,2]])]&,facetMinMax],1];
If[Length[possibleFacets]>0,
possibleFacets=Extract[tow,possibleFacets];

vertexFacetPos=Position[Round[possibleFacets[[All,All,{1,2}]],10^(-6.)],Round[ray,10^(-6.)]];
If[Length[vertexFacetPos]>0,
vertexFacets=Extract[possibleFacets,vertexFacetPos[[All,{1}]]];
vertexes=DeleteDuplicates[Extract[possibleFacets,vertexFacetPos],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];

vertexFacetPos=Map[Position[vertexFacets,#]&,vertexes][[All,All,{1}]];
vertexFacets=Map[Extract[vertexFacets,#]&,vertexFacetPos];

positiveFacets=Position[MapThread[vertexCheck[#1,#2,"z"]&,{vertexFacets,vertexes}],True];
positiveFacets=Extract[vertexFacets,positiveFacets][[All,1]];
,
positiveFacets={};
];

positiveFacetsPos=Position[facetCheck[possibleFacets,ray,"z"],True];
If[Length[positiveFacetsPos]>0,
If[Length[positiveFacets]>0,
positiveFacets=positiveFacets~Join~Extract[possibleFacets,positiveFacetsPos];
,
positiveFacets=Extract[possibleFacets,positiveFacetsPos];
];
];

zRange=DeleteDuplicates[Map[facetIntersect[#,ray,"z"]&,positiveFacets],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];
If[(EvenQ[Length[zRange]]&&Length[zRange]>0),
zRange=Partition[Sort[zRange],2];
zPos=Flatten[Table[Map[First[Position[towCenters[[3]],#]]&,Select[towCenters[[3]],(zRange[[i,1]]<#<zRange[[i,2]])&]]
,{i,Length[zRange]}]];
voxelsOut=Map[{#,x,y}&,zPos];
,
voxelsOut={};
];
,
voxelsOut={};
];
voxelsOut
,{x,Length[towCenters[[1]]]},{y,Length[towCenters[[2]]]}],2];
arrayOut=ReplacePart[towArray,positivePos->1];

(*xz plane*)
positivePos=Flatten[ParallelTable[
ray={towCenters[[1,x]],towCenters[[3,z]]};

possibleFacets=Position[Map[Boole[(#[[1,1]]<=ray[[1]]<=#[[1,2]]&&#[[3,1]]<=ray[[2]]<=#[[3,2]])]&,facetMinMax],1];
If[Length[possibleFacets]>0,
possibleFacets=Extract[tow,possibleFacets];

vertexFacetPos=Position[Round[possibleFacets[[All,All,{1,3}]],10^(-6.)],Round[ray,10^(-6.)]];
If[Length[vertexFacetPos]>0,
vertexFacets=Extract[possibleFacets,vertexFacetPos[[All,{1}]]];
vertexes=DeleteDuplicates[Extract[possibleFacets,vertexFacetPos],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];

vertexFacetPos=Map[Position[vertexFacets,#]&,vertexes][[All,All,{1}]];
vertexFacets=Map[Extract[vertexFacets,#]&,vertexFacetPos];

positiveFacets=Position[MapThread[vertexCheck[#1,#2,"y"]&,{vertexFacets,vertexes}],True];
positiveFacets=Extract[vertexFacets,positiveFacets][[All,1]];
,
positiveFacets={};
];

positiveFacetsPos=Position[facetCheck[possibleFacets,ray,"y"],True];
If[Length[positiveFacetsPos]>0,
If[Length[positiveFacets]>0,
positiveFacets=positiveFacets~Join~Extract[possibleFacets,positiveFacetsPos];
,
positiveFacets=Extract[possibleFacets,positiveFacetsPos];
];
];

yRange=DeleteDuplicates[Map[facetIntersect[#,ray,"y"]&,positiveFacets],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];
If[(EvenQ[Length[yRange]]&&Length[yRange]>0),
yRange=Partition[Sort[yRange],2];
yPos=Flatten[Table[Map[First[Position[towCenters[[2]],#]]&,Select[towCenters[[2]],(yRange[[i,1]]<#<yRange[[i,2]])&]]
,{i,Length[yRange]}]];
voxelsOut=Map[{z,x,#}&,yPos];
,
voxelsOut={};
];
,
voxelsOut={};
];
voxelsOut
,{x,Length[towCenters[[1]]]},{z,Length[towCenters[[3]]]}],2];
arrayOut=arrayOut+ReplacePart[towArray,positivePos->1];

(*yz plane*)
positivePos=Flatten[ParallelTable[
ray={towCenters[[2,y]],towCenters[[3,z]]};

possibleFacets=Position[Map[Boole[(#[[2,1]]<=ray[[1]]<=#[[2,2]]&&#[[3,1]]<=ray[[2]]<=#[[3,2]])]&,facetMinMax],1];
If[Length[possibleFacets]>0,
possibleFacets=Extract[tow,possibleFacets];

vertexFacetPos=Position[Round[possibleFacets[[All,All,{2,3}]],10^(-6.)],Round[ray,10^(-6.)]];
If[Length[vertexFacetPos]>0,
vertexFacets=Extract[possibleFacets,vertexFacetPos[[All,{1}]]];
vertexes=DeleteDuplicates[Extract[possibleFacets,vertexFacetPos],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];

vertexFacetPos=Map[Position[vertexFacets,#]&,vertexes][[All,All,{1}]];
vertexFacets=Map[Extract[vertexFacets,#]&,vertexFacetPos];

positiveFacets=Position[MapThread[vertexCheck[#1,#2,"x"]&,{vertexFacets,vertexes}],True];
positiveFacets=Extract[vertexFacets,positiveFacets][[All,1]];
,
positiveFacets={};
];

positiveFacetsPos=Position[facetCheck[possibleFacets,ray,"x"],True];
If[Length[positiveFacetsPos]>0,
If[Length[positiveFacets]>0,
positiveFacets=positiveFacets~Join~Extract[possibleFacets,positiveFacetsPos];
,
positiveFacets=Extract[possibleFacets,positiveFacetsPos];
];
];

xRange=DeleteDuplicates[Map[facetIntersect[#,ray,"x"]&,positiveFacets],Round[#1,10^(-6.)]==Round[#2,10^(-6.)]&];
If[(EvenQ[Length[xRange]]&&Length[xRange]>0),
xRange=Partition[Sort[xRange],2];
xPos=Flatten[Table[Map[First[Position[towCenters[[1]],#]]&,Select[towCenters[[1]],(xRange[[i,1]]<#<xRange[[i,2]])&]]
,{i,Length[xRange]}]];
voxelsOut=Map[{z,#,y}&,xPos];
,
voxelsOut={};
];
,
voxelsOut={};
];
voxelsOut
,{y,Length[towCenters[[2]]]},{z,Length[towCenters[[3]]]}],2];
arrayOut=arrayOut+ReplacePart[towArray,positivePos->1];

arrayOut=ReplaceAll[arrayOut,{1->0,2->t,3->t}];
logicArray[[towArrayPos[[3]],towArrayPos[[1]],towArrayPos[[2]]]]=logicArray[[towArrayPos[[3]],towArrayPos[[1]],towArrayPos[[2]]]]+arrayOut;

Print["Voxels In Tow "<>ToString[t]<>" = "<>ToString[Length[Position[Flatten[arrayOut],t]]]]
,{t,Length[tows]}];

{logicArray,meshRanges,elementCenters}
];


meshData[{gridIn_,nodeVectors_}]:=Module[
{dx,dy,dz,nodesPos,nodes,nodeList,elementOrder,partPos,elementSets,elementPos,elementList,nodeGrid,xSet,ySet,zSet,faceSets},

{dx,dy,dz}=Map[Length[#]&,nodeVectors];

nodesPos=Flatten[Array[{#2,#3,#1}&,{dz,dx,dy},1],2];
nodes=nodesPos;
nodes[[All,1]]=nodes[[All,1]]/.MapIndexed[First[#2]->#1&,nodeVectors[[1]]];
nodes[[All,2]]=nodes[[All,2]]/.MapIndexed[First[#2]->#1&,nodeVectors[[2]]];
nodes[[All,3]]=nodes[[All,3]]/.MapIndexed[First[#2]->#1&,nodeVectors[[3]]];
nodeList=Join[Partition[Range[Length[nodes]],1],nodes,2];

partPos=Flatten[gridIn,2];
partPos=Transpose[{Range[Length[partPos]],partPos}];
elementSets=GatherBy[SortBy[partPos,(#[[2]])&],(#[[2]])&][[All,All,1]];

elementPos=Flatten[Array[{#2,#3,#1}&,Dimensions[gridIn]],2];
elementOrder={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
elementList=Map[ConstantArray[#,8]+elementOrder&,elementPos];
elementList=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,elementList,{2}];
elementList=Join[Partition[Range[Length[elementList]],1],elementList,2];

nodeGrid=Array[{#2,#3,#1}&,{dz,dx,dy}];
xSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,1,All]],1],Flatten[nodeGrid[[All,-1,All]],1]},{2}];
ySet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,All,1]],1],Flatten[nodeGrid[[All,All,-1]],1]},{2}];
zSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[1,All,All]],1],Flatten[nodeGrid[[-1,All,All]],1]},{2}];
faceSets={xSet,ySet,zSet};

{nodeList,elementList,elementSets,faceSets}
];


voxelCenters[{gridIn_,voxelVectors_}]:=Module[
{elementCenters},

elementCenters=Flatten[MapIndexed[#2[[{2,3,1}]]~Join~{#1}&,gridIn,{3}],2];
elementCenters[[All,1]]=elementCenters[[All,1]]/.MapIndexed[First[#2]->#1&,voxelVectors[[1]]];
elementCenters[[All,2]]=elementCenters[[All,2]]/.MapIndexed[First[#2]->#1&,voxelVectors[[2]]];
elementCenters[[All,3]]=elementCenters[[All,3]]/.MapIndexed[First[#2]->#1&,voxelVectors[[3]]];


GatherBy[SortBy[elementCenters,(#[[4]])&],(#[[4]])&]
];


(* ::Subsection::Closed:: *)
(*Calculate Mechanical Properties*)


towProperties[{E1_,G1_,\[Nu]1_},{{Ea2_,Et2_},{Ga2_,Gt2_},{va2_,vt2_},c2_}]:=Module[
{Et11,Et22,Gt12,Gt23,x,\[Nu]t12,\[Nu]t21,\[Nu]t23,AxialTowProperties,TransverseTowProperties,MatrixProperties,c1,k1,k2,Ga1,Gt1,Ea1,Et1,va1,vt1,eta1,eta2,Gtr,Achr,Bchr,Cchr,sols,Eac,vac,vac2,Gac,Gtc,vtc,Etc,kc,mc},

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

{AxialTowProperties,MatrixProperties}

];


rebarProps[elemDims_,Vf_]:=Module[
{rebarA,rebarT,rebarS},

rebarA=DeleteDuplicates[{elemDims[[2]]*elemDims[[3]]*Vf,elemDims[[1]]*elemDims[[3]]*Vf}];
rebarT=DeleteDuplicates[{elemDims[[2]]*Vf,elemDims[[1]]*Vf}];
rebarS=rebarA/rebarT;

If[Length[rebarS]==1,
{First[rebarA],First[rebarS]}
,
{rebarA,rebarS}
]

];


(* ::Subsection::Closed:: *)
(*Create Abaqus Files*)


Options[localOrient]={
NumberofTows->Automatic
};


localOrient[lociPath_,nodeListIn_,elementListIn_,towElemSets_,OptionsPattern[]]:=Module[
{towLoci,nodeList,elementList,xTows,xTowLoci,xTowInterps,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,elemDims,xTowOrient,interp,elemCenters,elemNumbers,delta,yTowOrient},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];
nodeList=nodeListIn[[All,2;;]];
elementList=elementListIn[[All,2;;]];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]/2];
];
];

xTowLoci=towLoci[[1;;xTows,All,1;;3]];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowLoci=towLoci[[xTows+1;;,All,1;;3]];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets=towElemSets[[1;;xTows]];
xTowElemCenters=Map[Map[Mean[nodeList[[#]]]&,Extract[elementList,Partition[#,1]]]&,xTowElemSets];
yTowElemSets=towElemSets[[xTows+1;;]];
yTowElemCenters=Map[Map[Mean[nodeList[[#]]]&,Extract[elementList,Partition[#,1]]]&,yTowElemSets];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[Flatten[{#2,{Chop[{delta,0,interp'[#1]*delta}],Chop[{0.,-1.,0.}-Projection[{0.,-1.,0.},Chop[{delta,0,interp'[#1]*delta}]]]}}]&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}];

delta=elemDims[[2]]/2.;
yTowOrient=Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[Flatten[{#2,{Chop[{0,delta,interp'[#1]*delta}],Chop[{1.,0.,0.}-Projection[{1.,0.,0.},Chop[{0,delta,interp'[#1]*delta}]]]}}]&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}];

SortBy[Flatten[xTowOrient~Join~yTowOrient,1],First]
];


Options[rebarOrient]={
NumberofTows->Automatic,Tilt->"Rebar"
};


rebarOrient[lociPath_,nodeListIn_,elementListIn_,towElemSets_,OptionsPattern[]]:=Module[
{towLoci,nodeList,elementList,xTows,xTowElems,yTowElems,xTowNodes,yTowNodes,nodeNearest,xTowElemList,yTowElemList,xTowSets,yTowSets,xSet,ySet,zSet,newFaceSets,xTowLoci,xTowInterps,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,elemDims,xTowOrient,interp,elemCenters,elemNumbers,delta,xTowAngleSets,xTowAngles,xEdgeFracs,xTowFracs,yTowOrient,yTowAngleSets,yTowAngles,yEdgeFracs,yTowFracs,isoAngle},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];
nodeList=nodeListIn[[All,2;;]];
elementList=elementListIn[[All,2;;]];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]/2];
];
];

xTowLoci=towLoci[[;;xTows,All,1;;3]];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowLoci=towLoci[[xTows+1;;,All,1;;3]];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets=towElemSets[[1;;xTows]];
xTowElemCenters=Map[Map[Mean[nodeList[[#]]]&,Extract[elementList,Partition[#,1]]]&,xTowElemSets];
yTowElemSets=towElemSets[[xTows+1;;]];
yTowElemCenters=Map[Map[Mean[nodeList[[#]]]&,Extract[elementList,Partition[#,1]]]&,yTowElemSets];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=GatherBy[SortBy[Flatten[Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,0.1]}&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}],1],Last],Last];
xTowAngleSets=xTowOrient[[All,All,1]];
xTowAngles=Union[Flatten[xTowOrient[[All,All,2]]]];

delta=elemDims[[2]]/2.;
yTowOrient=GatherBy[SortBy[Flatten[Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,0.1]}&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}],1],Last],Last];
yTowAngleSets=yTowOrient[[All,All,1]];
yTowAngles=Union[Flatten[yTowOrient[[All,All,2]]]];

If[StringMatchQ[OptionValue[Tilt],"r*",IgnoreCase->True],
isoAngle[beta_]:=ArcTan[elemDims[[1]]/elemDims[[3]]Tan[beta Degree]]/Degree;
xTowAngles=Map[isoAngle[#]&,xTowAngles];
isoAngle[beta_]:=ArcTan[elemDims[[2]]/elemDims[[3]]Tan[beta Degree]]/Degree;
yTowAngles=Map[isoAngle[#]&,yTowAngles];

xTowOrient={xTowAngleSets,xTowAngles};
yTowOrient={yTowAngleSets,yTowAngles};
,
xEdgeFracs[theta_]:=Switch[(Tan[theta Degree]elemDims[[1]]/2.)<elemDims[[3]]/2.&&(Tan[theta Degree]elemDims[[1]]/2.)>-elemDims[[3]]/2.,True,{Round[(elemDims[[3]]/2.-Tan[theta Degree]elemDims[[1]]/2.)/elemDims[[3]],10^(-6.)],0,Round[(elemDims[[3]]/2.-Tan[theta Degree]elemDims[[1]]/2.)/elemDims[[3]],10^(-6.)],0},False,{0,Round[(elemDims[[1]]/2.+Tan[(90-theta )Degree]elemDims[[3]]/2.)/elemDims[[1]],10^(-6.)],0,Round[(elemDims[[1]]/2.+Tan[(90-theta ) Degree]elemDims[[3]]/2.)/elemDims[[1]],10^(-6.)]}];
xTowFracs=Map[xEdgeFracs[#]&,xTowAngles];
xTowOrient={xTowAngleSets,xTowFracs};

yEdgeFracs[theta_]:=Switch[(Tan[theta Degree]elemDims[[2]]/2.)<elemDims[[3]]/2.&&(Tan[theta Degree]elemDims[[2]]/2.)>-elemDims[[3]]/2.,True,{0,Round[(elemDims[[3]]/2.+Tan[theta Degree]elemDims[[2]]/2.)/elemDims[[3]],10^(-6.)],0,Round[(elemDims[[3]]/2.+Tan[theta Degree]elemDims[[2]]/2.)/elemDims[[3]],10^(-6.)]},False,{Round[(elemDims[[2]]/2.-Tan[(90-theta )Degree]elemDims[[3]]/2.)/elemDims[[2]],10^(-6.)],0,Round[(elemDims[[2]]/2.-Tan[(90-theta ) Degree]elemDims[[3]]/2.)/elemDims[[2]],10^(-6.)],0}];
yTowFracs=Map[yEdgeFracs[#]&,yTowAngles];
yTowOrient={yTowAngleSets,yTowFracs};
];

{xTowOrient,yTowOrient}
];


Options[rebarLayers]={
NumberofTows->Automatic,RoundTo->0.1
};


rebarLayers[lociPath_,nodeListIn_,elementListIn_,elemSets_,OptionsPattern[]]:=Module[
{towLoci,nodeList,elementList,xTows,xTowElems,yTowElems,xTowNodes,yTowNodes,nodeNearest,xTowElemList,yTowElemList,xTowSets,yTowSets,xSet,ySet,zSet,newFaceSets,xTowLoci,xTowInterps,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,xTowOrient,interp,elemCenters,elemNumbers,delta,xTowOrientSets,xTowOrients,yTowOrient,yTowOrientSets,yTowOrients,elemDims},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];
nodeList=nodeListIn[[All,2;;]];
elementList=elementListIn[[All,2;;]];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]/2];
];
];

xTowElems=Map[Extract[nodeList,Partition[#,1]]&,Map[Extract[elementList,Partition[#,1]]&,elemSets[[2;;1+xTows]]],{2}];
yTowElems=Map[Extract[nodeList,Partition[#,1]]&,Map[Extract[elementList,Partition[#,1]]&,elemSets[[2+xTows;;]]],{2}];

xTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{1,3}]])&]&,xTowElems,{2}],{3}][[All,All,{1,2,4,3}]];
yTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{2,3}]])&]&,yTowElems,{2}],{3}][[All,All,{1,2,4,3}]];

xTowNodes=DeleteDuplicates[Flatten[xTowElems,2]];
yTowNodes=DeleteDuplicates[Flatten[yTowElems,2]];

nodeNearest=Nearest[xTowNodes->Range[Length[xTowNodes]]];
xTowElemList=Map[First[nodeNearest[#]]&,Flatten[xTowElems,1],{2}];
nodeNearest=Nearest[yTowNodes->Range[Length[yTowNodes]]];
yTowElemList=Map[First[nodeNearest[#]]&,Flatten[yTowElems,1],{2}];

xTowNodes=Join[Partition[Range[Length[xTowNodes]],1],xTowNodes,2];
yTowNodes=Join[Partition[Range[Length[yTowNodes]],1],yTowNodes,2];

xTowElemList=Join[Partition[Range[Length[xTowElemList]],1],xTowElemList,2];
yTowElemList=Join[Partition[Range[Length[yTowElemList]],1],yTowElemList,2];

xTowSets=MapThread[(Range[Length[#1]]+#2)&,{xTowElems,Accumulate[{0}~Join~Map[Length[#]&,xTowElems]][[;;-2]]}];
yTowSets=MapThread[(Range[Length[#1]]+#2)&,{yTowElems,Accumulate[{0}~Join~Map[Length[#]&,yTowElems]][[;;-2]]}];

xTowLoci=towLoci[[;;xTows,All,1;;3]];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowLoci=towLoci[[xTows+1;;,All,1;;3]];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets= xTowSets;
xTowElemCenters=Map[Mean[#]&,xTowElems,{2}];
yTowElemSets=yTowSets;
yTowElemCenters=Map[Mean[#]&,yTowElems,{2}];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=GatherBy[SortBy[Flatten[Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}],1],Last],Last];
xTowOrientSets=xTowOrient[[All,All,1]];
xTowOrients=Union[Flatten[xTowOrient[[All,All,2;;]]]];
xTowOrients=Map[Flatten[{RotationMatrix[# Degree,{0.,-1.,0.}].{1.,0.,0.},RotationMatrix[# Degree,{0.,-1.,0.}].{0.,0.,1.}}]&,xTowOrients];

delta=elemDims[[2]]/2.;
yTowOrient=GatherBy[SortBy[Flatten[Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}],1],Last],Last];
yTowOrientSets=yTowOrient[[All,All,1]];
yTowOrients=Union[Flatten[yTowOrient[[All,All,2;;]]]];
yTowOrients=Map[Flatten[{RotationMatrix[# Degree,{1.,0.,0.}].{0.,1.,0.},RotationMatrix[# Degree,{1.,0.,0.}].{0.,0.,1.}}]&,yTowOrients];

{{nodeListIn,{xTowNodes,yTowNodes}},{elementListIn,{xTowElemList,yTowElemList}},{elemSets,{xTowSets,yTowSets}},{{xTowOrientSets,xTowOrients},{yTowOrientSets,yTowOrients}}}
];


Options[twoPlyRebarLayer]={
NumberofTows->Automatic,RoundTo->0.1
};


twoPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,OptionsPattern[]]:=Module[
{towLoci,xTows,replace1,replace2,twoPlyMesh,twoPlyVectors,zStep,dz,zShift,dx,dy,nodesPos,nodes,nodeList,partPos,elementSets,elementPos,elementOrder,elements,elementList,nodeGrid,xSet,ySet,zSet,faceSets,xTowElems,yTowElems,xTowNodes,yTowNodes,nodeNearest,xTowElemList,yTowElemList,xTowSets,yTowSets,xTowLoci1,xTowLoci2,xTowLoci,xTowInterps,yTowLoci1,yTowLoci2,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,elemDims,delta,xTowOrient,interp,elemCenters,elemNumbers,xTowOrientSets,xTowOrients,yTowOrient,yTowOrientSets,yTowOrients},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]];
];
];

replace1=MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows/2+1,xTows],Range[xTows/2+1,xTows]+xTows/2}];
replace2=MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows],Range[xTows/2+1,xTows]~Join~(Range[xTows/2+1,xTows]+xTows)}];

twoPlyMesh=Reverse[Join[Reverse[gridIn]/.replace2,Reverse[gridIn,2]/.replace1]];

twoPlyVectors=nodeVectors;
zStep=Mean[Differences[twoPlyVectors[[3]]]];
dz=Length[twoPlyMesh];
twoPlyVectors[[3]]=Range[-zStep*dz/2,zStep*dz/2,zStep];
zShift=twoPlyVectors[[3,{1,-1}]]-nodeVectors[[3,{1,-1}]];

{dx,dy,dz}=Map[Length[#]&,twoPlyVectors];

nodesPos=Flatten[Array[{#2,#3,#1}&,{dz,dx,dy},1],2];
nodes=nodesPos;
nodes[[All,1]]=nodes[[All,1]]/.MapIndexed[First[#2]->#1&,twoPlyVectors[[1]]];
nodes[[All,2]]=nodes[[All,2]]/.MapIndexed[First[#2]->#1&,twoPlyVectors[[2]]];
nodes[[All,3]]=nodes[[All,3]]/.MapIndexed[First[#2]->#1&,twoPlyVectors[[3]]];
nodeList=Join[Partition[Range[Length[nodes]],1],nodes,2];

partPos=Flatten[twoPlyMesh,2];
partPos=Transpose[{Range[Length[partPos]],partPos}];
elementSets=GatherBy[SortBy[partPos,(#[[2]])&],(#[[2]])&][[All,All,1]];

elementPos=Flatten[Array[{#2,#3,#1}&,Dimensions[twoPlyMesh]],2];
elementOrder={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
elements=Map[ConstantArray[#,8]+elementOrder&,elementPos];
elements=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,elements,{2}];
elementList=Join[Partition[Range[Length[elements]],1],elements,2];

nodeGrid=Array[{#2,#3,#1}&,{dz,dx,dy}];
xSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,1,All]],1],Flatten[nodeGrid[[All,-1,All]],1]},{2}];
ySet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,All,1]],1],Flatten[nodeGrid[[All,All,-1]],1]},{2}];
zSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[1,All,All]],1],Flatten[nodeGrid[[-1,All,All]],1]},{2}];
faceSets={xSet,ySet,zSet};

xTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[2;;1+xTows]]],{2}];
yTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[2+xTows;;]]],{2}];

xTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{1,3}]])&]&,xTowElems,{2}],{3}][[All,All,{1,2,4,3}]];
yTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{2,3}]])&]&,yTowElems,{2}],{3}][[All,All,{1,2,4,3}]];

xTowNodes=DeleteDuplicates[Flatten[xTowElems,2]];
yTowNodes=DeleteDuplicates[Flatten[yTowElems,2]];

nodeNearest=Nearest[xTowNodes->Range[Length[xTowNodes]]];
xTowElemList=Map[First[nodeNearest[#]]&,Flatten[xTowElems,1],{2}];
nodeNearest=Nearest[yTowNodes->Range[Length[yTowNodes]]];
yTowElemList=Map[First[nodeNearest[#]]&,Flatten[yTowElems,1],{2}];

xTowNodes=Join[Partition[Range[Length[xTowNodes]],1],xTowNodes,2];
yTowNodes=Join[Partition[Range[Length[yTowNodes]],1],yTowNodes,2];

xTowElemList=Join[Partition[Range[Length[xTowElemList]],1],xTowElemList,2];
yTowElemList=Join[Partition[Range[Length[yTowElemList]],1],yTowElemList,2];

xTowSets=MapThread[(Range[Length[#1]]+#2)&,{xTowElems,Accumulate[{0}~Join~Map[Length[#]&,xTowElems]][[;;-2]]}];
yTowSets=MapThread[(Range[Length[#1]]+#2)&,{yTowElems,Accumulate[{0}~Join~Map[Length[#]&,yTowElems]][[;;-2]]}];

xTowLoci1=towLoci[[;;xTows/2,All,1;;3]];
xTowLoci2=towLoci[[;;xTows/2,All,1;;3]];
xTowLoci1[[All,All,3]]=-xTowLoci1[[All,All,3]]+zShift[[1]];
xTowLoci2[[All,All,3]]=xTowLoci2[[All,All,3]]+zShift[[2]];
yTowLoci1=towLoci[[xTows/2+1;;,All,1;;3]];
yTowLoci2=towLoci[[xTows/2+1;;,All,1;;3]];
yTowLoci1[[All,All,3]]=-yTowLoci1[[All,All,3]]+zShift[[1]];
yTowLoci2[[All,All,3]]=yTowLoci2[[All,All,3]]+zShift[[2]];

xTowLoci1[[All,All,1]]=-xTowLoci1[[All,All,1]]+Max[xTowLoci1[[All,All,1]]];
yTowLoci1[[All,All,1]]=-yTowLoci1[[All,All,1]]+Max[xTowLoci1[[All,All,1]]];

xTowLoci=Join[xTowLoci1,xTowLoci2];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowLoci=Join[yTowLoci1,yTowLoci2];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets=xTowSets;
xTowElemCenters=Map[Mean[#]&,xTowElems,{2}];
yTowElemSets=yTowSets;
yTowElemCenters=Map[Mean[#]&,yTowElems,{2}];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=GatherBy[SortBy[Flatten[Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}],1],Last],Last];
xTowOrientSets=xTowOrient[[All,All,1]];
xTowOrients=Union[Flatten[xTowOrient[[All,All,2;;]]]];
xTowOrients=Map[Flatten[{RotationMatrix[# Degree,{0.,-1.,0.}].{1.,0.,0.},RotationMatrix[# Degree,{0.,-1.,0.}].{0.,0.,1.}}]&,xTowOrients];

delta=elemDims[[2]]/2.;
yTowOrient=GatherBy[SortBy[Flatten[Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}],1],Last],Last];
yTowOrientSets=yTowOrient[[All,All,1]];
yTowOrients=Union[Flatten[yTowOrient[[All,All,2;;]]]];
yTowOrients=Map[Flatten[{RotationMatrix[# Degree,{1.,0.,0.}].{0.,1.,0.},RotationMatrix[# Degree,{1.,0.,0.}].{0.,0.,1.}}]&,yTowOrients];

{{nodeList,{xTowNodes,yTowNodes}},{elementList,{xTowElemList,yTowElemList}},{elementSets,{xTowSets,yTowSets}},faceSets,{{xTowOrientSets,xTowOrients},{yTowOrientSets,yTowOrients}}}

];


Options[fourPlyRebarLayer]={
NumberofTows->Automatic,RoundTo->0.1
};


fourPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,OptionsPattern[]]:=Module[
{towLoci,xTows,replace,twoPlyMesh1,replace2,twoPlyMesh2,fourPlyMesh,xTowPos,yTowPos,fourPlyVectors,zStep,dz,zShifts,dx,dy,nodesPos,nodes,nodeList,partPos,elementSets,elementPos,elementOrder,elements,elementList,nodeGrid,xSet,ySet,zSet,faceSets,xTowElems,yTowElems,xTowNodes,yTowNodes,nodeNearest,xTowElemList,yTowElemList,xTowSets,yTowSets,xPly,xPly1,xPly2,xPly3,xPly4,xTowLoci,xTowInterps,yPly,yPly1,yPly2,yPly3,yPly4,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,elemDims,delta,xTowOrient,interp,elemCenters,elemNumbers,xTowOrientSets,xTowOrients,yTowOrient,yTowOrientSets,yTowOrients},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]/2];
];
];

replace=MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows*2],Range[xTows*2]+xTows*2}];
twoPlyMesh1=Reverse[Join[Reverse[gridIn]/.replace,Reverse[gridIn,2]]];

replace2=MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows*4],Range[xTows*4]+xTows*4}];
twoPlyMesh2=Reverse[twoPlyMesh1,{2,3}]/.replace2;
fourPlyMesh=twoPlyMesh2~Join~twoPlyMesh1;

xTowPos=IntegerPart[Flatten[Table[Range[xTows]+(i-1)*2*xTows+1,{i,4}]]];
yTowPos=IntegerPart[Flatten[Table[Range[xTows+1,xTows*2]+(i-1)*2*xTows+1,{i,4}]]];

fourPlyVectors=nodeVectors;
zStep=Mean[Differences[fourPlyVectors[[3]]]];
dz=Length[fourPlyMesh];
fourPlyVectors[[3]]=Range[-zStep*dz/2,zStep*dz/2,zStep];
zShifts=(Reverse[-1*Range[1,4,2]]~Join~Range[1,4,2]*1/2)*First[Differences[nodeVectors[[3,{1,-1}]]]];

{dx,dy,dz}=Map[Length[#]&,fourPlyVectors];

nodesPos=Flatten[Array[{#2,#3,#1}&,{dz,dx,dy},1],2];
nodes=nodesPos;
nodes[[All,1]]=nodes[[All,1]]/.MapIndexed[First[#2]->#1&,fourPlyVectors[[1]]];
nodes[[All,2]]=nodes[[All,2]]/.MapIndexed[First[#2]->#1&,fourPlyVectors[[2]]];
nodes[[All,3]]=nodes[[All,3]]/.MapIndexed[First[#2]->#1&,fourPlyVectors[[3]]];
nodeList=Join[Partition[Range[Length[nodes]],1],nodes,2];

partPos=Flatten[fourPlyMesh,2];
partPos=Transpose[{Range[Length[partPos]],partPos}];
elementSets=GatherBy[SortBy[partPos,(#[[2]])&],(#[[2]])&][[All,All,1]];

elementPos=Flatten[Array[{#2,#3,#1}&,Dimensions[fourPlyMesh]],2];
elementOrder={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
elements=Map[ConstantArray[#,8]+elementOrder&,elementPos];
elements=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,elements,{2}];
elementList=Join[Partition[Range[Length[elements]],1],elements,2];

nodeGrid=Array[{#2,#3,#1}&,{dz,dx,dy}];
xSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,1,All]],1],Flatten[nodeGrid[[All,-1,All]],1]},{2}];
ySet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,All,1]],1],Flatten[nodeGrid[[All,All,-1]],1]},{2}];
zSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[1,All,All]],1],Flatten[nodeGrid[[-1,All,All]],1]},{2}];
faceSets={xSet,ySet,zSet};

xTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[xTowPos]]],{2}];
yTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[yTowPos]]],{2}];

xTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{1,3}]])&]&,xTowElems,{2}],{3}][[All,All,{1,2,4,3}]];
yTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{2,3}]])&]&,yTowElems,{2}],{3}][[All,All,{1,2,4,3}]];

xTowNodes=DeleteDuplicates[Flatten[xTowElems,2]];
yTowNodes=DeleteDuplicates[Flatten[yTowElems,2]];

nodeNearest=Nearest[xTowNodes->Range[Length[xTowNodes]]];
xTowElemList=Map[First[nodeNearest[#]]&,Flatten[xTowElems,1],{2}];
nodeNearest=Nearest[yTowNodes->Range[Length[yTowNodes]]];
yTowElemList=Map[First[nodeNearest[#]]&,Flatten[yTowElems,1],{2}];

xTowNodes=Join[Partition[Range[Length[xTowNodes]],1],xTowNodes,2];
yTowNodes=Join[Partition[Range[Length[yTowNodes]],1],yTowNodes,2];

xTowElemList=Join[Partition[Range[Length[xTowElemList]],1],xTowElemList,2];
yTowElemList=Join[Partition[Range[Length[yTowElemList]],1],yTowElemList,2];

xTowSets=MapThread[(Range[Length[#1]]+#2)&,{xTowElems,Accumulate[{0}~Join~Map[Length[#]&,xTowElems]][[;;-2]]}];
yTowSets=MapThread[(Range[Length[#1]]+#2)&,{yTowElems,Accumulate[{0}~Join~Map[Length[#]&,yTowElems]][[;;-2]]}];

xPly=towLoci[[;;xTows,All,1;;3]];
yPly=towLoci[[xTows+1;;,All,1;;3]];
xPly1=xPly;
xPly1[[All,All,3]]=-xPly1[[All,All,3]]+zShifts[[3]];
xPly1[[All,All,1]]=-xPly1[[All,All,1]]+Max[xPly[[All,All,1]]];
yPly1=yPly;
yPly1[[All,All,3]]=-yPly1[[All,All,3]]+zShifts[[3]];
yPly1[[All,All,1]]=-yPly1[[All,All,1]]+Max[xPly[[All,All,1]]];

xPly2=xPly;
xPly2[[All,All,3]]=xPly2[[All,All,3]]+zShifts[[4]];
yPly2=yPly;
yPly2[[All,All,3]]=yPly2[[All,All,3]]+zShifts[[4]];

xPly3=xPly;
xPly3[[All,All,3]]=-xPly3[[All,All,3]]+zShifts[[1]];
xPly3[[All,All,2]]=-xPly3[[All,All,2]]+Max[yPly[[All,All,2]]];
yPly3=yPly;
yPly3[[All,All,3]]=-yPly3[[All,All,3]]+zShifts[[1]];
yPly3[[All,All,2]]=-yPly3[[All,All,2]]+Max[yPly[[All,All,2]]];

xPly4=xPly;
xPly4[[All,All,3]]=xPly4[[All,All,3]]+zShifts[[2]];
xPly4[[All,All,2]]=-xPly4[[All,All,2]]+Max[yPly[[All,All,2]]];
xPly4[[All,All,1]]=-xPly4[[All,All,1]]+Max[xPly[[All,All,1]]];
yPly4=yPly;
yPly4[[All,All,3]]=yPly4[[All,All,3]]+zShifts[[2]];
yPly4[[All,All,2]]=-yPly4[[All,All,2]]+Max[yPly[[All,All,2]]];
yPly4[[All,All,1]]=-yPly4[[All,All,1]]+Max[xPly[[All,All,1]]];

xTowLoci=Join[xPly1,xPly2,xPly3,xPly4];
yTowLoci=Join[yPly1,yPly2,yPly3,yPly4];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets=xTowSets;
xTowElemCenters=Map[Mean[#]&,xTowElems,{2}];
yTowElemSets=yTowSets;
yTowElemCenters=Map[Mean[#]&,yTowElems,{2}];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=GatherBy[SortBy[Flatten[Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}],1],Last],Last];
xTowOrientSets=xTowOrient[[All,All,1]];
xTowOrients=Union[Flatten[xTowOrient[[All,All,2;;]]]];
xTowOrients=Map[Flatten[{RotationMatrix[# Degree,{0.,-1.,0.}].{1.,0.,0.},RotationMatrix[# Degree,{0.,-1.,0.}].{0.,0.,1.}}]&,xTowOrients];

delta=elemDims[[2]]/2.;
yTowOrient=GatherBy[SortBy[Flatten[Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}],1],Last],Last];
yTowOrientSets=yTowOrient[[All,All,1]];
yTowOrients=Union[Flatten[yTowOrient[[All,All,2;;]]]];
yTowOrients=Map[Flatten[{RotationMatrix[# Degree,{1.,0.,0.}].{0.,1.,0.},RotationMatrix[# Degree,{1.,0.,0.}].{0.,0.,1.}}]&,yTowOrients];

{{nodeList,{xTowNodes,yTowNodes}},{elementList,{xTowElemList,yTowElemList}},{elementSets,{xTowSets,yTowSets}},faceSets,{{xTowOrientSets,xTowOrients},{yTowOrientSets,yTowOrients}}}

];


Options[nPlyRebarLayer]={
NumberofTows->Automatic,RoundTo->0.1
};


nPlyRebarLayer[gridIn_,nodeVectors_,lociPath_,plies_,OptionsPattern[]]:=Module[
{towLoci,xTows,nPlies,replace,twoPlyMesh,nPlyMesh,xTowPos,yTowPos,nPlyVectors,zStep,dz,zShifts,dx,dy,nodesPos,nodes,nodeList,partPos,elementSets,elementPos,elementOrder,elements,elementList,nodeGrid,xSet,ySet,zSet,faceSets,xTowElems,yTowElems,xTowNodes,yTowNodes,nodeNearest,xTowElemList,yTowElemList,xTowSets,yTowSets,xyTowLoci,xPly,xTowLoci,xTowInterps,yPly,yTowLoci,yTowInterps,xTowElemSets,xTowElemCenters,yTowElemSets,yTowElemCenters,elemDims,delta,xTowOrient,interp,elemCenters,elemNumbers,xTowOrientSets,xTowOrients,yTowOrient,yTowOrientSets,yTowOrients},

towLoci=Import[lociPath<>"allWEFTWARPsm_LociArea.xls"];

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
xTows=IntegerPart[Length[towLoci]/2];
];
];

If[OddQ[plies],
nPlies=IntegerPart[(plies-1)/2];
,
nPlies=IntegerPart[plies/2];
];

replace=MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows*2],Range[xTows*2]+xTows*2}];
twoPlyMesh=Reverse[Join[Reverse[gridIn]/.replace,Reverse[gridIn,2]]];

nPlyMesh=Join@@Table[twoPlyMesh/.MapThread[IntegerPart[#1]->IntegerPart[#2]&,{Range[xTows*4],Range[xTows*4]+4*xTows*(i-1)}],{i,nPlies}];

xTowPos=IntegerPart[Flatten[Table[Range[xTows]+(i-1)*2*xTows+1,{i,nPlies*2}]]];
yTowPos=IntegerPart[Flatten[Table[Range[xTows+1,xTows*2]+(i-1)*2*xTows+1,{i,nPlies*2}]]];

nPlyVectors=nodeVectors;
zStep=Mean[Differences[nPlyVectors[[3]]]];
dz=Length[nPlyMesh];
nPlyVectors[[3]]=Range[-zStep*dz/2,zStep*dz/2,zStep];
zShifts=(Reverse[-1*Range[1,nPlies*2,2]]~Join~Range[1,nPlies*2,2]*1/2)*First[Differences[nodeVectors[[3,{1,-1}]]]];

{dx,dy,dz}=Map[Length[#]&,nPlyVectors];

nodesPos=Flatten[Array[{#2,#3,#1}&,{dz,dx,dy},1],2];
nodes=nodesPos;
nodes[[All,1]]=nodes[[All,1]]/.MapIndexed[First[#2]->#1&,nPlyVectors[[1]]];
nodes[[All,2]]=nodes[[All,2]]/.MapIndexed[First[#2]->#1&,nPlyVectors[[2]]];
nodes[[All,3]]=nodes[[All,3]]/.MapIndexed[First[#2]->#1&,nPlyVectors[[3]]];
nodeList=Join[Partition[Range[Length[nodes]],1],nodes,2];

partPos=Flatten[nPlyMesh,2];
partPos=Transpose[{Range[Length[partPos]],partPos}];
elementSets=GatherBy[SortBy[partPos,(#[[2]])&],(#[[2]])&][[All,All,1]];

elementPos=Flatten[Array[{#2,#3,#1}&,Dimensions[nPlyMesh]],2];
elementOrder={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
elements=Map[ConstantArray[#,8]+elementOrder&,elementPos];
elements=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,elements,{2}];
elementList=Join[Partition[Range[Length[elements]],1],elements,2];

nodeGrid=Array[{#2,#3,#1}&,{dz,dx,dy}];
xSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,1,All]],1],Flatten[nodeGrid[[All,-1,All]],1]},{2}];
ySet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[All,All,1]],1],Flatten[nodeGrid[[All,All,-1]],1]},{2}];
zSet=Map[IntegerPart[(#[[1]]-1)dy+#[[2]]+(#[[3]]-1)dx dy]&,{Flatten[nodeGrid[[1,All,All]],1],Flatten[nodeGrid[[-1,All,All]],1]},{2}];
faceSets={xSet,ySet,zSet};

xTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[xTowPos]]],{2}];
yTowElems=Map[Extract[nodes,Partition[#,1]]&,Map[Extract[elements,Partition[#,1]]&,elementSets[[yTowPos]]],{2}];

xTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{1,3}]])&]&,xTowElems,{2}],{3}][[All,All,{1,2,4,3}]];
yTowElems=Map[Mean[#]&,Map[GatherBy[#,(#[[{2,3}]])&]&,yTowElems,{2}],{3}][[All,All,{1,2,4,3}]];

xTowNodes=DeleteDuplicates[Flatten[xTowElems,2]];
yTowNodes=DeleteDuplicates[Flatten[yTowElems,2]];

nodeNearest=Nearest[xTowNodes->Range[Length[xTowNodes]]];
xTowElemList=Map[First[nodeNearest[#]]&,Flatten[xTowElems,1],{2}];
nodeNearest=Nearest[yTowNodes->Range[Length[yTowNodes]]];
yTowElemList=Map[First[nodeNearest[#]]&,Flatten[yTowElems,1],{2}];

xTowNodes=Join[Partition[Range[Length[xTowNodes]],1],xTowNodes,2];
yTowNodes=Join[Partition[Range[Length[yTowNodes]],1],yTowNodes,2];

xTowElemList=Join[Partition[Range[Length[xTowElemList]],1],xTowElemList,2];
yTowElemList=Join[Partition[Range[Length[yTowElemList]],1],yTowElemList,2];

xTowSets=MapThread[(Range[Length[#1]]+#2)&,{xTowElems,Accumulate[{0}~Join~Map[Length[#]&,xTowElems]][[;;-2]]}];
yTowSets=MapThread[(Range[Length[#1]]+#2)&,{yTowElems,Accumulate[{0}~Join~Map[Length[#]&,yTowElems]][[;;-2]]}];

xyTowLoci=Table[xPly=towLoci[[;;xTows,All,1;;3]];
yPly=towLoci[[xTows+1;;,All,1;;3]];
If[OddQ[i],
xPly[[All,All,3]]=-xPly[[All,All,3]]+zShifts[[i]];
yPly[[All,All,3]]=-yPly[[All,All,3]]+zShifts[[i]];
xPly[[All,All,1]]=-xPly[[All,All,1]]+Max[xPly[[All,All,1]]];
yPly[[All,All,1]]=-yPly[[All,All,1]]+Max[xPly[[All,All,1]]];
,
xPly[[All,All,3]]=xPly[[All,All,3]]+zShifts[[i]];
yPly[[All,All,3]]=yPly[[All,All,3]]+zShifts[[i]];
];
{xPly,yPly}
,{i,nPlies*2}];
xTowLoci=Flatten[xyTowLoci[[All,1]],1];
yTowLoci=Flatten[xyTowLoci[[All,2]],1];
xTowInterps=Map[Interpolation[#[[All,{1,3}]]]&,xTowLoci];
yTowInterps=Map[Interpolation[#[[All,{2,3}]]]&,yTowLoci];

xTowElemSets=xTowSets;
xTowElemCenters=Map[Mean[#]&,xTowElems,{2}];
yTowElemSets=yTowSets;
yTowElemCenters=Map[Mean[#]&,yTowElems,{2}];

elemDims=Map[First[DeleteDuplicates[Differences[Union[#]],(Round[#1,10^(-6.)]==Round[#2,10^(-6.)])&]]&,Transpose[Flatten[xTowElemCenters~Join~yTowElemCenters,1]]];

delta=elemDims[[1]]/2.;
xTowOrient=GatherBy[SortBy[Flatten[Table[
interp[x_]:=xTowInterps[[i]][x];
elemCenters=xTowElemCenters[[i]];
elemNumbers=xTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,1]],elemNumbers}]

,{i,Length[xTowLoci]}],1],Last],Last];
xTowOrientSets=xTowOrient[[All,All,1]];
xTowOrients=Union[Flatten[xTowOrient[[All,All,2;;]]]];
xTowOrients=Map[Flatten[{RotationMatrix[# Degree,{0.,-1.,0.}].{1.,0.,0.},RotationMatrix[# Degree,{0.,-1.,0.}].{0.,0.,1.}}]&,xTowOrients];

delta=elemDims[[2]]/2.;
yTowOrient=GatherBy[SortBy[Flatten[Table[
interp[y_]:=yTowInterps[[i]][y];
elemCenters=yTowElemCenters[[i]];
elemNumbers=yTowElemSets[[i]];

MapThread[{#2,Round[ArcTan[delta,interp'[#1]*delta]/Degree,OptionValue[RoundTo]]}&,{elemCenters[[All,2]],elemNumbers}]

,{i,Length[yTowLoci]}],1],Last],Last];
yTowOrientSets=yTowOrient[[All,All,1]];
yTowOrients=Union[Flatten[yTowOrient[[All,All,2;;]]]];
yTowOrients=Map[Flatten[{RotationMatrix[# Degree,{1.,0.,0.}].{0.,1.,0.},RotationMatrix[# Degree,{1.,0.,0.}].{0.,0.,1.}}]&,yTowOrients];

{{nodeList,{xTowNodes,yTowNodes}},{elementList,{xTowElemList,yTowElemList}},{elementSets,{xTowSets,yTowSets}},faceSets,{{xTowOrientSets,xTowOrients},{yTowOrientSets,yTowOrients}}}

];


periodicFaceSets[faceSetsIn_,nodeListIn_]:=Module[
{periodicSets,xPlusNodes,yMinusNodes,yPlusNodes,zMinusNodes,zPlusNodes,badNodes},

periodicSets=faceSetsIn;

xPlusNodes=Extract[nodeListIn[[All,2;;]],Partition[periodicSets[[1,2]],1]];
yMinusNodes=Extract[nodeListIn[[All,2;;]],Partition[periodicSets[[2,1]],1]];
yPlusNodes=Extract[nodeListIn[[All,2;;]],Partition[periodicSets[[2,2]],1]];
zMinusNodes=Extract[nodeListIn[[All,2;;]],Partition[periodicSets[[3,1]],1]];
zPlusNodes=Extract[nodeListIn[[All,2;;]],Partition[periodicSets[[3,2]],1]];

badNodes=Map[First[Position[yMinusNodes,#]]&,Intersection[xPlusNodes,yMinusNodes]];
periodicSets[[2,1]]=Delete[periodicSets[[2,1]],badNodes];
badNodes=Map[First[Position[yPlusNodes,#]]&,Intersection[xPlusNodes,yPlusNodes]];
periodicSets[[2,2]]=Delete[periodicSets[[2,2]],badNodes];
badNodes=DeleteDuplicates[Map[First[Position[zMinusNodes,#]]&,Intersection[xPlusNodes,zMinusNodes]]~Join~Map[First[Position[zMinusNodes,#]]&,Intersection[yPlusNodes,zMinusNodes]]];
periodicSets[[3,1]]=Delete[periodicSets[[3,1]],badNodes];
badNodes=DeleteDuplicates[Map[First[Position[zPlusNodes,#]]&,Intersection[xPlusNodes,zPlusNodes]]~Join~Map[First[Position[zPlusNodes,#]]&,Intersection[yPlusNodes,zPlusNodes]]];
periodicSets[[3,2]]=Delete[periodicSets[[3,2]],badNodes];

{faceSetsIn,periodicSets}

];


Options[exportAbaqusFiles]={
NumberofTows->Automatic
};


exportAbaqusFiles[pathOut_,{nodeList_,elementList_,elementSets_,faceSets_,localOrient_},OptionsPattern[]]:=Module[
{xTows,continuumSets,rebarSets,xOrient,yOrient},

If[Length[OptionValue[NumberofTows]]==2,
xTows=IntegerPart[OptionValue[NumberofTows][[1]]];
,
If[NumberQ[OptionValue[NumberofTows]],
xTows=IntegerPart[OptionValue[NumberofTows]/2];
,
If[Length[elementSets]==2,
xTows=IntegerPart[Length[elementSets[[1,2;;]]]/2];
,
xTows=IntegerPart[Length[elementSets[[2;;]]]/2];
];
];
];

If[SameQ[DirectoryQ[pathOut],False],CreateDirectory[pathOut]];
If[SameQ[DirectoryQ[pathOut<>"\\InpFiles"],False],CreateDirectory[pathOut<>"\\InpFiles"]];
If[SameQ[DirectoryQ[pathOut<>"\\InpFiles\\sets"],False],CreateDirectory[pathOut<>"\\InpFiles\\sets"]];

(*Nodes*)
If[Length[nodeList]==2,
Export[pathOut<>"\\InpFiles\\ContinuumNodes.csv",nodeList[[1]]];

Export[pathOut<>"\\InpFiles\\xTowsNodes.csv",nodeList[[2,1]]];
Export[pathOut<>"\\InpFiles\\yTowsNodes.csv",nodeList[[2,2]]];
,
Export[pathOut<>"\\InpFiles\\nodes.csv",nodeList];
];

(*Elements*)
If[Length[elementList]==2,
Export[pathOut<>"\\InpFiles\\ContinuumElements.csv",elementList[[1]]];

Export[pathOut<>"\\InpFiles\\xTowsElements.csv",elementList[[2,1]]];
Export[pathOut<>"\\InpFiles\\yTowsElements.csv",elementList[[2,2]]];
,
Export[pathOut<>"\\InpFiles\\elements.csv",elementList];
];

(*Element Sets*)
If[Length[elementSets]==2,
continuumSets=elementSets[[1]];
rebarSets=elementSets[[2]];

Export[pathOut<>"\\InpFiles\\sets\\elem_continuumSet.csv",Flatten[continuumSets]];
Export[pathOut<>"\\InpFiles\\sets\\elem_continuumTowsSet.csv",Flatten[continuumSets[[2;;]]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_xContinuumTowsSet.csv",Flatten[continuumSets[[2;;1+xTows]]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_yContinuumTowsSet.csv",Flatten[continuumSets[[2+xTows;;]]]];
MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_xContinuumTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,continuumSets[[2;;1+xTows]]];
MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_yContinuumTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,continuumSets[[2+xTows;;]]];

Export[pathOut<>"\\InpFiles\\sets\\elem_xTowsSet.csv",Flatten[rebarSets[[1]]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_yTowsSet.csv",Flatten[rebarSets[[2]]]];
MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_xTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,rebarSets[[1]]];
MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_yTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,rebarSets[[2]]];
,
Export[pathOut<>"\\InpFiles\\sets\\elem_matrixSet.csv",elementSets[[1]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_towsSet.csv",Flatten[elementSets[[2;;]]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_xTowsSet.csv",Flatten[elementSets[[2;;1+xTows]]]];
Export[pathOut<>"\\InpFiles\\sets\\elem_yTowsSet.csv",Flatten[elementSets[[2+xTows;;]]]];

MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_xTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,elementSets[[2;;1+xTows]]];
MapIndexed[Export[pathOut<>"\\InpFiles\\sets\\elem_yTowSet"<>ToString[NumberForm[First[#2],{2,0},NumberPadding->{"0",""}]]<>"csv",#1]&,elementSets[[2+xTows;;]]];
];

(*Face Sets*)
If[Length[faceSets]==2,
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_xSet"<>#2<>".csv",#1]&,{faceSets[[1,1]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_ySet"<>#2<>".csv",#1]&,{faceSets[[1,2]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_zSet"<>#2<>".csv",#1]&,{faceSets[[1,3]],{"Minus","Plus"}}];

MapThread[Export[pathOut<>"\\InpFiles\\sets\\periodic_xSet"<>#2<>".csv",#1]&,{faceSets[[2,1]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\periodic_ySet"<>#2<>".csv",#1]&,{faceSets[[2,2]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\periodic_zSet"<>#2<>".csv",#1]&,{faceSets[[2,3]],{"Minus","Plus"}}];
,
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_xSet"<>#2<>".csv",#1]&,{faceSets[[1]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_ySet"<>#2<>".csv",#1]&,{faceSets[[2]],{"Minus","Plus"}}];
MapThread[Export[pathOut<>"\\InpFiles\\sets\\node_zSet"<>#2<>".csv",#1]&,{faceSets[[3]],{"Minus","Plus"}}];
];
If[Length[nodeList]==2,
Export[pathOut<>"\\InpFiles\\sets\\node_cornerSet.csv",First[Flatten[Position[nodeList[[1]],Map[Min[#]&,Transpose[nodeList[[1]]]]]]]];
,
Export[pathOut<>"\\InpFiles\\sets\\node_cornerSet.csv",First[Flatten[Position[nodeList,Map[Min[#]&,Transpose[nodeList]]]]]];
];

(*Orientation*)
If[Length[localOrient]==2,
xOrient=localOrient[[1]];
If[Length[xOrient]==2,
If[SameQ[DirectoryQ[pathOut<>"\\InpFiles\\Orient"],False],CreateDirectory[pathOut<>"\\InpFiles\\Orient"]];
MapIndexed[Export[pathOut<>"\\InpFiles\\Orient\\elem_xTowsOrientSet"<>ToString[First[#2]]<>".csv",#1]&,xOrient[[1]]];
Export[pathOut<>"\\InpFiles\\Orient\\xTowsOrient.MAT",xOrient[[2]]];
,
Export[pathOut<>"\\InpFiles\\xTowsOrientation.csv",Prepend[xOrient,{Null,1.,0.,0.,0.,1.,0.}]];
];
yOrient=localOrient[[2]];

If[Length[yOrient]==2,
If[SameQ[DirectoryQ[pathOut<>"\\InpFiles\\Orient"],False],CreateDirectory[pathOut<>"\\InpFiles\\Orient"]];
MapIndexed[Export[pathOut<>"\\InpFiles\\Orient\\elem_yTowsOrientSet"<>ToString[First[#2]]<>".csv",#1]&,yOrient[[1]]];
Export[pathOut<>"\\InpFiles\\Orient\\yTowsOrient.MAT",yOrient[[2]]];
,
Export[pathOut<>"\\InpFiles\\yTowsOrientation.csv",Prepend[yOrient,{Null,1.,0.,0.,0.,1.,0.}]];
];
,
Export[pathOut<>"\\InpFiles\\LocalOrientation.csv",Prepend[localOrient,{Null,1.,0.,0.,0.,1.,0.}]];
];

Print["Files Exported to "<>pathOut]
];


(* ::Subsection::Closed:: *)
(*Rebar INP*)


rebarInp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_]:=Module[
{xTows,yTows,fiberProp,xRebarA,yRebarA,xRebarS,yRebarS,elemT,xOrient,yOrient,constants,inp},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

fiberProp=rebarProp[[1]];
If[ListQ[rebarProp[[2,1]]],
{xRebarA,yRebarA}=rebarProp[[2,1]];
,
{xRebarA,yRebarA}={rebarProp[[2,1]],rebarProp[[2,1]]};
];
If[ListQ[rebarProp[[2,2]]],
{xRebarS,yRebarS}=rebarProp[[2,2]];
,
{xRebarS,yRebarS}={rebarProp[[2,2]],rebarProp[[2,2]]};
];
xOrient=rebarProp[[3,1]];
yOrient=rebarProp[[3,2]];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Continuum",
"*Node, input=InpFiles\\nodes.csv", 
"*Element, type=C3D8R, elset=elem_Continuum",
"*include, input=InpFiles\\elements.csv",
"**",
"*Elset,elset=elem_xTows",
"*include, input=InpFiles\\sets\\elem_xTowsSet.csv",
"**",
"*Elset,elset=elem_yTows",
"*include, input=InpFiles\\sets\\elem_yTowsSet.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_xTowsOrientSet"<>ToString[x],
"*include, input=InpFiles\\Orient\\elem_xTowsOrientSet"<>ToString[x]<>".csv"
}
,{x,Length[xOrient]}]]
,
Flatten[Table[{
"**",
"*Elset, elset=elem_yTowsOrientSet"<>ToString[y],
"*include, input=InpFiles\\Orient\\elem_yTowsOrientSet"<>ToString[y]<>".csv"
}
,{y,Length[yOrient]}]]
,
{
"**",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Continuum, material=MATRIX",
"**",
"*Rebar, Element=Continuum, Material=FIBERS, Geometry=Isoparametric, Name=xTowsRebar"
}
,
Flatten[MapIndexed[{
"elem_xTowsOrientSet"<>ToString[First[#2]]<>", "<>ToString[xRebarA]<>", "<>ToString[xRebarS]<>", "<>ToString[Minus[N[#1]]]<>", 0.5, 1, 3"
}&,xOrient]]
,
{
"**",
"*Rebar, Element=Continuum, Material=FIBERS, Geometry=Isoparametric, Name=yTowsRebar"
}
,
Flatten[MapIndexed[{
"elem_yTowsOrientSet"<>ToString[First[#2]]<>", "<>ToString[yRebarA]<>", "<>ToString[yRebarS]<>", "<>ToString[N[#1]]<>", 0.5, 2, 3"}&
,yOrient]]
,
{
"**",
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Continuum-1, part=Continuum",
"*End Instance",
"**",
"*Nset, nset=xMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*End Assembly",
"**",
"** MATERIALS",
"**",
"*Material, name=MATRIX",
"*Elastic",
StringTrim[ToString[N[continuumProp]],("{"|"}")],
"**",
"*Material, name=FIBERS",
"*Elastic",
StringTrim[ToString[N[fiberProp]],("{"|"}")],
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: PinCorner Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"corner, ENCASTRE",
"**",
"** Name: YSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"yMinus, YSYMM",
"**",
"** Name: ZSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"zMinus, ZSYMM",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=100000, NLGEOM=YES",
"*Static, stabilize",
"0.0001, 1., 0.00001, 0.1",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"yPlus, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: F-Output-1",
"**",
"*Output, field, variable=PRESELECT, time interval = 0.100000",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


(* ::Subsection::Closed:: *)
(*RebarLayer INP*)


rebarLayerInp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},disp_]:=Module[
{xTows,yTows,fiberProp,xRebarA,yRebarA,xRebarS,yRebarS,elemT,xOrient,yOrient,constants,inp},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

fiberProp=rebarProp[[1]];
If[ListQ[rebarProp[[2,1]]],
{xRebarA,yRebarA}=rebarProp[[2,1]];
,
{xRebarA,yRebarA}={rebarProp[[2,1]],rebarProp[[2,1]]};
];
If[ListQ[rebarProp[[2,2]]],
{xRebarS,yRebarS}=rebarProp[[2,2]];
,
{xRebarS,yRebarS}={rebarProp[[2,2]],rebarProp[[2,2]]};
];
xOrient=rebarProp[[3,1]];
yOrient=rebarProp[[3,2]];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Continuum",
"*Node, nset=nodes_Continuum",
"*include, input=InpFiles\\ContinuumNodes.csv", 
"*Element, type=C3D8R, elset=elem_Continuum",
"*include, input=InpFiles\\ContinuumElements.csv",
"**",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Continuum, material=MATRIX",
"*EndPart",
"**",
"*Part, name=xTows",
"*Node, input=InpFiles\\xTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_xTows",
"*include, input=InpFiles\\xTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_xTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_xTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[xOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=xTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,xOrient]]
,
Flatten[Table[{
"**",
"** Section: xTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_xTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=xTowOrient_"<>ToString[i],
"xTowsRebar_"<>ToString[i]<>", "<>ToString[xRebarA]<>", "<>ToString[xRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[xOrient]}]]
,
{
"*EndPart",
"**",
"*Part, name=yTows",
"*Node, input=InpFiles\\yTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_yTows",
"*include, input=InpFiles\\yTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_yTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_yTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[yOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=yTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,yOrient]]
,
Flatten[Table[{
"**",
"** Section: yTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_yTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=yTowOrient_"<>ToString[i],
"yTowsRebar_"<>ToString[i]<>", "<>ToString[yRebarA]<>", "<>ToString[yRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[yOrient]}]]
,
{
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Continuum-1, part=Continuum",
"*End Instance",
"**",
"*Instance, name=xTows-1, part=xTows",
"*End Instance",
"**",
"*Instance, name=yTows-1, part=yTows",
"*End Instance",
"**",
"*Nset, nset=xMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*Elset, elset=elem_xContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_xContinuumTowsSet.csv",
"**",
"*Elset, elset=elem_yContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_yContinuumTowsSet.csv",
"**",
"*Embedded Element, host elset=Continuum-1.elem_Continuum",
"xTows-1.elem_xTows, yTows-1.elem_yTows",
"*End Assembly",
"**",
"** MATERIALS",
"**",
"*Material, name=MATRIX",
"*Elastic",
StringTrim[ToString[N[continuumProp]],("{"|"}")],
"**",
"*Material, name=FIBERS",
"*Elastic",
StringTrim[ToString[N[fiberProp]],("{"|"}")],
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: PinCorner Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"corner, ENCASTRE",
"**",
"** Name: YSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"yMinus, YSYMM",
"**",
"** Name: ZSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"zMinus, ZSYMM",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=100000, NLGEOM=YES",
"*Static, stabilize",
"0.0001, 1., 0.00001, 0.1",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"yPlus, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: Continuum",
"**",
"*Output, field, time interval=0.05",
"*Node Output, nset=Continuum-1.nodes_Continuum",
"RF, U",
"*Element Output, elset=Continuum-1.elem_Continuum, directions=YES",
"NE, LE S",
"**",
"** FIELD OUTPUT: xTowRebar",
"**",
"*Element Output, elset=xTows-1.elem_xTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** FIELD OUTPUT: yTowRebar",
"**",
"*Element Output, elset=yTows-1.elem_yTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


rebarLayerP2Inp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},modelDims_,disp_]:=Module[
{xTows,yTows,fiberProp,xRebarA,yRebarA,xRebarS,yRebarS,elemT,xOrient,yOrient,constants,inp},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

fiberProp=rebarProp[[1]];
If[ListQ[rebarProp[[2,1]]],
{xRebarA,yRebarA}=rebarProp[[2,1]];
,
{xRebarA,yRebarA}={rebarProp[[2,1]],rebarProp[[2,1]]};
];
If[ListQ[rebarProp[[2,2]]],
{xRebarS,yRebarS}=rebarProp[[2,2]];
,
{xRebarS,yRebarS}={rebarProp[[2,2]],rebarProp[[2,2]]};
];
xOrient=rebarProp[[3,1]];
yOrient=rebarProp[[3,2]];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Fict",
"*Node",
"1, 1.000000, 1.000000, 1.000000,",
"2, 1.000000, 1.000000, 1.000000,",
"3, 1.000000, 1.000000, 1.000000,",
"*End Part",
"**",
"*Part, name=Continuum",
"*Node, nset=nodes_Continuum",
"*include, input=InpFiles\\ContinuumNodes.csv", 
"*Element, type=C3D8R, elset=elem_Continuum",
"*include, input=InpFiles\\ContinuumElements.csv",
"**",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Continuum, material=MATRIX",
"*EndPart",
"**",
"*Part, name=xTows",
"*Node, input=InpFiles\\xTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_xTows",
"*include, input=InpFiles\\xTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_xTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_xTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[xOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=xTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,xOrient]]
,
Flatten[Table[{
"**",
"** Section: xTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_xTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=xTowOrient_"<>ToString[i],
"xTowsRebar_"<>ToString[i]<>", "<>ToString[xRebarA]<>", "<>ToString[xRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[xOrient]}]]
,
{
"*EndPart",
"**",
"*Part, name=yTows",
"*Node, input=InpFiles\\yTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_yTows",
"*include, input=InpFiles\\yTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_yTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_yTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[yOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=yTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,yOrient]]
,
Flatten[Table[{
"**",
"** Section: yTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_yTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=yTowOrient_"<>ToString[i],
"yTowsRebar_"<>ToString[i]<>", "<>ToString[yRebarA]<>", "<>ToString[yRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[yOrient]}]]
,
{
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Fict-1, part=Fict",
"*End Instance",
"**",
"*Instance, name=Continuum-1, part=Continuum",
"*End Instance",
"**",
"*Instance, name=xTows-1, part=xTows",
"*End Instance",
"**",
"*Instance, name=yTows-1, part=yTows",
"*End Instance",
"**",
"*Nset, nset=f1, instance=Fict-1",
"1",
"**",
"*Nset, nset=f2, instance=Fict-1",
"2",
"**",
"*Nset, nset=xMinus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*Elset, elset=elem_xContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_xContinuumTowsSet.csv",
"**",
"*Elset, elset=elem_yContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_yContinuumTowsSet.csv",
"**",
"*Embedded Element, host elset=Continuum-1.elem_Continuum",
"xTows-1.elem_xTows, yTows-1.elem_yTows",
"**",
"*Equation",
"3",
"xPlus, 1, 1.000000, xMinus, 1, -1.000000, f1, 1, -1.000000,",
"*Equation",
"3",
"xPlus, 2, 1.000000, xMinus, 2, -1.000000, f1, 2, -1.000000,",
"*Equation",
"3",
"xPlus, 3, 1.000000, xMinus, 3, -1.000000, f1, 3, -1.000000,",
"**",
"*Equation",
"3",
"yPlus, 1, 1.000000, yMinus, 1, -1.000000, f2, 1, -1.000000,",
"*Equation",
"3",
"yPlus, 2, 1.000000, yMinus, 2, -1.000000, f2, 2, -1.000000,",
"*Equation",
"3",
"yPlus, 3, 1.000000, yMinus, 3, -1.000000, f2, 3, -1.000000,",
"**",
"*Equation",
"2",
"f1, 2, "<>ToString[modelDims[[2]]]<>", f2, 1, "<>ToString[-1.*modelDims[[1]]]<>",",
"**",
"*End Assembly",
"**",
"** MATERIALS",
"**",
"*Material, name=MATRIX",
"*Elastic",
StringTrim[ToString[N[continuumProp]],("{"|"}")],
"**",
"*Material, name=FIBERS",
"*Elastic",
StringTrim[ToString[N[fiberProp]],("{"|"}")],
"**",
"** BOUNDARY CONDITIONS",
"**",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=1000000, NLGEOM=YES",
"*Static, stabilize",
"0.0001, 1., 0.0000001, 0.1",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"f2, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: Continuum",
"**",
"*Output, field, time interval=0.05",
"*Node Output, nset=Continuum-1.nodes_Continuum",
"RF, U",
"*Element Output, elset=Continuum-1.elem_Continuum, directions=YES",
"NE, LE, S",
"**",
"** FIELD OUTPUT: Fict",
"*Node Output, nset=f2",
"RF, U",
"**",
"** FIELD OUTPUT: xTowRebar",
"**",
"*Element Output, elset=xTows-1.elem_xTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** FIELD OUTPUT: yTowRebar",
"**",
"*Element Output, elset=yTows-1.elem_yTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


rebarLayerP3Inp[pathOut_,modelName_,tows_,{rebarProp_,continuumProp_},modelDims_,disp_]:=Module[
{xTows,yTows,fiberProp,xRebarA,yRebarA,xRebarS,yRebarS,elemT,xOrient,yOrient,constants,inp},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

fiberProp=rebarProp[[1]];
If[ListQ[rebarProp[[2,1]]],
{xRebarA,yRebarA}=rebarProp[[2,1]];
,
{xRebarA,yRebarA}={rebarProp[[2,1]],rebarProp[[2,1]]};
];
If[ListQ[rebarProp[[2,2]]],
{xRebarS,yRebarS}=rebarProp[[2,2]];
,
{xRebarS,yRebarS}={rebarProp[[2,2]],rebarProp[[2,2]]};
];
xOrient=rebarProp[[3,1]];
yOrient=rebarProp[[3,2]];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Fict",
"*Node",
"1, 1.000000, 1.000000, 1.000000,",
"2, 1.000000, 1.000000, 1.000000,",
"3, 1.000000, 1.000000, 1.000000,",
"*End Part",
"**",
"*Part, name=Continuum",
"*Node, nset=nodes_Continuum",
"*include, input=InpFiles\\ContinuumNodes.csv", 
"*Element, type=C3D8R, elset=elem_Continuum",
"*include, input=InpFiles\\ContinuumElements.csv",
"**",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Continuum, material=MATRIX",
"*EndPart",
"**",
"*Part, name=xTows",
"*Node, input=InpFiles\\xTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_xTows",
"*include, input=InpFiles\\xTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_xTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_xTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[xOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=xTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,xOrient]]
,
Flatten[Table[{
"**",
"** Section: xTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_xTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=xTowOrient_"<>ToString[i],
"xTowsRebar_"<>ToString[i]<>", "<>ToString[xRebarA]<>", "<>ToString[xRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[xOrient]}]]
,
{
"*EndPart",
"**",
"*Part, name=yTows",
"*Node, input=InpFiles\\yTowsNodes.csv", 
"*Element, type=SFM3D4R, elset=elem_yTows",
"*include, input=InpFiles\\yTowsElements.csv"
}
,
Flatten[Table[{
"**",
"*Elset, elset=elem_yTowsOrientSet"<>ToString[i],
"*include, input=InpFiles\\Orient\\elem_yTowsOrientSet"<>ToString[i]<>".csv"
}
,{i,Length[yOrient]}]]
,
Flatten[MapIndexed[{
"**",
"*Orientation, System=Rectangular, Name=yTowOrient_"<>ToString[First[#2]],
StringTrim[ToString[N[#1]],("{"|"}")],
"3, 0.0"
}&
,yOrient]]
,
Flatten[Table[{
"**",
"** Section: yTowsRebar_"<>ToString[i],
"*Surface Section, elset=elem_yTowsOrientSet"<>ToString[i],
"*Rebar Layer, Orientation=yTowOrient_"<>ToString[i],
"yTowsRebar_"<>ToString[i]<>", "<>ToString[yRebarA]<>", "<>ToString[yRebarS]<>", , FIBERS, 0., 1"
}
,{i,Length[yOrient]}]]
,
{
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Fict-1, part=Fict",
"*End Instance",
"**",
"*Instance, name=Continuum-1, part=Continuum",
"*End Instance",
"**",
"*Instance, name=xTows-1, part=xTows",
"*End Instance",
"**",
"*Instance, name=yTows-1, part=yTows",
"*End Instance",
"**",
"*Nset, nset=f1, instance=Fict-1",
"1",
"**",
"*Nset, nset=f2, instance=Fict-1",
"2",
"**",
"*Nset, nset=f3, instance=Fict-1",
"3",
"**",
"*Nset, nset=xMinus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=Continuum-1, unsorted",
"*include, input=InpFiles\\sets\\periodic_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=Continuum-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*Elset, elset=elem_xContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_xContinuumTowsSet.csv",
"**",
"*Elset, elset=elem_yContinuumTows, instance=Continuum-1",
"*include, input=InpFiles\\sets\\elem_yContinuumTowsSet.csv",
"**",
"*Embedded Element, host elset=Continuum-1.elem_Continuum",
"xTows-1.elem_xTows, yTows-1.elem_yTows",
"**",
"*Equation",
"3",
"xPlus, 1, 1.000000, xMinus, 1, -1.000000, f1, 1, -1.000000,",
"*Equation",
"3",
"xPlus, 2, 1.000000, xMinus, 2, -1.000000, f1, 2, -1.000000,",
"*Equation",
"3",
"xPlus, 3, 1.000000, xMinus, 3, -1.000000, f1, 3, -1.000000,",
"**",
"*Equation",
"3",
"yPlus, 1, 1.000000, yMinus, 1, -1.000000, f2, 1, -1.000000,",
"*Equation",
"3",
"yPlus, 2, 1.000000, yMinus, 2, -1.000000, f2, 2, -1.000000,",
"*Equation",
"3",
"yPlus, 3, 1.000000, yMinus, 3, -1.000000, f2, 3, -1.000000,",
"**",
"*Equation",
"3",
"zPlus, 1, 1.000000, zMinus, 1, -1.000000, f3, 1, -1.000000,",
"*Equation",
"3",
"zPlus, 2, 1.000000, zMinus, 2, -1.000000, f3, 2, -1.000000,",
"*Equation",
"3",
"zPlus, 3, 1.000000, zMinus, 3, -1.000000, f3, 3, -1.000000,",
"**",
"*Equation",
"2",
"f1, 2, "<>ToString[modelDims[[2]]]<>", f2, 1, "<>ToString[-1.*modelDims[[1]]]<>",",
"*Equation",
"2",
"f1, 3, "<>ToString[modelDims[[3]]]<>", f3, 1, "<>ToString[-1.*modelDims[[1]]]<>",",
"*Equation",
"2",
"f2, 3, "<>ToString[modelDims[[3]]]<>", f3, 2, "<>ToString[-1.*modelDims[[2]]]<>",",
"**",
"*End Assembly",
"**",
"** MATERIALS",
"**",
"*Material, name=MATRIX",
"*Elastic",
StringTrim[ToString[N[continuumProp]],("{"|"}")],
"**",
"*Material, name=FIBERS",
"*Elastic",
StringTrim[ToString[N[fiberProp]],("{"|"}")],
"**",
"** BOUNDARY CONDITIONS",
"**",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=1000000, NLGEOM=YES",
"*Static, stabilize",
"0.0001, 1., 0.0000001, 0.1",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"f2, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: Continuum",
"**",
"*Output, field, time interval=0.05",
"*Node Output, nset=Continuum-1.nodes_Continuum",
"RF, U",
"*Element Output, elset=Continuum-1.elem_Continuum, directions=YES",
"NE, LE, S",
"**",
"** FIELD OUTPUT: Fict",
"*Node Output, nset=f2",
"RF, U",
"**",
"** FIELD OUTPUT: xTowRebar",
"**",
"*Element Output, elset=xTows-1.elem_xTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** FIELD OUTPUT: yTowRebar",
"**",
"*Element Output, elset=yTows-1.elem_yTows, rebar, directions=YES",
"NE, LE, RBANG, RBFOR, RBROT, S",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


(* ::Subsection::Closed:: *)
(*MultiPart INP*)


elasticInp[pathOut_,modelName_,tows_,{towProp_,matrixProp_},disp_]:=Module[
{xTows,yTows,inp},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Part-1",
"*Node, input=InpFiles\\nodes.csv", 
"*Element, type=C3D8R, input=InpFiles\\elements.csv",
"**",
"*Elset,elset=elem_Matrix",
"*include, input=InpFiles\\sets\\elem_matrixSet.csv",
"**",
"*Elset,elset=elem_Tows",
"*include, input=InpFiles\\sets\\elem_towsSet.csv",
"**",
(*"*Elset,elset=elem_xTows",
"*include, input=InpFiles\\sets\\elem_xTowsSet.csv",
"**",
"*Elset,elset=elem_yTows",
"*include, input=InpFiles\\sets\\elem_yTowsSet.csv",*)
"**",
"*Distribution, name=dist, location=ELEMENT, table=tab_localori, input=InpFiles\\LocalOrientation.csv",
"*Orientation, name=LocalOri, definition=COORDINATES",
"dist",
"3, 0.",
"**",
"** Section: Section-Tows",
"*Solid Section, elset=elem_Tows, orientation=LocalOri, material=TOWS",
",",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Matrix, material=MATRIX",
",",
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Part-1-1, part=Part-1",
"*End Instance",
"**"
}
,
Flatten[Table[{"*Elset, elset=elem_xTow"<>ToString[t]<>", instance=Part-1-1","*include, input=InpFiles\\sets\\elem_xTowSet"<>ToString[t]<>".csv","**"}
,{t,xTows}]]
,
Flatten[Table[{"*Elset, elset=elem_yTow"<>ToString[t]<>", instance=Part-1-1","*include, input=InpFiles\\sets\\elem_yTowSet"<>ToString[t]<>".csv","**"}
,{t,yTows}]]
,
{
"*Nset, nset=xMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*End Assembly",
"**",
"*Distribution Table, name=tab_localori",
"COORD3D, COORD3D",
"**",
"** MATERIALS",
"**",
"*Material, name=TOWS",
"*Elastic, type=ENGINEERING CONSTANTS",
StringTrim[ToString[N[towProp[[1;;8]]]],("{"|"}")],
ToString[N[towProp[[9]]]]<>",",
"*Material, name=MATRIX",
"*Elastic",
StringTrim[ToString[N[matrixProp]],("{"|"}")],
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: PinCorner Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"corner, ENCASTRE",
"**",
"** Name: YSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"yMinus, YSYMM",
"**",
"** Name: ZSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"zMinus, ZSYMM",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=10000",
"*Static, stabilize",
"1., 1., 0.001, 1.",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"yPlus, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: F-Output-1",
"**",
"*Output, field, variable=PRESELECT",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


elasticUMATInp[pathOut_,modelName_,tows_,{towProp_,myMatProp_},disp_]:=Module[
{xTows,yTows,inp,num,matrixSS,myMat,constants},

If[Length[tows]==0,
{xTows,yTows}=IntegerPart[{tows/2,tows/2}];
,
{xTows,yTows}=IntegerPart[tows];
];

num=myMatProp[[2]];
matrixSS=myMatProp[[1]];
myMat=Join[{num},{N[Length[matrixSS]]},Flatten[matrixSS]];

If[SameQ[Mod[Length[myMat],8],0],
myMat=Partition[myMat,8];
,
myMat=Partition[myMat[[1;;Floor[Length[myMat],8]]],8]~Join~{myMat[[Floor[Length[myMat],8]+1;;]]};
];
constants=Length[Flatten[myMat,1]];

myMat=Map[StringReplace[ToString[#],{"{"->"","}"->","}]&,myMat];
myMat[[-1]]=StringDrop[myMat[[-1]],-1];

inp=Join[{
"*Heading",
"**Job name: "<>modelName<>" name: Model-1",
"**",
"** PARTS",
"**",
"*Part, name=Part-1",
"*Node, input=InpFiles\\nodes.csv", 
"*Element, type=C3D8R, input=InpFiles\\elements.csv",
"**",
"*Elset,elset=elem_Matrix",
"*include, input=InpFiles\\sets\\elem_matrixSet.csv",
"**",
"*Elset,elset=elem_Tows",
"*include, input=InpFiles\\sets\\elem_towsSet.csv",
(*"**",
"*Elset,elset=elem_xTows",
"*include, input=InpFiles\\sets\\elem_xTowsSet.csv",
"**",
"*Elset,elset=elem_yTows",
"*include, input=InpFiles\\sets\\elem_yTowsSet.csv",*)
"**",
"*Distribution, name=dist, location=ELEMENT, table=tab_localori, input=InpFiles\\LocalOrientation.csv",
"*Orientation, name=LocalOri, definition=COORDINATES",
"dist",
"3, 0.",
"**",
"** Section: Section-Tows",
"*Solid Section, elset=elem_Tows, orientation=LocalOri, material=TOWS",
",",
"** Section: Section-Matrix",
"*Solid Section, elset=elem_Matrix, material=MyMat",
",",
"*EndPart",
"**",
"** ASSEMBLY",
"**",
"*Assembly, name=Assembly",
"**",
"*Instance, name=Part-1-1, part=Part-1",
"*End Instance",
"**"
}
,
Flatten[Table[{"*Elset, elset=elem_xTow"<>ToString[t]<>", instance=Part-1-1","*include, input=InpFiles\\sets\\elem_xTowSet"<>ToString[t]<>".csv","**"}
,{t,xTows}]]
,
Flatten[Table[{"*Elset, elset=elem_yTow"<>ToString[t]<>", instance=Part-1-1","*include, input=InpFiles\\sets\\elem_yTowSet"<>ToString[t]<>".csv","**"}
,{t,yTows}]]
,
{
"*Nset, nset=xMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_xSetMinus.csv",
"**",
"*Nset, nset=xPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_xSetPlus.csv",
"**",
"*Nset, nset=yMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_ySetMinus.csv",
"**",
"*Nset, nset=yPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_ySetPlus.csv",
"**",
"*Nset, nset=zMinus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_zSetMinus.csv",
"**",
"*Nset, nset=zPlus, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_zSetPlus.csv",
"**",
"*Nset, nset=corner, instance=PART-1-1",
"*include, input=InpFiles\\sets\\node_cornerSet.csv",
"**",
"*End Assembly",
"**",
"*Distribution Table, name=tab_localori",
"COORD3D, COORD3D",
"**",
"** MATERIALS",
"**",
"*Material, name=MyMat",
"*User material, CONSTANTS = "<>ToString[constants]
}
~Join~
myMat
~Join~
{
"** ----------------------------------------------------------------",
"**",
"*Material, name=TOWS",
"*Elastic, type=ENGINEERING CONSTANTS",
StringTrim[ToString[N[towProp[[1;;8]]]],("{"|"}")],
ToString[N[towProp[[9]]]]<>",",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: PinCorner Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"corner, ENCASTRE",
"**",
"** Name: YSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"yMinus, YSYMM",
"**",
"** Name: ZSymmetry Type: Symmetry/Antisymmetry/Encastre",
"*Boundary",
"zMinus, ZSYMM",
"** ----------------------------------------------------------------",
"**",
"** STEP: Loading",
"**",
"*Step, name=Loading, inc=100000",
"*Static, stabilize",
"0.0001, 1., 0.00001, 0.1",
"**",
"** BOUNDARY CONDITIONS",
"**",
"** Name: YDisp Type: Displacement/Rotation",
"*Boundary",
"yPlus, 2, 2, "<>ToString[disp],
"**",
"** OUPUT REQUESTS",
"**",
"*Restart, write, frequency=0",
"**",
"** FIELD OUTPUT: F-Output-1",
"**",
"*Output, field, variable=PRESELECT, time interval = 0.100000",
"**",
"** HISTORY OUTPUT: H-Output-1",
"**",
"*Output, history, variable=PRESELECT",
"*End Step"
}
];

(*inp File*)
Export[pathOut<>"\\"<>modelName<>".txt",inp];
If[Length[FileNames[pathOut<>"\\"<>modelName<>".inp"]]==1,
DeleteFile[pathOut<>"\\"<>modelName<>".inp"];
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
,
RenameFile[pathOut<>"\\"<>modelName<>".txt",pathOut<>"\\"<>modelName<>".inp"]
]
];


(* ::Subsection:: *)
(*End*)


End[ ]
Protect @@ Complement[Names["Voxelization`*"],
{
}];
EndPackage[ ]
