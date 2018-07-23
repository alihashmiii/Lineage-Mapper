(* ::Package:: *)

(* Mathematica Package *)

(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: vertextrack *)
(* :Context: vertextrack` *)
(* :Author: Ali Hashmi *)
(* :Date: 2018-07-22 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: 11.3 *)
(* :Copyright: (c) 2018 Ali Hashmi *)
(* :Keywords: *)
(* :Discussion: vertex tracking for lineage mapper. Applicable for epithelial tissues *)

BeginPackage["vertextrack`"]


trackVertices::usage = "trackVertices[images_,segments_] takes the sequence of binarized images and tracked labeled matrices and outputs a list of vertices
corresponding to the different cells";

extractCellVertices::usage = "extractCellVertices[track_,ind_] takes a list of tracked vertices and a cell index or a list of cell indices to filter for cell vertices";

visualizeCellVertices::usage = "visualizeCellVertices[imageStack_, trackVector_, ind_: All] utilizes the imagestack and the output of trackVertices[ ] to generate a visual
for tracked vertices. All vertices \"All\" or particular vertices can be defined as {1,2 ...}";


(*
Clear@associateVertices;
associateVertices[img_,seg_,DilationThresh_:1,maskDilation_:2]:=With[{dim =Reverse@ImageDimensions@img},
Module[{pts,segDil,morph,a,maskPixel,extractedMasks,elem,sArray,members,vertices,nearest},
pts = PixelValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}],1];
segDil=seg~Dilation~DilationThresh;
morph=MorphologicalComponents[ReplacePixelValue[ConstantImage[0,ImageDimensions@img],pts\[Rule]1],CornerNeighbors\[Rule]False];
{a,maskPixel}=Transpose[Values@ComponentMeasurements[morph,{"Area","Mask"}]];

If[MemberQ[a, 2.],
{extractedMasks,maskPixel} = Through[{Apply[Extract],Apply[Delete]}[{maskPixel,Position[a,2.]}]];
sArray = Map[(elem = #["NonzeroPositions"];
SparseArray[#\[Rule] 1,dim]&/@elem)&, extractedMasks];
maskPixel = Join[maskPixel, Level[sArray,{-1}]];
];
members = Block[{spArray,elems},
spArray=SparseArray[ImageData@Dilation[Image[#],maskDilation]]*segDil;
elems= Round@Union@spArray["NonzeroValues"]
]&/@maskPixel;
vertices = Cases[Thread[Round@members\[Rule] pts],HoldPattern[x:{__}/;Length@x\[GreaterEqual]3\[Rule] _]];
nearest = Nearest@Reverse[vertices,2];
KeyMap[Union@*Flatten]@GroupBy[
MapAt[Sort,(#\[Rule] nearest[#,{2,2}]&/@Values[vertices]),{All,2}],
Last\[Rule]First,Round@*N@*Mean]//Normal
]
];
*)


Clear@associateVertices;
associateVertices[img_,seg_,DilationThresh_:1,maskDilation_:2]:= With[{dim =Reverse@ImageDimensions@img },
Module[{pts,segDil,morph,mask,members,vertices,nearest},
pts = PixelValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}], 1];
segDil=seg~Dilation~DilationThresh;
members=Block[{spArray,elems},
elems=SparseArray[{First@dim-Last@#,First@#}-> 1,dim];
 spArray=SparseArray[ImageData@Dilation[Image@elems,maskDilation]]*segDil;
Round@Union@spArray["NonzeroValues"]
]&/@pts;
vertices=Cases[Thread[Round@members-> pts],HoldPattern[x:{__}/;Length@x >= 3 -> _]];
nearest=Nearest@Reverse[vertices, 2];
KeyMap[Union@*Flatten]@GroupBy[
MapAt[Sort,(#-> nearest[#,{2, 2}]&/@Values[vertices]),{All,2}],
Last->First,N@*Mean]//Normal 
]
];


Clear@trackVertices;
trackVertices[images_,segments_]:= ParallelTable[
associateVertices[images[[i]], segments[[i]]],{i,1,Length@images}];


Clear@subextractions;
subextractions[extractions_,ind_]:= FilterRules[extractions,{OrderlessPatternSequence[___,ind,___]}-> _];


Clear@extractCellVertices;
extractCellVertices[track_,ind_]:= First@*Last@Reap[
Catch@Scan[Block[{x},
x = Cases[#,PatternSequence[OrderlessPatternSequence[{___,Alternatives@@ind,___}]-> _]];
If[Length@x>0,Sow@x,Throw@"termination"]
]&,track]
];


Clear@visualizeCellVertices;
visualizeCellVertices[imageStack_,trackVector_,ind_: All]:=Block[{elems},
Switch[ind, All,
MapThread[HighlightImage,{Take[imageStack,Length@trackVector],Values@trackVector}]
, _ ,
(elems= Values@extractCellVertices[trackVector,ind];
MapThread[HighlightImage,{Take[imageStack,Length@elems],elems}]
)]//ListAnimate
];


(* `Private` *)
Begin["`Private`"]
End[] 


EndPackage[]
