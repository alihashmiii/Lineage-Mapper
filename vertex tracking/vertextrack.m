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


trackVertices::usage = "trackVertices[images_,segments_] takes the sequence of binarized images and tracked labeled matrices and
outputs a list of vertices corresponding to the different cells";

extractCellVertices::usage = "extractCellVertices[track_,ind_] takes a list of tracked vertices and a cell index or a list of cell
indices to filter for cell vertices";

visualizeCellVertices::usage = "visualizeCellVertices[imageStack_, trackVector_, ind_: All] utilizes the imagestack and the output
of trackVertices[] to generate a visual for tracked vertices. All vertices \"All\" or particular vertices can be defined as {1,2 ...}";

(*
Options[associateVertices]={"watershed"-> True,"dilSeg"-> 1};
associateVertices[img_,segt_,maskDil_:2,OptionsPattern[]]:= With[{dim =Reverse@ImageDimensions@img,watershed=OptionValue["watershed"],
dilSeg=OptionValue["dilSeg"]},
Module[{pts,members,vertices,nearest,segDil=segt},
pts = PixelValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}], 1];
If[!OptionValue["watershed"],segDil=Dilation[segDil,dilSeg]];
members=Block[{spArray,elems},
elems=SparseArray[{First@dim-Last@#,First@#}->1,dim];
spArray=SparseArray[ImageData@Dilation[Image@elems,maskDil]]*segDil;
Round@Union@spArray["NonzeroValues"]
]&/@pts;
vertices=Cases[Thread[Round@members-> pts],HoldPattern[pattern:{__}/;Length@pattern >= 3 -> _]];
nearest=Nearest@Reverse[vertices, 2];
KeyMap[Union@*Flatten]@GroupBy[
MapAt[Sort,(#-> nearest[#,{2, 2}]&/@Values[vertices]),{All,2}],
Last->First,N@*Mean]//Normal
]
];
*)


Options[associateVertices]={"watershed"-> True,"dilSeg"-> 1,"stringentCheck"-> True};
associateVertices[img_,segt_,maskDil_:2,OptionsPattern[]]:= With[{dim =Reverse@ImageDimensions@img,watershed=OptionValue["watershed"],
dilSeg=OptionValue["dilSeg"],stringentQ=OptionValue["stringentCheck"]},
Module[{pts,members,vertices,nearest,segC=segt,vertexset,likelymergers,imagegraph,imggraphweight,imggraphpts,vertexpairs,
posVertexMergers,meanVertices,Fn},
pts = PixelValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}], 1]; (* finding branch points *)

If[!OptionValue["watershed"],segC=Dilation[segC,dilSeg]]; (* dilate if connected components was used for segmentation *)

members = Block[{spArray,elems},
 elems = SparseArray[{First@dim-Last@#,First@#}->1,dim];
 spArray = SparseArray[ImageData@Dilation[Image@elems,maskDil]]*segC;
 Round@Union@spArray["NonzeroValues"]
 ]&/@pts;
 
vertices = Cases[Thread[Round@members-> pts],HoldPattern[pattern:{__}/;Length@pattern >= 3 -> _]];
(* finding vertices with three or more neighbouring cells *)

nearest = Nearest[Reverse[vertices, 2],DistanceFunction->ManhattanDistance]; (* nearest func for candidate vertices *)
Fn = GroupBy[MapAt[Sort,(#-> nearest[#,{All,2}]&/@Values[vertices]),{All,2}],Last->First,#]&;

Which[Not@stringentQ,
 (* merge if candidate vertices are 2 manhattan blocks away. Not a stringent check for merging *)
 KeyMap[Union@*Flatten]@Fn[N@*Mean]//Normal,
 stringentQ,
 (* a better check is to see the pixels separating the vertices are less than 3 blocks *)
 vertexset = Fn[Identity];
 (* candidates for merging*)
 likelymergers = Cases[Normal[vertexset],PatternSequence[{{__Integer}..}-> i:{__List}/;Length[i]>= 2]];
 (*defining graph properties of the image *)
 imagegraph = MorphologicalGraph@MorphologicalTransform[img,{"Fill"}];
 imggraphweight = AssociationThread[(EdgeList[imagegraph]/.UndirectedEdge->List )-> PropertyValue[imagegraph,EdgeWeight]];
 imggraphpts = Nearest@Reverse[Thread[VertexList[imagegraph]-> PropertyValue[imagegraph,VertexCoordinates]],2];
 (* corresponding nodes on the graph *)
 vertexpairs = Union@*Flatten@*imggraphpts/@(Values[likelymergers]);
 (* find pairs < than 3 edgeweights away, take a mean of vertices and update the association with mean position *)
 posVertexMergers = Position[Thread[Lookup[imggraphweight,vertexpairs]<3],True];
 If[posVertexMergers != {},
  meanVertices=MapAt[N@*Mean,likelymergers,Thread[{Flatten@posVertexMergers,2}]];
  Scan[(vertexset[#[[1]]]=#[[2]])&,meanVertices]
  ];
  KeyMap[Union@*Flatten]@vertexset//Normal]
  ]
]; 


trackVertices[image_Image,segment_]:= associateVertices[image, segment];
trackVertices[images_,segments_]:= ParallelTable[associateVertices[images[[i]], segments[[i]]],{i,1,Length@images}];


subextractions[extractions_,ind_]:= FilterRules[extractions,{OrderlessPatternSequence[___,ind,___]}-> _];


extractCellVertices[track_,ind_]:= First@*Last@Reap[
Catch@Scan[Block[{tv},
tv = Cases[#,PatternSequence[OrderlessPatternSequence[{___,Alternatives@@ind,___}]-> _]];
If[Length@tv>0,Sow@tv,Throw@"termination"]
]&,track]
];


visualizeCellVertices[image_Image,trackVector_]:=HighlightImage[image,{Green,Values@trackVector}];
visualizeCellVertices[imageStack_,trackVector_,ind_: All]:=Block[{elems},
Switch[ind, All,
 MapThread[HighlightImage[#1,{Green,#2}]&,{Take[imageStack,Length@trackVector],Values@trackVector}]
 , _ ,
 (elems= Values@extractCellVertices[trackVector,ind];
 MapThread[HighlightImage[#1,{Green,#2}]&,{Take[imageStack,Length@elems],elems}]
 )]//ListAnimate
];


(* `Private` *)
Begin["`Private`"]
End[]


EndPackage[]
