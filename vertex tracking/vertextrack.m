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

LaunchKernels[];

trackVertices::usage = "trackVertices[images_,segments_] takes the sequence of binarized images and tracked labeled matrices and
outputs a list of vertices corresponding to the different cells";

extractCellVertices::usage = "extractCellVertices[track_,ind_] takes a list of tracked vertices and a cell index or a list of cell
indices to filter for cell vertices";

visualizeCellVertices::usage = "visualizeCellVertices[imageStack_, trackVector_, ind_: All] utilizes the imagestack and the output
of trackVertices[] to generate a visual for tracked vertices. All vertices \"All\" or particular vertices can be defined as {1,2 ...}";


Options[associateVertices]= {"stringentCheck"-> True};
associateVertices[img_,segt_,maskDil_:2,OptionsPattern[]]:= With[{dim =Reverse@ImageDimensions@img,
stringentQ=OptionValue["stringentCheck"]},
Module[{pts,members,vertices,nearest,vertexset,likelymergers,imagegraph,imggraphweight,imggraphpts,vertexpairs,
posVertexMergers,meanVertices,Fn},

pts = PixelValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}], 1]; (* finding branch points *)
members = Block[{spArray,elems},
 elems = SparseArray[{First@dim-Last@#,First@#}->1,dim];
 spArray = SparseArray[ImageData@Dilation[Image@elems,maskDil]]*segt;
 Round@Union@spArray["NonzeroValues"]
 ]&/@pts;
 
(*pts = ImageValuePositions[MorphologicalTransform[img,{"Fill","SkeletonBranchPoints"}], 1]; (* finding branch points *)
members = Block[{elems},
 elems = Dilation[ReplaceImageValue[ConstantImage[0,Reverse@dim],#->1],1];
 DeleteCases[Union@Flatten@ImageData[elems*Image[segt]],0.]
 ]&/@pts;
*)
 
vertices = Cases[Thread[Round@members-> pts],HoldPattern[pattern:{__}/;Length@pattern >= 3 -> _]];
(* finding vertices with three or more neighbouring cells *)

nearest = Nearest[Reverse[vertices, 2]]; (* nearest func for candidate vertices *)
Fn = GroupBy[MapAt[Sort,(#-> nearest[#,{All,2}]&/@Values[vertices]),{All,2}],Last->First,#]&;

Which[Not@stringentQ,
 (* merge if candidate vertices are 2 euclidean distance away. Not a stringent check for merging *)
 KeyMap[Union@*Flatten]@Fn[List@*N@*Mean]//Normal,
 stringentQ,
 (* a better check is to see the non-zero pixels separating the vertices are less than 3 away *)
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
 posVertexMergers = Position[Thread[Lookup[imggraphweight,vertexpairs] < 3],True];
 If[posVertexMergers != {},
  meanVertices=MapAt[List@*N@*Mean,likelymergers,Thread[{Flatten@posVertexMergers,2}]];
  Scan[(vertexset[#[[1]]]=#[[2]])&,meanVertices]
  ];
  KeyMap[Union@*Flatten]@vertexset//Normal]
  ]
]; 


trackVertices[image_Image,segment_]:= associateVertices[image, segment];
DistributeDefinitions[associateVertices];
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
