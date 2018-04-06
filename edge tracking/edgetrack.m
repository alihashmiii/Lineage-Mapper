(* ::Package:: *)

(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: edgetrack *)
(* :Context: edgetrack` *)
(* :Author: Ali Hashmi *)
(* :Date: 2018-01-15 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 Ali Hashmi *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["edgetrack`"]

trackedEdgeMask::usage = "generate mask of an edge between two cells";
plotEdge::usage = "plot of the tracked edge on the colorized segmented matrix";


Begin["`Private`"] (* `Private` *)


(* ::Subsection:: *)
(*Edge Tracking*)


segmentEdgesHelper[img_]:= Module[{imagetemp},
  imagetemp = ImageFilter[If[#[[3,3]]==1 && Total[#[[2;;-2,2;;-2]],2]==3,1,0]&,img,2]; (* create edge component matrix *)
  MorphologicalComponents[imagetemp]
];


segmentEdges[labeledMat_]:= Module[{img},
  img = Binarize[Thinning@Closing[MorphologicalPerimeter[Colorize[labeledMat]],2]]; (* create binarized edge mask *)
  {img,Sow@segmentEdgesHelper[img]}
];


(* creates edge -> cell linkage i.e. edge label -> cell(s) bordering it *)
edgeAssoc[labeledMat_]:= Module[{cellM,edgeMat,img,maxCell,totalMat,smalledges,edgetocell,position,
relabeledgetocell,edgetocellfinal},
  cellM = labeledMat;
  {img,edgeMat} = segmentEdges@labeledMat; (* edge-mask and EdgeComponentMatrix*)
  cellM = (ImageData[ColorNegate[img]] cellM)//Dilation[Erosion[#,1],1]&; (*cell component matrix *)
  maxCell = Max[cellM]; (* maximum cell label *)
  edgeMat = Map[If[# != 0,# + maxCell, #]&, edgeMat, {2}]; (* relabel the edge matrix *)
  totalMat = edgeMat + cellM; (* sum of edge and cell component matrix *)
  smalledges = Keys@ComponentMeasurements[edgeMat, #Count <= 3 &]; (* small edges <=3 pixels to be excluded later *)
  edgetocell = Select[ComponentMeasurements[totalMat,"ExteriorNeighbors"],(First@# > maxCell &)]; (* edge -> cell linkage *)
  position = Flatten[Position[Keys@edgetocell, Alternatives@@smalledges]]; (*position of smaller edges in the edge-cell linkage list*)
  relabeledgetocell = Thread[Range@Length[#]-> #]&@Values[edgetocell]; (*relabel edges -> cell indices *)
  edgetocellfinal = DeleteCases[relabeledgetocell,PatternSequence[Alternatives@@position-> _]]; (* delete small edges *)
  {position, edgetocellfinal} (* indices for tiny edges, edge-to-cell linkage *)
];


SetAttributes[edgeAssociation, {HoldFirst}];
edgeAssociation[segments_] := edgeAssociation[segments] = Reap[edgeAssoc/@segments];


(* ::Subsubsection:: *)
(*track shape/mask of an edge between shared cells*)


trackEdgebetweenCells[edgeCellLink_,cell1_,cell2_]:=Module[{sharededges,position},
  sharededges = Replace[Flatten[Rest/@edgeCellLink, 1],
    HoldPattern[_ -> {_Integer}] :> Sequence[], {2}];
  position = Position[sharededges, PatternSequence[HoldPattern[_ -> {cell1,cell2}] | HoldPattern[_ -> {cell2,cell1}]]];
  Thread[Map[First]@position-> Keys@Extract[sharededges,position]]
];


SetAttributes[trackedEdgeMask,{HoldFirst}];
Options[trackedEdgeMask]={"property"-> "Shape"};
trackedEdgeMask[segstacks_, cell1_, cell2_, OptionsPattern[]]:= With[{edgeAssignMat = Composition[First,Last]@edgeAssociation[Unevaluated@segstacks]},
    Module[{edgeToCellLinkage = First@edgeAssociation[segstacks]},
    First@Last@Reap[
      Module[{tracked = trackEdgebetweenCells[edgeToCellLinkage,cell1,cell2],frame,elem,
        optional = OptionValue["property"]},
        Do[
          {frame,elem} = {Keys[i],Values[i]};
          Sow[ComponentMeasurements[edgeAssignMat[[frame]],optional][[elem,2]]],
          {i, tracked}]
          ]
        ]
      ]
    ];


(* ::Subsubsection:: *)
(*visualize edge track*)


(*
(* highlight edge on the colorized cell component matrix *)
SetAttributes[plotEdge,{HoldFirst}];
plotEdge[segstacks_, cell1_, cell2_] := First@Last@Reap@With[{edgeAssignMat = Composition[First,Last]@edgeAssociation[Unevaluated@segstacks],
edgeCellLink = First@edgeAssociation[segstacks]},
Module[{tracked,frame, elem, edgematrix, labelcellMat},
   tracked = trackEdgebetweenCells[edgeCellLink,cell1,cell2];
   Do[
   {frame,elem}={Keys[i],Values[i]};
   {edgematrix,labelcellMat} = {Part[edgeAssignMat,frame],Part[segstacks,frame]};
   Sow[Colorize[labelcellMat] ~ HighlightImage ~ Image@ComponentMeasurements[edgematrix,"Mask"][[elem,2]]],
   {i,tracked}];
   ]
];
*)


(* highlight edge on the colorized cell component matrix *)
SetAttributes[plotEdge,{HoldFirst}];
plotEdge[segstacks_, cell1_, cell2_] := With[{edgeAssignMat = Composition[First,Last]@edgeAssociation[Unevaluated@segstacks],
edgeCellLink = First@edgeAssociation[segstacks]},
Module[{tracked,frame, elem, edgematrix, labelcellMat},
   tracked = trackEdgebetweenCells[edgeCellLink,cell1,cell2];
  MapThread[
   ({frame, elem}={#1, #2};
   {edgematrix,labelcellMat} = {Part[edgeAssignMat,frame],Part[segstacks,frame]};
   Colorize[labelcellMat]~HighlightImage~Image@ComponentMeasurements[edgematrix,"Mask"][[elem,2]])&,
   {Keys@tracked,Values@tracked}]
   ]
];


End[] 

EndPackage[]
