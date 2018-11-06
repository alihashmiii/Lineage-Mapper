(* ::Package:: *)

(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: segmentStack *)
(* :Context: segmentStack` *)
(* :Author: Ali Hashmi *)
(* :Date: 2018-01-18 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2018 Ali Hashmi *)
(* :Keywords: *)
(* :Discussion: *)

BeginPackage["segmentStack`"]
(* Exported symbols added here with SymbolName::usage *)


(* ::Subsection:: *)
(*Segment Image*)


segmentImage[binarizedMask_?ImageQ,opt:"ConnectedComponents"|"Watershed":"Watershed",threshCellsize_:20000]:= Module[{seg,areas,
indexMaxarea,maxArea,indsmallareas={},$ind},
seg= Switch[opt,"ConnectedComponents",
(* assuming we input 0 as foreground and 1 as background. ConnectedComponents is a more general segmentation
framework *)
MorphologicalComponents[ColorNegate@binarizedMask, CornerNeighbors->False],
"Watershed",
(* for epithelial cells *)
WatershedComponents[binarizedMask, CornerNeighbors->False]
];
areas=ComponentMeasurements[seg,"Area"];
{indexMaxarea,maxArea}=First@MaximalBy[areas,Last]/.Rule-> List;
indsmallareas = Keys@Cases[areas,HoldPattern[_-> 1.]];
If[maxArea >= threshCellsize||indsmallareas!= {},
$ind={indexMaxarea}~Join~indsmallareas;
seg=ArrayComponents[seg,Length@areas,Thread[$ind->0]]
];
seg
];


Begin["`Private`"]

End[] (* `Private` *)

EndPackage[]
