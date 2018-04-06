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


segmentImage[binarizedMask_?ImageQ,threshCellsize_:20000]:= Module[{seg,areas,indexMaxarea,maxArea},
  seg=MorphologicalComponents@*ColorNegate@Dilation[binarizedMask,1];
  areas=ComponentMeasurements[seg,"Area"];
  {indexMaxarea,maxArea}=First@MaximalBy[areas,Last]/.Rule-> List;
  If[maxArea >= threshCellsize,ArrayComponents[seg,Length@areas,indexMaxarea->0],seg]~Dilation~1
];


(* segmentImage[binarizedMask_?ImageQ, thresh_:1000]:= Module[{seg, mincount},
  seg = DeleteBorderComponents@MorphologicalComponents[ColorNegate@binarizedMask, CornerNeighbors->False];
  mincount = Min@Values@ComponentMeasurements[seg,"Count"];
  ArrayComponents@If[mincount <= thresh, DeleteSmallComponents[seg, thresh], seg]
];*)


Begin["`Private`"]

End[] (* `Private` *)

EndPackage[]
