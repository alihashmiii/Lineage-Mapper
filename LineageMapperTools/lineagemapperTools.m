(* ::Package:: *)

(* ::Section:: *)
(*LineageMapperTools*)


BeginPackage["LineageMapperTools`"];


neighboursQuery::usage = "given a labeled matrix the function can be used to determine either the \"Neighbours\" or the \"NeighbourCount\".
option -> Neighbor || NeighborCount";


lineageTable::usage = "construct a table of parents and daughters along with the frame number prior division. For the function provide the
 second argument/result from LineageMapper";


centroidMap::usage = "centroidMap[segments] takes in the segmented stacks to plot the displacement of either the entire cell population or a selection of cells
specified by the user that survive some x duration (integer). The option for example is provided as \"cutoffLength\" -> 3 ";


birthDeathFrameRecord::usage = "birthDeathFrameRecord[segments_] takes the segmented stacks and for each cell index outputs the birth and death frame"


lineageTree::usage = "given the linkages (2nd arg) of the result from LineageMapper` plot a Lineage Tree of the cell population";


cellTrackMov::usage = "cellTrackMov[segments_, index_] generates a movie for a single cell index for its entire lifespan.
To have a cropped movie use \"cropped\" -> True. The option is set to False by default";


cellMesh::usage = "cellMesh[segments_,index_] generates a geometric mesh for a cell index either embedded in the image coordinates or cropped.
To have the mesh embedded in the coordinate system use \"cropped\" -> False. The option is True by default";


cellExtract::usage = "cellExtract[segments_, index_] extracts specific cell(s) from the segmented image-stack";


(* ::Subsection:: *)
(*Functions*)


Begin["`Private`"];


Options[neighboursQuery]={"option" -> "Neighbors"};
neighboursQuery[labeledMat_, OptionsPattern[]]:= Dataset@<|ComponentMeasurements[labeledMat~Dilation~1, OptionValue@"option"]|>;


lineageTable[linkage_] := Dataset@DeleteMissing[
    MapIndexed[If[# =!= Null,
       <|"parentFrame" -> First[#2], "Parent" -> List /@ (Keys@#), "Daughters" -> Values@#|>,
       Missing[]] &, linkage]
       ];


Options[centroidMap] = {"cutoffLength" -> All};
centroidMap[segments_, OptionsPattern[]] := 
 Module[{cells, values, labels, option = OptionValue@"cutoffLength"},
  cells = GroupBy[Flatten[ComponentMeasurements[#, "Centroid"] & /@ segments, 1], Keys -> Values];
  {labels, values} = Switch[option, All, {Keys@#, Values@#} &@cells,
    _Integer, {Keys@#, Values@#}&@Normal[Select[Dataset@cells, Length@# > option &]]];
  Show[Graphics[{Thick, Dashed, Line /@ values}, Frame -> True,
    FrameStyle -> Directive[{Black, Thick, 12}], 
    PlotLabel -> Style["Centroid Displacements",
      {Bold, FontFamily -> "Consolas" , FontSize -> 18}]], 
   ListPlot[MapThread[Labeled[#1, #2, Mean@#] &, {values,Style[#, Bold, FontSize -> 11] & /@ labels}]],
   ImageSize -> Full]
   ];


birthDeathFrameRecord[segments_] := Dataset@GroupBy[
    Flatten@MapIndexed[Thread[First@#2 -> #1] &,
      Values@ComponentMeasurements[#, "Label"] & /@ segments],
    Last -> First, MinMax];


lineageTree[linkages_] := Block[{func, arrowend},
  arrowend = Graphics[Line[{{{-1, 1/2}, {0, 0}, {-1, -1/2}}}]];
  func[Null] = Sequence[];
  func[<|patt : HoldPattern[_ -> {_, _}] ... |>]:= Apply[Thread[Rule[#, #2]] &, {patt}, {1}];
  TreePlot[Flatten[func /@ linkages], Left,
   EdgeRenderingFunction -> ({Arrowheads[{{Automatic, Automatic, arrowend}}], Arrow[#1, 0.1]} &), 
   VertexRenderingFunction -> ({Lighter@RandomColor[], EdgeForm[None], Disk[#, .1], 
   Style[Text[#2, #1], {Bold, Black, FontSize -> 12}]} &)
   ]
  ];


singleCellExtract[segments_, index_Integer] := With[{ind = index},
   Block[{maskRep, elem},
    Map[Image]@First[(Reap@Scan[
          (maskRep = <|ComponentMeasurements[#, "Mask"]|>;
            elem = maskRep[ind];
            If[MatchQ[elem, _SparseArray], Sow[elem]]) &, segments])[[2]]
      ]
    ]
 ];


(*cellExtract[segments_,ind_] := cellExtract[segments, {ind}]*)
cellExtract[segments_,ind_Integer] := singleCellExtract[segments, ind];
cellExtract[segments_,ind_] := With[{dim = Composition[Reverse,Dimensions,First][segments]},
First@*Last@Reap@Block[{maskassoc,elem},
Scan[(maskassoc= <|ComponentMeasurements[#,"Mask"]|>;
elem=maskassoc[#]&/@ind;
Sow@Switch[Count[elem,_SparseArray],
0,ConstantImage[0,dim],
Except[0],Image[Total@Cases[elem,x_SparseArray:> x]]]
)&,segments]
	]
];


Options[cellTrackMov] = {"cropped" -> False, "duration" -> 1};
cellTrackMov[segments_, index_, OptionsPattern[]] := ListAnimate[
   If[! OptionValue["cropped"], #, ImageCrop /@ #]&@cellExtract[segments, index],
    DefaultDuration -> OptionValue["duration"]
   ];


Options[cellMesh] = {"cropped" -> True};
cellMesh[segments_, index_Integer, OptionsPattern[]] := Block[{plotopt},
   plotopt = Switch[OptionValue@"cropped", True, Automatic,
      False, Thread[{0, Reverse@*Dimensions@*First@segments}]
     ];
   Map[ImageMesh[#, PlotRange -> plotopt] &, segments~singleCellExtract~index]
   ];


End[];


EndPackage[];