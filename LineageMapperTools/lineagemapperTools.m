(* ::Package:: *)

(* ::Section:: *)
(*LineageMapperTools*)


BeginPackage["LineageMapperTools`"];


neighboursQuery::usage = "given a labeled matrix the function can be used to determine either the \"Neighbours\" or the
\"NeighbourCount\". option -> Neighbor || NeighborCount";


lineageTable::usage = "construct a table of parents and daughters along with the frame number prior division. For the function
provide the second argument/result from LineageMapper";


centroidMap::usage = "centroidMap[segments] takes in the segmented stacks to plot the displacement of either the entire cell
population or a selection of cells specified by the user that survive some x duration (integer). The option for example is
provided as \"cutoffLength\" -> 3 ";


birthDeathFrameRecord::usage = "birthDeathFrameRecord[segments_] takes the segmented stacks and for each cell index outputs
the birth and death frame"


lineageTree::usage = "given the linkages (2nd arg) of the result from LineageMapper` plot a Lineage Tree of the cell population";


cellTrackMov::usage = "cellTrackMov[segments_, index_] generates a movie for a single cell index for its entire lifespan.
To have a cropped movie use \"cropped\" -> True. The option is set to False by default";


cellMesh::usage = "cellMesh[segments_,index_] generates a geometric mesh for a cell index either embedded in the image
coordinates or cropped. To have the mesh embedded in the coordinate system use \"cropped\" -> False. The option is True by default";


cellExtract::usage = "cellExtract[segments_, index_] extracts specific cell(s) from the segmented image-stack";

confidenceIndex::usage = "confidenceIndex[seg_,mincelllife_:32,dilationfact_:2] generates a confidence index for the tracked cells"

apoptoticCells::usage = "apoptoticCells[seg_, minlifeThresh_:32, celldeathdelta_:10] yields possible apoptotic cells from cell-size and
centroid information"

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


fusionTree[linkages_]:= Block[{normlinks = Normal@Last@linkages,color},
 color = (#-> RandomColor[]&)/@Cases[normlinks,_Integer,{3,4}];
 Graph[Flatten@Map[Thread]@Flatten[Reverse[#,2]&/@normlinks, 1], VertexLabels->"Name",VertexStyle->color,VertexSize-> 0.5]
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

gWrapper[seg_,arg_,prop_,func_]:= GroupBy[
 Flatten[ComponentMeasurements[arg,prop]&/@seg],Keys-> Values,func
 ];

SetAttributes[confidenceIndex, HoldAll];
confidenceIndex[seg_,mincelllife_:32,dilationfact_:2]:= Block[{a,b,c},
 a =gWrapper[seg, #, "AdjacentBorderCount",If[#, 0, 1]&@*MemberQ[Except[0]]];
 b = Function[1/(#+1)]@gWrapper[seg,Unevaluated@Dilation[#,dilationfact],"Neighbors",Composition[Length,Union,Flatten]]; 
 c = gWrapper[seg, #, "Label",If[#>=mincelllife, 1, 0]&@*Length];
 Dataset@ReverseSort[(a + b + c + 1.0)/4]
]

apoptoticCells[seg_,minlifeThresh_:32,celldeathdelta_:10]:=Module[{sizeD, centD, cellspan ,cellspankeys, mean,sizediff, 
highestframenum,pos,filtered,frame,meandistcomp,lookup},
 sizeD = gWrapper[seg,#,"Count",Identity];
 centD =  gWrapper[seg,#,"Centroid",Identity];
 cellspan = Cases[Normal@Map[Length,sizeD], HoldPattern[_ -> x_]/;x >= minlifeThresh];
 cellspankeys = Keys@cellspan;
 lookup = Lookup[sizeD,cellspankeys];
 mean = 0.1*Mean[{#}]&@@@lookup;
 sizediff = BlockMap[EuclideanDistance[Sequence@@#]&,#,2,1]&/@lookup;
 highestframenum = MapThread[Position[#,x_/;x<#2]/.{}-> {{}}&,{sizediff,mean}]/.{x:{__Integer}..}:>First@*Last@{x} +1;
 pos = Position[Thread[(Values@cellspan - highestframenum +1) > minlifeThresh],True];
 filtered = FilterRules[centD,Extract[cellspankeys,pos]];
 frame = Extract[highestframenum,pos];
 meandistcomp = MapThread[First@#1 ->Mean@BlockMap[EuclideanDistance[Sequence@@#]&,(Last@#1 )[[#2;;]],2,1]& ,{filtered,frame}];
 Dataset@Keys@Cases[meandistcomp,HoldPattern[_ -> x_] /;x <= celldeathdelta,\[Infinity]]
];

End[];


EndPackage[];
