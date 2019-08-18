(* ::Package:: *)

(* ::Section:: *)
(*Header*)


(* :Title: RQA *)
(* :Context: RQA` *)
(* :Author: Flip Phillips *)
(* :Summary: This package provides various things and stuff to Mathematica. *)
(* :Package Version: 0.3.0 *)
(* :Mathematica Version: 10.0+ *)
(* :Copyright: Copyright 2014-2019, Flip Phillips, All Rights Reserved.  *)
(* :History: *)
(* :Keywords: dynsys, dynamical systems,  *)
(* :Limitations:  *)
(* :Discussion:  *)


(* ::Text:: *)
(*Some special version notes:*)


(* ::Text:: *)
(*Version 0.1 - Based on a ton of papers that were not the Webber originals, but instead the interpretative work of others. *)
(*Version 0.2 - Based on the paper in the Dynsys NSF book.  Webber Jr, C. L., & Zbilut, J. P. (2005). Recurrence quantification analysis of nonlinear dynamical systems. Tutorials in contemporary nonlinear methods for the behavioral sciences, 94, 26-94. https://www.nsf.gov/pubs/2005/nsf05057/nmbs/chap2.pdf*)
(*Version 0.3 - Based on Webber Jr, C. L., & Marwan, N. (2015). Recurrence quantification analysis. Theory and Best Practices. Springer.*)


(* ::Section:: *)
(*Set up the package context, including any imports*)


BeginPackage["RQA`"]


(* ::Subsection:: *)
(*Unprotect exported symbols*)


Unprotect[{RQAEmbed,RQADistanceMap,RQARecurrenceMap,RQANeighbors,RQANearestNeighbors,RQANeighborDistances,RQAEstimateDimensionality,RQAEstimateLag,
	RQARecurrence,RQADeterminism,RQALaminarity,RQATrappingTime,RQATrend,RQAEntropy,RQADmax,RQAChaos,RQAVmax,
	RQAMakeTimeSeries,RQATimeSeriesEpochs,              
	RQAVerticalLineLengths,RQADiagonalLineLengths,
	RQANLines,RQANRecurrentPoints,RQANVerticalLines,RQANDiagonalLines,
	RQANVerticalRecurrentPoints,RQANDiagonalRecurrentPoints}]


(* ::Subsection:: *)
(*Usage Messages*)


(* ::Subsubsection::Closed:: *)
(*Package*)


RQA::usage="RQA is a package which provides a suite of basic Recurrance Quantification Analysis functions to Mathematica."


(* ::Subsubsection::Closed:: *)
(*Building the map and embeddings*)


RQAEmbeddedLastTime::usage="RQAEmbeddedLastTime[ts,d,tau] returns the last valid time for an embedding in dimension d for lag tau."


RQAEmbed::usage="RQAEmbed[ts,d,tau] embeds time series ts into d dimensions with a lag of tau."


RQADistanceMap::usage="RQADistanceMap[ts] computes a distance map of all samples to all other samples in time series ts. Note that ts is assumed to be regularly sampled."


RQARecurrenceMap::usage="RQARecurrenceMap[ts,r] computes a recurrence map from ts, using a recurrence raduis r."


RQANeighbors::usage="RQANeighbors[ts] computes a neighborhood for each sample in ts. Number of neighbors and distance functions are available as options."


RQANearestNeighbors::usage="RQANearestNeighbors[ts] computes a singular nearest neighbor for each value in ts."


RQANeighborDistances::usage="RQANeighborDistances[ts] computes the distance to each nearest neighbor for each value in ts."


RQAEstimateDimensionality::usage="RQAEstimateDimensionality[ts,tau] estimates the proper embedding dimension for a given lag tau using the false-neighbor algorithm."


RQAEstimateLag::usage="RQAEstimateLag[ts] estimates an appropriate lag / tau for a given time series using the decorrelation method."


(* ::Subsubsection::Closed:: *)
(*Structure helpers*)


RQAVerticalLineLengths::usage="RQAVerticalLineLengths[rm,tau] "


RQADiagonalLineLengths::usage="RQADiagonalLineLengths[rm,tau] "


RQANLines::usage="RQANLines[rm,tau] "


RQANRecurrentPoints::usage="RQANRecurrentPoints[rm,tau] "


RQANVerticalRecurrentPoints::usage="RQANVerticalRecurrentPoints[rm,tau] "


RQANDiagonalRecurrentPoints::usage="RQANDiagonalRecurrentPoints[rm,tau] "


(* ::Subsubsection::Closed:: *)
(*Estimating the parameters*)


RQARecurrence::usage="RQARecurrence[rm] calculates the recurrence for map rm."


RQADeterminism::usage="RQADeterminism[rm,l:2] calculates the determinism for map rm. Optionally accepts a minimum length run, defaults to 2."


RQALaminarity::usage="RQALaminarity[rm,l:2] calculates the laminarity for map rm. Optionally accepts a minimum length run, defaults to 2."


RQATrappingTime::usage="RQATrappingTime[rm,l:2] calculates the trapping time for map rm. Optionally accepts a minimum length run, defaults to 2."


RQATrend::usage="RQATrend[rm] calculates the trend for map rm."


RQAEntropy::usage="RQAEntropy[rm] calculates entropy for map rm. Optionally accepts a minimum length run, defaults to 2."


RQAChaos::usage="RQAChaos[rm] calculates chaos for map rm. Optionally accepts a minimum length run, defaults to 2."


RQADmax::usage="RQADmax[rm] calculates Dmax for map rm. Optionally accepts a minimum length run, defaults to 2."


RQAVmax::usage="RQAVmax[rm] calculates Vmax for map rm. Optionally accepts a minimum length run, defaults to 2."


(* ::Text:: *)
(*Utility Functions*)


RQAMakeTimeSeries::usage="RQAMakeTimeSeries[data,dt,t0] creates a regularly sampled time series from data, using dt as a time-step and, optionally t0 for starting time."


RQATimeSeriesEpochs::usage="RQATimeSeriesEpochs[ts,wid,overlap] creates a list of time series of a given width and overlap."


(* ::Subsection:: *)
(*Begin the private context*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Unprotect any system functions for which rules will be defined*)


(* ::Subsection:: *)
(*Definition of auxiliary functions and local (static) variables*)


(* ::Text:: *)
(*These are inline / adjacent as I don't like this particular format.*)


(* ::Subsection:: *)
(*Error messages for the exported objects*)


(* ::Text:: *)
(*These are inline / adjacent as I don't like this particular format.*)


(* ::Subsection:: *)
(*Definition of the exported functions*)


(* ::Subsubsection:: *)
(*TimeSeries functions*)


RQAMakeTimeSeries[data_,dt_,t0_:0]:=
  TimeSeries[
    TemporalData[data,{t0,Automatic,dt},
	  ValueDimensions->Length[Dimensions[data]]],
	  ResamplingMethod->Automatic]


RQATimeSeriesEpochs[ts_,wid_,overlap_]:=
  Table[TimeSeriesWindow[ts,{t0,t0+wid}],{t0,ts["FirstTime"],ts["LastTime"]-wid,wid-overlap}]


(* ::Subsubsection:: *)
(*Embedding*)


(* ::Text:: *)
(*This takes a time series ts, an embedding dimension D, and a time-delay \[Tau] and returns a TimeSeries. Note that it selects the time values such that there is no extrapolation.*)


embedHelper[f_,ts_,dim_,\[Tau]_]:=
	{#,Quiet@Table[f[#+(d*\[Tau])],{d,0,dim-1}]}&/@ts


RQAEmbeddedLastTime[ts_,d_,\[Tau]_]:=(ts["LastTime"]-(d-1)\[Tau])


Options[RQAEmbed]={
	"Truncate"->True,ResamplingMethod->Automatic,"Periodic"->False};


RQAEmbed::truncated="Resulting truncated series length would be < 0 for this d and \[Tau].";


RQAEmbed[ts_,d_,\[Tau]_,OptionsPattern[]]:=Module[{vals,validts,truncTime},
	validts=If[OptionValue["Truncate"],
		truncTime=RQAEmbeddedLastTime[ts,d,\[Tau]];
		If[truncTime<=0,Message[RQAEmbed::truncated];Return[Null],
			Map[Select[#,(#<=truncTime&)]&,ts["TimeList"]]],
		ts["TimeList"]];
	  
	vals=MapThread[embedHelper[#1,#2,d,\[Tau]]&,
		{ts["PathFunction",All],validts}];
	TimeSeries[TemporalData[vals],ResamplingMethod->OptionValue[ResamplingMethod]]]


(* ::Subsubsection:: *)
(*Map generation*)


(* ::Text:: *)
(*Helpers*)


brutalHeavisideTheta[x_]:=0/;PossibleZeroQ[x]
brutalHeavisideTheta[x_]:=HeavisideTheta[x]


(* ::Text:: *)
(*Compute a distance map (DM). *)


Options[RQADistanceMap]={
	DistanceFunction->SquaredEuclideanDistance,
	"Rescale"->False,"RescaleFunction"->Automatic};


RQADistanceMap::irreg="Time series is not regularly sampled.";


RQADistanceMap[ts_,OptionsPattern[]]:=Module[{df,rescale,map,rf},
	df=OptionValue[DistanceFunction];
	rescale=OptionValue["Rescale"];
	rf=OptionValue["RescaleFunction"];
	rf=If[rf===Automatic,Rescale,rf];
	
	(* handle this in the future, for now, this just throws an error. 
		we could have an option that tells what to do otherwise. *)
	If[!RegularlySampledQ[ts],Message[RQADistanceMap::irreg]];
	
	map=Outer[df,ts["Values"],ts["Values"],1];
	If[rescale,rf[map],map]]


(* ::Text:: *)
(*Create the RM using a simple, automatic threshold. As of version 0.2.1, this creates a non-binary recurrence map. Distances are included.*)


Options[RQARecurrenceMap]={
	DistanceFunction->SquaredEuclideanDistance,
	"Rescale"->False,"RescaleFunction"->Automatic}


(* ::Text:: *)
(*Two forms- if a time series is provided, compute the distance map and throw it away. *)


RQARecurrenceMap[ts_TemporalData,radius_:Automatic,opts:OptionsPattern[]]:=Module[{dm,r,map},
	dm=RQADistanceMap[ts,opts];
	r=If[radius===Automatic,Mean[Flatten[dm]],radius];
	SparseArray[Map[brutalHeavisideTheta[r-#]&,dm,{2}]*dm]]


(* ::Text:: *)
(*But, if a distance map is provided, make a different thing.*)


RQARecurrenceMap[dm_,radius_:Automatic,opts:OptionsPattern[]]:=Module[{r},
	r=If[radius===Automatic,Mean[Flatten[dm]],radius];
	SparseArray[Map[brutalHeavisideTheta[r-#]&,dm,{2}]*dm]] /;SquareMatrixQ[dm]


(* ::Subsubsection::Closed:: *)
(*Parameter estimation - \[Tau]*)


(* ::Text:: *)
(*For some reason, the Autocorrelation and CorrelationFunction in MMa 10/11 don't use the times of the time series, but instead seem to just use the raw "Values" so we need to figure out the time scale/samples/etc. First zero crossing of the autocorrelation*)


RQAEstimateLag[ts_,\[Tau]max_:Automatic]:=Module[{f,dt,smax,sols,t},
	dt=Mean[Differences[ts["Times"]]]; (* could just do First, or Median... *)
	smax=Round[If[\[Tau]max===Automatic,Length[ts["Values"]]/2.0,\[Tau]max/dt]];
	f=CorrelationFunction[TimeSeriesResample[ts,{1,Length[ts["Values"]],1}],{1,smax}];
	dt*t/.Quiet[FindRoot[f[t],{t,2,1,smax}]]]


(* ::Subsubsection::Closed:: *)
(*Parameter estimation - dimensionality*)


(* ::Text:: *)
(*Compute the nearest-neighbor of the time series in the current embedding dimension. usually you do this in the 1-D case, then see what happens as you increase D... but you don't have to.*)


(* ::Text:: *)
(*I should change this to find the point that is closest that isn't the point in question, more formally. Right now, it just finds the two closest points, and selects the second one since, the first one shouldn't be itself. But I think this might depend on the implementation of Nearest. *)


Options[RQANeighbors]={DistanceFunction->Automatic,"Length"->5}


(* ::Code:: *)
(*RQANeighbors[ts_,OptionsPattern[]]:=Module[{nf},*)
(*	nf=Nearest[ts["Values"],DistanceFunction->OptionValue[DistanceFunction]];*)
(*	Flatten[FirstPosition[ts["Values"],nf[#,2][[2]]]&/@ts["Values"]]]*)


(* ::Text:: *)
(*This version is a mess*)


RQANeighbors[ts_,OptionsPattern[]]:=Module[{nf,d,n},
	d=ts["Values"];
	n=OptionValue["Length"];
	nf=Nearest[d,DistanceFunction->OptionValue[DistanceFunction]];
	
	MapThread[
		If[Length[#1]==1,#1,DeleteCases[#1,#2]]&,
		{First[Flatten[Map[Position[d,#]&,(nf[#,n]&/@d),{2}],{3,4}]],
		 Range[Length[d]]}]]


RQANearestNeighbors[ts_,opts:OptionsPattern[]]:=First/@RQANeighbors[ts,opts]


(* ::Text:: *)
(*Given a ts, compute the distance from its nearest neighbor. As D increases, for some points, these neighbors will be the same or decrease, but, for others, that distance will increase, they are false neighbors that we will clean out in the next step.*)


Options[RQANeighborDistances]={DistanceFunction->SquaredEuclideanDistance,"Neighbors"->None};


(* ::Code:: *)
(*RQANeighborDistances[ts_,neighbors_:None,OptionsPattern[]]:=Module[{n,df},*)
(*	df=OptionValue[DistanceFunction];*)
(*	n=If[neighbors===None,RQANearestNeighbors[ts,DistanceFunction->df],neighbors];*)
(*	MapThread[OptionValue[DistanceFunction],{ts["Values"],ts["Values"][[n]]}]]*)


(* ::Text:: *)
(*The above might not be working anymore since the new embedding results in reduction of the length of the time series. Maybe I shouldn't be doing that?*)


RQANeighborDistances[ts_,opts:OptionsPattern[]]:=Module[{n,df},
	df=OptionValue[DistanceFunction];
	n=If[OptionValue["Neighbors"]===None,RQANearestNeighbors[ts,DistanceFunction->df],OptionValue["Neighbors"]];
	MapThread[OptionValue[DistanceFunction],{ts["Values"],ts["Values"][[n]]}]]


(* ::Text:: *)
(*Find 'false' neighbors as per  Kennel et al 1992 and 2002. Basically:*)
(*	\[Bullet] Find nearest neighbors in dimension d*)
(*	\[Bullet] Compute distances*)
(*	\[Bullet] Embed in d+1*)
(*	\[Bullet] Does distance increase? If so, false neighbor *)


(* ::Text:: *)
(*Total distance change-*)


neighborDistanceChange::uneq="Incompatable dimensions."


Options[neighborDistanceChange]={DistanceFunction->SquaredEuclideanDistance}


neighborDistanceChange[d1_,d2_,OptionsPattern[]]:=
	MapThread[OptionValue[DistanceFunction],{d1,d2}]


Options[falseNeighborP]={DistanceFunction->SquaredEuclideanDistance}


falseNeighborP[d1_,d2_,thresh_:Automatic,opts:OptionsPattern[]]:=Module[{d,t},
	d=neighborDistanceChange[Take[d1,Length[d2]],d2,opts];
	t=If[thresh===Automatic,Mean[d1],thresh];
	N[Count[d,x_/;x>t]/Length[d2]]
]


RQAEstimateDimensionality::dmax="Search exceeded maximum dimensionality: `1`"


Options[RQAEstimateDimensionality]={"Method"->Automatic,"Threshold"->Automatic}


RQAEstimateDimensionality[ts_,\[Tau]_,opts:OptionsPattern[]]:=Module[
	{pFalse,dsubE=2,pThresh,dMax=10,
	 to,
	 da,db,embed},

    (* options *)
	to=OptionValue["Threshold"];
	pThresh=If[to===Automatic||!NumberQ[to],0.01,to];

	(* bootstrap while loop *)
	da=RQANeighborDistances[ts];
	embed=RQAEmbed[ts,dsubE,\[Tau]];
	If[embed===Null,Return[1]];
	db=RQANeighborDistances[embed];
	pFalse=falseNeighborP[da,db,pThresh];Print[pFalse];
	
	(* iterate until we fall below pThresh false neighbors *)
	While[pFalse>pThresh && dsubE<=dMax,
		dsubE=dsubE+1;
		da=db;
			
		embed=RQAEmbed[ts,dsubE,\[Tau]];
		If[embed===Null,Return[1]];
		db=RQANeighborDistances[embed];
		pFalse=falseNeighborP[da,db,pThresh];Print[pFalse];
		];
			
	If[dsubE>=dMax,
	  Message[RQAEstimateDimensionality::dmax,dMax];1,
	  dsubE]
	]		
	


(* ::Subsubsection:: *)
(*Metric Helpers*)


(* ::Text:: *)
(*There are a lot of things that, in v0.1 and 0.2, were upper-triangular related. That makes sense for symmetric systems. (See issues where this might not be the case for certain, asymmetric neighborhoods) We're getting away from that so we don't have to generalize. Also, it turns out that some recurrence maps include the LOI. In this case, why bother, I think.*)


(* ::Text:: *)
(*big-N*)


size[rm_]:=First[Dimensions[rm]]


nWithLOI[rm_]:=size[rm]^2


nNoLOI[rm_]:=(nWithLOI[rm]-size[rm])


nUpperTriangular[rm_]:=(size[rm](size[rm]-1)/2)


(* ::Text:: *)
(*thresholders - there are probably better ways to do this, but for now.*)


(* ::Text:: *)
(*Segment counters - these don't care about triangularity*)


runLength[list_]:={First[#],Length[#]}&/@Split[list]


segmentLengths[m_,tLen_:2]:=Module[{rl,selFun,n},
	rl=runLength/@Normal[Unitize[m]];
	selFun=Select[#[[1]]==1\[And]#[[2]]>=tLen&];
	(selFun/@rl)/.{{}->{0},{1,n_}->n}]


(* ::Text:: *)
(*Segment makers*)


(* ::Code:: *)
(*verticalSegments[rm_]:= Transpose[UpperTriangularize[rm,1]]/;SquareMatrixQ[rm]*)


verticalSegments[rm_]:= Transpose[rm]/;SquareMatrixQ[rm]


(* ::Code:: *)
(*diagonalSegments[rm_]:= Diagonal[UpperTriangularize[rm,1],#]&/@Range[1,size[rm]-1] /;SquareMatrixQ[rm]*)


diagonalSegments[rm_]:= Diagonal[rm,#]&/@Range[-(size[rm]-1),size[rm]-1]/;SquareMatrixQ[rm]


(* ::Text:: *)
(*public stuff - these are questionable.*)


RQAVerticalLineLengths[rm_,thresh_:2]:=
	Select[Flatten[segmentLengths[verticalSegments[rm],thresh]],#>0&] /; SquareMatrixQ[rm]


RQADiagonalLineLengths[rm_,thresh_:2]:=
	Select[Flatten[segmentLengths[diagonalSegments[rm],thresh]],#>0&] /; SquareMatrixQ[rm]


RQANVerticalLines[rm_,thresh_:2]:=Length[RQAVerticalLineLengths[rm,thresh]]


RQANDiagonalLines[rm_,thresh_:2]:=Length[RQADiagonalLineLengths[rm,thresh]]


(* ::Text:: *)
(*Just wrong*)


(* ::Code:: *)
(*RQANLines[rm_,thresh_:2]:=RQANVerticalLines[rm,thresh]+RQANDiagonalLines[rm,thresh]*)


RQANVerticalRecurrentPoints[rm_,thresh_:2]:=Total[RQAVerticalLineLengths[rm,thresh]]


RQANDiagonalRecurrentPoints[rm_,thresh_:2]:=Total[RQADiagonalLineLengths[rm,thresh]]


(* ::Text:: *)
(*This was just, well, wrong.*)


(* ::Code:: *)
(*RQANRecurrentPoints[rm_,thresh_:2]:=RQANVerticalRecurrentPoints[rm,thresh]+RQANDiagonalRecurrentPoints[rm,thresh]*)


RQANRecurrentPoints[rm_SparseArray]:=Length[rm["NonzeroValues"]]


RQANRecurrentPoints[rm_]:=Count[rm,x_/;!PossibleZeroQ[x],{2}]


(* ::Subsubsection:: *)
(*Main estimators*)


RQARecurrence[rm_SparseArray]:=
	RQANRecurrentPoints[rm]/nNoLOI[rm]/;SquareMatrixQ[rm]


RQADeterminism[rm_,thresh_:2]:=
	(RQANDiagonalRecurrentPoints[rm,thresh]/RQANRecurrentPoints[rm]) /; SquareMatrixQ[rm]


RQALaminarity[rm_,thresh_:2]:=
	(RQANVerticalRecurrentPoints[rm,thresh]/RQANRecurrentPoints[rm,thresh]) /; SquareMatrixQ[rm]


RQATrappingTime[rm_,thresh_:2]:=
	(Mean[RQAVerticalLineLengths[rm,thresh]]) /; SquareMatrixQ[rm]


(* ::Text:: *)
(*This is broken*)


RQATrend[rm_]:=Module[{ut,w,d,lm,x},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];

	d=Transpose[
		{Range[1,w-1],Accumulate[Total[Diagonal[ut,#]]&/@
			Range[w-1]]/Accumulate[Range[w-1,1,-1]]}];

	lm=LinearModelFit[{1,100}#&/@d,x,x];
	1000*lm["ParameterTableEntries"][[2,1]]
]/;SquareMatrixQ[rm]


RQAEntropy[rm_,thresh_:2]:=
	(Entropy[2,RQADiagonalLineLengths[rm,thresh]]) /; SquareMatrixQ[rm]


RQAChaos[rm_,thresh_:2]:=
	(1/Max[RQADiagonalLineLengths[rm,thresh]]) /; SquareMatrixQ[rm]


(* ::Text:: *)
(*Diagonal max*)


RQADmax[rm_,thresh_:2]:= Max[RQADiagonalLineLengths[rm,thresh]] /; SquareMatrixQ[rm]


(* ::Text:: *)
(*Vertical max*)


RQAVmax[rm_,thresh_:2]:= Max[RQAVerticalLineLengths[rm,thresh]] /; SquareMatrixQ[rm]


(* ::Subsection:: *)
(*Rules for the system functions*)


(* ::Subsection:: *)
(*Restore protection of system functions*)


(* ::Subsection:: *)
(*End the private context*)


End[]


(* ::Subsection:: *)
(*Protect exported symbols*)


Protect[{RQAEmbed,RQADistanceMap,RQARecurrenceMap,RQANeighbors,RQANearestNeighbors,RQANeighborDistances,RQAEstimateDimensionality,RQAEstimateLag,
	RQARecurrence,RQADeterminism,RQALaminarity,RQATrappingTime,RQATrend,RQAEntropy,RQADmax,RQAChaos,RQAVmax,
	RQAMakeTimeSeries,RQATimeSeriesEpochs,
	RQAVerticalLineLengths,RQADiagonalLineLengths,
	RQANLines,RQANRecurrentPoints,RQANVerticalLines,RQANDiagonalLines,
	RQANVerticalRecurrentPoints,RQANDiagonalRecurrentPoints}]


(* ::Subsection:: *)
(*End the package context*)


EndPackage[]
