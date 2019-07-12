(* ::Package:: *)

(* ::Section:: *)
(*Header*)


(* :Title: RQA *)
(* :Context: RQA` *)
(* :Author: Flip Phillips *)
(* :Summary: This package provides various things and stuff to Mathematica. *)
(* :Package Version: 0.1.1 *)
(* :Mathematica Version: 10.0+ *)
(* :Copyright: Copyright 1988-2019, Flip Phillips, All Rights Reserved.  *)
(* :History: *)
(* :Keywords: packages, path, keywords *)
(* :Limitations:  *)
(* :Discussion:  *)


(* ::Section:: *)
(*Set up the package context, including any imports*)


BeginPackage["RQA`"]


(* ::Subsection:: *)
(*Usage Messages*)


(* ::Text:: *)
(*Package*)


RQA::usage="RQA is a package which provides RQA."


(* ::Text:: *)
(*Building the map*)


RQAEmbeddedLastTime::usage="RQAEmbeddedLastTime[ts,d,t] returns the last valid time for an embedding in dimension d for lag t."


RQAEmbed::usage="RQAEmbed[ts,dims,tau] but why."


RQADistanceMap::usage=""


RQARecurrenceMap::usage=""


RQANeighbors::usage=""


RQANearestNeighbors::usage=""


RQANeighborDistances::usage=""


RQAEstimateDimensionality::usage=""


RQAEstimateLag::usage=""


(* ::Text:: *)
(*Estimating the parameters*)


RQARecurrence::usage=""


RQADeterminism::usage=""


RQALaminarity::usage=""


RQATrappingTime::usage=""


RQATrend::usage=""


RQAEntropy::usage=""


RQAChaos::usage=""


RQADmax::usage=""


RQAVmax::usage=""


(* ::Text:: *)
(*Utility Functions*)


RQAMakeTimeSeries::usage=""


RQATimeSeriesEpochs::usage=""


(* ::Subsection:: *)
(*Begin the private context*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Unprotect any functions for which rules will be defined*)


Unprotect[{RQAEmbed,RQADistanceMap,RQARecurrenceMap,RQANeighbors,RQANearestNeighbors,RQANeighborDistances,RQAEstimateDimensionality,RQAEstimateLag,
	RQARecurrence,RQADeterminism,RQALaminarity,RQATrappingTime,RQATrend,RQAEntropy,RQADmax,RQAChaos,RQAVmax,
	RQAMakeTimeSeries,RQATimeSeriesEpochs}]


(* ::Subsection:: *)
(*Definition of auxiliary functions and local (static) variables*)


(* ::Subsection:: *)
(*Error messages for the exported objects*)


(* ::Subsection:: *)
(*Definition of the exported functions*)


(* ::Subsubsection:: *)
(*Utility functions*)


RQAMakeTimeSeries[data_,dt_,t0_:0]:=
	TimeSeries[
		TemporalData[data,{t0,Automatic,dt},
			ValueDimensions->Length[Dimensions[data]]],
		ResamplingMethod->Automatic]


RQATimeSeriesEpochs[ts_,wid_,overlap_]:=
Table[TimeSeriesWindow[ts,{t0,t0+wid}],{t0,ts["FirstTime"],ts["LastTime"]-wid,overlap}]


(* ::Subsubsection:: *)
(*Embedding and map generation*)


(* ::Text:: *)
(*This takes a time series ts, an embedding dimension D, and a time-delay \[Tau] and returns a TimeSeries. Note that it selects the time values such that there is no extrapolation.*)


embedHelper[f_,ts_,dim_,\[Tau]_]:=
	{#,Quiet@Table[f[#+(d*\[Tau])],{d,0,dim-1}]}&/@ts


RQAEmbeddedLastTime[ts_,d_,\[Tau]_]:=(ts["LastTime"]-(d-1)\[Tau])


Options[RQAEmbed]={
	"Truncate"->True,ResamplingMethod->Automatic};


RQAEmbed::truncated="Resulting truncated series length would be < 0 for this d and \[Tau].";


(* ::Code:: *)
(*RQAEmbed[ts_,d_,\[Tau]_,OptionsPattern[]]:=Module[{vals,validts},*)
(*	validts=If[OptionValue["Truncate"],*)
(*	  Map[*)
(*		Select[#,(#<=(ts["LastTime"]-\[Tau])&)]&,ts["TimeList"]],*)
(*	  ts["TimeList"]];*)
(*	  *)
(*	vals=MapThread[embedHelper[#1,#2,d,\[Tau]]&,*)
(*		{Take[ts["PathFunction",All],Length[validts]],validts}];*)
(*	TimeSeries[TemporalData[vals],ResamplingMethod->OptionValue[ResamplingMethod]]]*)


RQAEmbed[ts_,d_,\[Tau]_,OptionsPattern[]]:=Module[{vals,validts,truncTime},
	validts=If[OptionValue["Truncate"],
		truncTime=RQAEmbeddedLastTime[ts,d,\[Tau]];
		If[truncTime<=0,Message[RQAEmbed::truncated];Null,
			Map[Select[#,(#<=truncTime&)]&,ts["TimeList"]]],
		ts["TimeList"]];
	  
	vals=MapThread[embedHelper[#1,#2,d,\[Tau]]&,
		{ts["PathFunction",All],validts}];
	TimeSeries[TemporalData[vals],ResamplingMethod->OptionValue[ResamplingMethod]]]


(* ::Text:: *)
(*Compute a distance map (DM). Note that I should change the Rescale option to support other sorts of methods other than the max method.*)


Options[RQADistanceMap]={
	DistanceFunction->SquaredEuclideanDistance,
	"Rescale"->False,"RescaleFunction"->Automatic}


RQADistanceMap[ts_,OptionsPattern[]]:=Module[{df,rescale,map,rf},
	df=OptionValue[DistanceFunction];
	rescale=OptionValue["Rescale"];
	rf=OptionValue["RescaleFunction"];
	rf=If[rf===Automatic,Rescale,rf];
	
	map=Outer[df,ts["Values"],ts["Values"],1];
	If[rescale,rf[map],map]]


(* ::Text:: *)
(*Create the RM using a simple, automatic threshold. Again, I should be more clever about this.*)


Options[RQARecurrenceMap]={DistanceFunction->SquaredEuclideanDistance,"Rescale"->False}


RQARecurrenceMap[ts_,radius_:Automatic,opts:OptionsPattern[]]:=Module[{dm,r},
	dm=RQADistanceMap[ts,opts];
	r=If[radius===Automatic,Mean[Flatten[dm]],radius];
	Map[HeavisideTheta[r-#]&,dm,{2}]]


(* ::Subsubsection:: *)
(*Parameter estimation - \[Tau]*)


(* ::Text:: *)
(*For some reason, the Autocorrelation and CorrelationFunction in MMa 10/11 don't use the times of the time series, but instead seem to just use the raw "Values" so we need to figure out the time scale/samples/etc. First zero crossing of the autocorrelation*)


RQAEstimateLag[ts_,\[Tau]max_:Automatic]:=Module[{f,dt,smax,sols,t},
	dt=Mean[Differences[ts["Times"]]]; (* could just do First, or Median... *)
	smax=Round[If[\[Tau]max===Automatic,Length[ts["Values"]]/2.0,\[Tau]max/dt]];
	f=CorrelationFunction[TimeSeriesResample[ts,{1,Length[ts["Values"]],1}],{1,smax}];
	dt*t/.Quiet[FindRoot[f[t],{t,2,1,smax}]]]


(* ::Subsubsection:: *)
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
	{pFalse=1.0,dsubE=1,pThresh=0.001,dMax=10,to,
	 da,db,embeda,embedb,
	 res},

	(* bootstrap *)
	embedb=ts;
	db=RQANeighborDistances[embedb];
	
	to=OptionValue["Threshold"];
	pThresh=If[to===Automatic||!NumberQ[to],0.01,to];
	
	(* iterate until we fall below pThresh false neighbors *)
	res=Reap[
		While[pFalse>pThresh&&dsubE<dMax,
			dsubE=dsubE+1;
			da=db;embeda=embedb;
			
			embedb=RQAEmbed[ts,dsubE,\[Tau]];
			db=RQANeighborDistances[embedb];
			
			pFalse=falseNeighborP[da,db];
			Sow[{dsubE,pFalse}]]];
			
	If[dsubE>=dMax,Message[RQAEstimateDimensionality::dmax,dMax,0],		
	res[[2,-1,-1,1]]]]


(* ::Subsubsection:: *)
(*Metrics*)


RQARecurrence[rm_]:=Module[{ut,w,nMax,nRP},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];
	nMax=(w (w-1)/2);
	nRP=Total[Flatten[ut]];
	nRP/nMax
]/;SquareMatrixQ[rm]


runLength[list_]:={First[#],Length[#]}&/@Split[list]


segmentLengths[m_,tLen_:2]:=Module[{rl,selFun,n},
	rl=runLength/@m;
	selFun=Select[#[[1]]==1\[And]#[[2]]>=tLen&];
	(selFun/@rl)/.{{}->{0},{1,n_}->n}]


RQADeterminism[rm_,thresh_:2]:=Module[{ut,w,diags,segLens,lens,nRP},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];

	diags=Diagonal[ut,#]&/@Range[1,w-1];
	segLens=segmentLengths[diags,thresh];
	lens=Total/@segLens;
	nRP=Total[Flatten[ut]];
	Total[lens]/nRP
]/;SquareMatrixQ[rm]


RQALaminarity[rm_,thresh_:2]:=Module[{ut,verts,segLens,lens,nRP},
	ut=UpperTriangularize[rm,1];

	verts=Transpose[ut];
	segLens=segmentLengths[verts,thresh];
	lens=Total/@segLens;
	nRP=Total[Flatten[ut]];
	Total[lens]/nRP
]/;SquareMatrixQ[rm]


RQATrappingTime[rm_,thresh_:2]:=Module[{ut,verts,segLens,lens,nRP},
	ut=UpperTriangularize[rm,1];

	verts=Transpose[ut];
	segLens=Select[Flatten[segmentLengths[verts,thresh]],#>=thresh&];
	Mean[segLens]
]/;SquareMatrixQ[rm]


(* ::Text:: *)
(*I have resisted the 'percentage' stuff above, mainly because I don't really like thinking that way, but- this one I have done the gain/percentage thing because the stability is only really a thing when it is \[PlusMinus]5 units, according to their paper. Plus, the arbitrary rescaling is already done w/ the '1000' so why not go all the way to 100,000?*)


RQATrend[rm_]:=Module[{ut,w,d,lm,x},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];

	d=Transpose[
		{Range[1,w-1],Accumulate[Total[Diagonal[ut,#]]&/@
			Range[w-1]]/Accumulate[Range[w-1,1,-1]]}];

	lm=LinearModelFit[{1,100}#&/@d,x,x];
	1000*lm["ParameterTableEntries"][[2,1]]
]/;SquareMatrixQ[rm]


RQAEntropy[rm_,thresh_:2]:=Module[{ut,w,diags,segLens,lens,nRP},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];

	diags=Diagonal[ut,#]&/@Range[1,w-1];
	segLens=Select[Flatten[segmentLengths[diags,thresh]],#>0&];
	Entropy[2,segLens]
	]/;SquareMatrixQ[rm]


RQAChaos[rm_,thresh_:2]:=Module[{ut,w,diags,segLens,lens,nRP},
	ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];

	diags=Diagonal[ut,#]&/@Range[1,w-1];
	segLens=Select[Flatten[segmentLengths[diags,thresh]],#>0&];
	1/Max[segLens]
]/;SquareMatrixQ[rm]


RQADmax[rm_,thresh_:2]:=Module[{ut,w,diags,segLens,lens,nRP},ut=UpperTriangularize[rm,1];
	w=First[Dimensions[rm]];
	diags=Diagonal[ut,#]&/@Range[1,w-1];
	segLens=segmentLengths[diags,thresh];
	lens=Total/@segLens;
	Max[segLens]
	]/;SquareMatrixQ[rm]


RQAVmax[rm_,thresh_:2]:=Module[{ut,verts,segLens,lens,nRP},ut=UpperTriangularize[rm,1];
	verts=Transpose[ut];
	segLens=segmentLengths[verts,thresh];
	lens=Total/@segLens;
	Max[segLens]
]/;SquareMatrixQ[rm]


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
	RQAMakeTimeSeries,RQATimeSeriesEpochs}]


(* ::Subsection:: *)
(*End the package context*)


EndPackage[]
