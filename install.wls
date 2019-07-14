#!/usr/bin/env wolframscript
(* ::Package:: *)

packaegName="RQA";


dir = FileNameJoin[{$UserBaseDirectory, "Applications", packaegName}];
Print["Installing into "<>dir];


(* ::Text:: *)
(*OK, this is a shell script and won't work for windows installing, but, I don't use windows much. It just makes a symbolic link from here to make things a little easier to debug.*)


res = Run["ln -sf ./RQA "<>dir];


Print["Done: "<>ToString[res]]
