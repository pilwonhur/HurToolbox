(* ::Package:: *)

(* 	
	Hur - software for modeling and analysis of multibody systems    
	Copyright 2017, Pilwon Hur    

	Department of Mechanical Engineering 	
	Texas A&M University   	

	Hur Toolbox    
	Version 0.3, compatible with Mathematica 11.0  	

	Revision date: October 7, 2017 
*)  

BeginPackage["HurToolbox`"];

(* Usage statements *)
HurInitialize::usage="This procedure resets all global variables.";


HurInitialize := ( 	
	HurGlobalRF = {"n"}; HurGlobalDCM = List[RotationMatrix[0, {0, 0, 1}]];
)  




$VERSION$ = "0.3";
$EMAIL$ = "pilwonhur@tamu.edu";
Print["Hur Toolbox for modeling and analysis of multibody systems ", $VERSION$, "."];
Print["Copyright 2017 Pilwon Hur"];
Print["Department of Mechanical Engineering"];
Print["Texas A&M University"];
Print["All rights reserved."];
Print["Email questions, comments, or concerns to ", $EMAIL$, "."];





HurDefineRF::usage="HurDefineRF[rf_] defines a reference frame. Ex) HurDefineRF[b]: to define a reference frame B.";
HurGetNumGlobalRF::usage="HurGetNumGlobalRF[] returns the number of reference frames defined globally. Ex) n=HurGetNumGlobalRF[]";
HurGetIndexGlobalRF::usage="HurGetIndexGlobalRF[rf_] returns the index assigned for the reference frame. Ex) n=HurGetIndexGlobalRF[b]";
HurDefineDCM::usage="HurDefineDCM[rf_, dcm_] defines a direction cosine matrix of the reference frame with respect to the world reference frame. Ex) HurDefineDCM[b, RotationMatrix[theta[t],{0,0,1}]]";
HurGetAngularVel::usage="HurGetAngularVel[rf1_,rf2_] returns the angular velocity of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It uses skew symmetric matrix Sw such that Sw=R_dot x R'. Ex) w=HurGetAngularVel[b,d]";
HurGetAngularAcc::usage="HurGetAngularAcc[rf1_,rf2_] returns the angular acceleration of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It directly differentiate the angular velocity. Ex) w=HurGetAngularAcc[b,d]";
HurGetRelativeDCM::usage="HurGetRelativeDCM[rf1_, rf2_] returns the direction cosine matrix between rf1 and rf2. It is also equivalent to the rotation matrix of rf2 with respect to rf1. Ex) dcm=HurGetRelativeDCM[b,c]";
HurUnifyTriadPool::usage="HurUnifyTriadPool[v_, rf_] returns the list of triads to convert other triads with respect to the same reference frame. This is used only for the internal computation purpose. Ex) c1_in_b=HurUnifyTriadPool[c1, b]";
HurUnifyTriads::usage="HurUnifyTriads[v_, rf_] is used to represent the given vector with respect to rf. Ex) HurUnifyTriads[v, b]";
HurUnifyTriadsCoord::usage="HurUnifyTriadsCoord[v_, rf_] is used to display the coordinates of the given vector with respect to rf. Ex) HurUnifyTriadsCoord[v, b]";
HurCoordTriads::usage="HurCoordTriads[coord_] is used to represent the coorindate into the vector form. The reference frame of the coordinate is provided in the 4th element of the cooridnate list. Ex) HurCoordTriads[coord]";
HurCross::usage="HurCross[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a vector with respect to rf. Ex) HurCross[v1, v2, b]";
HurCrossCoord::usage="HurCrossCoord[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a coordinate with respect to rf. Ex) HurCross[v1, v2, b]";
HurVectorDiff::usage="HurVectorDiff[v_,rf1_,rf2_] is used to perform the vector differentiation. Vector v will be differentiated with respect to rf1. This function uses the vector differentiation formula for different RFs (i.e., rf2).";
HurAppendRF2Coord::usage="HurAppendRF2Coord[coord_, rf_]";

Begin["`Private`"];

HurDefineRF[rf_] := (AppendTo[HurGlobalRF, ToString[rf]]; 
  AppendTo[HurGlobalDCM, RotationMatrix[0, {0, 0, 1}]];)

HurGetNumGlobalRF[] := Length[HurGlobalRF];

HurGetIndexGlobalRF[rf_] := Position[HurGlobalRF, rf][[1]][[1]];

HurDefineDCM[rf_, dcm_] := (HurGlobalDCM[[HurGetIndexGlobalRF[ToString[rf]]]] = dcm;)

HurGetAngularVel[rf1_,rf2_] := (Rot5=HurGlobalDCM[[HurGetIndexGlobalRF[ToString[rf1]]]];
  Sw=D[Rot5, Global`t].Transpose[Rot5]//Simplify;
  ww=(-Sw[[2]][[3]])*ToExpression["n1"]+(Sw[[1]][[3]])*ToExpression["n2"]+(-Sw[[1]][[2]])*ToExpression["n3"]//FullSimplify;
  www=HurUnifyTriads[ww,rf2];
  Collect[www,ToExpression["{" <> ToString[rf2] <> "1," <> ToString[rf2] <> "2," <> ToString[rf2] <> "3}"]])

(*HurGetAngularAcc[rf1_,rf2_] := HurVectorDiff[HurGetAngularVel[rf1,rf2],"n",rf2];
*)

HurGetAngularAcc[rf1_,rf2_] := HurUnifyTriads[HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[HurGetAngularVel[rf1,"n"],"n"][[1;;3]],Global`t],"n"]],rf2]

HurUnifyTriadPool[v_, rf_] := (Rot1 = HurGlobalDCM[[HurGetIndexGlobalRF[StringTake[v, StringLength[v] - 1]]]]; 
  Rot2 = HurGlobalDCM[[HurGetIndexGlobalRF[rf]]]; 
  Rot3 = Transpose[Rot2].Rot1 // FullSimplify; 
  (Rot3[[1]][[ToExpression[StringTake[v, {StringLength[v], StringLength[v]}]]]])*ToExpression[ToString[rf] <> "1"]+(Rot3[[2]][[ToExpression[StringTake[v, {StringLength[v], StringLength[v]}]]]])* ToExpression[ToString[rf] <> "2"]+(Rot3[[3]][[ToExpression[StringTake[v, {StringLength[v], StringLength[v]}]]]])*ToExpression[ToString[rf] <> "3"])

HurGetRelativeDCM[rf1_, rf2_] := (Rot1=HurGlobalDCM[[HurGetIndexGlobalRF[ToString[rf1]]]];
	Rot2=HurGlobalDCM[[HurGetIndexGlobalRF[ToString[rf2]]]];
	Transpose[Rot1].Rot2// FullSimplify
	)

HurUnifyTriads[v_, rf_] := (temp=v /. 
   Flatten[Table[
     ToExpression[HurGlobalRF[[i]] <> ToString[j]] \[Rule]  
       HurUnifyTriadPool[HurGlobalRF[[i]] <> ToString[j], 
        ToString[rf]], {i, HurGetNumGlobalRF[]}, {j, 3}]]//FullSimplify;
   Collect[temp,ToExpression["{" <> ToString[rf] <> "1," <> ToString[rf] <> "2," <> ToString[rf] <> "3}"]])

HurUnifyTriadsCoord[v_, rf_] := (coord=Table[D[HurUnifyTriads[v, rf], 
   ToExpression[ToString[rf] <> ToString[i]]], {i, 3}]//FullSimplify;
	HurAppendRF2Coord[coord,rf])

HurCoordTriads[v_] := (v[[1]])*ToExpression[ToString[v[[4]]]<>"1"]+(v[[2]])*ToExpression[ToString[v[[4]]]<>"2"]+(v[[3]])*ToExpression[ToString[v[[4]]]<>"3"]

HurCross[v1_, v2_, rf_] := (coord = Cross[HurUnifyTriadsCoord[v1, rf][[1;;3]], HurUnifyTriadsCoord[v2, rf][[1;;3]]]//FullSimplify; 
  (coord[[1]])*ToExpression[ToString[rf] <> "1"]+(coord[[2]])*ToExpression[ToString[rf] <> "2"]+(coord[[3]])*ToExpression[ToString[rf] <> "3"])

HurCrossCoord[v1_, v2_, rf_] := 
 (temp=Cross[HurUnifyTriadsCoord[v1, rf][[1;;3]], HurUnifyTriadsCoord[v2, rf][[1;;3]]]//FullSimplify;
 	HurAppendRF2Coord[temp,rf]
 	)

HurAppendRF2Coord[coord_, rf_] := Flatten[List[coord,ToString[rf]]]

(*
HurVectorDiff[v_,rf1_,rf2_] := (df2dvdt=HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[v,rf2][[1;;3]],Global`t],rf2]];
  www=HurUnifyTriads[HurGetAngularVel[rf2]-HurGetAngularVel[rf1],rf2]//FullSimplify;
  wcrossv=HurCross[www,v,rf2]//FullSimplify;
  df1dvdt=HurUnifyTriadsCoord[df2dvdt+wcrossv,rf2]//FullSimplify;
  HurCoordTriads[df1dvdt]
  )
*)

HurVectorDiff[v_,rf1_,rf2_] := (df2dvdt=D[HurUnifyTriadsCoord[v,rf2][[1;;3]],Global`t];
  www=HurUnifyTriads[HurGetAngularVel[rf2,rf2]-HurGetAngularVel[rf1,rf2],rf2]//FullSimplify;
  wcrossv=HurCrossCoord[www,v,rf2]//FullSimplify;
  df1dvdt=df2dvdt+wcrossv[[1;;3]]//FullSimplify;
  HurCoordTriads[HurAppendRF2Coord[df1dvdt,rf2]]
  )

  
End[];

EndPackage[]

(*
Version History
0.1
Added the basic functions for vector analysis and reference frame manipulations.

0.2
Handled time derivative problem when t was treated as priviate variable. Forced to treat as global variable.
https://mathematica.stackexchange.com/questions/114769/derivative-from-my-package-function-returns-0
Added a function to compute the angular velocity directly from the rotation matrix.
Added a Vector Differentiation Function

0.3
Added a function to compute the DCM between two reference frames






  *)