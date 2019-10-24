(* ::Package:: *)

(* 	
	Hur - software for modeling and analysis of multibody systems    
	Copyright 2019, Pilwon Hur    

	Department of Mechanical Engineering 	
	Texas A&M University   	

	Hur Toolbox    
	Version 0.7, compatible with Mathematica 11.0  	

	Revision date: September 7, 2019 
*)  

BeginPackage["HurToolbox`"];

(* Usage statements *)
HurInitialize::usage="This procedure resets all global variables.";

$VERSION$ = "0.7.5";
$EMAIL$ = "pilwonhur@tamu.edu";
Print["Hur Toolbox for modeling and analysis of multibody systems ", $VERSION$, ". \nCopyright 2019 Pilwon Hur\nDepartment of Mechanical Engineering\nTexas A&M University\nAll rights reserved.\nEmail questions, comments, or concerns to ", $EMAIL$, "."];


HurDefineRF::usage="HurDefineRF[rf__] defines a reference frame. Ex) HurDefineRF[b]: to define a reference frame B. HurDefineRF[a,b]: to define reference frames a, and b at the same time.";
HurGetNumGlobalRF::usage="HurGetNumGlobalRF[] returns the number of reference frames defined globally. Ex) n=HurGetNumGlobalRF[]";
HurGetIndexGlobalRF::usage="HurGetIndexGlobalRF[rf_] returns the index assigned for the reference frame. Ex) n=HurGetIndexGlobalRF[b]";
HurDefineDCM::usage="HurDefineDCM[rf_, dcm_] defines a direction cosine matrix of the reference frame with respect to the world reference frame. Ex) HurDefineDCM[b, RotationMatrix[theta[t]+phi[t],{0,0,1}]]";
HurDefineDCMRelative::usage="HurDefineDCMRelative[rf1_, rf2_, dcm_] defines a direction cosine matrix of the reference frame rf1 with respect to rf2. Ex) HurDefineDCMRelative[b, a, RotationMatrix[phi[t],{0,0,1}]]";
HurDefineMass::usage="HurDefineMass[rf_, m_] defines the mass of the reference frame. Ex) HurDefineMass[b, m]";
HurDefineInertia::usage="HurDefineInertia[rf_, II_] defines the inertia matrix of the reference frame about the center of mass in terms of body reference frame. Ex) HurDefineInertia[b, {Ixx,Ixy,Ixz,Iyy,Iyz,Izz}]";
HurGetAngularVel::usage="HurGetAngularVel[rf1_,rf2_] returns the angular velocity of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It uses skew symmetric matrix Sw such that Sw=R_dot x R'. Ex) w=HurGetAngularVel[b,d]";
HurGetAngularAcc::usage="HurGetAngularAcc[rf1_,rf2_] returns the angular acceleration of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It directly differentiate the angular velocity. Ex) w=HurGetAngularAcc[b,d]";
HurGetRelativeDCM::usage="HurGetRelativeDCM[rf1_, rf2_] returns the direction cosine matrix between rf1 and rf2. It is also equivalent to the rotation matrix of rf1 with respect to rf2. Ex) dcm=HurGetRelativeDCM[b,c]";
HurUnifyTriadPool::usage="HurUnifyTriadPool[rf1_, rf2_] returns triads for rf1 expressed in terms of rf2. In other words, it returns the x,y,z axes of rf1 expressed in terms of rf2. Note that when rf1 is defined via HurDefineRF and HurDefineDCM, rf1 is expressed in terms of Newtonian reference frame. This is used only for the internal computation purpose. It is also equivalent to HurGetRelativeDCM[rf1_, rf2_]. Ex) HurUnifyTriadPool[a, b]";
HurUnifyTriads::usage="HurUnifyTriads[v_, rf_] is used to represent the given vector with respect to rf. Ex) HurUnifyTriads[v, b]";
HurUnifyTriadsCoord::usage="HurUnifyTriadsCoord[v_, rf_] is used to display the coordinates of the given vector with respect to rf. Ex) HurUnifyTriadsCoord[v, b]";
HurCoordTriads::usage="HurCoordTriads[coord_] is used to represent the coordinate into the vector form. The reference frame of the coordinate is provided in the 4th element of the cooridnate list. Ex) HurCoordTriads[coord]";
HurCross::usage="HurCross[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a vector with respect to rf. Ex) HurCross[v1, v2, b]";
HurDot::usage="HurDot[v1_, v2_] is used to perform the dot product of the two vectors. v1 and v2 can be in different RFs. Also, HurDot does not require RF information since the output is a scalar. Ex) HurDot[v1, v2]";
HurNorm::usage="HurNorm[v_] returns the norm of the vector. Vector v can have triads of mixed RFs."
HurCrossCoord::usage="HurCrossCoord[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a coordinate with respect to rf. It is equivalent to the following: HurUnifyTriadsCoord[ HurCross[v1_, v2_, rf_], rf_] Ex) HurCrossCoord[v1, v2, b]";
HurVectorDiff::usage="HurVectorDiff[v_,rf1_,rf2_] is used to perform the vector differentiation. Vector v will be differentiated with respect to rf1. This function uses the vector differentiation formula for different RFs (i.e., rf2). Note that if a vector is differentiated w.r.t. time, then it will be the same as HurVectorDiff[v_,rf1_,n].";
HurAppendRF2Coord::usage="HurAppendRF2Coord[coord_, rf_] explicitly specify the RF information to the given coordinate (i.e., 3 numbers) without RF. Regardless of the size of coord_, HurAppendRF2Coord[coord_, rf_] takes the first 3 components of coord_ and attach rf_ to it. Usually, this function is used for internal usage.";
HurDefineForces::usage="HurDefineForces[rf_, force_, r_]";
HurDefineMoments::usage="HurDefineMoments[rf_, moment_]";
HurNEResultantForce::usage="HurNEResultantForce[rf1_, rf2_]";
HurNEResultantMoment::usage="HurNEResultantMoment[rf1_, rf2_]";
HurGetAngularMomentum::usage="HurGetAngularMomentum[rf1_, rf2_]";
HurGetLinearMomentum::usage="HurGetLinearMomentum[rf1_, rf2_]";
HurGetLinearCOMVel::usage="HurGetLinearCOMVel[rf1_,rf2_] computes the linear velocity of the COM of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. Ex) w=HurGetLinearCOMVel[a,b]";
HurGetLinearCOMAcc::usage="HurGetLinearCOMAcc[rf1_,rf2_] computes the linear acceleration of the COM of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. Ex) w=HurGetAngularAcc[a,b]";
HurGetAngularVel::usage="HurGetAngularVel[rf1_,rf2_] computes the angular velocity of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. It uses skew symmetric matrix Sw such that Sw=R_dot x R'. Ex) w=HurGetAngularVel[b,d]";
HurGetAngularAcc::usage="HurGetLinearCOMAcc[rf1_,rf2_] computes the angular acceleration of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. It directly differentiate the angular velocity. Ex) w=HurGetAngularAcc[b,d]";
HurKinematics::usage="HurKinematics[] computes the linear velocity and acceleration of COM of all RFs, and angular velocity and acceleration of all RFs. HurKinematics[] uses information from HurGlobalCOMPos. In other words, HurDefineCOMPos[rf_,v_] needs to be run beforehand. Once done correctly, the following variables will be set: HurGlobalCOMPos, HurGlobalCOMVel, HurGlobalCOMAcc, HurGlobalAngularVel, HurGlobalAngularAcc.";
HurDefineCOMPos::usage="HurDefineCOMPos[rf_,v_] defines COM position (v) of the specified RF. Ex) HurDefineCOMPos[ a , 5 a1+2 a2 ]";
HurSetCOMVel::usage="HurSetCOMVel[rf_,v_] sets velocity of the COM of the specified RF. In the usual case, it is not recommended to use this function to forcefully define the COM velocity. If you know what you are doing for sure, then use it. It may cause inconsistency, though. Mostly, this function is used for internal usage by HurKinematics[].";
HurSetCOMAcc::usage="HurSetCOMAcc[rf_,acc_] sets acceleration of the COM of the specified RF. In the usual case, it is not recommended to use this function to forcefully define the COM acceleration. If you know what you are doing for sure, then use it. It may cause inconsistency, though. Mostly, this function is used for internal usage by HurKinematics[].";
HurSetAngularVel::usage="HurSetAngularVel[rf_,w_] sets the angular velocity of the specified RF. In the usual case, it is not recommended to use this function to forcefully define the angular velocity of the RF. If you know what you are doing for sure, then use it. It may cause inconsistency, though. Mostly, this function is used for internal usage by HurKinematics[].";
HurSetAngularAcc::usage="HurSetAngularAcc[rf_,alpha_] sets the angular acceleration of the specified RF. In the usual case, it is not recommended to use this function to forcefully define the angular acceleration of the RF. If you know what you are doing for sure, then use it. It may cause inconsistency, though. Mostly, this function is used for internal usage by HurKinematics[].";
HurGetNEEquation::usage="HurGetNEEquation[rf1_, rf2]";
HurGetNEForm::usage="HurGetNEForm[rf1_, rf2_]";
HurGetKineticE::usage="HurGetKineticE[rf__]";
HurGetPotentialE::usage="HurGetKineticE[rf__]";
HurGetLagrangian::usage="HurGetLagrangian[rf__]";
HurDefineVertical::usage="HurDefineVertical[v_]";
HurGetELEquation::usage="HurGetELEquation[rf_,gc_]";
HurGetMMatrix::usage="HurGetMMatrix[]";
HurGetCMatrix::usage="HurGetCMatrix[]";
HurGetGVector::usage="HurGetGVector[]";
HurNEEquation::usage="HurNEEquation[]";
HurDefineVariableList::usage="HurDefineVariableList[f__]";
HurSolveNEInverse::usage="HurSolveNEInverse[]";
HurNEForward::usage="HurNEForward[]";
HurMakeSymmetricMatrix::usage="HurMakeSymmetricMatrix[list_]";
HurDefineGeneralizedCoordinates::usage="HurDefineGeneralizedCoordinates[f__]";
HurELEquation::usage="HurELEquation[]";
HurDefineConstraints::usage="HurDefineConstraints[con__]";
HurDefineConstrainedJacobian::usage="HurDefineConstrainedJacobian[]";
HurConstrainedELEquation::usage="HurConstrainedELEquation[]";
HurConstrainedELInverse::usage="HurConstrainedELInverse[]";
HurELInverse::usage="HurELInverse[]";
HurDefineLambda::usage="HurDefineLambda[]";
HurDefineNonConservativeForces::usage="HurDefineNonConservativeForces[f__]"; 
HurDefineOtherPotentialE::usage="HurDefineOtherPotentialE[rf_, pe_]"; 
HurGetInertiaTensor::usage="HurGetInertiaTensor[rf_]"; 
HurProductMatVec::usage="HurProductMatVec[mat_,vec_,rf_]"; 
HurDumpSaveData::usage="HurDumpSaveData[filename__]";
HurSaveData::usage="HurSaveData[filename__]";
HurLoadData::usage="HurLoadData[filename_]";


HurInitialize[] := (  
  HurGlobalRF = {Global`n};
  HurGlobalDCM = List[RotationMatrix[0, {0, 0, 1}]]; HurGlobalMass={1}; HurGlobalInertia={{1,0,0,1,0,1}};
  HurGlobalForce = {}; HurGlobalMoment = {}; HurGlobalCOMPos = {0}; HurGlobalCOMVel = {0}; HurGlobalCOMAcc = {0};
  HurGlobalAngularVel = {0}; HurGlobalAngularAcc = {0};
  HurGlobalLinearMomentum = {0}; HurGlobalAngularMomentum = {0}; HurGlobalVertical = {0};
  HurGlobalNEEquation = {{}}; HurGlobalVariableList = {};
  HurGlobalKineticE = {0}; HurGlobalPotentialE = {0}; HurGlobalLagrangian = {0};
  HurGlobalELEquation = {0}; HurGlobalGeneralizedCoordinates = {};
  HurGlobalMMatrix = {}; HurGlobalCMatrix = {}; HurGlobalGVector = {};
  HurGlobalConstrainedJacobian = {0}; HurGlobalConstraints = {};
  HurGlobalLambda = {}; HurGlobalGeneralizedConstrainingForce = {};
  HurGlobalConstrainedELEquation = {0}; HurGlobalConstrainedModified = {0};
  HurGlobalNonConservativeForces = {0}; HurGlobalOtherPotentialE = {0};
  HurGlobalListTriads = {{n1,n2,n3}}; 
  HurGlobalTriadsConversion = {{n1->n1,n2->n2,n3->n3}};
  (*HurGlobalTriadsConversion = {0};*)
  HurGlobalSimplify = True;
)

Begin["`Private`"];

HurDefineRF[rf__] := (rfs=List[rf];narg=Length[rfs];
  Do[
      If[
        HurGetIndexGlobalRF[ rfs[[i]] ]===0
        ,
        AppendTo[HurGlobalRF, rfs[[i]] ]; 

        AppendTo[HurGlobalListTriads, Table[ Symbol[ ToString[ rfs[[i]] ] <> ToString[j] ], {j,3}] ]; 

        AppendTo[HurGlobalDCM, RotationMatrix[0, {0, 0, 1}]];
        HurDefineDCM[ rfs[[i]] ,RotationMatrix[0, {0, 0, 1}] ];
        AppendTo[HurGlobalMass, 0];
        AppendTo[HurGlobalInertia, {0,0,0,0,0,0}];
        AppendTo[HurGlobalAngularVel, 0];
        AppendTo[HurGlobalAngularAcc, 0];
        AppendTo[HurGlobalCOMPos, 0];
        AppendTo[HurGlobalCOMVel, 0];
        AppendTo[HurGlobalCOMAcc, 0];
        AppendTo[HurGlobalLinearMomentum, 0];
        AppendTo[HurGlobalAngularMomentum, 0];
        AppendTo[HurGlobalNEEquation, {}];
        AppendTo[HurGlobalKineticE, 0];
        AppendTo[HurGlobalPotentialE, 0];
        AppendTo[HurGlobalLagrangian, 0];
        AppendTo[HurGlobalOtherPotentialE, 0];
        ,
        Null
      ];
    , 
    {i,narg}
    ];

(*
  nrf=Length[HurGlobalRF];
  HurGlobalListTriads = Table[
    Symbol[ ToString[ HurGlobalRF[[i]] ] <> ToString[j] ]
    ,
    {i,nrf},{j,3}
    ]
*)
  )


HurGetNumGlobalRF[] := Length[HurGlobalRF];

HurGetIndexGlobalRF[rf_] := 
  If[
      Flatten[Position[HurGlobalRF, rf]]==={},
      0,
      Flatten[Position[HurGlobalRF, rf]][[1]]
    ]

HurDefineDCM[rf_, dcm_] := (HurDefineDCMRelative[rf, HurGlobalRF[[1]], dcm];
  )

HurDefineDCMRelative[rf1_, rf2_, dcm_] := (HurGlobalDCM[[ HurGetIndexGlobalRF[ rf1 ] ]] = HurGlobalDCM[[ HurGetIndexGlobalRF[ rf2 ] ]].dcm;
  nrf=Length[HurGlobalRF];

  HurGlobalTriadsConversion=Table[0,{i,nrf}];

  Do[
    HurGlobalTriadsConversion[[k]] = Table[
      Rot=HurUnifyTriadPool[  HurGlobalRF[[ i ]] , HurGlobalRF[[ k ]] ];
      Symbol[ ToString[ HurGlobalRF[[i]] ] <> ToString[j] ] \[Rule]  
      Dot[ Rot[[;;,j]],HurGlobalListTriads[[k]] ]
      ,
      {i,nrf},{j,3}
      ]
    ,
    {k,nrf}
    ]

  )

(* Kinematics *)

HurUnifyTriadPool[rf1_, rf2_] := (Rot1 = HurGlobalDCM[[ HurGetIndexGlobalRF[ rf1 ] ]]; Rot2 = HurGlobalDCM[[ HurGetIndexGlobalRF[rf2]] ]; 
  Rot3 = Transpose[Rot2].Rot1;
  If[HurGlobalSimplify, Simplify[Rot3], Rot3]
  )

HurGetAngularVel[rf1_,rf2_] := (Rot5=HurGlobalDCM[[ HurGetIndexGlobalRF[rf1] ]];
  Sw=D[Rot5, Global`t].Transpose[Rot5];
  ww=(-Sw[[2]][[3]])*HurGlobalListTriads[[1,1]]+(Sw[[1]][[3]])*HurGlobalListTriads[[1,2]]+(-Sw[[1]][[2]])*HurGlobalListTriads[[1,3]];
  ww1=If[HurGlobalSimplify, Simplify[ww], ww];
  www=HurUnifyTriads[ww1,rf2];
  HurSetAngularVel[rf1,www];
  www
  )

HurGetAngularAcc[rf1_,rf2_] := (
  alpha=HurUnifyTriads[HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[HurGetAngularVel[rf1,HurGlobalRF[[1]] ],HurGlobalRF[[1]] ][[1;;3]],Global`t],HurGlobalRF[[1]] ]],rf2];
  alpha1=If[HurGlobalSimplify, Simplify[alpha], alpha];
  HurSetAngularAcc[rf1,alpha1];
  alpha1
  )

HurSetAngularVel[rf_,w_] := (
  HurGlobalAngularVel[[ HurGetIndexGlobalRF[ rf ] ]] = w;
  )

HurSetAngularAcc[rf_,alpha_] := (
  HurGlobalAngularAcc[[HurGetIndexGlobalRF[ rf ] ]] = alpha;
  )

HurDefineCOMPos[rf_,v_] := (HurGlobalCOMPos[[HurGetIndexGlobalRF[ rf ] ]] = v;)

HurSetCOMVel[rf_,v_] := (HurGlobalCOMVel[[HurGetIndexGlobalRF[ rf ] ]] = v;)

HurSetCOMAcc[rf_,acc_] := (HurGlobalCOMAcc[[HurGetIndexGlobalRF[ rf ] ]] = acc;)

HurGetLinearCOMVel[rf1_,rf2_] := (
  v=HurUnifyTriads[ HurCoordTriads[ HurAppendRF2Coord[D[HurUnifyTriadsCoord[ HurGlobalCOMPos[[ HurGetIndexGlobalRF[ rf1 ] ]] ,HurGlobalRF[[1]] ][[1;;3]],Global`t],HurGlobalRF[[1]] ]],rf2];
  v1=If[HurGlobalSimplify, Simplify[v], v];
  HurSetCOMVel[rf1,v1];
  v1
  )

HurGetLinearCOMAcc[rf1_,rf2_] := (
  acc=HurUnifyTriads[HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[ HurGlobalCOMPos[[ HurGetIndexGlobalRF[ rf1 ] ]] ,HurGlobalRF[[1]] ][[1;;3]],Global`t,Global`t],HurGlobalRF[[1]] ]],rf2];
  acc1=If[HurGlobalSimplify, Simplify[acc], acc];
  HurSetCOMAcc[rf1,acc1];
  acc1
  )

HurGetRelativeDCM[rf1_, rf2_] := (Rot1=HurGlobalDCM[[ HurGetIndexGlobalRF[rf1] ]];
	Rot2=HurGlobalDCM[[HurGetIndexGlobalRF[ rf2 ] ]];
  Rot3=Transpose[Rot2].Rot1;
	If[HurGlobalSimplify, Simplify[Rot3], Rot3]
	)

HurUnifyTriads[v_, rf_] := (temp=v /. Flatten[ HurGlobalTriadsConversion[[ HurGetIndexGlobalRF[ rf ] ]] ];
   Collect[ temp,HurGlobalListTriads[[ HurGetIndexGlobalRF[ rf ] ]] ]
   )

HurUnifyTriadsCoord[v_, rf_] := (coord=Table[ D[ HurUnifyTriads[v, rf], 
   HurGlobalListTriads[[ HurGetIndexGlobalRF[ rf ],i ]] ] , {i, 3}];
  coord1=If[HurGlobalSimplify, Simplify[coord], coord];
	HurAppendRF2Coord[coord1,rf]
  )

HurCoordTriads[v_] := Dot[ v[[1;;3]],HurGlobalListTriads[[ HurGetIndexGlobalRF[ v[[4]] ] ]] ]    

HurCross[v1_, v2_, rf_] := (coord = Cross[HurUnifyTriadsCoord[v1, rf][[1;;3]], HurUnifyTriadsCoord[v2, rf][[1;;3]]]; 
  coord1=If[HurGlobalSimplify, Simplify[coord], coord];
  HurCoordTriads[ Flatten[ List[coord1,rf] ] ]
  )  

HurDot[v1_, v2_] := (tempD=Dot[HurUnifyTriadsCoord[v1, HurGlobalRF[[1]] ][[1;;3]], HurUnifyTriadsCoord[v2, HurGlobalRF[[1]] ][[1;;3]]];
  If[HurGlobalSimplify, Simplify[tempD], tempD]
  )

HurCrossCoord[v1_, v2_, rf_] := 
 (temp=Cross[HurUnifyTriadsCoord[v1, rf][[1;;3]], HurUnifyTriadsCoord[v2, rf][[1;;3]]];
  temp1=If[HurGlobalSimplify, Simplify[temp], temp];
 	HurAppendRF2Coord[temp1,rf]
 	)

HurNorm[v_] := (
  Norm[ HurUnifyTriadsCoord[ v, HurGlobalRF[[1]] ] [[1;;3]] ]
  )

HurAppendRF2Coord[coord_, rf_] := Flatten[ List[coord[[1 ;; 3]],rf] ]

(*
HurVectorDiff[v_,rf1_,rf2_] := (df2dvdt=HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[v,rf2][[1;;3]],Global`t],rf2]];
  www=HurUnifyTriads[HurGetAngularVel[rf2]-HurGetAngularVel[rf1],rf2]//Simplify;
  wcrossv=HurCross[www,v,rf2]//Simplify;
  df1dvdt=HurUnifyTriadsCoord[df2dvdt+wcrossv,rf2]//Simplify;
  HurCoordTriads[df1dvdt]
  )
*)

HurVectorDiff[v_, rf1_, rf2_] := (df2dvdt=D[HurUnifyTriadsCoord[v,rf2][[1;;3]],Global`t];
  www=HurUnifyTriads[HurGetAngularVel[rf2,rf2]-HurGetAngularVel[rf1,rf2],rf2];
  wcrossv=HurCrossCoord[www,v,rf2];
  df1dvdt=df2dvdt+wcrossv[[1;;3]];
  df1dvdt1=If[HurGlobalSimplify, Simplify[df1dvdt], df1dvdt];
  HurCoordTriads[HurAppendRF2Coord[df1dvdt1,rf2]]
  )

HurKinematics[] := (nrfs=HurGetNumGlobalRF[];
  If[
    nrfs===1
    ,
    Print["There is only one reference frame!"]
    ,
    Do[
      HurGetAngularVel[ HurGlobalRF[[i]] , HurGlobalRF[[i]] ]
      HurGetAngularAcc[ HurGlobalRF[[i]] , HurGlobalRF[[i]] ]
      HurGetLinearCOMVel[ HurGlobalRF[[i]] , HurGlobalRF[[i]] ]
      HurGetLinearCOMAcc[ HurGlobalRF[[i]] , HurGlobalRF[[i]] ]
      ,
      {i,2,nrfs}
      ]
    ]
  )









(* Newton-Euler Mechanics *)

HurProductMatVec[mat_,vec_,rf_] := (ww=HurUnifyTriadsCoord[vec,rf];
  HurUnifyTriads[HurCoordTriads[Flatten[{mat.{ww[[1]],ww[[2]],ww[[3]]},ww[[4]]}]],rf]
  )



HurGetInertiaTensor[rf_] := (IItemp=HurGlobalInertia[[ HurGetIndexGlobalRF[rf] ]];
  II={{IItemp[[1]],IItemp[[2]],IItemp[[3]]},{IItemp[[2]],IItemp[[4]],IItemp[[5]]},{IItemp[[3]],IItemp[[5]],IItemp[[6]]}};
  II)

HurGetAngularMomentum[rf1_, rf2_] := (II=HurGetInertiaTensor[rf1];
  ww=HurUnifyTriadsCoord[HurGetAngularVel[rf1,rf1],rf1];
  H=HurUnifyTriads[HurCoordTriads[Flatten[{II.{ww[[1]],ww[[2]],ww[[3]]},ww[[4]]}]],rf2];
  H1=If[HurGlobalSimplify, Simplify[H], H];
  HurGlobalAngularMomentum[[HurGetIndexGlobalRF[ rf1 ] ]] = H1;
  H1
  )


HurGetLinearMomentum[rf1_, rf2_] := (
  v=HurGlobalCOMVel[[HurGetIndexGlobalRF[ rf1 ] ]];
  m=HurGlobalMass[[HurGetIndexGlobalRF[ rf1] ]];
  ll=HurUnifyTriads[m*v,rf2];
  ll1=If[HurGlobalSimplify, Simplify[ll], ll];
  HurGlobalLinearMomentum[[HurGetIndexGlobalRF[ rf1 ] ]] = ll1;
  ll1
  )

HurDefineVertical[v_] := (HurGlobalVertical[[1]]=v;)

HurDefineMass[rf_, m_] := (HurGlobalMass[[ HurGetIndexGlobalRF[ rf ] ]] = m;)

HurDefineInertia[rf_, II_] := (HurGlobalInertia[[HurGetIndexGlobalRF[ rf ] ]] = II;)

HurDefineForces[rf_, force_, r_] := (AppendTo[HurGlobalForce,{rf,force,r}];)

HurDefineMoments[rf_, moment_] := (AppendTo[HurGlobalMoment,{rf,moment}];)

HurNEResultantForce[rf1_, rf2_] := HurUnifyTriads[Total[Table[HurGlobalForce[[i]][[2]],{i,Flatten[Position[HurGlobalForce[[;;,1]],rf1]] } ] ],rf2]

HurNEResultantMoment[rf1_, rf2_] := (Mom=Total[Table[HurGlobalMoment[[i]][[2]],{i,Flatten[Position[HurGlobalMoment[[;;,1]],rf1] ] } ] ];
  Force=Total[Table[HurCross[HurGlobalForce[[i]][[3]],HurGlobalForce[[i]][[2]],rf2],{i,Flatten[Position[HurGlobalForce[[;;,1]],rf1] ] } ] ];
  HurUnifyTriads[Mom+Force,rf2]
  )

HurGetNEEquation[rf1_, rf2_] := (
  eq1=HurNEResultantForce[rf1,rf2]-HurVectorDiff[HurGlobalLinearMomentum[[HurGetIndexGlobalRF[ rf1 ] ]],HurGlobalRF[[1]] ,rf2];
  eq2=HurNEResultantMoment[rf1,rf2]-HurVectorDiff[HurGlobalAngularMomentum[[HurGetIndexGlobalRF[ rf1 ] ]],HurGlobalRF[[1]] ,rf2];
  eq1temp=HurUnifyTriadsCoord[eq1,rf2];
  eq2temp=HurUnifyTriadsCoord[eq2,rf2];
  {eq1temp[[1]]==0,eq1temp[[2]]==0,eq1temp[[3]]==0,eq2temp[[1]]==0,eq2temp[[2]]==0,eq2temp[[3]]==0}
  )

HurGetNEForm[rf1_, rf2_] := (
  eq1=HurNEResultantForce[rf1,rf2]-HurVectorDiff[HurGlobalLinearMomentum[[HurGetIndexGlobalRF[ rf1 ] ]],HurGlobalRF[[1]] ,rf2];
  eq2=HurNEResultantMoment[rf1,rf2]-HurVectorDiff[HurGlobalAngularMomentum[[HurGetIndexGlobalRF[ rf1 ] ]],HurGlobalRF[[1]] ,rf2];
  eq1temp=HurUnifyTriadsCoord[eq1,rf2];
  eq2temp=HurUnifyTriadsCoord[eq2,rf2];
  {eq1temp[[1]],eq1temp[[2]],eq1temp[[3]],eq2temp[[1]],eq2temp[[2]],eq2temp[[3]]}
  )

HurNEEquation[] := (nrfs=HurGetNumGlobalRF[];
  If[
    nrfs===1
    ,
    Print["There is only one reference frame!"]
    ,
    Do[
      If[
        Position[ HurGlobalForce[[;; , 1]], HurGlobalRF[[i]] ] === {} && Position[HurGlobalMoment[[;; , 1]], HurGlobalRF[[i]] ] === {}
        ,
        Null
        ,
        HurGetLinearMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
        HurGetAngularMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
        HurGlobalNEEquation[[i]]=HurGetNEEquation[HurGlobalRF[[i]],HurGlobalRF[[i]]]        
        
        (*
          Print[HurGetNEEquation[HurGlobalRF[[i]],HurGlobalRF[[i]]]]
          AppendTo[ HurGlobalNEEquation[[i]], HurGetNEEquation[HurGlobalRF[[i]],HurGlobalRF[[i]]] ]
          HurGlobalNEEquation[[i]]=Flatten[ HurGlobalNEEquation[[i]] ]
        *)
        ]
      ,
      {i,2,nrfs}
      ];
    (*HurGlobalNEEquation[[3]]=temp[[2]];*)
    HurGlobalNEEquation
    ]
  )


HurNEForward[] := (
  tempEq=Flatten[HurGlobalNEEquation];
  neq=Length[tempEq];

  temp=Flatten[
      Table[
      First[ tempEq[[i]] ]
      ,
      {i, neq}
      ]
    ];
  AMatrix=Grad[temp, HurGlobalVariableList];
  bvector=temp/. Table[HurGlobalVariableList[[i]] -> 0, {i, Length[HurGlobalVariableList]}];
  List[AMatrix,-bvector]
  )
(*
HurNEForward[] := (index=Position[HurGlobalNEEquation, {__}, 1];
  temp=Flatten[
      Table[
      HurGetNEForm[ HurGlobalRF[[i]] , HurGlobalRF[[i]] ]        
      ,
      {i, Flatten[index]}
      ]
    ];
  AMatrix=Grad[temp, HurGlobalVariableList];
  bvector=temp/. Table[HurGlobalVariableList[[i]] -> 0, {i, Length[HurGlobalVariableList]}];
  List[AMatrix,-bvector]
  )
*)



HurDefineVariableList[f__] := (
  HurGlobalVariableList=List[f];
  )

HurSolveNEInverse[] := (temp=Solve[Flatten[HurGlobalNEEquation], HurGlobalVariableList];
  temp1=If[HurGlobalSimplify, Simplify[temp], temp];
  temp1
  )

HurMakeSymmetricMatrix[list_] := {{list[[1]], list[[2]], list[[3]] }, {list[[2]], list[[4]], list[[5]] },{list[[3]],list[[5]],list[[6]]}}












(* Euler-Lagrange Mechanics *)

HurDefineOtherPotentialE[rf_, pe_] := (
  HurGlobalOtherPotentialE[[ HurGetIndexGlobalRF[ rf ] ]] = pe;
  Total[HurGlobalOtherPotentialE]
  )


HurDefineNonConservativeForces[f__] := (
  HurGlobalNonConservativeForces=Flatten[ List[f] ];
  HurGlobalNonConservativeForces
  )

HurGetKineticE[rf__] := (rfs=Flatten[List[rf]];narg=Length[rfs];
  Do[
    temp=1/2*HurGlobalMass[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]] * HurDot[ HurGlobalCOMVel[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]], HurGlobalCOMVel[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]] ]+1/2*HurDot[ HurGlobalAngularVel[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]],HurGlobalAngularMomentum[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]] ];
    temp1=If[HurGlobalSimplify, Simplify[temp], temp];
    HurGlobalKineticE[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]] = temp1;
    ,
    {i,narg}
    ];
  Total[HurGlobalKineticE]
  )  

HurGetPotentialE[rf__] := (rfs=Flatten[List[rf]];narg=Length[rfs];
  Do[
    temp=Global`g*HurGlobalMass[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]]*HurDot[HurGlobalCOMPos[[HurGetIndexGlobalRF[ rfs[[i]] ] ]], HurGlobalVertical[[1]] ];
    temp1=If[HurGlobalSimplify, Simplify[temp], temp];
    HurGlobalPotentialE[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]] = temp1;
    ,
    {i,narg}
    ];
  Total[HurGlobalPotentialE]
  )

HurGetLagrangian[rf__] := (HurGetKineticE[rf];HurGetPotentialE[rf];HurGlobalLagrangian=HurGlobalKineticE-HurGlobalPotentialE-HurGlobalOtherPotentialE;
  Total[HurGlobalLagrangian]
  )

HurGetELEquation[gc__] := (gcs=Flatten[ List[gc] ];ngcs=Length[gcs];
  L=Total[HurGlobalLagrangian];
  
  Do[ 
    temp=HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gcs[[i]] ][[1]][[1]] ]] = D[ D[ L, D[gcs[[i]], Global`t] ], Global`t ] - D[ L , gcs[[i]] ];
    If[HurGlobalSimplify, Simplify[temp], temp]
    , 
    {i,ngcs} 
    ];
  HurGlobalELEquation
  )

(*
HurGetELEquation[rf__] := (HurGetLagrangian[rf];
  gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];ngcs=Length[gcs];rfs=Flatten[List[rf]];nrfs=Length[rfs];
  L=Total[
    Table[ 
      HurGlobalLagrangian[[ HurGetIndexGlobalRF[ rfs[[i]] ] ]]
      , 
      {i,nrfs} 
    ]
  ];
  HurGlobalELEquation=Table[ 
    temp=D[ D[ L, D[gcs[[i]], Global`t] ], Global`t ] - D[ L , gcs[[i]] ];
    If[HurGlobalSimplify, Simplify[temp], temp]
    , 
    {i,ngcs} 
    ];
  HurGlobalELEquation
  )
*)

HurELEquation[] := (nrfs=HurGetNumGlobalRF[];
  If[
    nrfs===1
    ,
    Print["There is only one reference frame!"]
    ,
    HurGetLagrangian[ HurGlobalRF[[2 ;; nrfs]] ];
    HurGetELEquation[ HurGlobalGeneralizedCoordinates ];
    HurGetMMatrix[];
    HurGetCMatrix[];
    HurGetGVector[];
    HurGlobalELEquation
    ]
  )

HurDefineLambda[] := (
  m=Length[HurGlobalConstraints];
  HurGlobalLambda=Table[ Symbol["lambda" <> ToString[i]] ,{i,m}];
  HurGlobalLambda
  )

HurConstrainedELEquation[] := (
  HurELEquation[];
  HurDefineConstrainedJacobian[];
  HurDefineLambda[];
  m=Length[HurGlobalLambda];
  HurGlobalGeneralizedConstrainingForce=Transpose[HurGlobalConstrainedJacobian].HurGlobalLambda;
  HurGlobalConstrainedELEquation=HurGlobalELEquation-HurGlobalGeneralizedConstrainingForce;
  HurGlobalConstrainedELEquation
  )


HurELInverse[] := (
  n=Length[HurGlobalELEquation];
  tempequations=Flatten[ Table[HurGlobalELEquation[[i]]==HurGlobalNonConservativeForces[[i]],{i,n} ] ];
  tempvariables=Flatten[ Table[ D[ HurGlobalGeneralizedCoordinates[[k]],Global`t,Global`t ], {k,n} ] ];
  temp=Solve[tempequations,tempvariables];
  If[HurGlobalSimplify, Simplify[temp], temp]
  )



HurConstrainedELInverse[] := (
  n=Length[HurGlobalConstrainedELEquation];
  m=Length[HurGlobalConstraints];
  (* GCdots=Table[ D[ HurGlobalGeneralizedCoordinates[[k]] , Global`t ] , {k,n} ]; *)
  GCddots=Table[ D[ HurGlobalGeneralizedCoordinates[[k]] , Global`t , Global`t ] , {k,n} ];
  Zeros=Table[0, {i, n}];
  HurGlobalConstrainedModified=Table[0, {i, m}];

  Do[
    temp=Grad[ HurGlobalConstraints[[i]],GCddots ];
    If[
        temp===Zeros
        ,
        tempp=D[HurGlobalConstraints[[i]],Global`t];
        temp=Grad[ tempp , GCddots ];
        If[
            temp===Zeros
            ,
            tempp=D[tempp,Global`t];    (* original constrains contain no dots. holonomic constraints *)
            HurGlobalConstrainedModified[[i]]=tempp;
            ,            
            HurGlobalConstrainedModified[[i]]=tempp;  (* original constrains contain single dots *)
          ]
        ,
        HurGlobalConstrainedModified[[i]]=HurGlobalConstraints[[i]];  (* original constrains contain double dots *)
      ]
    ,
    {i,m}
    ];

  tempequations=Flatten[ List[ 
    Table[HurGlobalConstrainedELEquation[[i]]==HurGlobalNonConservativeForces[[i]],{i,n} ]
    ,
    Table[HurGlobalConstrainedModified[[i]]==0,{i,m} ]
    ] 
  ];

  tempvariables=Flatten[ List[ Table[ D[ HurGlobalGeneralizedCoordinates[[k]],Global`t,Global`t  ], {k,n} ]   ,HurGlobalLambda ] ];
  tempk=Solve[tempequations,tempvariables];
  If[HurGlobalSimplify, Simplify[tempk], tempk]
  )






HurGetMMatrix[] := (gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];narg=Length[gcs];
  HurGlobalMMatrix=
  Table[ 
    temp=D[ HurGlobalELEquation[[i]], D[gcs[[j]], Global`t, Global`t] ];
    If[HurGlobalSimplify, Simplify[temp], temp]
    ,{i, narg},{j, narg}];
  HurGlobalMMatrix
  )

HurGetCMatrix[] := (gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];narg=Length[gcs];
  HurGlobalCMatrix=Table[
    Total[
      Table[
        1/2 (D[HurGlobalMMatrix[[k, j]], gcs[[i]]] + D[ HurGlobalMMatrix[[k, i]], gcs[[j]]] - D[ HurGlobalMMatrix[[i, j]], gcs[[k]]]) D[gcs[[i]],Global`t]
        , 
        {i, narg}
      ] 
    ]
    , 
    {k, narg}, {j, narg}
    ];
  HurGlobalCMatrix
  )

HurGetGVector[] := (gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];narg=Length[gcs];
  HurGlobalGVector=Table[
    temp=D[ Total[HurGlobalPotentialE], gcs[[i]] ];
    If[HurGlobalSimplify, Simplify[temp], temp]
      ,
      {i,narg}
    ];
  HurGlobalGVector
  )


HurDefineConstraints[con__] := (cons=Flatten[ List[ con ] ];ncon=Length[ cons ];
  Do[
      AppendTo[ HurGlobalConstraints,cons[[i]] ];
      ,
      {i,ncon}
    ];
  HurGlobalConstraints
  )


HurDefineConstrainedJacobian[] := (m=Length[HurGlobalConstraints];n=Length[HurGlobalGeneralizedCoordinates];
  GCdots=Table[ D[ HurGlobalGeneralizedCoordinates[[k]] , Global`t ] , {k,n} ];
  Zeros=Table[0, {i, n}];
  HurGlobalConstrainedJacobian=Table[0,{i,m},{k,n}];

  Do[
    temp=Grad[ HurGlobalConstraints[[i]],GCdots ];
    If[
      temp === Zeros
      ,
      tempp=D[HurGlobalConstraints[[i]],Global`t];
      temp=Grad[ tempp , GCdots ];
      HurGlobalConstrainedJacobian[[i]]=temp;
      ,
      HurGlobalConstrainedJacobian[[i]]=temp;
      ]
    ,
    {i,m}
    ];
  HurGlobalConstrainedJacobian
  )



(*
HurGetGMatrix[eq_, gc_] := (gcs=Flatten[ List[gc] ];narg=Length[gcs];
  Table[eq[[i]], {i, Length[eq]}] /. 
    Table[D[gcs[[i]], Global`t] -> 0, {i, narg}] /. 
   Table[D[D[gcs[[i]], Global`t], Global`t] -> 0, {i, narg}]// Simplify
   )
*)

HurDefineGeneralizedCoordinates[gc__] := (gcs=Flatten[ List[ gc ] ]; ngcs=Length[gcs];
  HurGlobalGeneralizedCoordinates=gcs;
  HurGlobalELEquation=Table[0,{i,ngcs}];
  HurGlobalNonConservativeForces=Table[0,{i,ngcs}];
  HurGlobalConstrainedELEquation=Table[0,{i,ngcs}];
  )


HurSaveData[filename__] := (filenames=Flatten[ List[ filename ] ];nargs=Length[filenames];
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalMass","HurGlobalInertia","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify"};
  If[
      nargs===1
      ,
      var=tempvar1;
      ,
      var1=Table[filenames[[i]],{i,2,nargs}];
      var=Join[tempvar1,var1];
    ];
  Save[filenames[[1]], Evaluate[ Table[ var[[j]] , {j,Length[var]} ]] ];
(*  Do[
      Save[filenames[[1]], Evaluate[ var[[j]] ] ];
      ,
      {j,Length[var]}
    ]
*)
  )

HurDumpSaveData[filename__] := (filenames=Flatten[ List[ filename ] ];nargs=Length[filenames];
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalMass","HurGlobalInertia","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify"};
  If[
      nargs===1
      ,
      var=tempvar1;
      ,
      var1=Table[filenames[[i]],{i,2,nargs}];
      var=Join[tempvar1,var1];
    ];
  DumpSave[filenames[[1]], Evaluate[ Table[ var[[j]] , {j,Length[var]} ]] ];
  )
  
HurLoadData[filename_] := (
  Get[filename];
  )

End[];

EndPackage[];

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

0.4
Added NE Equation
Added Automatic Procedures

0.5
Added EL Equation
Added Automatic Procedures

0.6
Added Saving and Loading all variables (both internal and user-defined)

  *)