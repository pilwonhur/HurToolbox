(* ::Package:: *)

(* 	
	Hur - software for modeling and analysis of multibody systems    
	Copyright 2019, Pilwon Hur    

	Department of Mechanical Engineering 	
	Texas A&M University   	

	Hur Toolbox    
	Version 0.7, compatible with Mathematica 10.0, 11.0 and 12.0	

	Revision date: September 7, 2019 


  Revision
  0.9.5
  DH parameter handling
  HurDHTable[], HurGlobalDHTable
  HurDHInertia[], HurGlobalDHInertia
  HurDHKinematics[] 
  HurDHELEquation[]

  0.9.4
  HurDefineDCM and HurDefineDCMRelative now also accept angle and axis.
  In other words, you can do HurDefineDCM[a,q1[t],{0,0,1}].
  This is the same as HurDefineDCM[a,RotationMatrix[q1[t],{0,0,1}]].
  If HurDefineDCM[a,q1[t],{0,0,1}] is used, angular velocity can be expressed simple. 
  HurGlobalAngularVelRel and HurGlobalAngularVelAbs are added.
  You need to run HurGetAngularVel1[]

  0.9.3
  HurSimplifyVariablesTimed[], HurFullSimplifyVariablesTimed[] are added.
  HurGlobalConstrainedELEquation is all replaced by HurGlobalELEquation

  0.9.2
  Please add ";" at the end of each line.
  "Set::write: Tag Times in 0 {} is Protected." was fixed by adding ";" at the 
  PE=Total[HurGlobalPotentialE];

  0.9.1
  Kane

  0.9.0
  revised HurUnifyTriads[]. It now uses HurCoordTriads

  0.8.9
  Added the explanation for HurDumpSaveData, HurSaveData, and HurLoadData
  
  0.8.8
  HurGetGVector[] is updated to iinclude spring energy and RayleighDissipative energy
 
  0.8.7
  HurGetJacobian[] is added

  0.8.6
  HurConstrainedNEInverse[] is added
  HurTurnOnSimplify[] is added
  HurTurnOffSimplify[] is added

  0.8.5 
  HurResetConstraints[] is added

  0.8.4
  HurSetELEquation[] is added
  HurSetLagrangian[] is added
  HurGetELEquationFromLagrangian[] is added

  0.8.3
  HurMakeSymmetricMatrix[] is modified to accept a list of 6 elements or simply 6 separate elements.
  HurMatrixVectorProduct[] is added to facilitate angular momentum calculation. H=Iw

  0.8.2
  HurNEInverse[] is added. HurNEInverse[] is the same as HurSolveNEInverse[]. HurNEInverse[] is simply matching HurELInverse[].

  0.8.1
  HurCoordTriads: accepts one or two inputs.
  Angular Momentum and Linear Momentum are automatically computed for EL equations.
  HurDefineVariableList for EL equations. HurDefineVariableList was already availalbe for NE equation.
  HurResetVariableList is added.
  HurConstrainedELInverse[] is updated.
  Nonconservativeforces are added into the EL equation directly. It used to be shown only in the inverse equation.

  0.8.0
  HurNorm: accepts one or two inputs.
  HurNormSquare: accepts one or two inputs.

  0.7.9

  0.6
  Added Saving and Loading all variables (both internal and user-defined)

  0.5
  Added EL Equation
  Added Automatic Procedures

  0.4
  Added NE Equation
  Added Automatic Procedures

  0.3
  Added a function to compute the DCM between two reference frames

  0.2
  Handled time derivative problem when t was treated as priviate variable. Forced to treat as global variable.
  https://mathematica.stackexchange.com/questions/114769/derivative-from-my-package-function-returns-0
  Added a function to compute the angular velocity directly from the rotation matrix.
  Added a Vector Differentiation Function

  0.1
  Added the basic functions for vector analysis and reference frame manipulations.
*)  

BeginPackage["HurToolbox`"];

(* Usage statements *)
HurInitialize::usage="This procedure resets all global variables.";

$VERSION$ = "0.9.5";
$EMAIL$ = "pilwonhur@tamu.edu";
Print["Hur Toolbox for modeling and analysis of multibody systems ", $VERSION$, ". \nCopyright 2019 Pilwon Hur\nDepartment of Mechanical Engineering\nTexas A&M University\nAll rights reserved.\nEmail questions, comments, or concerns to ", $EMAIL$, "."];


HurDefineRF::usage="HurDefineRF[rf__] defines a reference frame. Ex) HurDefineRF[b]: to define a reference frame B. HurDefineRF[a,b]: to define reference frames a, and b at the same time.";
HurGetNumGlobalRF::usage="HurGetNumGlobalRF[] returns the number of reference frames defined globally. Ex) n=HurGetNumGlobalRF[]";
HurGetIndexGlobalRF::usage="HurGetIndexGlobalRF[rf_] returns the index assigned for the reference frame. Ex) n=HurGetIndexGlobalRF[b]";
HurTurnOnSimplify::usage="HurTurnOnSimplify[] sets HurGlobalSimplify=True. Then every symbolic expression will be simplified. This may delay the computation.";
HurTurnOffSimplify::usage="HurTurnOffSimplify[] sets HurGlobalSimplify=False. Then every symbolic expression will be shown without simplification. You may try this if the computation time is significantly delayed.";
HurDefineDCM::usage="HurDefineDCM[rf_, dcm_] defines a direction cosine matrix of the reference frame with respect to the world reference frame. Ex) HurDefineDCM[b, RotationMatrix[theta[t]+phi[t],{0,0,1}]]";
HurDefineDCMRelative::usage="HurDefineDCMRelative[rf1_, rf2_, dcm_] defines a direction cosine matrix of the reference frame rf1 with respect to rf2. Ex) HurDefineDCMRelative[b, a, RotationMatrix[phi[t],{0,0,1}]]";
HurDefineMass::usage="HurDefineMass[rf_, m_] defines the mass of the reference frame. Ex) HurDefineMass[b, m]";
HurDefineInertia::usage="HurDefineInertia[rf_, II_] defines the inertia matrix of the reference frame about the center of mass in terms of body reference frame. Ex) HurDefineInertia[b, {Ixx,Ixy,Ixz,Iyy,Iyz,Izz}]";
HurGetAngularVel::usage="HurGetAngularVel[rf1_,rf2_] returns the angular velocity of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It uses skew symmetric matrix Sw such that Sw=R_dot x R'. Ex) w=HurGetAngularVel[b,d]";
HurGetAngularVel1::usage="HurGetAngularVel1[] returns the angular velocity of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It uses skew symmetric matrix Sw such that Sw=R_dot x R'. Ex) w=HurGetAngularVel[b,d]";
HurGetAngularAcc::usage="HurGetAngularAcc[rf1_,rf2_] returns the angular acceleration of the reference frame rf1 with respect to the world reference frame n using the triads of rf2. It directly differentiate the angular velocity. Ex) w=HurGetAngularAcc[b,d]";
HurGetRelativeDCM::usage="HurGetRelativeDCM[rf1_, rf2_] returns the direction cosine matrix between rf1 and rf2. It is also equivalent to the rotation matrix of rf1 with respect to rf2. Ex) dcm=HurGetRelativeDCM[b,c]";
HurUnifyTriadPool::usage="HurUnifyTriadPool[rf1_, rf2_] returns triads for rf1 expressed in terms of rf2. In other words, it returns the x,y,z axes of rf1 expressed in terms of rf2. Note that when rf1 is defined via HurDefineRF and HurDefineDCM, rf1 is expressed in terms of Newtonian reference frame. This is used only for the internal computation purpose. It is also equivalent to HurGetRelativeDCM[rf1_, rf2_]. Ex) HurUnifyTriadPool[a, b]";
HurUnifyTriads::usage="HurUnifyTriads[v_, rf_] is used to represent the given vector with respect to rf. Ex) HurUnifyTriads[v, b]";
HurUnifyTriadsCoord::usage="HurUnifyTriadsCoord[v_, rf_] is used to display the coordinates of the given vector with respect to rf. Ex) HurUnifyTriadsCoord[v, b]";
HurCoordTriads::usage="HurCoordTriads[coord__] is used to represent the coordinate into the vector form. HurCoordTriads[coord_]: coord contains the 4 elements. The first 3 are the vector coordinates and the 4th element is RF. Ex) If v=x n1 + y n2 + z n3, then coord={x,y,z,n}. The following is also possible. HurCoordTriads[coord_,RF]: If v=x n1 + y n2 + z n3, then coord={x,y,z}, and RF=n";
HurCross::usage="HurCross[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a vector with respect to rf. Ex) HurCross[v1, v2, b]";
HurDot::usage="HurDot[v1_, v2_] is used to perform the dot product of the two vectors. v1 and v2 can be in different RFs. Also, HurDot does not require RF information since the output is a scalar. Ex) HurDot[v1, v2]";
HurNorm::usage="HurNorm[v__] returns the norm of the vector. HurNorm[v] first convert v into n, then returns norm of v. HurNorm[v,b] converts v into b and returns norm of v. If your vector is already expressed in b, then please use HurNorm[v,b]. Vector v can have triads of mixed RFs. Also, note that both will return the same norm. You may need to run Simplify"
HurNormSqure::usage="HurNormSqure[v_] returns the norm squared of the vector. Vector v can have triads of mixed RFs."
HurCrossCoord::usage="HurCrossCoord[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a coordinate with respect to rf. It is equivalent to the following: HurUnifyTriadsCoord[ HurCross[v1_, v2_, rf_], rf_] Ex) HurCrossCoord[v1, v2, b]";
HurVectorDiff::usage="HurVectorDiff[v_,rf1_,rf2_] is used to perform the vector differentiation. Vector v will be differentiated with respect to rf1. This function uses the vector differentiation formula for different RFs (i.e., rf2). Note that if a vector is differentiated w.r.t. time, then it will be the same as HurVectorDiff[v_,n,rf2_].";
HurAppendRF2Coord::usage="HurAppendRF2Coord[coord_, rf_] explicitly specify the RF information to the given coordinate (i.e., 3 numbers) without RF. Regardless of the size of coord_, HurAppendRF2Coord[coord_, rf_] takes the first 3 components of coord_ and attach rf_ to it. Usually, this function is used for internal usage.";
HurDefineForces::usage="HurDefineForces[rf_, force_, r_] assigns the force to the RF rf at the position r. r vector is relative to the COM of RF rf.";
HurResetForces::usage="HurResetForces[] resets all forces."
HurDefineMoments::usage="HurDefineMoments[rf_, moment_] assigns the moment to the RF rf";
HurResetMoments::usage="HurResetMoments[] resets all moments."
HurResetConstraints::usage="HurResetConstraints[] resets all constraints."
HurNEResultantForce::usage="HurNEResultantForce[rf1_, rf2_] returns the resultant force for the rf1 expressed with triads of rf2.";
HurNEResultantMoment::usage="HurNEResultantMoment[rf1_, rf2_] returns the resultant moment for the rf1 expressed with triads of rf2.";
HurGetAngularMomentum::usage="HurGetAngularMomentum[rf1_, rf2_] returns the angular momentum of rf1 expressed with triads of rf2.";
HurGetLinearMomentum::usage="HurGetLinearMomentum[rf1_, rf2_] returns the linear momentum of rf1 expressed with triads of rf2.";
HurGetLinearCOMVel::usage="HurGetLinearCOMVel[rf1_,rf2_] computes the linear velocity of the COM of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. Ex) w=HurGetLinearCOMVel[a,b]";
HurGetLinearCOMAcc::usage="HurGetLinearCOMAcc[rf1_,rf2_] computes the linear acceleration of the COM of the reference frame rf1 with respect to the world reference frame n, and then express it in terms of rf2. Ex) w=HurGetAngularAcc[a,b]";
HurGetLinearVelTwo::usage="HurGetLinearVelTwo[v_, w_, r_, rf1_, rf2_] computes the linear velocity of a point Q of RF A given the linear velocity (v_) of point P of RF A. w_ is Angular velocity of RF A, and r_ is the relative position vector from point P to point Q. rf1_ is the RF of the rigid body of interests (i.e., w, P and Q belong to) and rf2_ is for the unification of triads.";
HurGetLinearAccTwo::usage="HurGetLinearAccTwo[a_, alpha_, w_, r_, rf1_, rf2_] computes the linear acceleration of a point Q of RF A given the linear acceleration (a_) of point P of RF A. alpha_ is the Angular acceleration, of RF A, w_ is Angular velocity of RF A, and r_ is the relative position vector from point P to point Q. rf1_ is the RF of the rigid body of interests (i.e., w, P and Q belong to) and rf2_ is for the unification of triads.";
HurGetLinearVelOne::usage="HurGetLinearVelOne[vfixed_, vrel_, rf_] computes the linear velocity of a point P that is moving on a RF rf. vfixed_ is the velocity of the point P fixed on RF rf w.r.t. Newtonian RF. vrel_ is the relative velocity of point P w.r.t. RF rf.  rf_ is for the unification of triads.";
HurGetLinearAccOne::usage="HurGetLinearAccOne[afixed_, arel_, w_, vrel_, rf1_, rf2_] computes the linear acceleration of a point P that is moving on a RF rf. afixed_ is the acceleration of the point P fixed on RF rf1 w.r.t. Newtonian RF. arel_ is the relative acceleration of point P w.r.t. RF rf1. vrel_ is the relative velocity of point P w.r.t. RF rf1. w_ is the angular velocity of RF rf1, rf2_ is for the unification of triads.";
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
HurSetLagrangian::usage="HurSetLagrangian[lag_,rf_] can let you specify Lagrangian by the user, not via automatic procedure.";
HurDefineVertical::usage="HurDefineVertical[v_]";
HurGetELEquation::usage="HurGetELEquation[gc__]";
HurSetELEquation::usage="HurSetELEquation[eq_,gc_]";
HurGetJacobian::usage="HurGetJacobian[vec_,rf1_,rf2_] computes Jacobian matrix. vec_ can be either a position vector of a point of interest or velocity vector. rf1_ is the RF for angular velocity of your interest. rf2_ is expression of your vec_. rf2_ is usually n. However, it can also be a or b depending on your convenience."; 
HurGetMMatrix::usage="HurGetMMatrix[] returns the inertia matrix.";
HurGetMMatrix1::usage="HurGetMMatrix1[] returns the inertia matrix.";
HurGetCMatrix::usage="HurGetCMatrix[] returns the matrix for Coriolis and centrifugal forces";
HurGetGVector::usage="HurGetGVector[] returns the gravity vector. Note that it will also return other forces due to spring and Rayleigh Dissipative (viscous damping) forces.";
HurNEEquation::usage="HurNEEquation[]";
HurDefineVariableList::usage="HurDefineVariableList[f__]";
HurResetVariableList::usage="HurResetVariableList[] resets all variable list for Inverse of either NE or EL equations";
HurSolveNEInverse::usage="HurSolveNEInverse[] returns the the variables defined in HurGlobalVariableList. Variables will be either forces, torques, or double dots of generalized coordinates.";
HurNEInverse::usage="HurNEInverse[] is exactly the same as HurSolveNEInverse[]. HurNEInverse[] is to match with HurELInverse[]";
HurConstrainedNEInverse::usage="HurConstrainedNEInverse[] computes the inverse problems by including the constraints in addition to EOMs. Please make sure that constraints are defined via HurDefineConstraints[] and the variable lists are appropriately chosen.";
HurNEForward::usage="HurNEForward[]";
HurMakeSymmetricMatrix::usage="HurMakeSymmetricMatrix[list_]";
HurMatrixVectorProduct::usage="HurMatrixVectorProduct[mat_,vec_,rf_]"
HurMatrixVectorProductTriads::usage="HurMatrixVectorProductTriads[mat_,vec_,rf_]"
HurDefineGeneralizedCoordinates::usage="HurDefineGeneralizedCoordinates[f__]";
HurELEquation::usage="HurELEquation[] computes the EL equations via automatic procedure.";
HurELEquation1::usage="HurELEquation1[] computes the EL equations via automatic procedure.";
HurGetELEquationFromLagrangian::usage="HurGetELEquationFromLagrangian[] computes the EL equations from the user-provided Lagrangian. Make sure that you have HurGlobalLagrangian defined. Otherwise, please run HurSetLagrangian[lag, gc] first.";
HurDefineConstraints::usage="HurDefineConstraints[con__]";
HurDefineConstrainedJacobian::usage="HurDefineConstrainedJacobian[]";
HurConstrainedELEquation::usage="HurConstrainedELEquation[]";
HurConstrainedELInverse::usage="HurConstrainedELInverse[]";
HurELInverse::usage="HurELInverse[]";
HurDefineLambda::usage="HurDefineLambda[]";
HurDefineNonConservativeForces::usage="HurDefineNonConservativeForces[f__]"; 
HurDefineOtherPotentialE::usage="HurDefineOtherPotentialE[rf_, pe_] accepts the additional potential energy other than gravity in Lagrangian mechanics. Example includes spring. If intrinsic generalized coordinates (e.g., joint angles) are used, elastic energy will be uniquely assigned to a RF. However, if your GC's are extrinsic, then assignment of elastic energy to an RF is ambiguous. However, it doesn't matter since we simply need the total sum of the elastic energies. Therefore, simply assign the elastic energy to an any reasonable RF."; 
HurDefineRayleighDissipationE::usage="HurDefineRayleighDissipationE[rf_, de_] accepts the velocity-proportional frictional forces in Lagrangian mechanics. Example includes viscous damping friction."; 
HurGetInertiaTensor::usage="HurGetInertiaTensor[rf_]"; 
HurProductMatVec::usage="HurProductMatVec[mat_,vec_,rf_]"; 
HurDefineGeneralizedSpeedsConstraints::usage="HurDefineGeneralizedSpeedsConstraints[gs__] defines the generalized speeds for Kane's method"; 
HurDHTable::usage="HurDHTable[dh_]";
HurDHInertia::usage="HurDHInertia[data_]";
HurSimplifyVariablesTimed::usage="HurSimplifyVariablesTimed[var_,time_] simplifies all the elements of the provided var_ within the provided time_ for each element. If timed out, it will return its original expression"; 
HurFullSimplifyVariablesTimed::usage="HurFullSimplifyVariablesTimed[var_,time_] fully simplifies all the elements of the provided var_ within the provided time_ for each element. If timed out, it will return its original expression"; 
HurDumpSaveData::usage="HurDumpSaveData[filename__] Please use .mx for the extension of the filename. It will save variables in binary (unreadable) expression, is very fast to load (with large data). However, this binary data are platform-specific. If saved in Mac, it cannot be used in Windows or Linux.";
HurSaveData::usage="HurSaveData[filename__]. Please use .m for the extension of the filename. It will save variables in (readable) portable expression (or ascii format), is very slow to load (with large data). This data are platform-independent. You can use in any platforms.";
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
  HurGlobalNonConservativeForces = {0}; HurGlobalOtherPotentialE = {0}; HurGlobalRayleighDissipationE = {0};
  HurGlobalListTriads = {{n1,n2,n3}}; 
  HurGlobalTriadsConversion = {{n1->n1,n2->n2,n3->n3}};
  (*HurGlobalTriadsConversion = {0};*)
  HurGlobalSimplify = True; HurGlobalGeneralizedSpeedsConstraints = {};
  HurGlobalGeneralizedSpeeds = {}; HurGlobalKaneEquation = {0};
  HurGlobalTemp = {0}; 
  HurGlobalAngularVelAbs = {0}; HurGlobalAngularVelRel = {{0,0,0}};
  HurGlobalDHTable = {0}; HurGlobalDHInertia = {0}; HurGlobalDHOrigin = {0};
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
        AppendTo[HurGlobalAngularVelAbs, "NA"];
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
        AppendTo[HurGlobalRayleighDissipationE, 0];
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

HurTurnOnSimplify[] := (
  HurGlobalSimplify=True
  )

HurTurnOffSimplify[] := (
  HurGlobalSimplify=False
  )

HurGetNumGlobalRF[] := Length[HurGlobalRF];

HurGetIndexGlobalRF[rf_] := 
  If[
      Flatten[Position[HurGlobalRF, rf]]==={},
      0,
      Flatten[Position[HurGlobalRF, rf]][[1]]
    ]

HurDefineDCM[rf_, rot__] := (HurDefineDCMRelative[rf, HurGlobalRF[[1]], rot];
  )

HurDefineDCMRelative[rf1_, rf2_, rot__] := (
  nrf=Length[HurGlobalRF];

  rots=List[rot];narg=Length[rots];
  If[Length[rots]===1
    ,
    (* traditional way when rotation matrix is provided *)
    dcm=rots[[1]];
    ,
    (* avoid duplicate *)
    pos=Flatten[Position[HurGlobalAngularVelRel[[;; , 1]], rf1]];
    If[Length[pos]===0
      ,
      triad=HurGlobalListTriads[[ HurGetIndexGlobalRF[rf2] ]][[ Flatten[Position[rots[[2]],1]] ]][[1]];
      AppendTo[HurGlobalAngularVelRel, {rf1,rf2,D[rots[[1]],Global`t]*triad}];
      ,
      triad=HurGlobalListTriads[[ HurGetIndexGlobalRF[rf2] ]][[ Flatten[Position[rots[[2]],1]] ]][[1]];
      HurGlobalAngularVelRel[[ pos[[1]] ]]={rf1,rf2,D[rots[[1]],Global`t]*triad};
      ];  
    dcm=RotationMatrix[rots[[1]],rots[[2]]];
    ];
  HurGlobalDCM[[ HurGetIndexGlobalRF[ rf1 ] ]] = HurGlobalDCM[[ HurGetIndexGlobalRF[ rf2 ] ]].dcm;
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

HurResetAngularVelAbs[] :=(
  nang=Length[HurGlobalAngularVelAbs];
  Do[
    HurGlobalAngularVelAbs[[i]]="NA";
    ,
    {i,2,nang}
    ];
  )

HurGetAngularVel1[] := (
  nangvelrel=Length[HurGlobalAngularVelRel];
  nrfs=HurGetNumGlobalRF[];
  pos=Flatten[Position[HurGlobalAngularVelRel[[;; , 2]], HurGlobalRF[[1]]]];

  HurResetAngularVelAbs[];

  (* construct HurGlobalAngularVelAbs based on RF N first *)
  Do[
    index=HurGetIndexGlobalRF[ HurGlobalAngularVelRel[[pos[[1]],1 ]] ];
    HurGlobalAngularVelAbs[[ index ]] = HurGlobalAngularVelRel[[pos[[1]], 3 ]];
    ,
    {i,Length[pos]}
    ];

  (* construct HurGlobalAngularVelAbs based on HurGlobalAngularVelRel *)
  Do[
    If[HurGlobalAngularVelAbs[[i]]==="NA"
      ,
      pos1=Flatten[Position[HurGlobalAngularVelRel[[;; , 1]], HurGlobalRF[[i]] ]];    
      pos2=HurGetIndexGlobalRF[HurGlobalAngularVelRel[[pos1[[1]] , 2]] ];
      If[HurGlobalAngularVelAbs[[pos2]]==="NA",
        Null
        ,
        HurGlobalAngularVelAbs[[ i ]] = HurGlobalAngularVelAbs[[ pos2 ]] + HurGlobalAngularVelRel[[ pos1[[1]] , 3]] ;
        ]
      ,
      Null
      ]
    ,
    {i,2,nrfs}
    ]
  )

HurGetAngularAcc[rf1_,rf2_] := (
  alpha=HurUnifyTriads[HurCoordTriads[HurAppendRF2Coord[D[HurUnifyTriadsCoord[HurGetAngularVel[rf1,HurGlobalRF[[1]] ],HurGlobalRF[[1]] ][[1;;3]],Global`t],HurGlobalRF[[1]] ]],rf2];
  alpha1=If[HurGlobalSimplify, Simplify[alpha], alpha];
  HurSetAngularAcc[rf1,alpha1];
  alpha1
  )

HurGetLinearVelTwo[v_, w_, r_, rf1_, rf2_] := (
  HurUnifyTriads[v+HurCross[w,r,rf1],rf2] 
  )

HurGetLinearAccTwo[a_, alpha_, w_, r_, rf1_, rf2_] := (
  HurUnifyTriads[a+HurCross[alpha,r,rf1]+HurCross[w,HurCross[w,r,rf1],rf1],rf2] 
  )

HurGetLinearVelOne[vfixed_, vrel_, rf_] := (
  HurUnifyTriads[vfixed+vrel,rf] 
  )

HurGetLinearAccOne[afixed_, arel_, w_, vrel_, rf1_, rf2_] := (
  HurUnifyTriads[afixed+arel+2*HurCross[w,vrel,rf1],rf2] 
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

HurUnifyTriads[v_, rf_] := (
  (*
  temp=v /. Flatten[ HurGlobalTriadsConversion[[ HurGetIndexGlobalRF[ rf ] ]] ];
  Collect[ temp,HurGlobalListTriads[[ HurGetIndexGlobalRF[ rf ] ]] ]
  *)
  HurCoordTriads[ HurUnifyTriadsCoord[v, rf] ]
  )

HurUnifyTriadsCoord[v_, rf_] := (temp=v /. Flatten[ HurGlobalTriadsConversion[[ HurGetIndexGlobalRF[ rf ] ]] ];
  coord=Table[ D[ temp, HurGlobalListTriads[[ HurGetIndexGlobalRF[ rf ],i ]] ] , {i, 3}];
  coord1=If[HurGlobalSimplify, Simplify[coord], coord];
  HurAppendRF2Coord[coord1,rf]
  )

HurCoordTriads[v__] := (
  vs=List[v];narg=Length[vs];
  If[
      narg===1
      ,
      vec=Dot[ v[[1;;3]],HurGlobalListTriads[[ HurGetIndexGlobalRF[ v[[4]] ] ]] ];
      ,
      If[
        narg===2
        ,
        vec=Dot[ vs[[1]],HurGlobalListTriads[[ HurGetIndexGlobalRF[ vs[[2]] ] ]] ];
        ,
        Null
      ];
    ];
  vec
  )    

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

HurNormSqure[v__] := (
  vs=List[v];narg=Length[vs];
  If[
      narg===1
      ,
      vec=HurUnifyTriadsCoord[ v, HurGlobalRF[[1]] ] [[1;;3]];    
      ,
      If[
        narg===2
        ,
        vec=HurUnifyTriadsCoord[ vs[[1]], vs[[2]] ] [[1;;3]];    
        ,
        Null
      ];
    ];
  vec[[1]]^2+vec[[2]]^2+vec[[3]]^2
  )

HurNorm[v__] := (
  Sqrt[HurNormSqure[v]]
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

HurResetForces[] := (HurGlobalForce={})

HurDefineMoments[rf_, moment_] := (AppendTo[HurGlobalMoment,{rf,moment}];)

HurResetMoments[] := (HurGlobalMoment={})

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
        HurGlobalNEEquation[[i]]=HurGetNEForm[HurGlobalRF[[i]],HurGlobalRF[[i]]]        
        
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

HurResetVariableList[] := (HurGlobalVariableList={})

HurSolveNEInverse[] := (
  TempNEEquation=Flatten[HurGlobalNEEquation];
  n=Length[TempNEEquation];
  temp=Solve[Table[TempNEEquation[[i]]==0,{i,n}] , HurGlobalVariableList];
  temp1=If[HurGlobalSimplify, Simplify[temp], temp];
  temp1
  )

HurNEInverse[] := (
  TempNEEquation=Flatten[HurGlobalNEEquation];
  n=Length[TempNEEquation];
  m=Length[HurGlobalVariableList];

  (* temp=Solve[Table[TempNEEquation[[i]]==0,{i,n}] , HurGlobalVariableList] *)
  If[
      n===m
      ,
      temp=Solve[Table[TempNEEquation[[i]]==0,{i,n}] , HurGlobalVariableList];
      If[HurGlobalSimplify, Simplify[temp], temp]
      ,
      Print["Please make sure that you have " <> ToString[n] <> " variables in HurGlobalVariableList."]
    ]
  )

HurConstrainedNEInverse[] := (
  n=Length[Flatten[HurGlobalNEEquation]];
  m=Length[HurGlobalConstraints];
  k=Length[Flatten[HurGlobalGeneralizedCoordinates]];
  GCddots=Table[ D[ HurGlobalGeneralizedCoordinates[[i]] , Global`t , Global`t ] , {i,k} ];
  Zeros=Table[0, {i, k}];
  HurGlobalConstrainedModified=Table[0, {i, m}];

  Do[
    temp=Grad[ HurGlobalConstraints[[i]],GCddots ];
    If[
        temp===Zeros  (* if no ddots *)
        ,
        tempp=D[HurGlobalConstraints[[i]],Global`t];
        temp=Grad[ tempp , GCddots ];
        If[
            temp===Zeros   (* if no dots *)
            ,
            tempp=D[tempp,Global`t];    (* original constrains contain no dots. holonomic constraints *)
            HurGlobalConstrainedModified[[i]]=tempp;
            ,            
            HurGlobalConstrainedModified[[i]]=tempp;  (* original constrains contain single dots -> nonholonomic *)
          ]
        ,
        HurGlobalConstrainedModified[[i]]=HurGlobalConstraints[[i]];  (* original constrains contain double dots *)
      ]
    ,
    {i,m}
    ];

  TempNEEquation=Flatten[HurGlobalNEEquation];
  tempequations=Flatten[ List[ 
    Table[TempNEEquation[[i]]==0,{i,n} ]
    ,
    Table[HurGlobalConstrainedModified[[i]]==0,{i,m} ]
    ] 
  ];

  tempvariables=HurGlobalVariableList;
  
  tempk=Solve[tempequations,tempvariables];
  If[HurGlobalSimplify, Simplify[tempk], tempk]
  )

HurMakeSymmetricMatrix[lists__] := (
  list=Flatten[List[lists]];
  {{list[[1]], list[[2]], list[[3]] }, {list[[2]], list[[4]], list[[5]] },{list[[3]],list[[5]],list[[6]]}}
  )

HurMatrixVectorProduct[mat_,vec_,rf_] := mat.HurUnifyTriadsCoord[vec,rf][[1;;3]]

HurMatrixVectorProductTriads[mat_,vec_,rf_] := HurCoordTriads[HurMatrixVectorProduct[mat,vec,rf],rf]



(* Euler-Lagrange Mechanics *)

HurDefineOtherPotentialE[rf_, pe_] := (
  HurGlobalOtherPotentialE[[ HurGetIndexGlobalRF[ rf ] ]] = pe;
  Total[HurGlobalOtherPotentialE]
  )

HurDefineRayleighDissipationE[rf_, de_] := (
  HurGlobalRayleighDissipationE[[ HurGetIndexGlobalRF[ rf ] ]] = de;
  Total[HurGlobalRayleighDissipationE]
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

HurSetLagrangian[lag_,rf_] := (
    HurGlobalLagrangian[[ HurGetIndexGlobalRF[ rf ]  ]] = lag;
    Total[HurGlobalLagrangian]
  )

HurGetJacobian[vec_,rf1_,rf2_] := (ngcs=Length[HurGlobalGeneralizedCoordinates];
  (* check if vec_ includes q only or qdot. If q only, then it is a position vector. If qdot is included, then vec_ is velocity vector. *)
  GCdots=Table[ D[ HurGlobalGeneralizedCoordinates[[i]],Global`t], {i,ngcs}];
  Zeros=Table[0, {i, ngcs}];
  vec1=HurUnifyTriads[vec,rf2];
  vel1=HurVectorDiff[vec1, HurGlobalRF[[1]] ,rf2];
  omega=HurUnifyTriadsCoord[ HurGetAngularVel[rf1,rf2], rf2 ];
  vel=If[Grad[vec1,GCdots]===Zeros,HurUnifyTriadsCoord[ vel1, rf2 ],HurUnifyTriadsCoord[ vec1, rf2 ] ];
  List[
    Grad[ vel[[1]],GCdots ],
    Grad[ vel[[2]],GCdots ],
    Grad[ vel[[3]],GCdots ],
    Grad[ omega[[1]],GCdots ],
    Grad[ omega[[2]],GCdots ],
    Grad[ omega[[3]],GCdots ]
  ]
  )

HurSetELEquation[eq_,gc_] := (
    HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gc ] [[1]][[1]] ]] = eq;
    HurGlobalELEquation
  )

HurGetELEquationFromLagrangian[] := (gcs=HurGlobalGeneralizedCoordinates;ngcs=Length[gcs];
  L=Total[HurGlobalLagrangian];
  DE=Total[ HurGlobalRayleighDissipationE ];
  
  Do[ 
    temp=HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gcs[[i]] ][[1]][[1]] ]] = D[ D[ L, D[gcs[[i]], Global`t] ], Global`t ] - D[ L , gcs[[i]] ] + D[ DE, D[gcs[[i]], Global`t] ] - HurGlobalNonConservativeForces[[i]];
    If[HurGlobalSimplify, Simplify[temp], temp]
    , 
    {i,ngcs} 
    ];
  HurGlobalELEquation
  )

HurGetELEquation[gc__] := (gcs=Flatten[ List[gc] ];ngcs=Length[gcs];
  L=Total[HurGlobalLagrangian];
  DE=Total[ HurGlobalRayleighDissipationE ];
  
  Do[ 
    temp = D[ D[ L, D[gcs[[i]], Global`t] ], Global`t ] - D[ L , gcs[[i]] ] + D[ DE, D[gcs[[i]], Global`t] ] - HurGlobalNonConservativeForces[[i]];
    temp=If[HurGlobalSimplify, Simplify[temp], temp];
    HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gcs[[i]] ][[1]][[1]] ]] = temp;
    temp
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

HurELEquation[] := (
  nrfs=HurGetNumGlobalRF[];
  If[
    nrfs===1
    ,
    Print["There is only one reference frame!"]
    ,
    Do[
      If[
        HurGlobalMass[[i]] === 0 
        ,
        Null
        ,
        HurGetLinearMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
        HurGetAngularMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];        
        (*
          Print[HurGetNEEquation[HurGlobalRF[[i]],HurGlobalRF[[i]]]]
          AppendTo[ HurGlobalNEEquation[[i]], HurGetNEEquation[HurGlobalRF[[i]],HurGlobalRF[[i]]] ]
          HurGlobalNEEquation[[i]]=Flatten[ HurGlobalNEEquation[[i]] ]
        *)
        ]
      ,
      {i,2,nrfs}
      ];

    HurGetLagrangian[ HurGlobalRF[[2 ;; nrfs]] ];
    HurGetELEquation[ HurGlobalGeneralizedCoordinates ];
    
    HurGetMMatrix[];
    HurGetCMatrix[];
    HurGetGVector[];
    HurGlobalELEquation
    ]
  )

(* when purely mechanical systems without deformable *)
HurELEquation1[] := (
	HurGetMMatrix1[];
	HurGetCMatrix[];
	HurGetGVector[];

	gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];
	ngcs=Length[gcs];

	tempM=HurGlobalMMatrix.Transpose[List[D[gcs, Global`t, Global`t]]];
	tempC=HurGlobalCMatrix.Transpose[List[D[gcs, Global`t]]];
	tempG=HurGlobalGVector;
	tempQ=Transpose[List[HurGlobalNonConservativeForces]];
	tempTot=tempM+tempC+tempG-tempQ;
	Do[ 
		temp=If[HurGlobalSimplify, Simplify[tempTot[[i]]], tempTot[[i]]];
    	HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gcs[[i]] ][[1]][[1]] ]] = temp;
    	, 
    	{i,ngcs} 
    ];
    HurGlobalELEquation
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
  HurGlobalELEquation=HurGlobalELEquation-HurGlobalGeneralizedConstrainingForce;
  HurGlobalELEquation
  )

HurELInverse[] := (
  n=Length[HurGlobalELEquation];
  tempequations=Flatten[ Table[HurGlobalELEquation[[i]]==0,{i,n} ] ];
  If[
      Total[HurGlobalVariableList] === 0
      ,
      tempvariables=Flatten[ Table[ D[ HurGlobalGeneralizedCoordinates[[k]],Global`t,Global`t ], {k,n} ] ];
      ,
      tempvariables=HurGlobalVariableList;
    ];

  temp=Solve[tempequations,tempvariables];
  If[HurGlobalSimplify, Simplify[temp], temp]
  )



HurConstrainedELInverse[] := (
  n=Length[HurGlobalELEquation];
  m=Length[HurGlobalConstraints];
  (* GCdots=Table[ D[ HurGlobalGeneralizedCoordinates[[k]] , Global`t ] , {k,n} ]; *)
  GCddots=Table[ D[ HurGlobalGeneralizedCoordinates[[k]] , Global`t , Global`t ] , {k,n} ];
  Zeros=Table[0, {i, n}];
  HurGlobalConstrainedModified=Table[0, {i, m}];

  Do[
    temp=Grad[ HurGlobalConstraints[[i]],GCddots ];
    If[
        temp===Zeros  (* if no ddots *)
        ,
        tempp=D[HurGlobalConstraints[[i]],Global`t];
        temp=Grad[ tempp , GCddots ];
        If[
            temp===Zeros   (* if no dots *)
            ,
            tempp=D[tempp,Global`t];    (* original constrains contain no dots. holonomic constraints *)
            HurGlobalConstrainedModified[[i]]=tempp;
            ,            
            HurGlobalConstrainedModified[[i]]=tempp;  (* original constrains contain single dots -> nonholonomic *)
          ]
        ,
        HurGlobalConstrainedModified[[i]]=HurGlobalConstraints[[i]];  (* original constrains contain double dots *)
      ]
    ,
    {i,m}
    ];

  tempequations=Flatten[ List[ 
    Table[HurGlobalELEquation[[i]]==0,{i,n} ]
    ,
    Table[HurGlobalConstrainedModified[[i]]==0,{i,m} ]
    ] 
  ];

  If[
      Total[HurGlobalVariableList] === 0
      ,
      tempvariables=Flatten[ List[ Table[ D[ HurGlobalGeneralizedCoordinates[[k]],Global`t,Global`t  ], {k,n} ] , HurGlobalLambda ] ];
      ,
      tempvariables=HurGlobalVariableList;
    ];
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

HurGetMMatrix1[] := (
  gcs=Flatten[ List[HurGlobalGeneralizedCoordinates] ];narg=Length[gcs];
  tempM=
  	Table[ 
    Jac=HurGetJacobian[ HurGlobalCOMPos[[i+1]], HurGlobalRF[[i+1]], HurGlobalRF[[1]] ];
 	Jacv=Jac[[1;;3,;;]];
 	Jacw=Jac[[4;;6,;;]];
 	temp1=HurGlobalMass[[i+1]]*Transpose[Jacv].Jacv;
 	rot=HurGetRelativeDCM[ HurGlobalRF[[i+1]],HurGlobalRF[[1]] ];
 	II=HurGetInertiaTensor[ HurGlobalRF[[i+1]] ];
 	temp2=Transpose[Jacw].rot.II.Transpose[rot].Jacw;
    temp=temp1+temp2;
    If[HurGlobalSimplify, Simplify[temp], temp]
    ,{i, narg}];
  HurGlobalMMatrix=Total[tempM];
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
  PE=Total[HurGlobalPotentialE]+Total[HurGlobalOtherPotentialE];
  DE=Total[HurGlobalRayleighDissipationE];

  HurGlobalGVector=Table[
    temp=D[ PE , gcs[[i]] ] + D[ DE , D[gcs[[i]],Global`t] ] ; 
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

HurResetConstraints[] := (HurGlobalConstraints={})


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

(* Kane Method *)
HurDefineGeneralizedSpeedsConstraints[gs__] := (gss=Flatten[ List[ gs ] ]; ngss=Length[gss];
  HurGlobalGeneralizedSpeedsConstraints=gss;
  HurGlobalKaneEquation=Table[0,{i,ngss}];

  HurGlobalGeneralizedSpeeds=Table[ ToExpression["u" <> ToString[i] <> "[t]" ] ,{i,ngss} ];
  (* Global`t *)

  ngcs=Length[HurGlobalGeneralizedCoordinates]
  (* ngcs-ngss is the number of nonholonomic constraints *)


  )

(*
  identify all nonholonomic constraints
  reduce the generalized coordinates to remove the constraints
  identify RFs that have inertia
  convert vCom with generalized speeds
  convert Hc with generalized speeds


  *)


(* Robotics equations via DH parameters *)

HurDHTable[dh_] := (
	ndh=Length[dh];
	HurGlobalDHTable=dh;
  	(* define RF*)
  	Do[
  		HurDefineRF[ Symbol[ "rf" <> ToString[i] ]];
  		, 
  		{i,ndh}
  	];
  	(* define DCM *)
  	Do[
  		HurDefineDCMRelative[ HurGlobalRF[[ dh[[i,1]]+1 ]] , HurGlobalRF[[ dh[[i,1]] ]],  RotationMatrix[dh[[i,3]],{0,0,1}].RotationMatrix[dh[[i,6]],{1,0,0}] ];
  		, 
  		{i,ndh}
  	];
	)

HurDHInertia[data_] := (
	ndata=Length[data];
	dh=HurGlobalDHTable;
	HurGlobalDHOrigin=Table[0, Length[HurGlobalCOMPos]];
	If[Length[dh]===ndata
		,
		HurGlobalDHInertia=data;
		(* define COM *)
  		Do[
  			posOrgRel=HurGlobalListTriads[[i, 3]]*dh[[i,4]]+HurGlobalListTriads[[i+1, 1]]*dh[[i,5]];
  			HurGlobalDHOrigin[[i+1]]=HurGlobalDHOrigin[[i]]+HurUnifyTriads[posOrgRel,HurGlobalRF[[1]]];
  			posCOMRel=HurCoordTriads[data[[i,3]], HurGlobalRF[[i+1]] ];
  			HurDefineCOMPos[ HurGlobalRF[[ data[[i,1]]+1 ]] , HurGlobalDHOrigin[[i+1]]+posCOMRel ];
  			HurDefineMass[ HurGlobalRF[[ data[[i,1]]+1 ]],data[[i,2]] ];
  			HurDefineInertia[ HurGlobalRF[[ data[[i,1]]+1 ]],data[[i,4]] ]; 
  			, 
  			{i,ndata}
  		];
		,
		Print["The number of inertial information does not match with the number of links."];
	];
	)

HurSimplifyVariablesTimed[var_,time_] := (
  dims=Dimensions[var];
  totalnum = Product[dims[[i]], {i, Length[dims]}];
  vec = ArrayReshape[var, {totalnum}];
  
  Print["There are " , ToString[totalnum] , " components to be simplified."];
  curTimeConst=OptionValue[Simplify, TimeConstraint];
  If[
    curTimeConst<time
    ,
    Print["The current Maximum Time for Simplify is set to ", ToString[curTimeConst] , " seconds. However, you requested longer time (" , ToString[time] , " seconds). The Maximum Time will be temporarily set to " , ToString[time] , " seconds."];
    SetOptions[Simplify, TimeConstraint -> time];
    ,
    Null
    ]

  Do[
    temp = Timing[TimeConstrained[Simplify[ vec[[i]] ], time]];
    If[
      temp[[2]] === $Aborted
      ,
        Print[ToString[i], "/" , ToString[totalnum], " components could not be simplified in " , ToString[time] , " seconds." ];
      ,
        Print[ToString[i], "/" , ToString[totalnum], " components was successfully simplified in " , ToString[temp[[1]]], " seconds."];
        vec[[i]]=temp[[2]];
      ];
    ,
    {i,totalnum}
    ];
  SetOptions[Simplify, TimeConstraint -> curTimeConst];
  ArrayReshape[vec, dims]
  )

HurFullSimplifyVariablesTimed[var_,time_] := (
  dims=Dimensions[var];
  totalnum = Product[dims[[i]], {i, Length[dims]}];
  vec = ArrayReshape[var, {totalnum}];
  
  Print["There are " , ToString[totalnum] , " components to be fully simplified."];
  curTimeConst=OptionValue[FullSimplify, TimeConstraint];
  If[
    curTimeConst<time
    ,
    Print["The current Maximum Time for FullSimplify is set to ", ToString[curTimeConst] , " seconds. However, you requested longer time (" , ToString[time] , " seconds). The Maximum Time will be temporarily set to " , ToString[time] , " seconds."];
    SetOptions[FullSimplify, TimeConstraint -> time];
    ,
    Null
    ]

  Do[
    temp = Timing[TimeConstrained[FullSimplify[ vec[[i]] ], time]];
    If[
      temp[[2]] === $Aborted
      ,
        Print[ToString[i], "/" , ToString[totalnum], " components could not be fully simplified in " , ToString[time] , " seconds." ];
      ,
        Print[ToString[i], "/" , ToString[totalnum], " components was successfully fully simplified in " , ToString[temp[[1]]], " seconds."];
        vec[[i]]=temp[[2]];
      ];
    ,
    {i,totalnum}
    ];
  SetOptions[FullSimplify, TimeConstraint -> curTimeConst];
  ArrayReshape[vec, dims]
  )

HurSaveData[filename__] := (filenames=Flatten[ List[ filename ] ];nargs=Length[filenames];
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalMass","HurGlobalInertia","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalRayleighDissipationE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify"};
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
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalMass","HurGlobalInertia","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalRayleighDissipationE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify"};
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