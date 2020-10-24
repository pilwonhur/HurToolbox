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
  2.0.5
  Revised HurNEInverse[]. When 2D NE EOMs are derived, HurNEInverse[] does not work since HurToolbox assumes 3D formulation.
  If 2D NE formulation is used, HurNEInverse[] drops all zero components in the HurGlobalNEEquation in the computation.
  Remember that appropriate HurDefineVariableList[] should be defined. That is, include forces and moments used in 2D.

  2.0.4
  Added HurKinematics[] in HurDHInertia[] function. 

  2.0.3
  Changed the order of arguments in HurDHJacobian as follows to be consistent with HurDHJacobian[rf6]
  HurDHJacobian[{x,y,z},rf6] -> HurDHJacobian[rf6,{x,y,z}]
  Updated HurDHTable[] to correct definition of rf number with dh[[i,1]]
  Added HurDHMMatrix[]. This is exactly the same as HurGetMMatrix1[]

  2.0.2
  Revised HurDHJacobian[data__] to handle the following cases.
  HurDHJacobian[rf6]: computes the Jacobian matrix for rf6 at its origin with respect to rf0.
  HurDHJacobian[{x,y,z},rf6]: computes the Jacobian matrix for rf6 at {x,y,z} point from the origin of rf6 with respect to rf0.

  2.0.1
  Added HurDHJacobian[rf_] that computes the Jacobian of the RF at its origin.
  Added HurDHRetrieveGeneralizedCoordinate[dh_] to retrieve the generalized coordinates from DH table.

  2.0.0
  Unified the expression for reference frames and triads.
  old: n, a, b, c, d, ... (this method has conflicts with other naming, e.g., g, n, m, ...)
  new: rf0, rf1, rf2, rf3, rf4, ... for reference frames
       i0, j0, k0, i1, j1, k1, i2, j2, k2, ... for triads



  1.2.1
  Corrected simple version conflict

  1.2.0
  There is a big jump here.
  The naming rule for reference frames is changed.
  Previously, "n", "a", "b", "c", "d", "e" are used.
  Also, 1, 2, 3 are used to indicate x, y, z directions.
  As the number of reference frames increases, it is not a good idea to use 
  "n", "a", "b", "c", "d", "e" as the names for reference frames.
  Usually, it will be better to use something like "rf1","rf2","rf3"...
  In this way, we can avoid conflict of using "g" as a reference frame and
  "g" as a gravity acceleration constant.
  Another issue is 1, 2, 3 for the indication of x, y, z.
  If "rf1", "rf2", "rf3" are used, then we cannot simply add 1, 2, 3 to indicate x, y, z.
  What is "rf11"? Is it x axis of "rf1" or the 11th reference frame?
  So, I decided to use x, y, z instead.

  The Newtonian reference frame is still following the same convention as "n".
  However, if we prefer to use different naming, for example, "rf0", then we can change it by
  using a function HurChangeNewtonianRF[]

  1.0.2
  Updated DH part to get the honogeneous transform.

  1.0.1
  Revised the initial mass, inertia, and HurGlobalInertiaDyadic
  That is, mass, inertia for RF n are set to zero.

  1.0.0
  Automatic Differentiation functions are added.
  HurAutoDiff[value,value_dot] is added 
  HurDefineAutoDiffFunc[fn_] is added 
  HurFunctionForm4Var
  HurGetFunctionForm
  HurGetGradAutoDiff
  HurGetGradAutoDiffValue  


  0.9.9
  HurToMatlab[] is added.
  Kane method;
  HurDot[] is modified to accept dot product between dyadic, vector, and coordindate form.
  HurDefineInertia[] defines not only inertia matrix but also inertia dyadic.

  0.9.8
  HurHamiltonEquation[] now includes all nonconservative forces.
  HurConstrainedHamiltonEquation[] is added.
  HurResetNonConservativeForces[] is added.
  HurToJulia[] is added
  HurToPython[] is added
  HurResetNonConservativeForces[] is updated.

  0.9.7
  The following functions for the Hamiltonian mechanics are added.
  HurDefineGeneralizedMomentumSymbol[]
  HurDefineGeneralizedMomentumExpression[]
  HurGetHamiltonian[]
  HurHamiltonEquation[]

  0.9.6
  The following two equations are modified to accept no-argument cases.
  In this way, user-defined Lagrangian and the corresponding EL equations can be computed.
  HurGetLagrangian[]
  HurGetELEquation[]
  HurELEquation[]
  
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

$VERSION$ = "2.0.5";
$EMAIL$ = "pilwonhur@tamu.edu";
Print["HurToolbox for modeling and analysis of multibody systems ", $VERSION$, ". \nHurToolbox mainly uses vector manipulation (vectors, dyadics).\nCoordinates and matrix representation of the dyadics are also available.\nAvailable methods: Newton-Euler Method, Euler-Lagrange Method, Hamiltonian Method, Kane's Method.\nCopyright 2019 Pilwon Hur\nDepartment of Mechanical Engineering\nTexas A&M University\nAll rights reserved.\nEmail questions, comments, or concerns to ", $EMAIL$, "."];


HurDefineRF::usage="HurDefineRF[rf__] defines a reference frame. Ex) HurDefineRF[b]: to define a reference frame B. HurDefineRF[a,b]: to define reference frames a, and b at the same time.";
HurChangeNewtonianRF::usage="HurChangeNewtonianRF[rf0_] replaces the symbol for the Newtonian reference frame from n to rf0."
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
HurGetRelativeDCM::usage="HurGetRelativeDCM[rf1_, rf2_] returns the direction cosine matrix of rf1 with respect to rf2. It is also equivalent to the rotation matrix of rf1 with respect to rf2.Ex) dcm=HurGetRelativeDCM[b,c]";
HurGetRelativeDCMDyadic::usage="HurGetRelativeDCMDyadic[rf1_, rf2_] returns the rotation dyadics for rf1 with respect to rf2. (vec in rf2)=(HurGetRelativeDCMDyadic[rf1_, rf2_]) dot (vec in rf1) ";
HurUnifyTriadPool::usage="HurUnifyTriadPool[rf1_, rf2_] returns triads for rf1 expressed in terms of rf2. In other words, it returns the x,y,z axes of rf1 expressed in terms of rf2. Note that when rf1 is defined via HurDefineRF and HurDefineDCM, rf1 is expressed in terms of Newtonian reference frame. This is used only for the internal computation purpose. It is also equivalent to HurGetRelativeDCM[rf1_, rf2_]. Ex) HurUnifyTriadPool[a, b]";
HurUnifyTriads::usage="HurUnifyTriads[v_, rf_] is used to represent the given vector with respect to rf. Ex) HurUnifyTriads[v, b]";
HurUnifyTriadsCoord::usage="HurUnifyTriadsCoord[v_, rf_] is used to display the coordinates of the given vector with respect to rf. Ex) HurUnifyTriadsCoord[v, b]";
HurCoordTriads::usage="HurCoordTriads[coord__] is used to represent the coordinate into the vector form. HurCoordTriads[coord_]: coord contains the 4 elements. The first 3 are the vector coordinates and the 4th element is RF. Ex) If v=x n1 + y n2 + z n3, then coord={x,y,z,n}. The following is also possible. HurCoordTriads[coord_,RF]: If v=x nx + y ny + z nz, then coord={x,y,z}, and RF=n";
HurCross::usage="HurCross[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a vector with respect to rf. Ex) HurCross[v1, v2, b]";
HurDot::usage="HurDot[v1_, v2_] is used to perform the dot product of the two vectors. v1 and v2 can be in different RFs. Also, HurDot does not require RF information since the output is a scalar. It computes dot product between dyadics, vectors, coordiates Ex) HurDot[v1, v2]";
HurDotSimple::usage="HurDotSimple[v1_, v2_] is the subroutine that HurDot uses. HurDotSimple[] is the leaf node. It computes dot product between vectors.";
HurFindFirstTriad::usage="HurFindFirstTriad[vec_] finds the first triad in the vec. For example, HurFindFirstTriad[r*bx+2*cx] will return the index of either bx or cx whatever that comes first in the HurGlobalListTriads.";
HurNorm::usage="HurNorm[v__] returns the norm of the vector. HurNorm[v] first convert v into n, then returns norm of v. HurNorm[v,b] converts v into b and returns norm of v. If your vector is already expressed in b, then please use HurNorm[v,b]. Vector v can have triads of mixed RFs. Also, note that both will return the same norm. You may need to run Simplify"
HurNormSqure::usage="HurNormSqure[v_] returns the norm squared of the vector. Vector v can have triads of mixed RFs."
HurCrossCoord::usage="HurCrossCoord[v1_, v2_, rf_] is used to perform the cross product of the two vectors. v1 and v2 can be in different RFs. The result of the cross product will be expressed as a coordinate with respect to rf. It is equivalent to the following: HurUnifyTriadsCoord[ HurCross[v1_, v2_, rf_], rf_] Ex) HurCrossCoord[v1, v2, b]";
HurVectorDiff::usage="HurVectorDiff[v_,rf1_,rf2_] is used to perform the vector differentiation. Vector v will be differentiated with respect to rf1. This function uses the vector differentiation formula for different RFs (i.e., rf2). Note that if a vector is differentiated w.r.t. time, then it will be the same as HurVectorDiff[v_,n,rf2_].";
HurList2Column::usage="HurList2Column[list__] converts list to a column vector so that matrix multiplication can be performed."
HurAppendRF2Coord::usage="HurAppendRF2Coord[coord_, rf_] explicitly specify the RF information to the given coordinate (i.e., 3 numbers) without RF. Regardless of the size of coord_, HurAppendRF2Coord[coord_, rf_] takes the first 3 components of coord_ and attach rf_ to it. Usually, this function is used for internal usage.";
HurGetCongruenceTransform::usage="HurGetCongruenceTransform[mat_,level_:1]";
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
HurGetLagrangian::usage="HurGetLagrangian[rf__]. If no arguments are inserted, HurGlobalRF[[2 ;; nrfs]] are assumed.";
HurSetLagrangian::usage="HurSetLagrangian[lag_,rf_] can let you specify Lagrangian by the user, not via automatic procedure.";
HurDefineVertical::usage="HurDefineVertical[v_]";
HurGetELEquation::usage="HurGetELEquation[gc__]";
HurSetELEquation::usage="HurSetELEquation[eq_,gc_]";
HurGetJacobian::usage="HurGetJacobian[vec_,rf1_,rf2_] computes Jacobian matrix. vec_ can be either a position vector of a point of interest or velocity vector. rf1_ is the RF for angular velocity of your interest. rf2_ is expression of your vec_. rf2_ is usually n. However, it can also be a or b depending on your convenience."; 
HurGetMMatrix::usage="HurGetMMatrix[] returns the inertia matrix.";
HurGetMMatrix1::usage="HurGetMMatrix1[] returns the inertia matrix.";
HurDHMMatrix::usage="HurDHMMatrix[] returns the inertia matrix based on DH table. Make sure that HurDHJacobian[], HurDHInertia[] are run beforehand. Once these two are run, then HurDHMMatrix[] is exactly the same as HurGetMMatrix1[]";
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
HurELEquation::usage="HurELEquation[] computes the EL equations via automatic procedure. It internally computes LinearMomentum, AngularMomentum, Lagrangian, and the subsequent EL Equations. From the EL Equations, M-C-G Matrices are computed. If you want to simply apply EL equations given the Lagrangian, please use HurGetELEquation[] instead.";
HurELEquation1::usage="HurELEquation1[] computes the EL equations via automatic procedure. Unlike HurELEquation[], HurELEquation1[] does not differentiate Lagrangian directly. HurELEquation1[] computes i) Inertia matrix directly from the structure of the Kinetic Energy, ii) C matrix from Christoffel symbols, and iii) G vector from (gravitational and spring) potential energies and rayleigh dissipation energy. From the M-C-G matrices, EL equations are reconstructed backward via robotic equation.";
HurGetELEquationFromLagrangian::usage="(obsolete) HurGetELEquationFromLagrangian[] computes the EL equations from the user-provided Lagrangian. Make sure that you have HurGlobalLagrangian defined. Otherwise, please run HurSetLagrangian[lag, gc] first. This function is equivalent to HurGetELEquation[]. Please use HurGetELEquation[] instead.";
HurDefineConstraints::usage="HurDefineConstraints[con__]";
HurDefineConstrainedJacobian::usage="HurDefineConstrainedJacobian[]";
HurConstrainedELEquation::usage="HurConstrainedELEquation[]";
HurConstrainedELInverse::usage="HurConstrainedELInverse[]";
HurELInverse::usage="HurELInverse[]";
HurDefineLambda::usage="HurDefineLambda[]";
HurDefineNonConservativeForces::usage="HurDefineNonConservativeForces[f__]"; 
HurResetNonConservativeForces::usage="HurResetNonConservativeForces[]"; 
HurDefineOtherPotentialE::usage="HurDefineOtherPotentialE[rf_, pe_] accepts the additional potential energy other than gravity in Lagrangian mechanics. Example includes spring. If intrinsic generalized coordinates (e.g., joint angles) are used, elastic energy will be uniquely assigned to a RF. However, if your GC's are extrinsic, then assignment of elastic energy to an RF is ambiguous. However, it doesn't matter since we simply need the total sum of the elastic energies. Therefore, simply assign the elastic energy to an any reasonable RF."; 
HurDefineRayleighDissipationE::usage="HurDefineRayleighDissipationE[rf_, de_] accepts the velocity-proportional frictional forces in Lagrangian mechanics. Example includes viscous damping friction. de is the energy due to dissipation energy. rf does not matter."; 
HurGetInertiaTensor::usage="HurGetInertiaTensor[rf_]"; 
HurProductMatVec::usage="HurProductMatVec[mat_,vec_,rf_]"; 
HurDefineGeneralizedMomentumSymbol::usage="HurDefineGeneralizedMomentumSymbol[] defines the symbols for the generalized momenta. Please make sure that generalied coordinates are defined beforehand. The numbers of generalized momenta and generalized coordinates are the same."; 
HurDefineGeneralizedMomentumExpression::usage="HurDefineGeneralizedMomentumExpression[] returns the expression of the generalized momentum as defined by the partial derivative of Lagrangian with first time derivative of the generalized coordinates."; 
HurGetHamiltonian::usage="HurGetHamiltonian[] returns the Hamiltonian in terms of both q and p."; 
HurHamiltonEquation::usage="HurHamiltonEquation[] returns the Hamilton canonical equations."; 
HurConstrainedHamiltonEquation::usage="HurConstrainedHamiltonEquation[]"; 
HurDefineGeneralizedSpeeds::usage="HurDefineGeneralizedSpeeds[gs__] defines the generalized speeds for Kane's method"; 
HurDefineDependentGeneralizedSpeeds::usage="HurDefineDependentGeneralizedSpeeds[gs__]";
HurGetSimpleNonholonomicGeneralizedSpeedsJacobian::usage="HurGetSimpleNonholonomicGeneralizedSpeedsJacobian[]"
HurDefineDyadic::usage="HurDefineDyadic[data__] defines a dyadic. Make sure that data__ includes 11 inputs: 9 values, and 2 reference frames. 9 values correspond to the outer product of the 2 reference frames.";
HurDefineUnitDyadic::usage="HurDefineUnitDyadic[rf_]";
HurTransposeDyadic::usage="HurTransposeDyadic[data__] transposes the provided dyadic.";
HurDHTable::usage="HurDHTable[dh_]";
HurDHInertia::usage="HurDHInertia[data_]";
HurGetHomogeneousTransformDH::usage="HurGetHomogeneousTransformDH[dh_]";
HurGetRelativeHomogeneousTransform::usage="HurGetRelativeHomogeneousTransform[dh_]";
HurInverseHomogeneousTransform::usage="HurInverseHomogeneousTransform[dh_]";
HurDHJacobian::usage="HurDHJacobian[data__] that computes the Jacobian. Two usages: i) HurDHJacobian[rf6] returns the Jacobian of rf6 at its origin. i) HurDHJacobian[rf6,{x,y,z}] returns the Jacobian of rf6 at the offset {x,y,z} from its origin.";
HurDHRetrieveGeneralizedCoordinate::usage="HurDHRetrieveGeneralizedCoordinate[dh_] retrieves the generalized coordinates from the DH table.";
HurDHELEquation::usage="HurDHELEquation[] computes the M,C,G. Once run, check HurGlobalMMatrix, HurGlobalCMatrix, HurGlobalGVector";
HurSimplifyVariablesTimed::usage="HurSimplifyVariablesTimed[var_,time_] simplifies all the elements of the provided var_ within the provided time_ for each element. If timed out, it will return its original expression"; 
HurFullSimplifyVariablesTimed::usage="HurFullSimplifyVariablesTimed[var_,time_] fully simplifies all the elements of the provided var_ within the provided time_ for each element. If timed out, it will return its original expression"; 
HurDumpSaveData::usage="HurDumpSaveData[filename__] Please use .mx for the extension of the filename. It will save variables in binary (unreadable) expression, is very fast to load (with large data). However, this binary data are platform-specific. If saved in Mac, it cannot be used in Windows or Linux.";
HurSaveData::usage="HurSaveData[filename__]. Please use .m for the extension of the filename. It will save variables in (readable) portable expression (or ascii format), is very slow to load (with large data). This data are platform-independent. You can use in any platforms.";
HurLoadData::usage="HurLoadData[filename_]";
HurToJulia::usage="HurToJulia[exp_]";
HurToPython::usage="HurToPython[exp_]";
HurToMatlab::usage="HurToMatlab[exp_]";
HurAutoDiff::usage="HurAutoDiff[value,value_dot]";
HurDefineAutoDiffFunc::usage="HurDefineAutoDiffFunc[fn_]";
HurFunctionForm4Var::usage="HurFunctionForm4Var[] provides the rule q[t]->q, q'[t]->qd, and q''[t]->qdd.";
HurConvertVar::usage="HurConvertVar[x__] converts {q[t], q'[t], q''[t]} to {q, qd, qdd} the the Rule from HurFunctionForm4Var[]";
HurGetFunctionForm::usage="HurGetFunctionForm[var_, fn_] defines a function to be used for any generic purposes. AutoDifferentiation will use this.";
HurGetGradAutoDiff::usage="HurGetGradAutoDiff[fn_, val__] performs the Automatic Differentiation of a function fn_ about the point val_. It will return only the gradient.";
HurGetGradAutoDiffValue::usage="HurGetGradAutoDiffValue[fn_, val__] is the same as HurGetGradAutoDiff[] except that HurGetGradAutoDiffValue returns the function value as well as the gradient.";


HurInitialize[data__:Global`n] := (
  HurGlobalCoordinateSystemsOption=1;
  datas=Flatten[List[data]];
  Switch[Length[datas],
   0, NewtonianRF=Global`n;,
   1, NewtonianRF=datas[[1]];,
   2, NewtonianRF=datas[[1]];If[datas[[2]]==="xyz",HurGlobalCoordinateSystemsOption=2,Null];
   ];

  HurGlobalRF = {NewtonianRF};
  HurGlobalDCM = List[RotationMatrix[0, {0, 0, 1}]]; HurGlobalMass={0}; HurGlobalInertia={{0,0,0,0,0,0}};
  HurGlobalInertiaDyadic={{0,0,0,0,0,0,0,0,0,NewtonianRF,NewtonianRF}};
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

  tempAxisString={"x", "y", "z"};
  If[HurGlobalCoordinateSystemsOption===1,
    HurGlobalListTriads = {{i0,j0,k0}};HurGlobalTriadsConversion = {{i0->i0,j0->j0,k0->k0}};
    ,
    HurGlobalListTriads = List[Table[Symbol[ToString[NewtonianRF]<>tempAxisString[[j]]],{j,3}]];
    HurGlobalTriadsConversion = List[Table[Symbol[ToString[NewtonianRF]<>tempAxisString[[j]]] \[Rule] Symbol[ToString[NewtonianRF]<>tempAxisString[[j]]],{j,3}]];
    ];
  
  (*HurGlobalTriadsConversion = {0};*)
  HurGlobalSimplify = True; HurGlobalGeneralizedSpeedsExpression = {};
  HurGlobalGeneralizedSpeedsSymbol = {}; HurGlobalKaneEquation = {0};
  HurGlobalTemp = {0}; HurGlobalDependentGeneralizedSpeeds = {0};
  HurGlobalGC2GS = {0}; HurGlobalConstraintsGS = {0};
  HurGlobalGSConstrainedJacobian = {0};
  HurGlobalAngularVelAbs = {0}; HurGlobalAngularVelRel = {{0,0,0}};
  HurGlobalDHTable = {0}; HurGlobalDHInertia = {0}; HurGlobalDHOrigin = {0};
  HurGlobalGeneralizedMomentumSymbol = {}; HurGlobalGeneralizedMomentumExpression = {};
  HurGlobalHamiltonian = {0}; HurGlobalHamiltonEquation = {{0},{0}};
  HurGlobalHomogeneousTransform = List[IdentityMatrix[4]];
)

Begin["`Private`"];

HurDefineRF[rf__] := (rfs=List[rf];narg=Length[rfs];
  tempAxisString={"x", "y", "z"};
  tempAxisijkString={"i", "j", "k"};
  Do[
      If[
        HurGetIndexGlobalRF[ rfs[[i]] ]===0
        ,
        AppendTo[HurGlobalRF, rfs[[i]] ]; 

        If[HurGlobalCoordinateSystemsOption===1,
          AppendTo[HurGlobalListTriads, Table[ Symbol[ tempAxisijkString[[j]] <> ToString[ HurGetIndexGlobalRF[ rfs[[i]] ]-1 ] ], {j,3}] ];
          ,
          AppendTo[HurGlobalListTriads, Table[ Symbol[ ToString[ rfs[[i]] ] <> tempAxisString[[j]] ], {j,3}] ];
          ];
        

        AppendTo[HurGlobalDCM, RotationMatrix[0, {0, 0, 1}]];
        HurDefineDCM[ rfs[[i]] ,RotationMatrix[0, {0, 0, 1}] ];
        AppendTo[HurGlobalHomogeneousTransform, IdentityMatrix[4]];
        AppendTo[HurGlobalMass, 0];
        AppendTo[HurGlobalInertia, {0,0,0,0,0,0}];
        AppendTo[HurGlobalInertiaDyadic, {0,0,0,0,0,0,0,0,0,rfs[[i]],rfs[[i]]}];
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

HurChangeNewtonianRF[rf0_] := (
  HurGlobalRF[[1]]=rf0;
  HurGlobalInertiaDyadic[[1]]={0,0,0,0,0,0,0,0,0,rf0,rf0};
  tempAxisString={"x", "y", "z"};
  HurGlobalListTriads[[1]] = Table[ Symbol[ ToString[rf0] <> tempAxisString[[j]] ], {j,3}]; 
  HurGlobalTriadsConversion[[1]] = Table[
    Symbol[ ToString[rf0] <> tempAxisString[[j]] ] \[Rule] Symbol[ ToString[rf0] <> tempAxisString[[j]] ]
    ,{j,3}
  ];
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

  tempAxisString={"x", "y", "z"};
  tempAxisijkString={"i", "j", "k"};
  Do[
    HurGlobalTriadsConversion[[k]] = Table[
      Rot=HurUnifyTriadPool[  HurGlobalRF[[ i ]] , HurGlobalRF[[ k ]] ];
      If[HurGlobalCoordinateSystemsOption===1,
        Symbol[ tempAxisijkString[[j]] <> ToString[ i-1 ] ] \[Rule]  
        Dot[ Rot[[;;,j]],HurGlobalListTriads[[k]] ]
        ,
        Symbol[ ToString[ HurGlobalRF[[i]] ] <> tempAxisString[[j]] ] \[Rule]  
        Dot[ Rot[[;;,j]],HurGlobalListTriads[[k]] ]
        ]
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
	Rot3=If[HurGlobalSimplify, Simplify[Rot3], Rot3];
  Rot3
	)

HurGetRelativeDCMDyadic[rf1_, rf2_] := (
  rot=HurGetRelativeDCM[rf1, rf2];
  HurDefineDyadic[rot, rf2, rf1]
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

HurFindFirstTriad[vec_] := (
  First@SparseArray[Coefficient[vec, HurGlobalListTriads]]["NonzeroPositions"]
  )

HurDotSimple[v1_, v2_] := (
  If[v1===0 || v2===0
    ,
    tempAns=0;
    ,
    If[MemberQ[Flatten[HurGlobalListTriads],v1]  (* simple vector. e.g. b1 *)
      ,
      tempPos=Flatten[Position[HurGlobalListTriads, v1]];
      tempAns=HurUnifyTriadsCoord[v2, HurGlobalRF[[ tempPos[[1]] ]]][[ tempPos[[2]] ]];
      ,
      If[MemberQ[Flatten[HurGlobalListTriads],v2]  (* simple vector. e.g. b1 *)
        ,
        tempPos=Flatten[Position[HurGlobalListTriads, v2]];
        tempAns=HurUnifyTriadsCoord[v1, HurGlobalRF[[ tempPos[[1]] ]]][[ tempPos[[2]] ]];
        ,
        tempPos=HurFindFirstTriad[v1];
        tempD = Dot[HurUnifyTriadsCoord[v1,HurGlobalRF[[ tempPos[[1]] ]]][[1 ;; 3]], 
        HurUnifyTriadsCoord[v2, HurGlobalRF[[ tempPos[[1]] ]]][[1 ;; 3]]]; 
        tempAns = If[HurGlobalSimplify, Simplify[tempD], tempD];
        ];
      ];
    ];
  tempAns
  )

HurDot[v1_, v2_] := (
  Switch[Length[Flatten[List[v1]]]*Length[Flatten[List[v2]]],
   1, tempAns=HurDotSimple[v1,v2];,
   4, If[Length[Flatten[List[v1]]] === 1, 
     tempAns = HurDot[v1, HurCoordTriads[v2]], 
     tempAns = HurDot[HurCoordTriads[v1], v2]];,
   11, If[Length[Flatten[List[v1]]] === 11,
     data = ArrayReshape[v1[[1 ;; 9]], {3, 3}];
     tempAns = Total[Flatten[Table[ If[data[[j,i]]===0,0,data[[j,i]]*
          HurDot[HurGlobalListTriads[[HurGetIndexGlobalRF[v1[[11]]],i]],v2]*
          HurGlobalListTriads[[HurGetIndexGlobalRF[v1[[10]]],j]] ],{j,3}, {i, 3}]]];,
     data = ArrayReshape[v2[[1 ;; 9]], {3, 3}];
     tempAns = Total[Flatten[Table[ If[data[[j,i]]===0,0,data[[j,i]]*
          HurDot[v1,HurGlobalListTriads[[HurGetIndexGlobalRF[v2[[10]]],j]]]*
          HurGlobalListTriads[[HurGetIndexGlobalRF[v2[[11]]],i]] ],{j,3},{i,3}]]];];,
   16, tempAns = HurDot[HurCoordTriads[v1],HurCoordTriads[v2]];,
   44, If[Length[Flatten[List[v1]]] === 11, 
     data = ArrayReshape[v1[[1 ;; 9]], {3, 3}];
     tempAns = Total[Flatten[Table[ If[data[[j,i]]===0,0,data[[j,i]]*
          HurDot[HurGlobalListTriads[[HurGetIndexGlobalRF[v1[[11]]],i]],HurCoordTriads[v2]]*
          HurGlobalListTriads[[HurGetIndexGlobalRF[v1[[10]]],j]] ],{j,3},{i, 3}]]];, 
     data = ArrayReshape[v2[[1 ;; 9]], {3, 3}];
     tempAns = Total[Flatten[Table[ If[data[[j,i]]===0,0,data[[j,i]]*
          HurDot[HurCoordTriads[v1],HurGlobalListTriads[[HurGetIndexGlobalRF[v2[[10]]],j]]]*
          HurGlobalListTriads[[HurGetIndexGlobalRF[v2[[11]]],i]] ],{j,3},{i,3}]]];];,
   121, tempAns=Flatten[List[Table[Total[Flatten[Table[
    v1[[i + 3*(k - 1)]]*v2[[3*j - 3 + l]]*
     HurDot[HurGlobalListTriads[[HurGetIndexGlobalRF[v1[[11]]], i]], 
      HurGlobalListTriads[[HurGetIndexGlobalRF[v2[[10]]], j]]], {i, 
     3}, {j, 3}]]], {k, 3}, {l, 3}],v1[[10]],v2[[11]]]]];
  tempAns = If[HurGlobalSimplify, Simplify[tempAns], tempAns];
  tempAns
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

HurList2Column[list__] := (
  ls=Flatten[List[list]];
  nls=Length[ls];
  ArrayReshape[ls, {nls, 1}]
  )


HurGetCongruenceTransform[mat_] := (
  dim=Dimensions[mat];
  n=dim[[1]];
  level=1;
  Tcumul=IdentityMatrix[n];
  mat1=mat;
  If[dim[[1]]!=dim[[2]],
    Print["The matrix must be square matrix!"];
    Tcumul=Null;
    ,
    For[level = 1, level < n, level++,
      T=IdentityMatrix[n];
      Do[
        T[[level,i+1]]=-mat1[[level,i+1]]/mat1[[level,level]];
        ,
        {i,level,n-1}
        ];
        mat1=Transpose[T].mat1.T;
        Tcumul=Tcumul.T;
      ];
    ];
  List[Tcumul,mat1]
  )




(* Newton-Euler Mechanics *)

HurProductMatVec[mat_,vec_,rf_] := (ww=HurUnifyTriadsCoord[vec,rf];
  HurUnifyTriads[HurCoordTriads[Flatten[{mat.{ww[[1]],ww[[2]],ww[[3]]},ww[[4]]}]],rf]
  )



HurGetInertiaTensor[rf_] := (IItemp=HurGlobalInertia[[ HurGetIndexGlobalRF[rf] ]];
  {{IItemp[[1]],IItemp[[2]],IItemp[[3]]},{IItemp[[2]],IItemp[[4]],IItemp[[5]]},{IItemp[[3]],IItemp[[5]],IItemp[[6]]}}
  )

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

HurDefineInertia[rf_, II_] := (HurGlobalInertia[[HurGetIndexGlobalRF[ rf ] ]] = II;
  HurGlobalInertiaDyadic[[HurGetIndexGlobalRF[ rf ] ]]={II[[1]],II[[2]],II[[3]],II[[2]],II[[4]],II[[5]],II[[3]],II[[5]],II[[6]],rf,rf };
  )

HurDefineForces[rf_, force_, r_] := (AppendTo[HurGlobalForce,{rf,force,r}];)

HurResetForces[] := (HurGlobalForce={};)

HurDefineMoments[rf_, moment_] := (AppendTo[HurGlobalMoment,{rf,moment}];)

HurResetMoments[] := (HurGlobalMoment={};)

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
  TempNEEquation=DeleteCases[TempNEEquation,0]
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

HurResetNonConservativeForces[] := (
  ngcs=Length[HurGlobalNonConservativeForces];
  HurGlobalNonConservativeForces=Table[0, {i, ngcs}];
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

HurGetLagrangian[rf__] := (
  nrfs=HurGetNumGlobalRF[];
  Do[
    If[
      HurGlobalMass[[i]] === 0 
      ,
      Null
      ,
      HurGetLinearMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
      HurGetAngularMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];        
      ];
    ,
    {i,2,nrfs}
    ];
  HurGetKineticE[rf];HurGetPotentialE[rf];HurGlobalLagrangian=HurGlobalKineticE-HurGlobalPotentialE-HurGlobalOtherPotentialE;
  Total[HurGlobalLagrangian]
  )
(* when no arguments are inserted *)
HurGetLagrangian[] := (
  nrfs=HurGetNumGlobalRF[];
  HurGetLagrangian[ HurGlobalRF[[2 ;; nrfs]] ]  
  (* 
  rf=HurGlobalRF[[2 ;; nrfs]];
  HurGetKineticE[rf];HurGetPotentialE[rf];HurGlobalLagrangian=HurGlobalKineticE-HurGlobalPotentialE-HurGlobalOtherPotentialE;
  Total[HurGlobalLagrangian]
  *)
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
  temp=List[
    Grad[ vel[[1]],GCdots ],
    Grad[ vel[[2]],GCdots ],
    Grad[ vel[[3]],GCdots ],
    Grad[ omega[[1]],GCdots ],
    Grad[ omega[[2]],GCdots ],
    Grad[ omega[[3]],GCdots ]
  ];
  If[HurGlobalSimplify, Simplify[temp], temp]
  )

HurSetELEquation[eq_,gc_] := (
    HurGlobalELEquation[[ Position[ HurGlobalGeneralizedCoordinates,gc ] [[1]][[1]] ]] = eq;
    HurGlobalELEquation
  )

(* obsolete *)
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

HurGetELEquation[gc__] := (
  gcs=Flatten[ List[gc] ];ngcs=Length[gcs];
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

(* use this when Lagrangian is modified by users *)
HurGetELEquation[] := (
  HurGetELEquation[HurGlobalGeneralizedCoordinates]
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
    (* the following Do[] part is inserted into HurGetLagrangian[]
    Do[
      If[
        HurGlobalMass[[i]] === 0 
        ,
        Null
        ,
        HurGetLinearMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
        HurGetAngularMomentum[HurGlobalRF[[i]],HurGlobalRF[[i]]];
        ];
      ,
      {i,2,nrfs}
      ];
    *)

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
  nrfs=HurGetNumGlobalRF[];
  tempM=
  	Table[ 
      Jac=HurGetJacobian[ HurGlobalCOMPos[[i]], HurGlobalRF[[i]], HurGlobalRF[[1]] ];
 	    Jacv=Jac[[1;;3,;;]];
 	    Jacw=Jac[[4;;6,;;]];
 	    temp1=HurGlobalMass[[i]]*Transpose[Jacv].Jacv;
 	    rot=HurGetRelativeDCM[ HurGlobalRF[[i]],HurGlobalRF[[1]] ];
 	    II=HurGetInertiaTensor[ HurGlobalRF[[i]] ];
 	    temp2=Transpose[Jacw].rot.II.Transpose[rot].Jacw;
      temp=temp1+temp2;
      If[HurGlobalSimplify, Simplify[temp], temp]
      ,{i,2,nrfs}
      ];
  HurGlobalMMatrix=Total[tempM];
  HurGlobalMMatrix
  )

(* HurGetMMatrix1[] exactly the same as HurDHMMatrix[] *)
HurDHMMatrix[] := (
  nrfs=HurGetNumGlobalRF[];
  tempM=
    Table[ 
      Jac=HurDHJacobian[HurGlobalRF[[i]], HurGlobalDHInertia[[i-1, 3]] ];
      (* Jac=HurGetJacobian[ HurGlobalCOMPos[[i]], HurGlobalRF[[i]], HurGlobalRF[[1]] ]; *) (* This is the same *)
      Jacv=Jac[[1;;3,;;]];
      Jacw=Jac[[4;;6,;;]];
      temp1=HurGlobalMass[[i]]*Transpose[Jacv].Jacv;
      rot=HurGetRelativeDCM[ HurGlobalRF[[i]],HurGlobalRF[[1]] ];
      II=HurGetInertiaTensor[ HurGlobalRF[[i]] ];
      temp2=Transpose[Jacw].rot.II.Transpose[rot].Jacw;
      temp=temp1+temp2;
      If[HurGlobalSimplify, Simplify[temp], temp]
      ,{i,2,nrfs}
      ];
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
  HurGlobalCMatrix=If[HurGlobalSimplify, Simplify[HurGlobalCMatrix], HurGlobalCMatrix]
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

HurResetConstraints[] := (HurGlobalConstraints={};)


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

(* Hamiltonian mechanics *)
HurDefineGeneralizedMomentumSymbol[] :=(
  ngcs=Length[HurGlobalGeneralizedCoordinates];
  HurGlobalGeneralizedMomentumSymbol=Table[ ToExpression["p" <> ToString[i] <> "[t]" ] ,{i,ngcs} ];
  )

HurDefineGeneralizedMomentumExpression[] :=(
  ngcs=Length[HurGlobalGeneralizedCoordinates];
  HurGetLagrangian[];
  L = Total[HurGlobalLagrangian];
  temp=Table[
    qdot=D[HurGlobalGeneralizedCoordinates[[i]],Global`t];
    D[L,qdot]
    ,
    {i,ngcs}
    ];
  HurGlobalGeneralizedMomentumExpression=If[HurGlobalSimplify, Simplify[temp], temp]
  )

HurGetHamiltonian[] :=(
  HurGetMMatrix1[];
  qdot=LinearSolve[HurGlobalMMatrix,HurGlobalGeneralizedMomentumSymbol];
  temp = 1/2 HurGlobalGeneralizedMomentumSymbol.qdot + Total[HurGlobalPotentialE];
  HurGlobalHamiltonian=If[HurGlobalSimplify, Simplify[temp], temp]
  )

HurHamiltonEquation[] :=(
  HurDefineGeneralizedMomentumSymbol[];
  HurDefineGeneralizedMomentumExpression[];
  HurGetHamiltonian[];

  temp1=Grad[HurGlobalHamiltonian, HurGlobalGeneralizedMomentumSymbol];
  HurGlobalHamiltonEquation[[1]]=If[HurGlobalSimplify, Simplify[temp1], temp1];
  temp2=-Grad[HurGlobalHamiltonian, HurGlobalGeneralizedCoordinates]+HurGlobalNonConservativeForces;
  HurGlobalHamiltonEquation[[2]]=If[HurGlobalSimplify, Simplify[temp2], temp2];
  HurGlobalHamiltonEquation
  )

HurConstrainedHamiltonEquation[] :=(
  HurHamiltonEquation[];
  HurDefineConstrainedJacobian[];
  HurDefineLambda[];
  m=Length[HurGlobalLambda];
  HurGlobalGeneralizedConstrainingForce=Transpose[HurGlobalConstrainedJacobian].HurGlobalLambda;
  HurGlobalHamiltonEquation[[2]]=HurGlobalHamiltonEquation[[2]]+HurGlobalGeneralizedConstrainingForce;
  HurGlobalHamiltonEquation
  )
  
(* Kane Method *)
HurDefineGeneralizedSpeeds[gs__] := (gss=Flatten[ List[ gs ] ]; ngss=Length[gss];
  HurGlobalGeneralizedSpeedsExpression=gss;
  HurGlobalKaneEquation=Table[0,{i,ngss}];

  HurGlobalGeneralizedSpeedsSymbol=Table[ ToExpression["u" <> ToString[i] <> "[t]" ] ,{i,ngss} ];
  (* Global`t *)

  ngcs=Length[HurGlobalGeneralizedCoordinates];
  (* ngcs-ngss is the number of nonholonomic constraints *)

  // HurDefineConstrainedJacobian[];
  eq=Table[HurGlobalGeneralizedSpeedsExpression[[i]]==HurGlobalGeneralizedSpeedsSymbol[[i]],{i,ngcs}];
  GCDot=D[HurGlobalGeneralizedCoordinates,Global`t];
  HurGlobalGC2GS=Flatten[Solve[eq,GCDot]];
  )

HurDefineDependentGeneralizedSpeeds[gs__] := (
  HurGlobalDependentGeneralizedSpeeds=Flatten[ List[ gs ] ];
  )

HurGetSimpleNonholonomicGeneralizedSpeedsJacobian[] := (
  ncon=Length[HurGlobalConstraints];
  ndgs=Length[HurGlobalDependentGeneralizedSpeeds];
  GSDepIndex=Flatten[Table[Position[HurGlobalGeneralizedSpeedsSymbol, HurGlobalDependentGeneralizedSpeeds[[i]]],{i,ndgs}],1];
  GSIndep=Delete[HurGlobalGeneralizedSpeedsSymbol,GSDepIndex];

  HurGlobalConstraintsGS=Flatten[HurGlobalConstraints/.HurGlobalGC2GS];
  eq=Table[HurGlobalConstraintsGS[[i]] == 0, {i, ncon}];
  ans=Flatten[Solve[eq,HurGlobalDependentGeneralizedSpeeds]];

  HurGlobalGSConstrainedJacobian=Table[
    temp=HurGlobalDependentGeneralizedSpeeds[[i]]/.ans;
    Grad[temp,GSIndep]
    ,
    {i,ndgs}
    ];
  HurGlobalGSConstrainedJacobian
  )

HurDefineDyadic[data__] := (
  datavec=Flatten[ List[ data ] ];
  ndata=Length[datavec];
  If[ndata===11,
    temp=datavec;
    ,
    Print["You entered wrong data. Returning zero dyadic."]
    temp={0,0,0,0,0,0,0,0,0,n,n};
    ];
  temp
  )

HurDefineUnitDyadic[rf_] := (
  HurDefineDyadic[1,0,0,0,1,0,0,0,1,rf,rf]
  )

HurTransposeDyadic[data__] := (
  datavec=Flatten[ List[ data ] ];
  ndata=Length[datavec];
  If[ndata===11,
    elem=Flatten[Transpose[ArrayReshape[datavec[[1;;9]],{3,3}]]];
    temp=Flatten[List[elem,datavec[[11]],datavec[[10]]]];
    ,
    Print["The argument you entered is not in the right format for dyadics. No changes happened."]
    temp=datavec;
    ];
  temp
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
      HurDefineRF[ Symbol[ "rf" <> ToString[ dh[[i,1]] ] ]];
      , 
      {i,ndh}
    ];
    (* define DCM *)
    Do[
      HurDefineDCMRelative[ HurGlobalRF[[ dh[[i,1]]+1 ]] , HurGlobalRF[[ dh[[i,1]] ]],  RotationMatrix[dh[[i,3]],{0,0,1}].RotationMatrix[dh[[i,6]],{1,0,0}] ];
      HurGetHomogeneousTransformDH[dh[[i]]];
      , 
      {i,ndh}
    ];
    HurDefineGeneralizedCoordinates[Table[HurDHRetrieveGeneralizedCoordinate[dh[[i]]], {i,Length[dh]}]];
  )

HurGetHomogeneousTransformDH[dh_] := (
  R1R2=HurGetRelativeDCM[ HurGlobalRF[[ dh[[1]]+1 ]], HurGlobalRF[[ dh[[1]] ]] ];
  R1=RotationMatrix[dh[[3]],{0,0,1}];
  TransVec=R1.{{dh[[5]]},{0},{dh[[4]]}};
  HurGlobalHomogeneousTransform[[ dh[[1]]+1 ]]=ArrayFlatten[{{R1R2, TransVec}, {{{0, 0, 0}}, 1}}];
  )

HurGetRelativeHomogeneousTransform[rf1_, rf2_] := (
  TempMat=IdentityMatrix[4];
  If[HurGetIndexGlobalRF[rf1]>=HurGetIndexGlobalRF[rf2]+1
    ,
    Do[
      TempMat=TempMat.HurGlobalHomogeneousTransform[[i]];
      , {i,HurGetIndexGlobalRF[rf2]+1,HurGetIndexGlobalRF[rf1]}
      ];
    ,
    Do[
      TempMat=TempMat.HurInverseHomogeneousTransform[HurGlobalHomogeneousTransform[[i]]];
      , {i,HurGetIndexGlobalRF[rf2],HurGetIndexGlobalRF[rf1]+1,-1}
      ];
    ];
  If[HurGlobalSimplify, Simplify[TempMat], TempMat]
  )

HurInverseHomogeneousTransform[mat_] := (
  R=mat[[1;;3,1;;3]];
  d=HurList2Column[mat[[1;;3,4]]];
  ArrayFlatten[{{Transpose[R], -Transpose[R].d}, {{{0, 0, 0}}, 1}}]
  )

HurDHJacobian[data__] := (
  datavec=List[ data ];
  ndata=Length[datavec];
  If[ndata===1,
    rf=datavec[[1]];
    T = HurGetRelativeHomogeneousTransform[rf, HurGlobalRF[[1]]];
    d=T[[1 ;; 3, 4]];
    vec=HurCoordTriads[Flatten[d], HurGlobalRF[[1]]];
    Jac=HurGetJacobian[vec,rf,HurGlobalRF[[1]]];
    ,
    offset=HurList2Column[ datavec[[2]] ];
    rf=datavec[[1]];
    T = HurGetRelativeHomogeneousTransform[rf, HurGlobalRF[[1]]];
    Rot=T[[1 ;; 3, 1 ;; 3]];
    d=T[[1 ;; 3, 4]]+Rot.offset;
    vec=HurCoordTriads[Flatten[d], HurGlobalRF[[1]]];
    Jac=HurGetJacobian[vec,rf,HurGlobalRF[[1]]];
  ];
  Jac
  )

HurDHRetrieveGeneralizedCoordinate[dh_] :=(
  var=Total[D[dh, Global`t ]];
  ToExpression[StringReplace[ToString[FullForm[var]], "[1]" -> "[0]"]]
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
  HurKinematics[];
  )

HurDHELEquation[] := (
    HurGetLagrangian[];
    HurDHMMatrix[];
    HurGetCMatrix[];
    HurGetGVector[];
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
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalCoordinateSystemsOption","HurGlobalMass","HurGlobalInertia","HurGlobalInertiaDyadic","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalRayleighDissipationE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify","HurGlobalGeneralizedSpeedsExpression","HurGlobalGeneralizedSpeedsSymbol","HurGlobalKaneEquation","HurGlobalTemp","HurGlobalDependentGeneralizedSpeeds","HurGlobalGC2GS","HurGlobalConstraintsGS","HurGlobalGSConstrainedJacobian","HurGlobalAngularVelAbs","HurGlobalAngularVelRel","HurGlobalDHTable","HurGlobalDHInertia","HurGlobalDHOrigin","HurGlobalGeneralizedMomentumSymbol","HurGlobalGeneralizedMomentumExpression","HurGlobalHamiltonian","HurGlobalHamiltonEquation","HurGlobalHomogeneousTransform"};
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
  tempvar1={"HurGlobalRF","HurGlobalDCM","HurGlobalCoordinateSystemsOption","HurGlobalMass","HurGlobalInertia","HurGlobalInertiaDyadic","HurGlobalForce","HurGlobalMoment","HurGlobalCOMPos","HurGlobalCOMVel","HurGlobalCOMAcc","HurGlobalAngularVel","HurGlobalAngularAcc","HurGlobalLinearMomentum","HurGlobalAngularMomentum","HurGlobalVertical","HurGlobalNEEquation","HurGlobalVariableList","HurGlobalKineticE","HurGlobalPotentialE","HurGlobalLagrangian","HurGlobalELEquation","HurGlobalGeneralizedCoordinates","HurGlobalMMatrix","HurGlobalCMatrix","HurGlobalGVector","HurGlobalConstrainedJacobian","HurGlobalConstraints","HurGlobalLambda","HurGlobalGeneralizedConstrainingForce","HurGlobalConstrainedELEquation","HurGlobalConstrainedModified","HurGlobalNonConservativeForces","HurGlobalOtherPotentialE","HurGlobalRayleighDissipationE","HurGlobalListTriads","HurGlobalTriadsConversion","HurGlobalSimplify","HurGlobalGeneralizedSpeedsExpression","HurGlobalGeneralizedSpeedsSymbol","HurGlobalKaneEquation","HurGlobalTemp","HurGlobalDependentGeneralizedSpeeds","HurGlobalGC2GS","HurGlobalConstraintsGS","HurGlobalGSConstrainedJacobian","HurGlobalAngularVelAbs","HurGlobalAngularVelRel","HurGlobalDHTable","HurGlobalDHInertia","HurGlobalDHOrigin","HurGlobalGeneralizedMomentumSymbol","HurGlobalGeneralizedMomentumExpression","HurGlobalHamiltonian","HurGlobalHamiltonEquation","HurGlobalHomogeneousTransform"};
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





JuliaForm[Rational[a_, b_]] := JuliaForm[a] <> "/" <> JuliaForm[b];
JuliaForm[Complex[a_, b_]] := "(" <> JuliaForm[a] <> "+" <> JuliaForm[b] <> "im" <> ")";
JuliaForm[Times[a_, b_]] := JuliaForm[a] <> "*" <> JuliaForm[b];
JuliaForm[Plus[a_, b_]] := "(" <> JuliaForm[a] <> "+" <> JuliaForm[b] <> ")";
(* Most of the commands starting with Capital letter will be lowered. 
  Ex: Sin[xxx] -> sin(xxx). Note that q1[t] -> q1 for simulation purpose. 
  In HurToolbox, only generalized coordinates and generalized momenta are 
  functions of time. For simulation purpose, the [t] from q1[t], p1[t] 
  will be removed. This treatment is only for the combination of HurToolbox 
  and other numerical packages (e.g., Julia). 
*)
JuliaForm[h_[args__]] := If[JuliaForm[args] === "t", JuliaForm[h], 
    ToLowerCase[JuliaForm[h]] <> "(" <> JuliaForm[args] <> ")"];
(* treat Exp -> exp and Sqrt -> sqrt. All other Power will be treated as ()^() *)
JuliaForm[Power[a_, b_]] := If[JuliaForm[b] === "1/2", 
    "sqrt" <> "(" <> JuliaForm[a] <> ")", 
    If[JuliaForm[a] === "exp(1)", 
      "(" <> "exp" <> "(" <> JuliaForm[b] <> ")", 
      "(" <> JuliaForm[a] <> ")" <> "^" <> "(" <> JuliaForm[b] <> ")"]];
JuliaForm[a_ListQ] := StringReplace[ToString[a], 
    {"}," -> "];", "{" -> "[", "}}" -> "]]", "," -> " "}];
JuliaForm[Derivative[a_][b_]] := Switch[a, 1, JuliaForm[b] <> "d", 
    2, JuliaForm[b] <> "dd", 3, JuliaForm[b] <> "ddd"];
(*
  JuliaForm[Sin]:="sin";
  JuliaForm[Cos]:="cos";
  JuliaForm[Tan]:="tan";
  JuliaForm[Csc]:="csc";
  JuliaForm[Sec]:="sec";
  JuliaForm[Cot]:="cot";
  JuliaForm[Sinh]:="sinh";
  JuliaForm[Cosh]:="cosh";
  JuliaForm[Tanh]:="tanh";
  JuliaForm[Csch]:="csch";
  JuliaForm[Sech]:="sech";
  JuliaForm[Coth]:="coth";
  JuliaForm[Exp]:="exp";
*)
(* Handling arrays *)
(*
JuliaForm[List[args__]] := "(" <> "[" <> Table[JuliaForm[{args}[[ii]]] <> " ", 
  {ii, 1, Length@{args}}] <> "]" <> ")";
*)
JuliaForm[List[args__]] := (Switch[Length[Dimensions[List[args]]], 
   1, 
   temp = "[" <> Table[JuliaForm[List[args][[i]]] <> " ", {i,Length[List[args]]}] <> "]",
   2, 
   temp = "[" <> Table[JuliaForm[List[args][[i]]] <> ";", {i, Length[List[args]]}] <> "]"
   ];
  temp = StringReplace[temp, {"[[" -> "[", "];[" -> ";", "];]" -> "]"}];
  StringReplace[temp, {" ]" -> "]"}])

(*Pi and E*)
JuliaForm[\[Pi]] = "pi";
JuliaForm[Pi] = "pi";
JuliaForm[E] = "exp(1)";
(*real numbers,engineering notation*)
JuliaForm[r_Real] := Block[{a = MantissaExponent[r]}, 
    If[r >= 0, ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]], 
     "(" <> ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]] <> ")"]];
(*Greek characters*)
greekrule = {"\[Alpha]" -> "alpha", "\[Beta]" -> "beta", 
    "\[Gamma]" -> "gamma", "\[Delta]" -> "delta", 
    "\[CurlyEpsilon]" -> "curlyepsilon", "\[Zeta]" -> "zeta", 
    "\[Eta]" -> "eta", "\[Theta]" -> "theta", "\[Iota]" -> "iota", 
    "\[Kappa]" -> "kappa", "\[Lambda]" -> "lambda", "\[Mu]" -> "mu", 
    "\[Nu]" -> "nu", "\[Xi]" -> "xi", "\[Omicron]" -> "omicron", 
    "\[Pi]" -> "pi", "\[Rho]" -> "rho", 
    "\[FinalSigma]" -> "finalsigma", "\[Sigma]" -> "sigma", 
    "\[Tau]" -> "tau", "\[Upsilon]" -> "upsilon", 
    "\[CurlyPhi]" -> "curlyphi", "\[Chi]" -> "chi", "\[Psi]" -> "psi",
     "\[Omega]" -> "omega", "\[CapitalAlpha]" -> "Alpha", 
    "\[CapitalBeta]" -> "Beta", "\[CapitalGamma]" -> "Gamma", 
    "\[CapitalDelta]" -> "Delta", 
    "\[CapitalEpsilon]" -> "CurlyEpsilon", "\[CapitalZeta]" -> "Zeta",
     "\[CapitalEta]" -> "Eta", "\[CapitalTheta]" -> "Theta", 
    "\[CapitalIota]" -> "Iota", "\[CapitalKappa]" -> "Kappa", 
    "\[CapitalLambda]" -> "Lambda", "\[CapitalMu]" -> "Mu", 
    "\[CapitalNu]" -> "Nu", "\[CapitalXi]" -> "Xi", 
    "\[CapitalOmicron]" -> "Omicron", "\[CapitalPi]" -> "Pi", 
    "\[CapitalRho]" -> "Rho", "\[CapitalSigma]" -> "Sigma", 
    "\[CapitalTau]" -> "Tau", "\[CapitalUpsilon]" -> "Upsilon", 
    "\[CapitalPhi]" -> "CurlyPhi", "\[CapitalChi]" -> "Chi", 
    "\[CapitalPsi]" -> "Psi", "\[CapitalOmega]" -> "Omega"};
(*Everything else*)
JuliaForm[allOther_] := 
   StringReplace[ToString[allOther, FortranForm], greekrule];

HurToJulia[exp_] := (
  CopyToClipboard[JuliaForm[exp]];
  JuliaForm[exp]
  )
  




ToPython[x_, numpyprefix_: "numpy"] := 
 Module[{expression = x, greekrule, PythonForm, 
   numpypre = numpyprefix, lp, rp, a, 
   b},(*FUNCTION TO CONVERT MATHEMATICA EXPRESSION TO NUMPY;
  ----------------------------------------------------;
  INPUT ARGUMENTS;
  x:your mathematica expression,it can be numbers,literals,
  complexes or lists;
  numpy\[LetterSpace]prefix:string defining your Numpy import prefix,
  e.g.:if your used "import numpy as np",
  your prefix should be the string "np" if your used "from numpy \
import *",your prefix should be the empty string "";
  OUTPUT;
  the Numpy python-ready expression (to be copied as a string);
  !The formatted expression will be copied ot your clipboard,
  ready to paste on Python!;
  ------------------------------------------------------;
  Not tested for every possible combination;use at your risk,
  by Gustavo Wiederhecker*)
  If[numpyprefix == "", sep = "", 
   sep = "."];(*if no prefix is included,
  the "." separator is not used*)lp = "( ";
  rp = " )";
  PythonForm[Rational[a_, b_]] := 
   PythonForm[a] <> "/" <> PythonForm[b];
  PythonForm[Complex[a_, b_]] := 
   "complex" <> lp <> PythonForm[a] <> "," <> PythonForm[b] <> rp;
  PythonForm[Times[a_, b_]] := PythonForm[a] <> " * " <> PythonForm[b];
  PythonForm[Plus[a_, b_]] := 
   lp <> PythonForm[a] <> " + " <> PythonForm[b] <> rp;
  PythonForm[h_[args__]] := 
   numpypre <> sep <> ToLowerCase[PythonForm[h]] <> lp <> 
    PythonForm[args] <> rp;
  PythonForm[Power[a_, b_]] := 
   lp <> PythonForm[a] <> rp <> "**" <> lp <> PythonForm[b] <> rp;
  PythonForm[a_ListQ] := 
   numpypre <> sep <> "array" <> 
    StringReplace[ToString[a], {"{" -> "[", "}" -> "]"}];
  PythonForm[Arg] = numpypre <> sep <> "angle";
  (*Some functions that are note defined in numpy*)
  
  PythonForm[Csc] := "1/" <> numpypre <> sep <> "sin";
  PythonForm[Sec] := "1/" <> numpypre <> sep <> "cos";
  PythonForm[Cot] := "1/" <> numpypre <> sep <> "tan";
  PythonForm[Csch] := "1/" <> numpypre <> sep <> "sinh";
  PythonForm[Sech] := "1/" <> numpypre <> sep <> "cosh";
  PythonForm[Coth] := "1/" <> numpypre <> sep <> "tanh";
  (*Handling arrays*)
  
  PythonForm[List[args__]] := 
   numpypre <> sep <> "array" <> lp <> "[" <> 
    Table[PythonForm[{args}[[ii]]] <> ",", {ii, 1, Length@{args}}] <> 
    "]" <> rp;
  (*Pi and E*)
  PythonForm[\[Pi]] = numpypre <> sep <> "pi";
  PythonForm[E] = numpypre <> sep <> "e";
  (*real numbers,engineering notation*)
  PythonForm[r_Real] := 
   Block[{a = MantissaExponent[r]}, 
    If[r >= 0, ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]], 
     "(" <> ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]] <> ")"]];
  (*Greek characters*)
  greekrule = {"\[Alpha]" -> "alpha", "\[Beta]" -> "beta", 
    "\[Gamma]" -> "gamma", "\[Delta]" -> "delta", 
    "\[CurlyEpsilon]" -> "curlyepsilon", "\[Zeta]" -> "zeta", 
    "\[Eta]" -> "eta", "\[Theta]" -> "theta", "\[Iota]" -> "iota", 
    "\[Kappa]" -> "kappa", "\[Lambda]" -> "lambda", "\[Mu]" -> "mu", 
    "\[Nu]" -> "nu", "\[Xi]" -> "xi", "\[Omicron]" -> "omicron", 
    "\[Pi]" -> "pi", "\[Rho]" -> "rho", 
    "\[FinalSigma]" -> "finalsigma", "\[Sigma]" -> "sigma", 
    "\[Tau]" -> "tau", "\[Upsilon]" -> "upsilon", 
    "\[CurlyPhi]" -> "curlyphi", "\[Chi]" -> "chi", "\[Psi]" -> "psi",
     "\[Omega]" -> "omega", "\[CapitalAlpha]" -> "Alpha", 
    "\[CapitalBeta]" -> "Beta", "\[CapitalGamma]" -> "Gamma", 
    "\[CapitalDelta]" -> "Delta", 
    "\[CapitalEpsilon]" -> "CurlyEpsilon", "\[CapitalZeta]" -> "Zeta",
     "\[CapitalEta]" -> "Eta", "\[CapitalTheta]" -> "Theta", 
    "\[CapitalIota]" -> "Iota", "\[CapitalKappa]" -> "Kappa", 
    "\[CapitalLambda]" -> "Lambda", "\[CapitalMu]" -> "Mu", 
    "\[CapitalNu]" -> "Nu", "\[CapitalXi]" -> "Xi", 
    "\[CapitalOmicron]" -> "Omicron", "\[CapitalPi]" -> "Pi", 
    "\[CapitalRho]" -> "Rho", "\[CapitalSigma]" -> "Sigma", 
    "\[CapitalTau]" -> "Tau", "\[CapitalUpsilon]" -> "Upsilon", 
    "\[CapitalPhi]" -> "CurlyPhi", "\[CapitalChi]" -> "Chi", 
    "\[CapitalPsi]" -> "Psi", "\[CapitalOmega]" -> "Omega"};
  (*Everything else*)
  PythonForm[allOther_] := 
   StringReplace[ToString[allOther, FortranForm], greekrule];
  (*Copy results to clipboard*)CopyToClipboard[PythonForm[expression]];
  PythonForm[expression]]


PythonForm[Rational[a_, b_]] := PythonForm[a] <> "/" <> PythonForm[b];
PythonForm[Complex[a_, b_]] := "complex" <> "(" <> PythonForm[a] <> "," <> PythonForm[b] <> ")";
PythonForm[Times[a_, b_]] := PythonForm[a] <> "*" <> PythonForm[b];
PythonForm[Plus[a_, b_]] := "(" <> PythonForm[a] <> "+" <> PythonForm[b] <> ")";
PythonForm[h_[args__]] := ToLowerCase[PythonForm[h]] <> "(" <> PythonForm[args] <> ")";
PythonForm[Power[a_, b_]] := "(" <> PythonForm[a] <> ")" <> "**" <> "(" <> PythonForm[b] <> ")";
PythonForm[a_ListQ] := StringReplace[ToString[a], {"{" -> "[", "}" -> "]"}];
PythonForm[Arg] = "angle";

(*Some functions that are note defined in numpy*)
PythonForm[Csc] := "1/" <> "sin";
PythonForm[Sec] := "1/" <> "cos";
PythonForm[Cot] := "1/" <> "tan";
PythonForm[Csch] := "1/" <> "sinh";
PythonForm[Sech] := "1/" <> "cosh";
PythonForm[Coth] := "1/" <> "tanh";

(*Handling arrays*)
PythonForm[List[args__]] := "[" <> Table[PythonForm[{args}[[ii]]] <> ",", {ii, 1, Length@{args}}] <> "]";

(*Pi and E*)
PythonForm[Pi] = "pi";
PythonForm[\[Pi]] = "pi";
PythonForm[E] = "e";

(*real numbers,engineering notation*)
PythonForm[r_Real] := Block[{a = MantissaExponent[r]}, 
    If[r >= 0, ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]], 
     "(" <> ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]] <> ")"]];

(*Everything else*)
PythonForm[allOther_] := StringReplace[ToString[allOther, FortranForm], greekrule];

HurToPython[exp_] := (
  CopyToClipboard[PythonForm[exp]];
  PythonForm[exp]
  )









MatlabForm[Rational[a_, b_]] := MatlabForm[a] <> "/" <> MatlabForm[b];
MatlabForm[Complex[a_, b_]] := "(" <> MatlabForm[a] <> "+" <> MatlabForm[b] <> "j" <> ")";
MatlabForm[Times[a_, b_]] := MatlabForm[a] <> "*" <> MatlabForm[b];
MatlabForm[Plus[a_, b_]] := "(" <> MatlabForm[a] <> "+" <> MatlabForm[b] <> ")";
(* Most of the commands starting with Capital letter will be lowered. 
  Ex: Sin[xxx] -> sin(xxx). Note that q1[t] -> q1 for simulation purpose. 
  In HurToolbox, only generalized coordinates and generalized momenta are 
  functions of time. For simulation purpose, the [t] from q1[t], p1[t] 
  will be removed. This treatment is only for the combination of HurToolbox 
  and other numerical packages (e.g., Julia, Matlab). 
*)

MatlabForm[h_[args__]] := If[MatlabForm[args] === "t", MatlabForm[h], 
    ToLowerCase[MatlabForm[h]] <> "(" <> MatlabForm[args] <> ")"];
(* treat Exp -> exp and Sqrt -> sqrt. All other Power will be treated as ()^() *)
MatlabForm[Power[a_, b_]] := If[MatlabForm[b] === "1/2", 
    "sqrt" <> "(" <> MatlabForm[a] <> ")", 
    If[MatlabForm[a] === "exp(1)", 
      "(" <> "exp" <> "(" <> MatlabForm[b] <> ")", 
      "(" <> MatlabForm[a] <> ")" <> "^" <> "(" <> MatlabForm[b] <> ")"]];
MatlabForm[a_ListQ] := StringReplace[ToString[a],{"{" -> "[", "}" -> "]"}];
MatlabForm[Derivative[a_][b_]] := Switch[a, 1, MatlabForm[b] <> "d", 
    2, MatlabForm[b] <> "dd", 3, MatlabForm[b] <> "ddd"];

(* Handling matrix *)

MatlabForm[List[args__]] := (Switch[ Length[Dimensions[List[args]]],
    1,
    temp="["<> Table[MatlabForm[List[args][[i]]]<>",",{i,Length[List[args]]}] <>"]",
    2,
    temp="["<> Table[MatlabForm[List[args][[i]]]<>";",{i,Length[List[args]]}] <>"]"
    ];
  temp=StringReplace[temp, {"[[" -> "[", "];[" -> ";", "];]" -> "]"}];
  StringReplace[temp, {",]" -> "]"}]
  )

(*
StringReplace[ToString[a], 
    {"},{" -> ";", "{{" -> "[", "}}" -> "]"}];

"(" <> "[" <> Table[MatlabForm[{args}[[ii]]] <> " ", 
  {ii, 1, Length@{args}}] <> "]" <> ")";
*)

(*Pi and E*)
MatlabForm[\[Pi]] = "pi";
MatlabForm[Pi] = "pi";
MatlabForm[E] = "exp(1)";
(*real numbers,engineering notation*)
MatlabForm[r_Real] := Block[{a = MantissaExponent[r]}, 
    If[r >= 0, ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]], 
     "(" <> ToString[N[a[[1]], 6]] <> "e" <> ToString[a[[2]]] <> ")"]];
(*Greek characters*)
greekrule = {"\[Alpha]" -> "alpha", "\[Beta]" -> "beta", 
    "\[Gamma]" -> "gamma", "\[Delta]" -> "delta", 
    "\[CurlyEpsilon]" -> "curlyepsilon", "\[Zeta]" -> "zeta", 
    "\[Eta]" -> "eta", "\[Theta]" -> "theta", "\[Iota]" -> "iota", 
    "\[Kappa]" -> "kappa", "\[Lambda]" -> "lambda", "\[Mu]" -> "mu", 
    "\[Nu]" -> "nu", "\[Xi]" -> "xi", "\[Omicron]" -> "omicron", 
    "\[Pi]" -> "pi", "\[Rho]" -> "rho", 
    "\[FinalSigma]" -> "finalsigma", "\[Sigma]" -> "sigma", 
    "\[Tau]" -> "tau", "\[Upsilon]" -> "upsilon", 
    "\[CurlyPhi]" -> "curlyphi", "\[Chi]" -> "chi", "\[Psi]" -> "psi",
     "\[Omega]" -> "omega", "\[CapitalAlpha]" -> "Alpha", 
    "\[CapitalBeta]" -> "Beta", "\[CapitalGamma]" -> "Gamma", 
    "\[CapitalDelta]" -> "Delta", 
    "\[CapitalEpsilon]" -> "CurlyEpsilon", "\[CapitalZeta]" -> "Zeta",
     "\[CapitalEta]" -> "Eta", "\[CapitalTheta]" -> "Theta", 
    "\[CapitalIota]" -> "Iota", "\[CapitalKappa]" -> "Kappa", 
    "\[CapitalLambda]" -> "Lambda", "\[CapitalMu]" -> "Mu", 
    "\[CapitalNu]" -> "Nu", "\[CapitalXi]" -> "Xi", 
    "\[CapitalOmicron]" -> "Omicron", "\[CapitalPi]" -> "Pi", 
    "\[CapitalRho]" -> "Rho", "\[CapitalSigma]" -> "Sigma", 
    "\[CapitalTau]" -> "Tau", "\[CapitalUpsilon]" -> "Upsilon", 
    "\[CapitalPhi]" -> "CurlyPhi", "\[CapitalChi]" -> "Chi", 
    "\[CapitalPsi]" -> "Psi", "\[CapitalOmega]" -> "Omega"};
(*Everything else*)
MatlabForm[allOther_] := 
   StringReplace[ToString[allOther, FortranForm], greekrule];

HurToMatlab[exp_] := (
  CopyToClipboard[MatlabForm[exp]];
  MatlabForm[exp]
  )




Clear[HurAutoDiff]
(*elementary diff rules*)
(*linear*)
HurAutoDiff /: HurAutoDiff[a_, adx_] + HurAutoDiff[b_, bdx_] := HurAutoDiff[a + b, adx + bdx]
HurAutoDiff /: HurAutoDiff[a_, adx_]*b_ := HurAutoDiff[a*b, adx*b]
(*product rule*)
HurAutoDiff /: HurAutoDiff[a_, adx_]*HurAutoDiff[b_, bdx_] := HurAutoDiff[a*b, adx*b + a*bdx]
(*constant rule cannot be defined. *)
(*HurAutoDiff /: b_:=HurAutoDiff[b,0]*)
HurAutoDiff /: HurAutoDiff[a_, adx_] + b_ := HurAutoDiff[a + b, adx] 
(*Power rule*)
HurAutoDiff /: HurAutoDiff[a_, adx_]^b_ := HurAutoDiff[a^b, a^(b - 1)*adx*b]
(*Exponent rule*)
HurAutoDiff /: a_^HurAutoDiff[b_, bdx_] := HurAutoDiff[a^b, a^b*Log[a]*bdx]
HurAutoDiff /: HurAutoDiff[a_, adx_]^HurAutoDiff[b_, bdx_] := HurAutoDiff[a^b, a^b*(b*adx/a+Log[a]*bdx)]
HurAutoDiff /: HurAutoDiff[a_, adx_]^(1/HurAutoDiff[b_, bdx_])_ := HurAutoDiff[a^(1/b), a^(1/b)*(adx/a/b-Log[a]*bdx/b^2)]
(*quocient rule*)
HurAutoDiff /: HurAutoDiff[a_, adx_]/HurAutoDiff[b_, bdx_] := HurAutoDiff[a/b, (adx*b - a*bdx)/b^2]
HurAutoDiff /: Sin[HurAutoDiff[a_, adx_]] := HurAutoDiff[Sin[a], adx*Cos[a]]
HurAutoDiff /: Cos[HurAutoDiff[a_, adx_]] := HurAutoDiff[Cos[a], -adx*Sin[a]]
HurAutoDiff /: Tan[HurAutoDiff[a_, adx_]] := HurAutoDiff[Tan[a], adx*Sec[a]^2]
HurAutoDiff /: Cot[HurAutoDiff[a_, adx_]] := HurAutoDiff[Cot[a], -adx*Csc[a]^2]
HurAutoDiff /: Sec[HurAutoDiff[a_, adx_]] := HurAutoDiff[Sec[a], adx*Sec[a]*Tan[a]]
HurAutoDiff /: Csc[HurAutoDiff[a_, adx_]] := HurAutoDiff[Csc[a], -adx*Cot[a]*Csc[a]]
HurAutoDiff /: Sinh[HurAutoDiff[a_, adx_]] := HurAutoDiff[Sinh[a], adx*Cosh[a]]
HurAutoDiff /: Cosh[HurAutoDiff[a_, adx_]] := HurAutoDiff[Cosh[a], adx*Sinh[a]]
HurAutoDiff /: ArcSin[HurAutoDiff[a_, adx_]] := HurAutoDiff[ArcSin[a], adx/Sqrt[1-a^2]]
HurAutoDiff /: ArcCos[HurAutoDiff[a_, adx_]] := HurAutoDiff[ArcCos[a], -(adx/Sqrt[1-a^2])]
HurAutoDiff /: ArcTan[HurAutoDiff[a_, adx_]] := HurAutoDiff[ArcTan[a], adx/(1+a^2)]
HurAutoDiff /: Exp[HurAutoDiff[a_, adx_]] := HurAutoDiff[E^a, adx*E^a]
HurAutoDiff /: Log[HurAutoDiff[a_, adx_]] := HurAutoDiff[Log[a], adx/a]
HurAutoDiff /: Log10[HurAutoDiff[a_, adx_]] := HurAutoDiff[Log[a]/Log[10], adx/(a*Log[10])]
HurAutoDiff /: Log2[HurAutoDiff[a_, adx_]] := HurAutoDiff[Log[a]/Log[2], adx/(a*Log[2])]

Clear[HurDefineAutoDiffFunc]
HurDefineAutoDiffFunc[fn_] := HurAutoDiff /: fn[HurAutoDiff[a_, adx_]] = HurAutoDiff[fn[a], D[fn[a + x*adx], x] /. x -> 0]



HurFunctionForm4Var[Derivative[a_][b_]] := Switch[a, 1, Symbol[ToString[b] <> "d"], 2, Symbol[ToString[b] <> "dd"], 3, Symbol[ToString[b] <> "ddd"]];
HurFunctionForm4Var[allOther_] := allOther;
HurConvertVar[x__] := (xs = Flatten[List[x]]; nxs = Length[xs];
  Table[
    (*Clear[HurFunctionForm4Var[Head[xs[[i]]]]];*)
    xs[[i]] -> HurFunctionForm4Var[Head[xs[[i]]]], {i, nxs}]
  )

HurGetFunctionForm[var_, fn_] := (var1 = Flatten[var];
  rep = HurConvertVar[var1];
  varconvt = var1 /. rep;
  func = fn /. rep;
  Function[Evaluate[varconvt], Evaluate[func]]
  )

HurGetGradAutoDiff[fn_, val__] := (
  nval = Length[val];
  Table[
   val1 = val;
   val1[[i]] = ToExpression["HurAutoDiff[" <> ToString[1] <> ",1]"];
   valSq = Flatten[List[val1]] /. List -> Sequence;
   fnval = fn[valSq];
   If[Head[fnval] === HurAutoDiff, fnval[[2]], 0]
   , {i, nval}]
  )

HurGetGradAutoDiffValue[fn_, val__] := (
  nval = Length[val];
  Table[
   val1 = val;
   val1[[i]] = ToExpression["HurAutoDiff[" <> ToString[1] <> ",1]"];
   valSq = Flatten[List[val1]] /. List -> Sequence;
   fn[valSq]
   , {i, nval}]
  )








End[];

EndPackage[];