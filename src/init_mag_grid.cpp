// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>

#include "../include/aether.h"

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_mag_grid(Planets planet, Inputs input, Report &report) {
  
  std::string function = "Grid::init_mag_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
  // This is just an example:
  
  Inputs::grid_input_struct grid_input = input.get_grid_inputs();
  
  int64_t iLon, iLat, iAlt;


  
  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lon1d(nLons);
  float dlon = (grid_input.lon_max - grid_input.lon_min) / (nLons-2*nGCs);
  for (iLon=0; iLon < nLons; iLon++)
    lon1d[iLon] = grid_input.lon_min + (iLon-nGCs+0.5) * dlon;
  
  for (iLat=0; iLat < nLats; iLat++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      magLon_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;
      magPhi_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;
    }
  }
  cout << "!!!!!!!!!!! HA1 "<< nLats<<endl;
  
  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lat1d(nLats);
  //cout << "!!!!!!!!!!! HA1a "<<endl;


  fvec lshell(nLats);

  float dlat = (grid_input.lat_max - grid_input.lat_min) / (nLats-2*nGCs);
  
  lshell(0) = 1/pow(cos(grid_input.lat_min),2.0);
  lshell(nLats-1) = 1/pow(cos(grid_input.lat_max),2.0);
  float dlshell = (lshell(nLats-1)-lshell(0))/nLats;
  
  for (iLat=1; iLat < nLats; iLat++){
    lshell[iLat] = lshell[iLat-1]+dlshell;
    //cout << "Lshell = " << lshell(iLat)<<endl;
    //lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;
  }
  
  for (iLon=0; iLon < nLons; iLon++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      
      //magP_scgc.subcube(iLon, 0, iAlt, iLon, nLats-1, iAlt) = lshell;
      for (iLat=0; iLat<nLats; iLat++){
	//cout << "L fill " << lshell[iLat] << endl;

	magP_scgc(iLon, iLat, iAlt) = lshell[iLat];
	//if (iLon==12 and iLat==10 and iAlt==1){
	  //cout << iLon<<" "<< iLat <<" "<< iAlt <<" P fill " << magP_scgc(iLon,iLat,iAlt) << endl;
	//}
	//cout << iLon<<" "<< iLat <<" "<< iAlt <<" P fill " << magP_scgc(12,10,1) << endl;
	//cout << "magP_scgc[iLon, iLat, iAlt] = lshell[iLat];"<<" "<<magP_scgc[iLon, iLat, iAlt] <<" " <<lshell[iLat] << endl;
      }
    }
  }

  //cout << magP_scgc(12, 10, 1)<<endl;
  
  
  
  // fill along the field 
  float qS;
  float qN;
  float Lshell;
  float Lon;
  float Gamma;
  double q[nZ];  

  // set the min alt and gamma factor for filling line
  float AltMin = grid_input.alt_min; 
  Gamma = 2.0;
  
  

  for (iLon=0; iLon < nLons; iLon++) {
    for (iLat=0; iLat < nLats; iLat++) {
      Lshell  = magP_scgc(iLon,iLat,1);
      Lon     = magPhi_scgc(iLon,iLat,1);
      
      // get the q value for N and S hemisphere
      cout << iLon << " " << iLat << " L= "<<Lshell<<" "<< magP_scgc(iLon,iLat,1)<<endl;
      auto Qvals = lshell_to_qn_qs(planet, Lshell, Lon, AltMin, report);
      qN = Qvals.first;
      qS = Qvals.second;

      //cout << qN << endl;
      //cout << qS << endl;

      //fill in the q array for this P and Phi
      fill_dipole_q_line(qN, qS, Gamma, nZ, Lshell, Lon, q, report);

      //copy this q array into
      for (iAlt=0;iAlt<nAlts; iAlt++){
	magQ_scgc(iLon, iLat, iAlt) = q[iAlt];
      }
    }
  }
  

  // fill x y z values
  float Llr[3], Xyz[3];
  int iX, iY, iZ;
  
  for (iX=0; iX<nX; iX++){
      for (iY=0; iY<nY; iY++){
	  for (iZ=0; iZ<nZ; iZ++){
	    
	    
	    // For given q and p we can now find cooresponding r and theta in
	    // dipole coord. Starty by numerically solving for r (normalized) 
	    // using equation 4 of Huba et al 2000 given by q^2r^4+1/p r -1 = 0
	    auto rtheta =
	      p_q_to_r_theta(magP_scgc(iX,iY,iZ), magQ_scgc(iX,iY,iZ));
	    
	    float r = rtheta.first;
	    float theta = rtheta.second;
  
	    //cout << "i, q " << i << "  " << q[i] << endl;
	    //cout << "i, x " << i << "  " << x[i] << endl;
	    //cout << "i, r, theta " << i << "  " << r[i]<<" "<<theta[i] << endl << endl;
	    
	    cout << iX <<" "<< iY << " " << iZ << " " << magPhi_scgc(iX,iY,iZ) << endl;
	    Llr[0] = magPhi_scgc(iX,iY,iZ);
	    Llr[1] = 0.5*cPI-theta;
	    Llr[2] = r;

	    //if (iZ==nZ/2){
	    //  cout << magP_scgc(iX,iY,iZ) <<endl;
	    //  cout << Llr[0]<<" "<<Llr[1]<<" "<<Llr[2] << endl;
	    //}
	    
	    transform_llr_to_xyz(Llr, Xyz);
	    magX_scgc(iX,iY,iZ)=Xyz[0];
	    magY_scgc(iX,iY,iZ)=Xyz[1];
	    magZ_scgc(iX,iY,iZ)=Xyz[2];
	  }
      }
  }


  // save 3D mag grid for examination
  std::fstream gridfile;
  gridfile.open ("grid3D.dat",ios::out);
  gridfile.precision(std::numeric_limits<long double>::digits10);
  
  // write header
  gridfile << "VARIABLES = \"X\", \"Y\", \"Z\" " << endl;
  gridfile << "Zone I = "<< nZ << ",J = " << nY << ",K = "<< nX
	   <<", DATAPACKING=POINT" << endl;
  
  // write grid data
  for (iX=0; iX<nX; iX++){
    for (iY=0; iY<nY; iY++){
      for (iZ=0; iZ<nZ; iZ++){
	
	gridfile << std::fixed << magX_scgc(iX,iY,iZ)
		 <<" "<< std::fixed << magY_scgc(iX,iY,iZ)
		 <<" "<< std::fixed << magZ_scgc(iX,iY,iZ)
		 << endl;
      }
    }
  }
  gridfile.close();


  // save 3D mag grid slice for examination
  std::fstream gridfileslice;
  gridfileslice.open ("grid_slice.dat",ios::out);

  // write header
  gridfileslice << "VARIABLES = \"X\", \"Y\", \"Z\" " << endl;
  gridfileslice << "Zone I = " << nZ << ",J = "<< nY
	   <<", DATAPACKING=POINT" << endl;
  
  // write grid data
  for (iY=0; iY<nY; iY++){
    for (iZ=0; iZ<nZ; iZ++){
      gridfileslice << magX_scgc(1,iY,iZ)
		    <<" "<< magY_scgc(1,iY,iZ)
		    <<" "<< magZ_scgc(1,iY,iZ)
		    << endl;
    }
  }
  
  gridfileslice.close();

  IsGeoGrid = 0;
  
  // Calculate the radius, etc:
  
  //  fill_grid_radius(planet, report);
  
  //  fill_grid_bfield(planet, input, report);
  
  
  // We want to now set up our p, phi, s coordinates as defined by Huba et al 2000
  // this cooresponds to x, y, and z in our grid object
  

//  float AltMin = grid_input.alt_min; 
//  float qS;
//  float qN;
//  float Lshell;
//  float Lon;
//  float Gamma;
//  
//  Lshell  = 4.0;
//  Lon     = 0.0;
//  Gamma = 2.0;
//  
//  auto Qvals = lshell_to_qn_qs(planet, Lshell, Lon, AltMin, report);
//
//  qN = Qvals.first;
//  qS = Qvals.second;
//
//  cout << qN << endl;
//  cout << qS << endl;
//
//  double q[nZ];
//  fill_dipole_q_line(qN, qS, Gamma, nZ, Lshell, Lon, q, report);
//  
  report.exit(function);
}

// ----------------------------------------------------------------------
// Routine to find q_N and q_S for a given L 
// 
// ----------------------------------------------------------------------
std::pair<float,float> Grid::lshell_to_qn_qs(Planets planet, float Lshell, float Lon, float AltMin, Report &report) {
  std::string function = "Grid::lshell_to_qn_qs";
  static int iFunction = -1;
  report.enter(function, iFunction);

  float qN,qS;
  
  float XyzDipoleLeft[3], XyzDipoleMid[3], XyzDipoleRight[3];
  float XyzGeoLeft[3], XyzGeoMid[3], XyzGeoRight[3];
  float rGeoLeft, rGeoMid, rGeoRight;
  float LlrDipoleLeft[3], LlrDipoleMid[3], LlrDipoleRight[3];
  float ThetaTilt, PhiTilt;
  float Lat, Radius, rMin;
  // Named dimension constants
  static int Lon_= 0, Lat_= 1, Radius_= 2;
 
  //bound vars for bisection search
  float ThetaRight, ThetaLeft, ThetaMid;
  float rDipoleLeft,rDipoleMid,rDipoleRight;
  
  //Stopping condition for bisection search
  float DeltaTheta;
  float Tolerance = 1e-4;
    
  // status vars for bisection search
  int iStatusLeft, iStatusRight, iStatusMid;
  // note we normalize Lshell by equatorial radius
  float RadiusEq = planet.get_radius(0.0);


  

  // loop for qN and qS
  for(int iQ = 0; iQ < 2; iQ++){

    if (iQ == 0){
      // set initial left, mid, right bounds for bisection search for qN
      ThetaRight = 0.5*cPI;
      ThetaLeft = 1.0*cDtoR;
      ThetaMid = 0.5*(ThetaRight+ThetaLeft);
    }else{
      // set initial left, mid, right bounds for bisection search for qS
      ThetaLeft = 0.5*cPI;
      ThetaRight = 179.0*cDtoR;
      ThetaMid = 0.5*(ThetaRight+ThetaLeft);
    }
      
    // Initial stopping condition stopping condition
    DeltaTheta = abs(ThetaLeft-ThetaRight);
    
    
    // start bisection search for qN
    while( DeltaTheta > Tolerance ) {
      
      // find rDipole that cooresponds to these Left,Mid,Right
      // ThetaDipole values 
      rDipoleLeft   = Lshell * pow(sin(ThetaLeft),2.0);
      rDipoleMid    = Lshell * pow(sin(ThetaMid),2.0);
      rDipoleRight  = Lshell * pow(sin(ThetaRight),2.0);
      
      // Compute XyzDipole for left, mid,right states
      LlrDipoleLeft[Lon_] = Lon;
      LlrDipoleLeft[Lat_] = 0.5*cPI-ThetaLeft;
      LlrDipoleLeft[Radius_] = rDipoleLeft;
      transform_llr_to_xyz(LlrDipoleLeft, XyzDipoleLeft);
      
      LlrDipoleMid[Lon_] = Lon;
      LlrDipoleMid[Lat_] = 0.5*cPI-ThetaMid;
      LlrDipoleMid[Radius_] = rDipoleMid;
      transform_llr_to_xyz(LlrDipoleMid, XyzDipoleMid);
      
      LlrDipoleRight[Lon_] = Lon;
      LlrDipoleRight[Lat_] = 0.5*cPI-ThetaRight;
      LlrDipoleRight[Radius_] = rDipoleRight;
      transform_llr_to_xyz(LlrDipoleRight, XyzDipoleRight);
      
      // Transform to XyzGeo and unnormalize
      convert_dipole_geo_xyz(planet, XyzDipoleLeft, XyzGeoLeft);
      convert_dipole_geo_xyz(planet, XyzDipoleMid, XyzGeoMid);
      convert_dipole_geo_xyz(planet, XyzDipoleRight, XyzGeoRight);
      
      //cout << "XyzGeoLeft[0]" << XyzGeoLeft[0] << endl;
      //cout << "XyzGeoLeft[1]" << XyzGeoLeft[1] << endl;
      //cout << "XyzGeoLeft[2]" << XyzGeoLeft[2] << endl;

      XyzGeoLeft[0]=XyzGeoLeft[0]*RadiusEq;
      XyzGeoLeft[1]=XyzGeoLeft[1]*RadiusEq;
      XyzGeoLeft[2]=XyzGeoLeft[2]*RadiusEq;

      abort;
      
      XyzGeoMid[0]=XyzGeoMid[0]*RadiusEq;
      XyzGeoMid[1]=XyzGeoMid[1]*RadiusEq;
      XyzGeoMid[2]=XyzGeoMid[2]*RadiusEq;
      
      XyzGeoRight[0]=XyzGeoRight[0]*RadiusEq;
      XyzGeoRight[1]=XyzGeoRight[1]*RadiusEq;
      XyzGeoRight[2]=XyzGeoRight[2]*RadiusEq;
      
      // Compute radius in geo coordinate for comparison to rmin
      rGeoLeft = 
	sqrt(pow(XyzGeoLeft[0],2)+pow(XyzGeoLeft[1],2)+pow(XyzGeoLeft[2],2));
      rGeoMid = 
	sqrt(pow(XyzGeoMid[0],2)+pow(XyzGeoMid[1],2)+pow(XyzGeoMid[2],2));
      rGeoRight =
	sqrt(pow(XyzGeoRight[0],2)+pow(XyzGeoRight[1],2)+pow(XyzGeoRight[2],2));
      
      // get rmin for given latitude. Radius is lat dependent in general.
      // also find status in (0) or out (1) of rMin
      Lat    = 0.5*cPI-acos(XyzGeoLeft[2]/rGeoLeft);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoLeft < rMin){
	iStatusLeft = 0;  
      }else{
	iStatusLeft = 1;  
      }
      
      Lat    = 0.5*cPI-acos(XyzGeoMid[2]/rGeoMid);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoMid < rMin){
	iStatusMid = 0;  
      }else{
	iStatusMid = 1;  
      }
      
      Lat    = 0.5*cPI-acos(XyzGeoRight[2]/rGeoRight);
      Radius = planet.get_radius(Lat);
      rMin   = Radius+AltMin;
      if (rGeoRight < rMin){
	iStatusRight = 0;  
      }else{
	iStatusRight = 1;  
      }
      
      // Use status values to update left, right and mid values of theta
      if (iStatusMid == 0) {
	if (iStatusRight == 1){
	  // Mid becomes left and right stays right
	  ThetaLeft = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}else{
	  // Mid becomes right and left stays left
	  ThetaRight = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}
      }else{
	if (iStatusRight == 0){
	  // Mid becomes left and right stays right
	  ThetaLeft = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}else{
	  // Mid becomes right and left stays left
	  ThetaRight = ThetaMid;
	  ThetaMid = 0.5*(ThetaLeft+ThetaRight);
	}
      }
      // Update stopping condition
      DeltaTheta = abs(ThetaLeft-ThetaRight);
    }
    
    //set the q value
    rDipoleMid    = Lshell * pow(sin(ThetaMid),2.0);
    if (iQ ==0){
      qN = pow(rDipoleMid,-2.0)*cos(ThetaMid);
      cout << "!!! For L = " << Lshell << endl;
      cout << "!!! qN = " << qN << endl;
      cout << "!!! ThetaMid = " << ThetaMid*cRtoD << endl;
    }else{
      qS = pow(rDipoleMid,-2.0)*cos(ThetaMid);
      cout << "!!! qS = " << qS << endl;
    }
  }

  report.exit(function);
  return {qN,qS};
  
}

// -----------------------------------------------------------------------
// Convert XyzDipole to XyzGeo
//  
// -----------------------------------------------------------------------

void Grid::convert_dipole_geo_xyz(Planets planet, float XyzDipole[3], float XyzGeo[3]) {
  float XyzRemoveShift[3];
  float XyzRemoveTilt[3];
  float XyzRemoveRot[3];

  // get planetary parameters, use radius at equator for Lshell reference
  float magnetic_pole_tilt = planet.get_dipole_tilt();
  float magnetic_pole_rotation = planet.get_dipole_rotation();
  float radius = planet.get_radius(0.0);

  
  // get the dipole shift, but normalize it to equatorial radius 
  float dipole_center[3];
  std::vector<float> temp_dipole_center = planet.get_dipole_center();
  transform_float_vector_to_array(temp_dipole_center, dipole_center);

  dipole_center[0]=dipole_center[0]/radius;
  dipole_center[1]=dipole_center[1]/radius;
  dipole_center[2]=dipole_center[2]/radius;

  // Remove Tilt
  transform_rot_y(XyzDipole, magnetic_pole_tilt, XyzRemoveTilt);

  // Remove Rot
  transform_rot_z(XyzRemoveTilt, magnetic_pole_rotation, XyzRemoveRot);

  // Remove Shift
  vector_add(XyzRemoveRot, dipole_center, XyzGeo);

//  cout << "XyzDipole[0]" << XyzDipole[0] << endl;
//  cout << "XyzDipole[1]" << XyzDipole[1] << endl;
//  cout << "XyzDipole[2]" << XyzDipole[2] << endl;
//
//  cout << "XyzGeo[0]" << XyzGeo[0] << endl;
//  cout << "XyzGeo[1]" << XyzGeo[1] << endl;
//  cout << "XyzGeo[2]" << XyzGeo[2] << endl;
  
}

// ----------------------------------------------------------------------
// Routine to fill in the q values for a particular L and lon
// using equations 7-8 from Huba et al 2000
// ----------------------------------------------------------------------
void Grid::fill_dipole_q_line(float qN, float qS, float Gamma, int nZ, float Lshell, float Lon, double *q,Report &report) {
  std::string function = "Grid::fill_dipole_q_line";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
  double x[nZ];
  double r[nZ];
  double theta[nZ];
  double Dx;
  float Llr[3], Xyz[3];
  int DoTestLine = 0;
  cout <<"test"<<endl;
  //open test file for writing the grid data for plotting
  std::fstream gridfile;
  if (DoTestLine==1){
    gridfile.open ("grid.dat",ios::out);
  }
  
  // set the Dx (c in equation 8 from Huba et al 2000)
  // Note this equation has a typo in it. The proper
  // version  would be found by defining the bounds of
  // x from equation 7, and then dividing that into
  // equal segments.
  //Dx = 2.0*(1.0-sinh(Gamma*qN))/((static_cast<float>(nZ)-1.0)*sinh(Gamma*qS));
  Dx = (sinh(Gamma*qS)-sinh(Gamma*qN))/((static_cast<float>(nZ)-1.0)*sinh(Gamma*qS));
  //Dx = 2.0/(static_cast<float>(nZ)-1.0);
  //Dx = (static_cast<double>(qN)-static_cast<double>(qS))/(static_cast<double>(nZ)-1.0);

  //cout << "Dx = " << Dx << endl;
  //cout << "Gamma = " << Gamma << endl;
  //cout << "nZ = " << nZ << endl;
  //cout << "qN = " << qN << endl;
  //cout << "qS = " << qS << endl;

  // set initial x[0] value using eq. 7  from Huba et al with qi=qN
  x[0] =  sinh(Gamma*qN)/sinh(Gamma*qS);
  q[0] =  asinh(x[0]*sinh(Gamma*qS))/Gamma;
    
  // fill x_i=x_(i-1)+Dx
  for(int i = 1; i < nZ; i++){
    x[i] = x[i-1]+Dx;
    q[i] = asinh(x[i]*sinh(Gamma*qS))/Gamma;

    // For given q and p we can now find cooresponding r and theta in
    // dipole coord. Starty by numerically solving for r (normalized) using
    // equation 4 of Huba et al 2000 given by q^2r^4+1/p r -1 = 0
     auto rtheta = p_q_to_r_theta(Lshell, q[i]);
     
     r[i] = rtheta.first;
     theta[i] = rtheta.second;
     
     //cout << "i, q " << i << "  " << q[i] << endl;
     //cout << "i, x " << i << "  " << x[i] << endl;
     //cout << "i, r, theta " << i << "  " << r[i]<<" "<<theta[i] << endl << endl;

     Llr[0] = Lon;
     Llr[1] = 0.5*cPI-theta[i];
     Llr[2] = r[i];
       
     transform_llr_to_xyz(Llr, Xyz);

     if (DoTestLine==1){
       gridfile << Xyz[0] <<" "<< Xyz[1] <<" "<< Xyz[2] << endl;

     }
  }

  if (DoTestLine==1){
    gridfile.close();
  }
  
  report.exit(function);
  return ;
}



// ----------------------------------------------------------------------
// Routine to convert p and q to r and theta. Appraoch is to first solve
// for r using eq 4 from Huba et al 2000. q^2*r^4+1/q*r-1=0
// This is solved numerically using Newton-Raphson (NR) technique.
// Once we know r we can reover theta from p=r*1/(sin theta)^2.
// note r here is normalized to planet radius. 
// 
// ----------------------------------------------------------------------
std::pair<float,float> Grid::p_q_to_r_theta(float p, float q) {
  //return quanties
  float r, theta;
  // function value and derivative for NR method
  float Func, dFunc;
  // tolerance for root finding
  float Tolerance = 0.00001;

  // initial guess for r
  r = 1.0;
  
  Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
  dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;

  // apply NR iterations to get 
  while( abs(Func/dFunc) > Tolerance ) { 
    Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
    dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;
    
    // in NR method r(i+1)=r(i)-f/f' for each iteration
    r=r-Func/dFunc;
  }
  
  // now that r is determined we can solve for theta
  //theta = asin(sqrt(r/p));
  theta = acos(q*pow(r,2.0));
  
  
  //cout << "for p,q = " << p <<" "<< q << endl;
  //cout << "  r     = " << r << endl;
  //cout << "  theta = " << theta << endl;
  //cout << endl;
  
  return {r,theta};
}
