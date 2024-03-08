// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include "../include/aether.h"


#define WLEVEL 0
#define QGRDTYPE 2
 
//  void NeutralSpeciesFillMagLine(){

//     int64_t nAlts = grid.get_nAlts();

//   for (int iAlt = 1; iAlt < nAlts; iAlt++) {
//     species[iSpecies].density_scgc.slice(iAlt) =
//       temperature_scgc.slice(iAlt - 1) /
//       temperature_scgc.slice(iAlt) %
//       species[iSpecies].density_scgc.slice(iAlt - 1) %
//       exp(-grid.dalt_lower_scgc.slice(iAlt) /
//           species[iSpecies].scale_height_scgc.slice(iAlt));
//   }


//  }

void NeutralSpeciesFillFullMagGrid(){

// cycle through the range of the grid
  
}




void HalfAndHalfLinrange(std::vector<double>& qLin, int size, double xL, double xR, 
                    double offSetFromMid){    
    assert(size % 2 == 0 && "size must be even in HalfAndHalfLinrange ");  
    std::vector<double> qLinPos;
    
    qLin = linspace(xL, -offSetFromMid, size/2); //#p_i from #q_s
    for (auto element : qLin) {
      std::cout << element << " ";
    }

    qLinPos = linspace(offSetFromMid, xR, size/2);
    
    qLin.insert(qLin.end(), qLinPos.begin(), qLinPos.end()); //symmetrical about 0; 0 is not included     
}

void GetDipoleXyz(double r_i, double ph_j, double q_i, double pLine, 
                 double &x, double &y, double &z, double &th){
      x = sqrt(pow(r_i,3)/pLine)*cos(ph_j);
      y = sqrt(pow(r_i,3)/pLine)*sin(ph_j);
      z = pow(r_i,3)*q_i;
      th = acos(pow(r_i,2)*q_i);
}

void QuartFunAndDeriv(double &f, double &dfdx, double rt, double  pMag, double q) { 
  double r=rt;
  f = pow(q,2)*pow(r,4) + r/pMag -1.0;
  dfdx = 4*pow(q,2)*pow(r,3) + 1/pMag;
  return ;
}


void SolveDipoleEq(double &xRt, double pMag, double qMag){

  double maxit =100, f,df,dx, dxold,x1,x2,xl,xh,fl,fh,dxL,dxH,dxN,tmp1,
         tolRt_ = 1e-3,
         tol_ = 1e-3;
        
  x1=1.; //x low
  x2=1.; //x high

  int maxExpandTrials=40;

  QuartFunAndDeriv(fl,df, x1,pMag,qMag);

  for (int i =0; i<=maxExpandTrials; i++){ //correct the f+/- interval
      QuartFunAndDeriv(fh,df,x2,pMag,qMag);
      // SHOW(x1);SHOW(x2);SHOW(fl);SHOW(fh); 
      if (fabs(fl) <= tol_){xRt=x1; return;};
      if (fabs(fh) <= tol_){xRt=x2; return;};
  
      if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)){      
        #if WLEVEL==1
        cout<<i<<":  Root must be bracketed in SolveDipoleEq=  "<<fl<<"  "<<fh<<endl;
        #endif
        x2*=1.3;            
      } else {
        break;
      }
    
  }
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
      throw std::runtime_error("Root must be bracketed in SolveDipoleEq");

  if (fl < 0.0){
      xl=x1;
      xh=x2;
  } else {
    xh=x1; //flip them!
    xl=x2; 
  }
  xRt=0.5*(x1+x2);   //guess
  dxold=fabs(x2-x1); //step
  dx=dxold;

  QuartFunAndDeriv(f,df,xRt,pMag,qMag);

  for(int i=1; i<=maxit; i++){
        
    dxL = (xRt - xl)*df-f;
    dxH = (xRt - xh)*df-f; //for out-of-range check
    dxN = fabs(f);

    if ((dxL*dxH>0) || fabs(2*dxN > df*dxold)){
      // out of range or not decreasing fast enough
      dx=0.4*(xh-xl);
      xRt = xl+dx;
      if (fabs(xl-xRt) <tol_){return;}
    } else { //NR step
      dxold=dx;      
      dx=f/df;
      tmp1=xRt;
      xRt = xRt - dx;
      if( fabs(tmp1-xRt) <tol_){return;} //root's not changing

    }
    if (fabs(dx) < tolRt_) return; //converging
    
    QuartFunAndDeriv(f,df,xRt,pMag,qMag);
    
    if (f < 0.0){
      xl=xRt;
      }
    else {
      xh=xRt;
    }
    if(i==maxit){
        throw std::runtime_error("max_iterations reached in 'SolveDipoleEq' ");  
    }
  }
  // SHOW("SolveDipoleEq");SHOW(x0);SHOW(pMag);SHOW(qMag);SHOW(x1);SHOW(df);SHOW(y0);
  // SHOW("SolveDipoleEq result:"); SHOW(r);SHOW(i);
  
}

void GetPsFromQs(double &ps, double qs, double &thetaOfPs){  
  //for a given q calculates corresponding p_s and l_s
  double psSmall = 1.e-3;
  ps = 1.0/(1.0-pow(qs,2)); //find corresponding p_s(q_s)
  thetaOfPs = asin( sqrt( 1./ std::max(ps,psSmall) ) );

  return;
}
      
void ProjectQsOnAllMagLinesAnsMask(Grid *grid, std::vector<double> q1d, 
                            int j_ph, int iOfQs, int lOfPs){
  assert(q1d.size() % 2 == 0 && "size of q1d must be even ");  

  int iMax=q1d.size();
  int iEq = iMax/2;
  
  double qs = q1d[iOfQs];


  if (qs<0){
    for (int i = 0; i < iEq; i++) 
        grid->magMask(j_ph, lOfPs, i) = (i<iOfQs) ? 0.0 : 1.0;
    
  } else if (qs>0){
  
    for (int i = iEq; i < iMax; i++) 
        grid->magMask(j_ph, lOfPs, i) = (i<iOfQs) ? 1.0 : 0.0;
  
  } else{
    throw("qs=0: not allowed in init_mag.cpp");
  }

  return;
}

void ScatPlotArmaArray(std::vector<double> q1d, std::vector<double> p1d, arma_cube magMask, int IdLon){
  //print three columns of data for a scatter plotp1d.size();
  
  char cwd1[100];
  std::string fdir(getcwd(cwd1, sizeof(cwd1)));
  auto fileName = fdir + "/" + "_2ddata.dat";   
  SHOW(fileName) 
  int len1 = q1d.size();
  int len2 = p1d.size();
  
    std::ofstream output_file(fileName);
    
    for (int j = 0; j < len2; j++) {        
      for (int i = 0; i < len1; i++) {                
            output_file << q1d[i] <<"\t"<< p1d[j]<<"\t"<<magMask(IdLon,j,i) <<std::endl;           
      }
    }
    output_file.close();    
  const std::string pyFileToExecute = fdir + "/" + "./_my_scat_plot_arma_array.py";  
  SHOW(pyFileToExecute)  
  system(pyFileToExecute.c_str());
}

void FillAllCoordsOnMagGrid(Grid *grid, Planets planet, std::vector<double> q1d, 
                                         std::vector<double> p1d, 
                                         std::vector<double> phiLon, int idOfLon){
  int Nlats=grid->get_nLats();
  int Nalts=grid->get_nAlts();  
  int nQmag=p1d[q1d.size()];
  precision_t radius0 = planet.get_radius(0.0);
  
  int nPmag=p1d.size();  
  double pj,qi,ri = 0.9;
  double x, y, z, th;  
  int isOnGrid = -1;

  double ph_j = phiLon[idOfLon];

  for (int j_p=0; j_p<Nlats; j_p++){
    pj = p1d[j_p];
    double psSmall = 1.e-3;  
    double thetaOfPs = asin( sqrt( 1./ std::max(pj,psSmall) ) );


    for(int i_q=0; i_q<Nalts; i_q++){
      
      isOnGrid = grid->magMask(idOfLon, j_p, i_q);
//    isOnGrid = 1;

      qi = q1d[i_q];
      if(isOnGrid==1) {
        SolveDipoleEq(ri, pj, qi); 
        } else {
          ri=0.0;  
        }
      

      GetDipoleXyz(ri, ph_j, qi, pj, x, y, z, th); 

      grid->magX_scgc(idOfLon, j_p, i_q) = x;
      grid->magY_scgc(idOfLon, j_p, i_q) = y;
      grid->magZ_scgc(idOfLon, j_p, i_q) = z;
      grid->magQ_scgc(idOfLon, j_p, i_q) = q1d[i_q]; 
      grid->magQ_scgc(idOfLon, j_p, i_q) = p1d[j_p]; 

      grid->geoLon_scgc(idOfLon,j_p,i_q) = ph_j;
      grid->geoLat_scgc(idOfLon,j_p,i_q) = thetaOfPs;
      grid->geoAlt_scgc(idOfLon,j_p,i_q) = ri*radius0;

      SHOW(grid->geoAlt_scgc(idOfLon,j_p,i_q))
      SHOW(idOfLon); SHOW(j_p); SHOW(i_q);

    }
  }
}


bool Grid::init_maggrid(Planets planet, Inputs input, Report &report) {
  // filling the mag grid great again

  std::string function = "Grid::init_maggrid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  double cPI = 3.141592653589793238;

  precision_t radius0 = planet.get_radius(0.0);

  double th_startDeg= 30;
  // double th_endDeg =  180.0 - th_startDeg;

  double th_endDeg = 90 - th_endDeg;
  
  const double deg2rad = cPI/180.0;

  double th_s = deg2rad*th_startDeg;
  double th_e = deg2rad*th_endDeg;
  
  int Nlons=this->nX;
  int Nlats=this->nY;
  int Nalts=this->nZ;

  SHOW(this->nY); SHOW(nLats)


  int Nq = Nalts,
  tPower =1,
  nPhi = Nlons,

  N_mLines = Nlats; //nFootPntsPerLon = N_mLines

  std::vector<double> t01; //parameter along the B-line
  std::vector<double> th_; //theta of b-field foot for a meridional cross-section 
  std::vector<double> ph_; 

  fvec LinesPfoot(N_mLines);
  fvec LinesThetaFoot(N_mLines);
  // fvec ph_(nPhi); 
  
  arma_cube LinesQq(nPhi, N_mLines, Nq); //2do: convert to 2d = no phi
  arma_cube LinesXx(nPhi, N_mLines, Nq);
  arma_cube LinesYy(nPhi, N_mLines, Nq);
  arma_cube LinesZz(nPhi, N_mLines, Nq);
  arma_cube LinesRr(nPhi, N_mLines, Nq);
  
  t01 = linspace(0.0, 1.0, Nq);
  ph_ = linspace(0.0, cPI, nPhi); //longitude coordinate ϕϕ = [0..π]
  
  fvec qqS(N_mLines,fill::zeros);  //q, p of the foot point
  fvec ppS(N_mLines,fill::zeros); 
  fvec xS(N_mLines,fill::zeros);  //positions of the foot points
  fvec zS(N_mLines,fill::zeros);

  th_ = linspace(th_s,th_e, N_mLines);

  // SPLOT(t01,t01,Nq)
  // if(N_mLines>1){
  // th_ = linspace(th_s,th_e, N_mLines);
  // } 
  // th_ = linspace(th_s,th_e, N_mLines);
  // SHOW(th_[0]); exit(10);

  // else if (N_mLines==1) 
  // {th_[0]=th_s;}

  double qs,qe,ps,th_i,ph_i,th_j,q_i,r_i;

  

  
         

#if QGRDTYPE==1
    for (int j_ph=1; j_ph<nPhi; j_ph++){
        double ph_j = ph_[j_ph];
        
        qe = 0.;
        for (int i_line=0; i_line<N_mLines; i_line++){         
         // get the base of lines for a merid plane
         th_i = th_[i_line];
         LinesXx(j_ph,i_line,0) = xS[i_line]=sin(th_i)*cos(ph_j); 
         LinesZz(j_ph,i_line,0) = zS[i_line]=cos(th_i);          
        //  SHOW(ph_j); SHOW(th_i ); SHOW(zS[i_line]); SHOW(xS[i_line]);
        //  SHOW(LinesZz(j_ph,i_line,0))
        }
// exit(10);

// SPLOT(xS,zS,N_mLines)

        for (int l = N_mLines-1; l>=0; l--){
        
          qs = LinesZz(j_ph,l,0); //for R=1;
          ps = 1/pow(LinesXx(j_ph,l,0),2);
          

          LinesRr(j_ph,l,0) = 1.0;
          LinesQq(j_ph,l,0) = LinesZz(j_ph,l,0);                    

          LinesXx(j_ph,l,0) = sqrt(1.0/ps)*cos(ph_j);
          LinesZz(j_ph,l,0) = LinesQq(j_ph,l,0);
        
           r_i =0.9; //start from the given radius

          for (int i_q=0; i_q < Nq; i_q++){
              
              LinesQq(j_ph, l, i_q) = q_i = qs + pow(t01[i_q],tPower) * (qe-qs);                                     
              SolveDipoleEq(r_i, ps, q_i);

SHOW(q_i)
              LinesRr(j_ph, l, i_q)=r_i;                      
              LinesXx(j_ph,l,i_q) = sqrt(pow(r_i,3)/ps)*cos(ph_i);
              LinesZz(j_ph,l,i_q) = pow(r_i,3)*q_i;
            
            // th_j = acos(pow(r_i,2)*q_i);            
            // magP_scgc(j_ph, l, i_q) = 
            // magQ_scgc(j_ph, l, i_q) = LinesQq(j_ph, l, i_q);
                               
            // this->geoLon_scgc(iX,iY,iZ) = magPhi_scgc(iX,iY,iZ);             
            // this->geoLat_scgc(j_ph,iY,iZ) = theta;
            // this->geoAlt_scgc(iX,iY,iZ) = r;

              // SHOW(r_i);SHOW(i_q);
              SHOW(LinesZz(j_ph, l, i_q));
        
        }   //for i_q over a line            
        
        std::sort(LongQ.begin(), LongQ.end());
            
          
        
        
        // exit(10);

          // SPLOT(LinesXx.tube(j_ph, l),LinesZz.tube(j_ph, l),Nq)



          exit(10);        

        }//for over separate lines



// IMPORTANT! the pq-grid is chosen to be orthogonal; from this and from the fact that this should respect the "dipolar" equation,
// several no-trivial consequences for the pq-grid follow

// There is only one *ordered* grid in magnetic case: p,q,ph: as it was said, it is not *fully* arbitrary: p is calculated to respect
// 1) dipole connectiviy and 2) the specific choice of the q grid spacing; see Figure.. 
// the dependence of the spacing of p on q comes from that on the surface of r_i=1, the qs(ps) = cos(th(ps)). 
// ps = 1/sin^2(th); The q-grid should be such as to contain the qs-s; if we would choose p completely independently the q-grid cannot be
// guaranteed to be everywhere orhogonal to p-grid;

// Notice: the geoLon_scgc, geoLat_scgc, geoAlt_scgc, magX,Y,Z etc are specifically for the maggrid. The most important connected Grid
// id the magP(i,j,k) and magQ(i,j,k) grids; everything else is dependent;
// that it for a given i,j,k => p,q => r_ijk, th_ijk, ph_ijk.

// in the ph-p-q grid the "lon"=ph, "lat"=p, "alt"=q

#elif QGRDTYPE==2
   int lOfPs;
   double pOfQs,thetaOfPs;

   if (nLats != nAlts){
    SHOW(nLats); SHOW(nAlts)
    throw std::runtime_error("It should be nLats /= nAlts if QGRDTYPE==2 in  Grid::init_maggrid");
   } 

    SHOW(nLons); SHOW(nLats); SHOW(nAlts);
    double x_ij, y_ij, z_ij, th_ij;
    qe = 0.;
    
    double qL = - 0.9,
      qR = fabs(qL),
      minAbsQ = 0.05;

    // construct a temporary q-grid
    std::vector<double> q1d;        
    
    // fill q vec with values
    HalfAndHalfLinrange(q1d, Nalts,qL, qR, minAbsQ);
    
    std::vector<double> p1d(q1d.size());  
    // SHOW(q1d.size())
    // SPLOT(q1d,q1d,Nlats); exit(10);
    
    // begin populating a 3D: ph,p,q -grid
    int nGhostZones = 4;
    int nLonsMax = nLons-nGhostZones;
    
    //phi -iteration:
    for (int j_lon=0; j_lon < nLonsMax; j_lon++){ 
      
      double ph_j = ph_[j_lon];
      
      // SHOW(j_lon); SHOW(nLonsMax); SHOW(ph_j); exit(10);

      this->magLon_scgc.subcube(j_lon, 0, 0, j_lon, nLats-1, nAlts-1).fill(ph_j);
      this->geoLon_scgc.subcube(j_lon, 0, 0, j_lon, nLats-1, nAlts-1).fill(ph_j);
      
      for(int i_q=0; i_q<Nalts; i_q++){
        double qs = q1d[i_q];
        lOfPs = i_q;

        GetPsFromQs(pOfQs,qs,thetaOfPs);                
        p1d[lOfPs] = pOfQs;
        //p-grid is not independent: defined by the qs spacing:
        this->magP_scgc.tube(j_lon, lOfPs).fill(pOfQs); 
        // SHOW(p1d[lOfPs]);SHOW(q1d[i_q]);SHOW(thetaOfPs);SHOW(" ")
        
        
        ProjectQsOnAllMagLinesAnsMask(this,q1d,j_lon,i_q,lOfPs);
      
      }

      FillAllCoordsOnMagGrid(this,planet,q1d,p1d,ph_,j_lon);


// PLOT(this->geoAlt_scgc.tube(j_lon,5))


      // ScatPlotArmaArray(q1d,p1d,magMask,j_lon);

      SHOW(j_lon);                         
      
    } //phi integration
#endif        
// turn the switch on! 
  this->set_IsMagGrid(1);
  
return(true);

} 


//             SHOW(this->magQ_scgc.n_rows)
//             SHOW(this->magQ_scgc.n_cols)
//             SHOW(this->magQ_scgc.n_slices)
    
            
//           for (int l_qs = 0; l_qs < Nlats; l_qs++){  // for (int l_qs = q1d.size()-1; l_qs >= 0; l_qs--){
            
//             double qs = q1d[l_qs];
//             double ps = 1.0/(1.0-pow(qs,2)); //find corresponding p_s(q_s)

//             SHOW(ps)
//             SHOW(magP_scgc.n_rows)
//             SHOW(magP_scgc.n_cols)
//             SHOW(magP_scgc.n_slices); 
//             SHOW(nAlts);
//             SHOW(q1d.size())
           
//             this->magP_scgc.tube(j_ph, l_qs).fill(ps); 
//             //p-grid is not independent: defined by the qs spacing                        

//             double th_i = acos(qs);  

//             std::cout << "\n qs,ps = "<< qs << "; " << ps << "; " << th_i/cPI*180<<endl;
                        
//             int delt_i = (qs>0) ? -1:1;            
            
//             int i_q = l_qs;
            
//             double  qNotCrossZero = q1d[i_q]*q1d[i_q+delt_i];
            
//             SHOW(i_q); SHOW(delt_i); SHOW(qNotCrossZero)
            
//             while (qNotCrossZero > 0){
//               // for simplicity, we populate q-grid as we go along a given magnetic line from bottom towards the end          
//               // that it the direction of this passage is important if we want to have a condistent q-grid

//               q_i = q1d[i_q];
              
//               SolveDipoleEq(r_i, ps, q_i);
//               GetDipoleXyz(r_i, ph_j, q_i, ps, x_ij, y_ij, z_ij, th_ij);              

                            
//               this->geoLat_scgc(j_ph, l_qs, i_q) = this->magLat_scgc(j_ph, l_qs, i_q) = th_ij;
//               this->geoAlt_scgc(j_ph, l_qs, i_q) = this->magAlt_scgc(j_ph, l_qs, i_q) = r_i*radius0;              
                                                                                    

//               std::cout <<" i_q,q,r_i:  " << i_q << "; " << q1d[i_q] << "; "<<r_i<<endl;
//               qNotCrossZero = q1d[i_q]*q1d[i_q+delt_i]; 
// SHOW(qNotCrossZero)
//               i_q = i_q + delt_i;              

//             } 
            

//             cout<<" end of the line -------- "<<endl;                  
//             // SPLOT(LinesXx.tube(j_ph, l),LinesZz.tube(j_ph, l),Nq)

//             // copy to the grid:
//             for(int i=0; i<nAlts; i++){
//               this->magQ_scgc(j_ph, l_qs, i) = q1d[i]; // SHOW(q1d[i])
//             }


//           }
                                    

// #elif QGRDTYPE==3
// double x_ij, y_ij, z_ij, th_ij, ph_j;
// double qL = - 0.9, qR = fabs(qL), minAbsQ = 0.05;
// construct a temporary q-grid
// std::vector<double> q1d;       
// HalfAndHalfLinrange(q1d, Nalts,qL, qR, minAbsQ);    
   
// begin populating a 3D: ph,p,q -grid
// for (int j_ph=0; j_ph<nLons; j_ph++){  
//   ph_j = ph_[j_ph];

//   magLon_scgc.subcube(j_ph, 0, 0, j_ph, nLats-1, nAlts-1).fill(ph_j);
//   this->geoLon_scgc.subcube(j_ph, 0, 0, j_ph, nLats-1, nAlts-1).fill(ph_j);
//   this->magLon_scgc.subcube(j_ph, 0, 0, j_ph, nLats-1, nAlts-1).fill(ph_j);

//   for (int l=0; l<Nlats; l++){
//     if(l>0){

//     } else {
//       xS[l]=sin(th_[l])*cos(ph_j);
//       zS[l]=cos(th_[l]);
//     }

//   }

// }    

// #else
//           throw std::runtime_error("unknown QGRDTYPE in  Grid::init_mgrid");
  
  
// {
// cout<<magPhi_scgc.n_rows<<"\n"<<magPhi_scgc.n_cols<<"\n"<<magPhi_scgc.n_slices<<"\n";
// }
    



// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------




bool Grid::init_mGrid(Planets planet, Inputs input, Report &report) {
  
  std::string function = "Grid::init_mag_grid";
  static int iFunction = -1;
  bool didWork = true;
  report.enter(function, iFunction);
  
// turn the switch on! 
  this->set_IsMagGrid(1);

  // This is just an example:  
  // AD ucomment below:
  Inputs::grid_input_struct grid_input = input.get_mgrid_inputs(); 
  
  
  

  int64_t iLon, iLat, iAlt;
  
  // SHOW(grid_input.dalt); exit(10);

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
      
      // AD: does below make sense?
      // magPhi_scgc.subcube(0, iLat, iAlt, nLons-1, iLat, iAlt) = lon1d;

    }
  }
  
  // cout << "from"<< function<< "nLats=  "<< nLats<<"\n"<<endl;

  magPhi_scgc = magLon_scgc;

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  
  fvec lat1d(nLats);
  fvec lshell(nLats);

  float dlat = (grid_input.lat_max - grid_input.lat_min) / (nLats-2*nGCs);
  
  lshell(0) = 1/pow(cos(grid_input.lat_min),2.0);
  
  lshell(nLats-1) = 1/pow(cos(grid_input.lat_max),2.0);

  float dlshell = (lshell(nLats-1)-lshell(0))/nLats;
  
  for (iLat=1; iLat < nLats; iLat++){

    lshell[iLat] = lshell[iLat-1]+dlshell;
        
    //<
    lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;
    //>
  }
  
  // SHOW(lshell);
  // SHOW(lat1d); exit(10);


  //<
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      this->magLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }
  //>


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
      // cout << iLon << " " << iLat << " L= "<<Lshell<<" "<< magP_scgc(iLon,iLat,1)<<endl;
      
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
        
        //cout << "iX+ magPhi: " << iX <<" "<< iY << " " << iZ << " " << magPhi_scgc(iX,iY,iZ) << endl;
        
        // Llr: lat, lon, rad

        Llr[0] = magPhi_scgc(iX,iY,iZ);
        Llr[1] = 0.5*cPI - theta;
        Llr[2] = r;

        //if (iZ==nZ/2){
        //  cout << magP_scgc(iX,iY,iZ) <<endl;
        //  cout << Llr[0]<<" "<<Llr[1]<<" "<<Llr[2] << endl;

        magX_scgc(iX,iY,iZ)=Xyz[0];
        magY_scgc(iX,iY,iZ)=Xyz[1];
        magZ_scgc(iX,iY,iZ)=Xyz[2];


        precision_t radius0 = planet.get_radius(0.0);
        // SHOW(radius0)
        // exit(10);

        //<

        // if (iX == 5 & iY == 5)
        //   cout << "lon, lat, alt: " << magPhi_scgc(iX,iY,iZ) << " "
        //   << theta << " " << r 
        //  << " " << magP_scgc(iX,iY,iZ) 
        //  << " " << magQ_scgc(iX,iY,iZ) << "\n";

        this->geoLon_scgc(iX,iY,iZ) = magPhi_scgc(iX,iY,iZ);
        this->geoLat_scgc(iX,iY,iZ) = theta;
        this->geoAlt_scgc(iX,iY,iZ) = r;

        this->geoAlt_scgc(iX,iY,iZ) *= radius0;


        // - grid_input.alt_min ?

        //>
        
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
  
  fill_grid_radius(planet);
  
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
  return didWork;
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
    
    // SHOW(DeltaTheta)
    // SHOW(Tolerance); abort();

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
      // cout << "!!! For L = " << Lshell << endl;
      // cout << "!!! qN = " << qN << endl;
      // cout << "!!! ThetaMid = " << ThetaMid*cRtoD << endl;
    }else{
      qS = pow(rDipoleMid,-2.0)*cos(ThetaMid);
      // cout << "!!! qS = " << qS << endl;
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
       gridfile << Xyz[0] <<" " 
                << Xyz[1] <<" "
                << Xyz[2] <<" "
                << Lshell <<  endl;

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
  r = 100.0;

  Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
  dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;
  
  // cout<< "p,q="<<p<<" "<<q << endl;
  // cout<< Func<<" "<<dFunc;
  // cout<<endl;

  int itr=0;
  int maxItr=100;

  // apply NR iterations to get 
  while( abs(Func/dFunc) > Tolerance) { 
    try {
      Func= pow(q,2.0) * pow(r,4.0) + 1.0/p*r-1;
      dFunc= 4.0*pow(q,2.0) * pow(r,3.0) + 1.0/p;
      
      // in NR method r(i+1)=r(i)-f/f' for each iteration
      
      r = r - Func/dFunc;

      if (++itr > maxItr){ throw(itr);}
    }
    catch (int itr){
        cout<<"WARN: exceeded max #iterations.. exiting ";
        exit(10);
    }
    // cout << r << " " << Func << " "<< dFunc << endl;
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

arma_vec get_r3_spacing(precision_t lat, precision_t rMin, 
                        precision_t rMax, int64_t nPts, int64_t nGcs) {

  precision_t rMaxReal = rMax;

  precision_t lShell = get_lshell(lat, rMin);
  if (lShell < rMaxReal) {
    rMaxReal = lShell;
    std::cout << "Limiting rMaxReal from " << rMax << " to " << rMaxReal << "\n";
  }
  precision_t rMin3 = pow(rMin, 1.0/3.0);
  precision_t rMax3 = pow(rMaxReal, 1.0/3.0);
  precision_t dr3 = (rMax3 - rMin3) / (nPts-nGcs*2);
  arma_vec r(nPts);
  for (int64_t iPt = 0; iPt < nPts; iPt++) {
    r(iPt) = pow(rMin3 + dr3 * (iPt - nGcs), 3);
  }
  return r;
}

// ----------------------------------------------------------------------
// Initialize the geographic grid.  At the moment, this is a simple
// Lon/Lat/Alt grid.  The grid structure is general enough that each
// of the lon, lat, and alt can be a function of the other variables.
// ----------------------------------------------------------------------

void Grid::init_dipole_grid(Quadtree quadtree, Planets planet, Inputs input, Report &report) {
  
  std::string function = "Grid::init_dipole_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);
  
// turn the switch on! 
  IsMagGrid = true;

  int64_t iLon, iLat, iAlt;

  // This is just an example:  
  report.print(3, "Getting mgrid_inputs inputs in dipole grid");

  Inputs::grid_input_struct grid_input = input.get_mgrid_inputs();

  report.print(3, "Setting inputs in dipole grid");
  precision_t min_apex = grid_input.min_apex;
  precision_t min_alt = grid_input.alt_min;
  precision_t max_alt = grid_input.alt_max;
  precision_t planetRadius = planet.get_radius(0.0);
  precision_t min_lshell = (min_apex + planetRadius)/planetRadius;
  precision_t min_r = (min_alt + planetRadius)/planetRadius;
  precision_t max_r = (max_alt + planetRadius)/planetRadius;
  precision_t min_lat = get_lat_from_r_and_lshell(min_r, min_lshell);
  precision_t stretch = (cPI/2 - min_lat) / (cPI/2);
  report.print(3, "Done setting inputs in dipole grid");

  // Get some coordinates and sizes in normalized coordinates:
  arma_vec lower_left_norm = quadtree.get_vect("LL");
  arma_vec size_right_norm = quadtree.get_vect("SR");
  arma_vec size_up_norm = quadtree.get_vect("SU");

  precision_t dlon = size_right_norm(0) * cPI / (nLons - 2 * nGCs);
  precision_t lon0 = lower_left_norm(0) * cPI;
  arma_vec lon1d(nLons);

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLon = 0; iLon < nLons; iLon++)
    lon1d(iLon) = lon0 + (iLon - nGCs + 0.5) * dlon;

  for (iLat = 0; iLat < nLats; iLat++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLon_scgc.subcube(0, iLat, iAlt, nLons - 1, iLat, iAlt) = lon1d;
  }

  geoLon_scgc = magLon_scgc;

  precision_t dlat = size_up_norm(1) * cPI / (nLats - 2 * nGCs);
  precision_t lat0 = lower_left_norm(1) * cPI;

  arma_vec lat1d(nLats);

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  for (iLat = 0; iLat < nLats; iLat++) {
    lat1d(iLat) = lat0 + (iLat - nGCs + 0.5) * dlat;
    std::cout << "Original : " <<  lat1d(iLat) << " ";
    if (lat1d(iLat) >= 0) {
      lat1d(iLat) = min_lat + lat1d(iLat) * stretch;
    } else {
      lat1d(iLat) = -min_lat + lat1d(iLat) * stretch;
    }
    std::cout << "Final : " <<  lat1d(iLat) << "\n";

  }
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      magLat_scgc.subcube(iLon, 0, iAlt, iLon, nLats - 1, iAlt) = lat1d;
  }

  arma_vec rNorm1d, lat1dAlong;
  arma_cube r3d(nLons, nLats, nAlts);

  precision_t lShell;

  for (iLat = 0; iLat < nLats; iLat++) {
    lat0 = lat1d(iLat);
    if (lat0 > cPI/2) lat0 = cPI - lat0;
    if (lat0 < -cPI/2) lat0 = -cPI - lat0;
    lShell = get_lshell(lat0, min_r);
    std::cout << "lShell : " << lat0 * cRtoD << " " << lShell << "\n";
    std::cout << "min_r : " << min_r << " " << max_r << " " << nAlts << " " << nGCs << "\n";
    rNorm1d = get_r3_spacing(lat0, min_r, max_r, nAlts, nGCs);
    lat1dAlong = get_lat_from_r_and_lshell(rNorm1d, lShell);
    if (lat0 < 0)
      lat1dAlong = -1.0 * lat1dAlong;
    for (iLon = 0; iLon < nLons; iLon++) {
      r3d.tube(iLon, iLat) = rNorm1d * planetRadius;
      magLat_scgc.tube(iLon, iLat) = lat1dAlong;
    }
  }

  geoLat_scgc = magLat_scgc;
  magAlt_scgc = r3d - planetRadius;

  std::vector<arma_cube> llr, xyz, xyzRot1, xyzRot2;
  llr.push_back(magLon_scgc);
  llr.push_back(magLat_scgc);
  llr.push_back(r3d);
  xyz = transform_llr_to_xyz_3d(llr);

  precision_t magnetic_pole_rotation = planet.get_dipole_rotation();
  precision_t magnetic_pole_tilt = planet.get_dipole_tilt();

  // Reverse our dipole rotations:
  xyzRot1 = rotate_around_y_3d(xyz, magnetic_pole_tilt);
  xyzRot2 = rotate_around_z_3d(xyzRot1, magnetic_pole_rotation);

  // transform back to lon, lat, radius:
  // from transform.cpp

  // AD uncomment below:AD 
  llr = transform_xyz_to_llr_3d(xyzRot2);

  

  geoLon_scgc = llr[0];
  geoLat_scgc = llr[1];
  geoAlt_scgc = llr[2] - planetRadius;

  std::cout << "magLon : " << magLon_scgc(12,10,5) * cRtoD << " "
              << magLat_scgc(12,10,5) * cRtoD << " "
              << magAlt_scgc(12,10,5) / 1000.0 << "\n";

  std::cout << "geoLon : " << geoLon_scgc(12,10,5) * cRtoD << " "
              << geoLat_scgc(12,10,5) * cRtoD << " "
              << geoAlt_scgc(12,10,5) / 1000.0 << "\n";

  report.exit(function);
  return;

}

