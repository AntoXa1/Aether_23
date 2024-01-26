// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
#include <iostream>
#include <filesystem>

#include "../include/aether.h"
#include <type_traits>

// #include <termios.h>
// #include <unistd.h>
// #define PAUS() \
//     (std::cout << "Press any key to continue...", \
//      std::cout.flush(), \
//      getch())

// (cd ../ && make clean && make -j)

// -----------------------------------------------------------------------------
// Main file for the Aether model.  This is needed when Aether is not used
// as a library in another code, such as the SWMF.
// -----------------------------------------------------------------------------

class DipoleLine{
public:
  int numElem;
  
  DipoleLine(int numElemIn,int tPow); 

// friend class Grid;

  std::vector<double> xx;
  std::vector<double> zz;
  std::vector<double> qq;
  std::vector<double> rr;
  std::vector<double> tt;
  double tPower;

};

DipoleLine::DipoleLine(int numElemIn,int tPow):tPower(tPow),numElem(numElemIn), 
                        xx(numElemIn),zz(numElemIn),qq(numElemIn),
                        rr(numElemIn),tt(numElemIn){
          tt = linspace(0.0, 1.0, numElem);
};


void WriteScatteredGridToFile(Grid grid){
//      WriteScatterPntsDataToFile
//      assumes pointwise data: x_col y_col f(x,y)_col

        char cwd[100];        
        const std::string name = "scatgrid.dat";

        std::string  fdir(getcwd(cwd, sizeof(cwd)));
        std::string fileName = fdir + "/" + "_" + name+"_dbg.dat";   


        std::cout << "Current working directory: " << fileName << std::endl;
        std::ofstream output_file(fileName);
        
        int Nlons=grid.get_nLons();
        int Nlats=grid.get_nLats();
        int Nalts=grid.get_nAlts();
        
        for (int j_ph=0; j_ph<Nlons; j_ph++){
          for (int i_lat = 0; i_lat < Nlats; i_lat++){
            for (int k_q = 0; k_q < Nlats; k_q++){

        // int n_rows=xS.n_rows;
        // for (int i = 0; i < n_rows; i++) {  
        //   output_file << xS(i)<< " " << zS(i) << std::endl;


            }// k-q-loop
          }// i_lat-loop         
        } // j-loop 
        output_file.close();
        
}

int main() {

// asdf


SHOW("test mgrid START")
  
  int iErr = 0;
  bool DidWork = true;
  Times time;
  Report report;

  // Define the function and report:
  std::string function = "test_mgrid";
  static int iFunction = -1;
  report.enter(function, iFunction);


cout<<"entering : "<< function <<endl;
// std::filesystem :: path cwd = std::filesystem::current_path();


// string :: path_to_rundir = ""
  // Create inputs (reading the input file):
  Inputs input(time, report);

  // Initialize the EUV system:
  // Euv euv(input, report);
// cout<<"passed .. \n";
  

  Quadtree quadtree(input, report);

  DidWork = init_parallel(input, quadtree, report);

  // Initialize the planet:
  Planets planet(input, report);

  // Initialize the indices (and read the files):
  Indices indices(input);
  iErr = read_and_store_indices(indices, input, report);

  // Initialize Geographic grid:
  Grid gGrid(input.get_nLonsGeo(),
       input.get_nLatsGeo(),
       input.get_nAltsGeo(), nGeoGhosts);

  if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");
  gGrid.init_geo_grid(quadtree, planet, input, report);
  
  gGrid.fill_grid(planet, report);

  gGrid.report_grid_boundaries();

  SHOW(gGrid.geoX_scgc(1,1,1)) 
  SHOW(gGrid.geoAlt_scgc(1,1,1))



  {
    int iAlt=0; int iLon=0, iLat =0;
    
    //  for (int i=1; i < gGrid.get_nAlts(); i++)
    //      SHOW( gGrid.radius_scgc(iLon,iLat, i) );          
  }

  Neutrals neutrals(gGrid, planet, time, indices, input, report);


  //------------Initialize Magnetic grid------------------
  // SHOW(input.get_nAltsMag())
  
  int  Nq =10;
  int tPower =1;

  DipoleLine mLine(Nq,tPower);

// SHOW(input.get_nLonsMag()); exit(10);

  Grid mGrid(input.get_nLonsMag(),
       input.get_nLatsMag(),
       input.get_nAltsMag(), nMagGhosts);

  // Anton's code:
  
  mGrid.initMagneticGrid(planet, input, report);
  cout << "initMagneticGrid done"<< endl;
  

  // Aaron's code:
  // mGrid.init_dipole_grid(quadtree, planet, input, report);

  // mGrid.init_mag_grid(planet, input, report);

  mGrid.fill_grid(planet, report);
  mGrid.fill_grid_radius(planet, report);

    int iAlt=0; int iLon=0, iLat =0;
    for (int i=1; i<mGrid.get_nAlts(); i++){
      SHOW(mGrid.radius_scgc(iLon,iLat, i) ); 
      SHOW(mGrid.geoAlt_scgc(iLon,iLat, i) );    
      SHOW(i);
    }
    
  

// WriteScatteredGridToFile(mGrid);

// cout<< "radius_scgc = " << mGrid.radius_scgc.subcube(iLon,iLat, 0, iLon,iLat,10 ) <<" "<<iAlt <<endl; 
// SHOW(gGrid.dalt_center_scgc.slice(iAlt))    
  
//-----------------------------------------------------  

  // iterate p,q; convert to r,theta,phi; 
  // p,q, is uniform, while rThPhi is non-uniform
  // feed and initialize m_neutrals with rThPhi

  // Initialize Neutrals on dipole grid:
  Neutrals mNeutrals(mGrid, planet, time, indices, input, report);
  cout<<" Initialize Neutrals on dipole grid: done .."<<endl;



// Initialize Ions on m-geographic grid:
  Ions mIons(mGrid, planet, input, report);
  cout<<"Initialize Ions on m-geographic grid: done .."<<endl;


// Once EUV, neutrals, and ions have been defined, pair cross sections

// Initialize the EUV system:
    
Euv euv(input, report);
cout<<"Initialize the EUV system: done .."<<endl;


  
if (!euv.is_ok())
      throw std::string("EUV initialization failed!");
    
euv.pair_euv(mNeutrals, mIons, report);



// Initialize Chemical scheme (including reading file):
Chemistry m_chemistry(mNeutrals, mIons, input, report);




// Read in the collision frequencies and other diffusion coefficients:
    read_collision_file(mNeutrals, mIons, input, report);
cout<<"read_collision_file done"<<"\n";


    // Initialize ion temperatures from neutral temperature
    mIons.init_ion_temperature(mNeutrals, mGrid, report);

cout<<"init_ion_temperature done"<<"\n";
// system("pause");


double dt_couple = time.get_end() - time.get_current();
time.increment_intermediate(dt_couple);

// AJR - added this stuff for one time-step:
  mGrid.calc_sza(planet, time, report);
  
  mNeutrals.calc_mass_density(report);
  
  mNeutrals.calc_specific_heat(report);

  time.calc_dt();

cout<<"calc_dt done"<<"\n";

  iErr = calc_euv(planet,
                  mGrid,
                  time,
                  euv,
                  mNeutrals,
                  mIons,
                  indices,
                  input,
                  report);
  

 

//cout<< " No calc_chemistry "<<endl;

m_chemistry.calc_chemistry(mNeutrals, mIons, time, mGrid, report); 

// advance chem

//iErr = output(neutrals,
//	      ions,
//	      gGrid,
//	      time,
//	      planet,
//	      input,
//	      report);
while (time.get_current() < time.get_end()) {

}

iErr = output(mNeutrals,
 	      mIons,
 	      mGrid,
 	      time,
 	      planet,
 	      input,
 	      report);
 
cout<<"----------" << "init mGrid.mag_grid is ok"<<endl;

// exit(10);
  report.exit(function);
  report.times();

  return iErr;
}
