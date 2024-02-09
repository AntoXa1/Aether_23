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
      
  int iErr = 0;
  bool didWork = true;
  Times time;
  Report report;

  // Define the function and report:
  std::string function = "test_mgrid";
  static int iFunction = -1;
  report.enter(function, iFunction);

  cout<<"entering : \n"<< function <<endl;
  // std::filesystem :: path cwd = std::filesystem::current_path();

  // string :: path_to_rundir = ""
    // Create inputs (reading the input file):
    input = Inputs(time);
    

    // Initialize the EUV system:
    // Euv euv(input, report);
  // cout<<"passed .. \n";
    
  Quadtree quadtree;
  if (!quadtree.is_ok()) throw std::string("quadtree initialization failed!");    
  SHOW("mgrid: did quadtree \n")

  didWork = init_parallel(quadtree);
  if (!didWork) throw std::string("init_parallel(quadtree) failed!");
// SHOW("mgrid: did until this place \n")


  // Initialize the planet:
  Planets planet;   
  if (!planet.is_ok()) throw std::string("planet initialization failed!"); 

  
  // Initialize the indices (and read the files):
  Indices indices;
  didWork = read_and_store_indices(indices);    

  // Initialize Geographic grid:
  Grid gGrid(input.get_nLonsGeo(),
        input.get_nLatsGeo(),
        input.get_nAltsGeo(),
        nGeoGhosts);
  didWork = gGrid.init_geo_grid(quadtree, planet);
  // MPI_Barrier(aether_comm);
  if (!didWork) throw std::string("init_geo_grid failed!");  

  // gGrid.fill_grid(planet);

  gGrid.report_grid_boundaries();

  SHOW(gGrid.geoX_scgc(1,1,1)) 
  SHOW(gGrid.geoAlt_scgc(1,1,1))

  // Initialize Neutrals on geographic grid:
  Neutrals neutrals(gGrid, planet, time, indices);

  // Initialize Ions on geographic grid:
  Ions ions(gGrid, planet);

  // Intitialize dipole magnetic grid
  Grid mGrid(input.get_nLonsMag(),
        input.get_nLatsMag(),
        input.get_nAltsMag(), nMagGhosts);

  didWork = mGrid.init_maggrid(planet, input, report);
  // MPI_Barrier(aether_comm);
  if (!didWork)throw std::string("init_geo_grid failed!");
  cout << "initMagneticGrid done"<< endl;

  // mGrid.fill_grid(planet);
  mGrid.fill_grid_radius(planet);

  int iAlt=0; int iLon=0, iLat =0;
  for (int i=1; i<mGrid.get_nAlts(); i++){
    SHOW(mGrid.radius_scgc(iLon,iLat, i));
    SHOW(mGrid.geoAlt_scgc(iLon,iLat, i));    
    SHOW(i);
  }

  // Initialize Neutrals on dipole grid:
    Neutrals mNeutrals(mGrid, planet, time, indices);  
    cout<<" Initialize Neutrals on dipole grid: done .."<<endl;

  // Initialize Ions on m-geographic grid:
    Ions mIons(mGrid, planet);
    cout<<"Initialize Ions on m-geographic grid: done .."<<endl;

  // Once EUV, neutrals, and ions have been defined, pair cross sections
  // Initialize the EUV system:    
  Euv euv;
  cout<<"Initialize the EUV system: done .."<<endl;  
  if (!euv.is_ok())
        throw std::string("EUV initialization failed!");    
  euv.pair_euv(mNeutrals, mIons);

  // Initialize Chemical scheme (including reading file):
  Chemistry m_chemistry(mNeutrals, mIons);

    
      
  // for (int i=1; i < gGrid.get_nAlts(); i++)
  // SHOW( gGrid.radius_scgc(iLon,iLat, i) );          
  // SHOW(input.get_nAltsMag())    
  // int  Nq =10;
  // int tPower =1;
  // DipoleLine mLine(Nq,tPower);

    // SHOW(input.get_nLonsMag()); exit(10);
    //   int iAlt=0; int iLon=0, iLat =0;
    //   for (int i=1; i<mGrid.get_nAlts(); i++){
    //     SHOW(mGrid.radius_scgc(iLon,iLat, i) ); 
    //     SHOW(mGrid.geoAlt_scgc(iLon,iLat, i) );    
    //     SHOW(i);
    // WriteScatteredGridToFile(mGrid);


  // Read in the collision frequencies and other diffusion coefficients:
  read_collision_file(mNeutrals, mIons);
  cout<<"read_collision_file done"<<"\n";


  // Initialize ion temperatures from neutral temperature
  mIons.init_ion_temperature(mNeutrals, mGrid);
  cout<<"init_ion_temperature done"<<"\n";
  // system("pause");


  double dt_couple = time.get_end() - time.get_current();
  time.increment_intermediate(dt_couple);

  // AJR - added this stuff for one time-step:
    mGrid.calc_sza(planet, time);
    
    mNeutrals.calc_mass_density();
    
    mNeutrals.calc_specific_heat();

    // time.calc_dt();
    // cout<<"calc_dt done"<<"\n";

    iErr = calc_euv(planet,
                    mGrid,
                    time,
                    euv,
                    mNeutrals,
                    mIons,
                    indices);
    

  

  //cout<< " No calc_chemistry "<<endl;
  // m_chemistry.calc_chemistry(mNeutrals, mIons, time, mGrid); 

  // advance chem
  while (time.get_current() < time.get_end()) {
    // didWork = advance(planet,
    // 		  gGrid,
    // 		  time,
    // 		  euv,
    // 		  neutrals,
    // 		  ions,
    // 		  mchemistry,
    // 		  electrodynamics,
    // 		  indices,
    // 		  logfile);

    if (!didWork) throw std::string("Error in advance!");  
        

  }
  //iErr = output(neutrals,
  //	      ions,
  //	      gGrid,
  //	      time,
  //	      planet,
  //	      input,
  //	      report);

  // iErr = output(mNeutrals,
  //  	      mIons,
  //  	      mGrid,
  //  	      time,
  //  	      planet);

  



  report.exit(function);
  report.times();

  cout<< "init mGrid.mag_grid is ok"<<endl;
  return iErr;
}
