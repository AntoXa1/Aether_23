// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md
#include <iostream>
#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Main file for the Aether model.  This is needed when Aether is not used
// as a library in another code, such as the SWMF.
// -----------------------------------------------------------------------------


int main() {

  int iErr = 0;
  bool DidWork = true;
  Times time;
  Report report;

  // Define the function and report:
  std::string function = "test_mgrid";
  static int iFunction = -1;
  report.enter(function, iFunction);


cout<<"entering : "<< function <<endl;

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



//   gGrid.init_geo_grid(planet, input, report);

    
  if (!quadtree.is_ok())
      throw std::string("quadtree initialization failed!");
  gGrid.init_geo_grid(quadtree, planet, input, report);
  
  gGrid.fill_grid(planet, report);

  // Initialize Magnetic grid:
  cout << "input.get_nAltsMag() " <<input.get_nAltsMag()<<endl;

  Grid mGrid(input.get_nLonsMag(),
       input.get_nLatsMag(),
       input.get_nAltsMag(), nMagGhosts);


  mGrid.init_mag_grid(planet, input, report);
  
  // iterate p,q; convert to r,theta,phi; 
  // p,q, is uniform, while rThPhi is non-uniform
  // feed and initialize m_neutrals with rThPhi

  // Initialize Neutrals on dipole grid:
  Neutrals m_neutrals(mGrid, planet, time, indices, input, report);

exit(10);
// Initialize Ions on m-geographic grid:
    Ions m_ions(mGrid, planet, input, report);

// Once EUV, neutrals, and ions have been defined, pair cross sections

// Initialize the EUV system:
    
Euv euv(input, report);
    
if (!euv.is_ok())
      throw std::string("EUV initialization failed!");
    
euv.pair_euv(m_neutrals, m_ions, report);

    // Initialize Chemical scheme (including reading file):
Chemistry m_chemistry(m_neutrals, m_ions, input, report);

// Read in the collision frequencies and other diffusion coefficients:
    read_collision_file(m_neutrals, m_ions, input, report);

    // Initialize ion temperatures from neutral temperature
    m_ions.init_ion_temperature(m_neutrals, mGrid, report);

double dt_couple = time.get_end() - time.get_current();
time.increment_intermediate(dt_couple);

m_chemistry.calc_chemistry(m_neutrals, m_ions, time, mGrid, report); 
// advance chem

cout<<"----------" << "init mGrid.mag_grid is ok"<<endl;


  report.exit(function);
  report.times();

  return iErr;
}
