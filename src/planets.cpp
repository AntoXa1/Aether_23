// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "../include/aether.h"

// -----------------------------------------------------------------------------
// Constructor (initiaze the class):
// -----------------------------------------------------------------------------

Planets::Planets(Inputs input, Report report) {
  IsOk = read_file(input, report);

  if (IsOk)
    IsOk = set_planet(input, report);

  if (IsOk)
    IsOk = read_planet_specific_file(input, report);

  IsOk = sync_across_all_procs(IsOk);
}

// -----------------------------------------------------------------------------
// Get the offset between the longitude and local time
// -----------------------------------------------------------------------------

precision_t Planets::get_longitude_offset(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.longitude_offset;
}

// -----------------------------------------------------------------------------
// Get the sin of the planet's declination angle
// -----------------------------------------------------------------------------

precision_t Planets::get_sin_dec(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.sin_dec;
}

// -----------------------------------------------------------------------------
// Get the cos of the planet's declination angle
// -----------------------------------------------------------------------------

precision_t Planets::get_cos_dec(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.cos_dec;
}

// -----------------------------------------------------------------------------
// Get the radius of the planet as a function of latitude (meters)
// -----------------------------------------------------------------------------

precision_t Planets::get_radius(precision_t latitude) {
  // Should modify this to allow an oblate spheriod, but not now.
  return planet.radius;
}

// -----------------------------------------------------------------------------
// Get the rotation of the dipole (in radians, rotation in longitude)
// -----------------------------------------------------------------------------

precision_t Planets::get_dipole_rotation() {
  return planet.dipole_rotation;
}

// -----------------------------------------------------------------------------
// Get the tilt of the dipole away from the rotation axis of the
// planet (radians)
// -----------------------------------------------------------------------------

precision_t Planets::get_dipole_tilt() {
  return planet.dipole_tilt;
}

// -----------------------------------------------------------------------------
// Get the strength of the dipole
// -----------------------------------------------------------------------------

precision_t Planets::get_dipole_strength() {
  return planet.dipole_strength;
}

// -----------------------------------------------------------------------------
// Get the location of the dipole center
// -----------------------------------------------------------------------------

std::vector<float> Planets::get_dipole_center() {
  return planet.dipole_center;
}

// -----------------------------------------------------------------------------
// Get the mu of the planet (mass * gravitational constant)
// -----------------------------------------------------------------------------

precision_t Planets::get_mu() {
  return planet.mu;
}

// -----------------------------------------------------------------------------
// Get the distance from the star to the planet (meters)
// -----------------------------------------------------------------------------

precision_t Planets::get_star_to_planet_dist(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.star_planet_distance;
}

// -----------------------------------------------------------------------------
// Get planetary orbit angle (radians)
// -----------------------------------------------------------------------------

precision_t Planets::get_orbit_angle(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.orbit_angle;
}

// -----------------------------------------------------------------------------
// Get planetary declination (radians)
// -----------------------------------------------------------------------------

precision_t Planets::get_declination(Times time) {
  int iErr = update(time);

  if (iErr > 0)
    std::cout << "Error in setting time!" << '\n';

  return planet.declination;;
}

// -----------------------------------------------------------------------------
// Update lots of information about the planet, such as the location
// declination, and offset to compute local time from longitude
// -----------------------------------------------------------------------------

int Planets::update(Times time) {

  int iErr = 0;

  // Update planetary stuff once a minute:
  if (time.get_current() - planet.update_time > 60.0) {

    planet.update_time = time.get_current();

    precision_t sma = planet.semimajoraxis +
                      planet.rates_semimajoraxis * time.get_orbittime();
    precision_t ecc = planet.eccentricity +
                      planet.rates_eccentricity * time.get_orbittime();
    precision_t inc = planet.inclination +
                      planet.rates_inclination * time.get_orbittime();
    precision_t meanlon = planet.meanlongitude +
                          planet.rates_meanlongitude * time.get_orbittime();
    precision_t node_long = planet.nodelongitude +
                            planet.rates_nodelongitude * time.get_orbittime();
    precision_t perilon = planet.perihelionlongitude +
                          planet.rates_perihelionlongitude * time.get_orbittime();

    // Compute argument of perihelion and mean anomaly
    precision_t arg_peri = perilon - node_long;

    // computation of M for Jupiter and out is supposed to be modified
    // by additional terms for the time interval 3000BC -
    // 30000AD... This probably doesn't matter.
    precision_t meananomaly = fmod(meanlon - perilon + 360.0, 360.0);

    if (meananomaly > 180.)
      meananomaly = meananomaly - 360.0;

    // Need to solve Kepler's equation by iterating
    precision_t d_ecc_anomaly = 10000.0;
    precision_t tol = 1.0e-6;
    precision_t ecc_deg = ecc * cRtoD;
    precision_t ecc_anomaly = meananomaly + ecc_deg * sin(meananomaly * cDtoR);
    int i = 0;
    precision_t dm = 0.0;

    while (abs(d_ecc_anomaly) > tol && i < 100) {
      dm = meananomaly - (ecc_anomaly - ecc_deg * sin(ecc_anomaly * cDtoR));
      d_ecc_anomaly = dm / (1 - ecc * cos(ecc_anomaly * cDtoR));
      ecc_anomaly = ecc_anomaly + d_ecc_anomaly;
      i++;
    }

    // Get heliocentric coordinates, TrueAnomaly and sunplanetdistance
    precision_t x_heliocentric = sma * (cos(ecc_anomaly * cDtoR) - ecc);
    precision_t y_heliocentric = sma * sqrt(1 - ecc * ecc) * sin(
                                   ecc_anomaly * cDtoR);

    precision_t true_anomaly = atan2(y_heliocentric, x_heliocentric) * cRtoD;

    planet.star_planet_distance = sqrt(x_heliocentric * x_heliocentric +
                                       y_heliocentric * y_heliocentric);

    // convert to J2000 coordinates with x-axis aligned with vernal equinox so
    // we can get solar longitude in the coorect system.  We don't need z.
    precision_t x_ecl =
      x_heliocentric *
      (cos(arg_peri * cDtoR) * cos(node_long * cDtoR) -
       sin(arg_peri * cDtoR) * sin(node_long * cDtoR) * cos(inc * cDtoR)) +
      y_heliocentric *
      (-sin(arg_peri * cDtoR) * cos(node_long * cDtoR) -
       cos(arg_peri * cDtoR) * sin(node_long * cDtoR) * cos(inc * cDtoR));

    precision_t y_ecl =
      x_heliocentric *
      (cos(arg_peri * cDtoR) * sin(node_long * cDtoR) +
       sin(arg_peri * cDtoR) * cos(node_long * cDtoR) * cos(inc * cDtoR)) +
      y_heliocentric *
      (-sin(arg_peri * cDtoR) * sin(node_long * cDtoR) +
       cos(arg_peri * cDtoR) * cos(node_long * cDtoR) * cos(inc * cDtoR));

    // Calculate orbit angle, aka Ls. In this CS, need the angle from -x axis.
    planet.orbit_angle = atan(y_ecl / x_ecl);

    if (x_ecl > 0)
      planet.orbit_angle = planet.orbit_angle + cPI;

    if (x_ecl < 0 && y_ecl > 0)
      planet.orbit_angle = planet.orbit_angle + cTWOPI;

    planet.declination = atan(tan(planet.planet_tilt) * sin(planet.orbit_angle));

    planet.sin_dec = sin(planet.declination);
    planet.cos_dec = cos(planet.declination);

    // This is take from :
    //   https://www.aa.quae.nl/en/reken/zonpositie.html (section 7)
    //   Basically, it identifies the longitude that is pointed along
    //   the sun-planet line at 00 UT on Jan 1 2000.  The rotation
    //   that is used on the webpage is to calculate the sidereal
    //   time, though, which we want the local time.  So, we use the
    //   length of a day instead of the rotation rate.  Ok, this is
    //   idea is correct, but the angles reported on this website at
    //   not correct.  They may be providing the angle for RAAN or
    //   something like that. It is not defining the longitude for
    //   midnight or noon. Ugh. (Fixed Earth, but no other planet.)
    double nEarthDaysSince2000 = (time.get_julian_day() - cJULIAN2000);
    double deg_per_day = 360.0 / (planet.length_of_day * cStoD);
    double deg_since_2000 = deg_per_day * nEarthDaysSince2000;
    double rotations = deg_since_2000 / 360.0;
    double left_over = (rotations - static_cast<int>(rotations)) * 360.0;
    // put into radians, so it is consistent with the rest of the code:
    planet.longitude_offset =
      fmod(planet.longitude_jb2000 + left_over, 360.0) * cDtoR;
  }

  return iErr;
}

// -----------------------------------------------------------------------------
// Set the specific planet that we are interested in simulating and
// move the appropriate data over to the planet structure
// -----------------------------------------------------------------------------

bool Planets::set_planet(Inputs input, Report report) {

  bool DidWork = true;
  int IsFound = 0;

  int iSize = planets.size();

  for (int i = 0; i < iSize; i++) {
    if (planets[i].name == input.get_planet()) {
      IsFound = 1;
      planet.name = planets[i].name;

      if (report.test_verbose(2))
        std::cout << "Planet set to : " << planet.name << "\n";

      planet.semimajoraxis = planets[i].semimajoraxis;
      planet.eccentricity = planets[i].eccentricity;
      planet.inclination = planets[i].inclination;
      planet.meanlongitude = planets[i].meanlongitude;
      planet.perihelionlongitude = planets[i].perihelionlongitude;
      planet.nodelongitude = planets[i].nodelongitude;
      planet.rates_semimajoraxis = planets[i].rates_semimajoraxis;
      planet.rates_eccentricity = planets[i].rates_eccentricity;
      planet.rates_inclination = planets[i].rates_inclination;
      planet.rates_meanlongitude = planets[i].rates_meanlongitude;
      planet.rates_perihelionlongitude = planets[i].rates_perihelionlongitude;
      planet.rates_nodelongitude = planets[i].rates_nodelongitude;
      planet.planet_tilt = planets[i].planet_tilt;

      double mu = cMSOL * cG;
      double a = planet.semimajoraxis * cAUtoM;
      double loy = cTWOPI * sqrt(a * a * a / mu);

      if (abs(loy - planets[i].length_of_year) / loy > 0.01) {
        std::cout << "Hmmm.... For some reason, the calculated length ";
        std::cout << "of year is different\n";
        std::cout << "Read in : "
                  << planet.length_of_year * cStoD
                  << " (days) \n";
        std::cout << "Calculated from SMA : "
                  << loy * cStoD << " (days)\n";
        std::cout << "Trusting SMA version!\n";
        DidWork = false;
      }

      planet.length_of_year = loy;
      planet.length_of_day = planets[i].length_of_day;
      planet.longitude_jb2000 = planets[i].longitude_jb2000;
      precision_t rotrate = cTWOPI /
                            planet.length_of_day * (1.0 + 1.0 / (loy * cStoD));
      planet.omega = rotrate;  // frequency (rad/s)
      planet.rotation_period = cTWOPI / rotrate; // (seconds)

      planet.mass = planets[i].mass;
      planet.mu = planets[i].mass * cG;
      planet.equator_radius = planets[i].equator_radius * 1000.0;  // km -> m
      planet.polar_radius = planets[i].polar_radius * 1000.0;  // km -> m
      // Looking at Earth and Saturn, it seems like the Volumetric
      // mean radius is at roughly 47 deg (cos(47)=0.68) Obviously an
      // approximation...
      planet.radius = 0.68 * planet.equator_radius + 0.32 * planet.polar_radius;

      if (report.test_verbose(2))
        std::cout << "Planet Radius set to : "
                  << planet.radius / 1000.0 << " (km)\n";

      planet.dipole_strength = planets[i].dipole_strength;
      planet.dipole_rotation = planets[i].dipole_rotation;
      planet.dipole_tilt = planets[i].dipole_tilt;

      for (int j = 0; j < 3; j++)
        planet.dipole_center[j] = planets[i].dipole_center[j];

      planet.update_time = -1e32;

      break;
    }
  }

  if (!IsFound) {
    std::cout << "Can't file planet " << input.get_planet()
              << " in planet file information!\n";
    DidWork = false;
  }

  return DidWork;
}

// -----------------------------------------------------------------------------
// Read the planetary characteristics file that contains all of the
// information for all of the different planets.
// -----------------------------------------------------------------------------

bool Planets::read_file(Inputs input, Report report) {

  planet_chars tmp;
  std::string line, col;
  std::ifstream myFile;
  bool DidWork = true;

  report.print(1, "Reading planetary file : " + input.get_planetary_file());

  myFile.open(input.get_planetary_file());

  if (!myFile.is_open()) {
    if (iProc == 0)
      std::cout << "Could not open planetary file : "
                << input.get_planetary_file() << "\n";

    DidWork = false;
  } else {

    if (myFile.good()) {

      std::vector<std::vector<std::string>> csv = read_csv(myFile);

      int nLines = csv.size();

      if (nLines <= 2) {
        report.print(0, "Number of lines in planetary CSV file too small");
        DidWork = false;
      } else {

        for (int iLine = 2; iLine < nLines; iLine++) {

          // Some final rows can have comments in them, so we want to
          // skip anything where the length of the string in column 2
          // is == 0:

          if (csv[iLine][1].length() > 0) {
            tmp.name = make_lower(csv[iLine][0]);
            tmp.semimajoraxis = stof(csv[iLine][1]);
            tmp.eccentricity = stof(csv[iLine][2]);
            tmp.inclination = stof(csv[iLine][3]);
            tmp.meanlongitude = stof(csv[iLine][4]);
            tmp.perihelionlongitude = stof(csv[iLine][5]);
            tmp.nodelongitude = stof(csv[iLine][6]);
            tmp.rates_semimajoraxis = stof(csv[iLine][7]);
            tmp.rates_eccentricity = stof(csv[iLine][8]);
            tmp.rates_inclination = stof(csv[iLine][9]);
            tmp.rates_meanlongitude = stof(csv[iLine][10]);
            tmp.rates_perihelionlongitude = stof(csv[iLine][11]);
            tmp.rates_nodelongitude = stof(csv[iLine][12]);
            tmp.length_of_day = stof(csv[iLine][13]) * cHtoS;
            tmp.length_of_year = stof(csv[iLine][14]) * cDtoS;
            tmp.longitude_jb2000 = stof(csv[iLine][15]);
            tmp.mass = stof(csv[iLine][16]);
            tmp.equator_radius = stof(csv[iLine][17]);
            tmp.polar_radius = stof(csv[iLine][18]);
            tmp.planet_tilt = stof(csv[iLine][19]) * cDtoR;
            tmp.dipole_strength = stof(csv[iLine][20]);
            tmp.dipole_rotation = stof(csv[iLine][21]) * cDtoR;
            tmp.dipole_tilt = stof(csv[iLine][22]) * cDtoR;

            for (int j = 0; j < 3; j++)
              tmp.dipole_center[j] = stof(csv[iLine][23 + j]);

            planets.push_back(tmp);
          }  // if length
        }  // for iLine
      }  // else nLines
    }  // if good file

    myFile.close();
  }  // else open file

  return DidWork;
}

// --------------------------------------------------------------------------
// check to see if class is ok
// --------------------------------------------------------------------------

bool Planets::is_ok() {
  return IsOk;
}

// --------------------------------------------------------------------------
// returns neutrals json for neutral density BCs
// --------------------------------------------------------------------------

json Planets::get_neutrals() {
  return neutrals;
}

// --------------------------------------------------------------------------
// returns neutral temperature json for neutral temperature BCs
// --------------------------------------------------------------------------

json Planets::get_temperatures() {
  return temperatures;
}

// --------------------------------------------------------------------------
// returns ions json for ion density initial conditions
// --------------------------------------------------------------------------

json Planets::get_ions() {
  return ions;
}

// -----------------------------------------------------------------------------
// Read in the planet specific file that describes the species
// -----------------------------------------------------------------------------

bool Planets::read_planet_specific_file(Inputs input, Report report) {

  bool DidWork = true;
  std::string hash;
  std::ifstream infile_ptr;

  std::string function = "Planets::read_planet_specific_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iDebug = 4;

  infile_ptr.open(input.get_planet_species_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open input file: "
              << input.get_planet_species_file() << "!!!\n";
    DidWork = false;
  } else {

    int IsDone = 0;

    while (!IsDone) {

      hash = find_next_hash(infile_ptr);

      if (report.test_verbose(iDebug))
        std::cout << "hash : -->" << hash << "<--\n";

      if (hash == "#neutrals") {
        report.print(iDebug, "Found #neutrals!");
        std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);
        neutrals = put_csv_in_json_w_name(lines);

        if (report.test_verbose(iDebug))
          std::cout << neutrals << "\n";
      } // #neutrals

      if (hash == "#temperature") {
        report.print(iDebug, "Found #temperatures!");
        std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);
        temperatures = put_csv_in_json_wo_name(lines);

        if (report.test_verbose(iDebug))
          std::cout << temperatures << "\n";
      }  // #temperature

      if (hash == "#ions") {
        report.print(iDebug, "Found #ions!");
        std::vector<std::vector<std::string>> lines = read_csv(infile_ptr);
        ions = put_csv_in_json_w_name(lines);

        if (report.test_verbose(iDebug))
          std::cout << ions << "\n";
      }  // #ions

      if (infile_ptr.eof())
        IsDone = 1;
    }   // while !IsDone

    infile_ptr.close();
  }  // else (isopen)

  report.exit(function);
  return DidWork;
}

