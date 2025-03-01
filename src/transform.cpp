// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <math.h>
#include <vector>

#include "../include/aether.h"

// -----------------------------------------------------------------------
// copy from c++ vector to c-native array
// -----------------------------------------------------------------------

void copy_vector_to_array(std::vector<float> vector_in,
                          int64_t nElements,
                          float *array_out) {

  for (int64_t i = 0; i < nElements; i++)
    array_out[i] = vector_in[i];
}

// -----------------------------------------------------------------------
// copy from armidillo cube to 3d c-native array
// -----------------------------------------------------------------------

void copy_cube_to_array(arma_cube cube_in,
                        float *array_out) {

  int64_t nX = cube_in.n_rows;
  int64_t nY = cube_in.n_cols;
  int64_t nZ = cube_in.n_slices;
  int64_t iX, iY, iZ, index;

  for (iX = 0; iX < nX; iX++) {
    for (iY = 0; iY < nY; iY++) {
      for (iZ = 0; iZ < nZ; iZ++) {
        index = iX * nY * nZ + iY * nZ + iZ;
        array_out[index] = cube_in(iX, iY, iZ);
      }
    }
  }
}



// -----------------------------------------------------------------------
// Transform Longitude (llr[0]), Latitude (llr[1]), Radius (llr[2]) to
// X, Y, Z
// Use armidillo cubes
// -----------------------------------------------------------------------

std::vector<arma_cube> transform_llr_to_xyz_3d(std::vector<arma_cube> llr) {
  std::vector<arma_cube> xyz;
  xyz.push_back(llr[2] % cos(llr[1]) % cos(llr[0]));
  xyz.push_back(llr[2] % cos(llr[1]) % sin(llr[0]));
  xyz.push_back(llr[2] % sin(llr[1]));
  return xyz;
}

// -----------------------------------------------------------------------
// Transform Longitude, Latitude, Radius to X, Y, Z
// -----------------------------------------------------------------------

void transform_llr_to_xyz(precision_t llr_in[3], precision_t xyz_out[3]) {
  // llr_in[0] = longitude (in radians)
  // llr_in[1] = latitude (in radians)
  // llr_in[2] = radius
  xyz_out[0] = llr_in[2] * cos(llr_in[1]) * cos(llr_in[0]);
  xyz_out[1] = llr_in[2] * cos(llr_in[1]) * sin(llr_in[0]);
  xyz_out[2] = llr_in[2] * sin(llr_in[1]);
}

// -----------------------------------------------------------------------
// Rotate 3D array (cube) around the z-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

std::vector<arma_cube> rotate_around_z_3d(std::vector<arma_cube> XYZ_in,
                                          precision_t angle) {

  arma_cube X = XYZ_in[0];
  arma_cube Y = XYZ_in[1];
  arma_cube Z = XYZ_in[2];
  std::vector<arma_cube> XYZ_out;

  precision_t ca = cos(angle);
  precision_t sa = sin(angle);

  XYZ_out.push_back( X * ca + Y * sa);
  XYZ_out.push_back(-X * sa + Y * ca);
  XYZ_out.push_back(Z);

  return XYZ_out;
}

// -----------------------------------------------------------------------
// Rotate 3D array (cube) around the y-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

std::vector<arma_cube> rotate_around_y_3d(std::vector<arma_cube> XYZ_in,
                                          precision_t angle) {

  arma_cube X = XYZ_in[0];
  arma_cube Y = XYZ_in[1];
  arma_cube Z = XYZ_in[2];
  std::vector<arma_cube> XYZ_out;

  precision_t ca = cos(angle);
  precision_t sa = sin(angle);

  XYZ_out.push_back(X * ca - Z * sa);
  XYZ_out.push_back(Y);
  XYZ_out.push_back(X * sa + Z * ca);

  return XYZ_out;
}

// -----------------------------------------------------------------------
// Rotate 3D array (cube) around the x-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

std::vector<arma_cube> rotate_around_x_3d(std::vector<arma_cube> XYZ_in,
                                          precision_t angle) {

  arma_cube X = XYZ_in[0];
  arma_cube Y = XYZ_in[1];
  arma_cube Z = XYZ_in[2];
  std::vector<arma_cube> XYZ_out;

  precision_t ca = cos(angle);
  precision_t sa = sin(angle);

  XYZ_out.push_back(X);
  XYZ_out.push_back( Y * ca + Z * sa);
  XYZ_out.push_back(-Y * sa + Z * ca);

  return XYZ_out;
}

// -----------------------------------------------------------------------
// Rotate around the z-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

void transform_rot_z(precision_t xyz_in[3], precision_t angle_in,
                     precision_t xyz_out[3]) {
  precision_t ca = cos(angle_in);
  precision_t sa = sin(angle_in);
  xyz_out[0] =  xyz_in[0] * ca + xyz_in[1] * sa;
  xyz_out[1] = -xyz_in[0] * sa + xyz_in[1] * ca;
  xyz_out[2] =  xyz_in[2];
}

// -----------------------------------------------------------------------
// Rotate around the y-axis
//  - Angle needs to be in radians!!!
// -----------------------------------------------------------------------

void transform_rot_y(precision_t xyz_in[3], precision_t angle_in,
                     precision_t xyz_out[3]) {
  precision_t ca = cos(angle_in);
  precision_t sa = sin(angle_in);
  xyz_out[0] = xyz_in[0] * ca - xyz_in[2] * sa;
  xyz_out[1] = xyz_in[1];
  xyz_out[2] = xyz_in[0] * sa + xyz_in[2] * ca;
}

// -----------------------------------------------------------------------
// Simply move data from a vector to a C-native array (float)
// -----------------------------------------------------------------------

void transform_float_vector_to_array(std::vector<float> input,
                                     precision_t output[3]) {
  for (int i = 0; i < 3; i++)
    output[i] = input[i];
}

// -----------------------------------------------------------------------
// Rotate a vector from XYZ to East, North, Vertical
// -----------------------------------------------------------------------

void transform_vector_xyz_to_env(precision_t xyz_in[3],
                                 precision_t lon,
                                 precision_t lat,
                                 precision_t env_out[3]) {

  env_out[2] =   xyz_in[0] * cos(lat) * cos(lon) +
                 xyz_in[1] * cos(lat) * sin(lon) + xyz_in[2] * sin(lat);
  env_out[1] = -(xyz_in[0] * sin(lat) * cos(lon) +
                 xyz_in[1] * sin(lat) * sin(lon) -
                 xyz_in[2] * cos(lat));
  env_out[0] = - xyz_in[0] * sin(lon) +
               xyz_in[1] * cos(lon);
}

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// These are not realy transform functions.  I need a better place to
// put them.
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Simple 3-element vector difference
// -----------------------------------------------------------------------

void vector_diff(precision_t vect_in_1[3],
                 precision_t vect_in_2[3],
                 precision_t vect_out[3]) {
  for (int i = 0; i < 3; i++)
    vect_out[i] = vect_in_1[i] - vect_in_2[i];
}

// -----------------------------------------------------------------------
// Simple 3-element vector addition
// -----------------------------------------------------------------------

void vector_add(float vect_in_1[3],
                 float vect_in_2[3],
                 float vect_out[3]) {
  for (int i = 0; i < 3; i++) vect_out[i] = vect_in_1[i] + vect_in_2[i];
}
