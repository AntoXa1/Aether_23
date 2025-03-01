// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "../include/aether.h"

// -------------------------------------------------------------------
// Read a line of the file, set to lower case, remove spaces until
// the end of the line or 2+ spaces are reached.
// -------------------------------------------------------------------

std::string make_lower(std::string instring) {
  std::string outstring = instring;
  int len = instring.length();

  for (int i = 0; i < len; i++)
    outstring[i] = tolower(outstring[i]);

  return outstring;
}

// -------------------------------------------------------------------
// Check to see if a string has 2+ spaces in it, and cut it off
// after that.
// -------------------------------------------------------------------

std::string strip_string_end(std::string instring) {
  int i;
  int iStart = -1;
  std::string outstring;
  int len = instring.length();

  for (i = 0; i < len; i++) {
    if (instring[i] == ' ') {
      if (i > 0 && iStart == i - 1)
        break;
      else
        iStart = i;
    }
  }

  if (iStart > -1 && i < len)
    outstring = instring.substr(0, iStart);
  else
    outstring = instring;

  return outstring;
}

// -------------------------------------------------------------------
// strip spaces out of the string.
// -------------------------------------------------------------------

std::string strip_spaces(std::string instring) {
  int i;
  std::string outstring;
  int len = instring.length();

  for (i = 0; i < len; i++) {
    if (instring[i] != ' ')
      outstring = outstring + instring[i];
  }

  return outstring;
}

// -------------------------------------------------------------------
// strip spaces out of the string.
// -------------------------------------------------------------------

std::string replace_spaces_with_commas(std::string instring) {
  int i;
  std::string outstring;
  int len = instring.length();
  i = 0;

  while (i < len) {
    if (instring[i] != ' ') {
      outstring = outstring + instring[i];
      i++;
    } else {
      outstring = outstring + ',';
      i++;

      while (i < len && instring[i] == ' ')
        i++;
    }
  }

  return outstring;
}

// -------------------------------------------------------------------
// Read a string, clean it up, and report error if length is 0
// -------------------------------------------------------------------

std::string read_string(std::ifstream &file_ptr, std::string hash) {

  std::string line = "";

  if (!file_ptr.is_open()) {
    std::cout << "File is not open (read_string)!\n";
    std::cout << "hash : " << hash << "\n";
  } else {
    getline(file_ptr, line);
    // Clearly these could all be combined, but I like separating them:
    line = strip_string_end(line);
    line = strip_spaces(line);

    if (line.length() < 1) {
      std::cout << "Issue in read_inputs!\n";
      std::cout << "Should be:\n";
      std::cout << hash << "\n";
      std::cout << "string variable with length > 0\n";
    }
  }

  return line;
}

// -------------------------------------------------------------------
// Read a string, clean it up, and convert it to an integer
// -------------------------------------------------------------------

int read_int(std::ifstream &file_ptr, std::string hash) {

  std::string line = "";
  int output = -1;

  if (!file_ptr.is_open()) {
    std::cout << "File is not open (read_int)!\n";
    std::cout << "hash : " << hash << "\n";
  } else {
    getline(file_ptr, line);
    line = strip_string_end(line);

    try {
      output = stoi(line);
    } catch (...) {
      std::cout << "Issue in read_integer!\n";
      std::cout << "In hash: ";
      std::cout << hash << "\n";
      std::cout << "Trying to read an integer, but got this: ";
      std::cout << line << "\n";
    }
  }

  return output;
}

// -------------------------------------------------------------------
// Read a string, clean it up, and convert it to an float
// -------------------------------------------------------------------

precision_t read_float(std::ifstream &file_ptr, std::string hash) {

  std::string line = "";
  precision_t output = -1;

  if (!file_ptr.is_open()) {
    std::cout << "File is not open (read_float)!\n";
    std::cout << "hash : " << hash << "\n";
  } else {
    getline(file_ptr, line);
    line = strip_string_end(line);

    try {
      output = stoi(line);
    } catch (...) {
      std::cout << "Issue in read_float!\n";
      std::cout << "In hash: ";
      std::cout << hash << "\n";
      std::cout << "Trying to read a float, but got this: ";
      std::cout << line << "\n";
    }
  }

  return output;
}

// -------------------------------------------------------------------
// Read the file until it gets to a # as the first character
// -------------------------------------------------------------------

std::string find_next_hash(std::ifstream &file_ptr) {

  std::string line;

  if (!file_ptr.is_open())
    std::cout << "File is not open (find_next_hash)!\n";

  else {
    while (getline(file_ptr, line)) {
      // Clearly these could all be combined, but I like separating them:
      line = strip_string_end(line);
      line = strip_spaces(line);
      line = make_lower(line);

      if (line[0] == '#')
        return line;
    }
  }

  line = "";
  return line;
}

// -------------------------------------------------------------------
// This takes a string with comma separated values and returns a
// vector of values.  Basically, it is one row of a CSV file.
// -------------------------------------------------------------------

std::vector<std::string> parse_csv_row_into_vector(std::string line) {

  std::vector<std::string> row;
  std::string col;
  std::stringstream ss(line);

  while (getline(ss, col, ','))
    row.push_back(col);

  return row;
}

// -------------------------------------------------------------------
// This can read a generic comma separated value file
// This is somewhat generic in that it can read from the current
// position to the first line that is blank, so it can be used
// multiple times per file.
// Returns a 2D array of strings.
// -------------------------------------------------------------------

std::vector<std::vector<std::string>> read_csv(std::ifstream &file_ptr) {

  std::vector<std::vector<std::string>> data;
  std::vector<std::string> row;
  std::string line, col;
  line = "  ";

  if (!file_ptr.is_open())
    std::cout << "File is not open (read_csv)!\n";

  else {
    // This assumes that the CSV file's layout is perfect - that the
    // number of columns is the same in each row.  If that is not the
    // case, then bad stuff happens.  I need to add more debugging
    // stuff in here.
    while (getline(file_ptr, line) && line.length() > 1) {
      line = strip_string_end(line);
      line = strip_spaces(line);
      row = parse_csv_row_into_vector(line);
      data.push_back(row);
    }
  }

  return data;
}

// -------------------------------------------------------------------
// This can read a generic SPACE separated value file
// This is somewhat generic in that it can read from the current
// position to the first line that is blank, so it can be used
// multiple times per file.
// Returns a 2D array of strings.
// -------------------------------------------------------------------

std::vector<std::vector<std::string>> read_ssv(std::ifstream &file_ptr) {

  std::vector<std::vector<std::string>> data;
  std::vector<std::string> row;
  std::string line, col;
  line = "  ";

  if (!file_ptr.is_open())
    std::cout << "File is not open (read_ssv)!\n";

  else {
    // This assumes that the SSV file's layout is perfect - that the
    // number of columns is the same in each row.  If that is not the
    // case, then bad stuff happens.  I need to add more debugging
    // stuff in here.
    while (getline(file_ptr, line) && line.length() > 1) {
      line = replace_spaces_with_commas(line);
      row = parse_csv_row_into_vector(line);
      data.push_back(row);
    }
  }

  return data;
}

// -------------------------------------------------------------------
// read in a time format (6 integers + 1 extra added as a zero pad)
// -------------------------------------------------------------------
std::vector<int> read_itime(std::ifstream &file_ptr, std::string hash) {

  int iErr = 0;
  std::vector<int> itime(7, 0);

  // Read in the series of numbers:
  std::vector<std::vector<std::string>> stimes = read_csv(file_ptr);

  // Interpret whether it is a column (first) or row (second):
  if (stimes.size() == 6) {  // column
    for (int i = 0; i < 6; i++)
      itime[i] = stoi(stimes[i][0]);
  } else {
    if (stimes.size() == 1 && stimes[0].size() >= 6) {  // row
      for (int i = 0; i < 6; i++)
        itime[i] = stoi(stimes[0][i]);
    } else
      iErr = 1;
  }

  if (iErr == 0)
    itime[6] = 0;

  else {
    itime[0] = -1;
    std::cout << "Issue in read_inputs!\n";
    std::cout << "Should be:\n";
    std::cout << hash << "\n";
    std::cout << "year     (int)\n";
    std::cout << "month    (int)\n";
    std::cout << "day      (int)\n";
    std::cout << "hour     (int)\n";
    std::cout << "min      (int)\n";
    std::cout << "sec      (int)\n";
    std::cout << "\n";
    std::cout << "or \n";
    std::cout << "\n";
    std::cout << hash << "\n";
    std::cout << "year, month, day, hour, minute, second\n";
  }

  return itime;
}


// -----------------------------------------------------------------------------
// Move CSV to json format (with name as first column)
// -----------------------------------------------------------------------------

json put_csv_in_json_w_name(std::vector<std::vector<std::string>>
                            csvLines) {

  bool DidWork = true;
  json output;

  // This assumes the first column is "name" and is a string
  json name;
  int64_t nRows = csvLines.size();

  for (int iRow = 1; iRow < nRows; iRow++)
    name.push_back(csvLines[iRow][0]);

  output["name"] =  name;

  // The rest of the lines can be in any order, really
  int64_t nCols = csvLines[0].size();

  for (int iCol = 1; iCol < nCols; iCol++) {
    json var;

    for (int iRow = 1; iRow < nRows; iRow++)
      var.push_back(stof(csvLines[iRow][iCol]));

    output[csvLines[0][iCol]] = var;
  }

  return output;
}

// -----------------------------------------------------------------------------
// Move CSV to json format (no name as first column)
// -----------------------------------------------------------------------------

json put_csv_in_json_wo_name(std::vector<std::vector<std::string>>
                             csvLines) {

  json output;

  // The rest of the lines can be in any order, really
  int64_t nRows = csvLines.size();
  int64_t nCols = csvLines[0].size();

  for (int iCol = 0; iCol < nCols; iCol++) {
    json var;

    for (int iRow = 1; iRow < nRows; iRow++)
      var.push_back(stof(csvLines[iRow][iCol]));

    output[csvLines[0][iCol]] = var;
  }

  return output;
}

