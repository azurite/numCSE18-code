#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <iomanip>

///
/// Writes the contents of the vector 'data' to the textfile filename
/// The output format should be load-able in MATLAB and Numpy using the
/// load-command.
///
template<typename T>
void writeToFile(const std::string& filename,
         const std::vector<T>& data) {
    
    std::ofstream file(filename.c_str());
    // Set highest possible precision, this way we are sure we are
    file << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    // Loop over vector and write output to file
    for(int i = 0; i < data.size(); ++i) {
    file << data[i] << " ";
    }
    file << std::endl;

    // File closes automatically at end of scope!
}
