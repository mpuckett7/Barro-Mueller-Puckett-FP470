/**
 * This file is not a part of the original OpenLB library. This was written to
 * give examples a way to change the lattice resolution without needing to
 * recompile the code. 
 * 
 * To use: ...
 */

#ifndef RESOLUTION_INPUT_HH
#define RESOLUTION_INPUT_HH
#include <iostream>
#include <cstdlib>
#include <unistd.h>

namespace olb {
namespace util {

int parseResolution(int argc, char *argv[], int defaultResolution) {
    int ch = 0;
    int resolution = defaultResolution;
    while ((ch = getopt (argc, argv, "r:")) != -1) {
        switch (ch) {
            case 'r':
                if ((resolution = static_cast<int>(std::strtol(optarg, NULL, 10))) > 0) {
                    return resolution; 
                }
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " -r <resolution>" << std::endl;
                exit(EXIT_FAILURE);
                break;
        }
    }
    return defaultResolution;
}
}
}
#endif