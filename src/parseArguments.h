#ifndef PARSE_ARGUMENTS_H_INCLUDED
#define PARSE_ARGUMENTS_H_INCLUDED

#include <vector>
#include <string>

struct parameters {
    float alpha; ///< regularization parameter
    int levels=0;
    bool segment;
    bool affine;
    bool rigid;
    std::vector<int> grid_spacing;
    std::vector<int> search_radius;
    std::vector<int> quantisation;
    std::string fixed_file;
    std::string moving_file;
    std::string output_stem;
    std::string moving_seg_file;
    std::string affine_file;
    std::string deformed_file;
};

char *realpathEx(const char *path, char *buff);

void parseCommandLine(parameters& args,int argc, char * const argv[]);

#endif // PARSE_ARGUMENTS_H_INCLUDED
