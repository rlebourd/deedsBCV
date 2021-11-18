#ifndef PARSE_ARGUMENTS_H_INCLUDED
#define PARSE_ARGUMENTS_H_INCLUDED

#include <vector>
#include <string>

struct parameters{
    float alpha; int levels=0; bool segment,affine,rigid;
    std::vector<int> grid_spacing; std::vector<int> search_radius;
    std::vector<int> quantisation;
    std::string fixed_file,moving_file,output_stem,moving_seg_file,affine_file,deformed_file;
};

char *realpathEx(const char *path, char *buff);

void parseCommandLine(parameters& args,int argc, char * const argv[]);

#endif // PARSE_ARGUMENTS_H_INCLUDED
