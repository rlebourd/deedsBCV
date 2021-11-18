#ifndef PARSE_ARGUMENTS_H_INCLUDED
#define PARSE_ARGUMENTS_H_INCLUDED

struct parameters; // forward declaration of parameters

char *realpathEx(const char *path, char *buff);

void parseCommandLine(parameters& args,int argc, char * const argv[]);

#endif // PARSE_ARGUMENTS_H_INCLUDED
