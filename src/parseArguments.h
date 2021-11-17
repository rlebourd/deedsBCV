#ifndef PARSE_ARGUMENTS_H_INCLUDED
#define

char *realpathEx(const char *path, char *buff);

void parseCommandLine(parameters& args,int argc, char * const argv[]);

#endif // PARSE_ARGUMENTS_H_INCLUDED
