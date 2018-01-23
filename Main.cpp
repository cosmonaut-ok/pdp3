#include <algorithm>
#include "Load_init_param.h"

using namespace std;

Particles_struct specie;
// #define BUILD_OPENCL

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}



int main(int argc, char **argv)
{
  char * filename;

  if (cmdOptionExists(argv, argv+argc, "-h"))
    {
      cerr << "USAGE:\n  pdp3 [ -f configfile_path ]\n";
      return 1;
    }
  
  if (cmdOptionExists(argv, argv+argc, "-f"))
    {
      filename = getCmdOption(argv, argv + argc, "-f");
      if (filename == NULL)
	{
	  cerr << "ERROR: configuration path is not specified" << endl;
	  return 1;
	}
    }
  else
    {
      filename = (char*)"parameters.xml";
    }

  Load_init_param init_param(filename);

  init_param.run();

  return 0;
}
