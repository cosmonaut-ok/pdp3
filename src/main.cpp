#include <algorithm>
#include "loadInitParam.h"

using namespace std;

char *get_cmd_option(char  **begin, char  **end, const std::string & option)
{
  char  **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    return *itr;
  return 0;
}

bool cmd_option_exists(char **begin, char **end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}



int main(int argc, char **argv)
{
  char  *filename;

  if (cmd_option_exists(argv, argv+argc, "-h"))
  {
    cerr << "USAGE:" << endl << "  pdp3 [ -f configfile_path ]" << endl;
    return 1;
  }

  if (cmd_option_exists(argv, argv+argc, "-f"))
  {
    filename = get_cmd_option(argv, argv + argc, "-f");
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

  LoadInitParam init_param(filename);

  init_param.run();

  return 0;
}
