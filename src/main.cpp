#include "config.h"
#include "loadInitParam.h"

using namespace std;

int main(int argc, char **argv)
{
  char * filename;

  if (lib::cmd_option_exists(argv, argv+argc, "-h"))
  {
    cerr << "USAGE:" << endl << "  pdp3 [ --version | -f path/to/parameters.xml ]" << endl;
    return 1;
  }

  if (lib::cmd_option_exists(argv, argv+argc, "--version"))
  {
    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
    return 0;
  }

  if (lib::cmd_option_exists(argv, argv+argc, "-f"))
  {
    filename = lib::get_cmd_option(argv, argv + argc, "-f");
    if (filename == NULL)
    {
      cerr << "ERROR: configuration path is not specified" << endl;
      return 1;
    }
  }
  else
    filename = (char*)"parameters.xml";

  cout << endl;
  cout << "Wellcome to " << PACKAGE_NAME << ", version: " << PACKAGE_VERSION << endl;
  cout << endl;
#ifdef EXPERIMENTAL
  cerr << "WARNING! Experimental features enabled" << endl;
#endif
#ifdef DEBUG
  cerr << "INFO! Running in debug mode" << endl;
#endif
#ifdef SINGLETHREAD
  cerr << "INFO! Running in single-thread mode" << endl;
#endif
#ifdef SPEEDUP
  cerr << "INFO! Running with fast-math option" << endl;
#endif
#ifdef PROFILER
  cerr << "INFO! Running with profiler" << endl;
#endif
#ifdef PUSHER_BORIS_CLASSIC
  cerr << "INFO! Using non-relativistig Boris particles pusher" << endl;
#endif
#ifdef PUSHER_BORIS_ADAPTIVE
  cerr << "INFO! Using adaptive Boris particles pusher" << endl;
#endif
#ifdef PUSHER_BORIS_RELATIVISTIC
  cerr << "INFO! Using fully-relativistic Boris particles pusher" << endl;
#endif
#ifdef PUSHER_VAY
  cerr << "INFO! Using fully-relativistic Vay particles pusher" << endl;
#endif
#ifdef TESTMODE
  cerr << "WARNING! Using TESTMODE with reproducible random numbers to compare result" << endl;
#endif

  LoadInitParam init_param(filename);

  init_param.run();

  return 0;
}
