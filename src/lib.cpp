#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include <algorithm>
#include <iomanip>
#include "constant.h"

// #include <ctime>
#include <math.h>
#include "lib.h"

// #include <fstream>
// #include <iostream>

// #include "particles.h"
// #include "eField.h"
// #include "hField.h"
// #include "pdp3Time.h"
// #include "triple.h"
// #include <math.h>
// #include <string.h>
// #include <cstdlib>
// #include "constant.h"

using namespace std;
using namespace constant;

namespace lib

{
	// C^2 define c^2 to decrease number of operations
	const double LIGHT_SPEED_POW_2 = pow (LIGHT_SPEED, 2);

	bool to_bool(string str)
	{
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
		std::istringstream is(str);
		bool b;
		is >> std::boolalpha >> b;
		return b;
	}

	double get_gamma (double velocity)
	{
		double gamma, beta;

		beta = pow(velocity, 2) / LIGHT_SPEED_POW_2;

		if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
		{
			cerr << "CRITICAL!(get_gamma): Lorentz factor aka gamma is complex. Can not continue." << endl
					 << "\tvelocity is: " << velocity << endl;
			exit(1);
		}

		gamma = pow(1.0 - beta, -0.5);

		if (isinf(gamma) == 1) { // avoid infinity values
			cerr << "WARNING(get_gamma): gamma (Lorenz factor) girects to infinity for velocity" << velocity << endl;
			return 1e100; // just return some very big value
		}
		else
		{
			return gamma;
		}
	}

	double get_gamma_inv (double velocity)
	{
		double gamma, beta;

		beta = pow(velocity, 2) / LIGHT_SPEED_POW_2;

		if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
		{
			cerr << "CRITICAL!(get_gamma_inv): Lorentz factor aka gamma is complex. Can not continue." << endl
					 << "\tvelocity is: " << velocity << endl;
			exit(1);
		}

		// gamma = pow(1.0 - beta, -0.5);
		gamma = pow(1.0 + beta, 0.5);

		if (isinf(gamma) == 1) { // avoid infinity values
			cerr << "WARNING(get_gamma_inv): gamma (Lorenz factor) girects to infinity for velocity" << velocity << endl;
			return 1e100; // just return some very big value
		}
		else
		{
			return gamma;
		}
	}
}



// gamma = pow(1.0 + beta, 0.5);
