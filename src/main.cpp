/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/********************************************************************************/

#include <Minia.hpp>
#include "build_info.hpp"

using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{

	if(argc > 1 && (   strcmp(argv[1],"--version")==0 || strcmp(argv[1],"-v")==0 || strcmp(argv[1],"-version")==0    )     ){
        std::cout << "Minia version " << STR_MINIA_VERSION << std::endl;
#ifdef GIT_SHA1
        std::cout << "git commit " << GIT_SHA1 << std::endl;
#endif
     	std::cout << "Using gatb-core version "<< System::info().getVersion() << std::endl;
     	std::cout << "OS: " << STR_MINIA_OPERATING_SYSTEM << std::endl;
        std::cout << "compiler: " << STR_MINIA_COMPILER << std::endl;
		return EXIT_SUCCESS;
	}


    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We execute dsk. */
        Minia().run (argc, argv);
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

