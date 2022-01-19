#include <iostream>
#include <cmath>

// Include general/common includes, utilities and initialisation
#include "../../../src/general/common_includes.h"
#include "../../../src/general/utilities.h"
#include "../../../src/general/initialise.h"

using namespace scicellxx;

// This is a basic project template
int main(int argc, char *argv[])
{
 // Initialise scicell++
 initialise_scicellxx();
 
 // Output for testing/validation
 std::ofstream output_test("output_test.dat", std::ios_base::out);

#ifdef SCICELLXX_PANIC_MODE
 // Output for testing/validation
 std::ofstream deb("deb.dat", std::ios_base::out);
#endif // #ifdef SCICELLXX_PANIC_MODE

 unsigned i = 0;
 Real x = 1.0e-3;

 std::cout << "Basic output" << std::endl;
 std::cout << "i: " << i << std::endl;
 std::cout << "x: " << x << std::endl;
 
 output_test << "Basic output to file" << std::endl;
 output_test << "i: " << i << std::endl;
 output_test << "x: " << x << std::endl;
 
#ifdef SCICELLXX_PANIC_MODE
 DEB(i);
 DEB(x);
 DEB2(i,x);
 DEB_TO_FILE(deb, i);
 DEB_TO_FILE(deb, x);
 DEB_TO_FILE2(deb, i, x);
#endif // #ifdef SCICELLXX_PANIC_MODE
 
 // Close the output for test
 output_test.close();
 
#ifdef SCICELLXX_PANIC_MODE
 deb.close();
#endif // #ifdef SCICELLXX_PANIC_MODE
 
 // Finalise scicell++
 finalise_scicellxx();
 
 return 0;
 
}
