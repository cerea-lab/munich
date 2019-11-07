
//////////////
// INCLUDES //

#define SELDONDATA_DEBUG_LEVEL_4

#include "AtmoData.hxx"
#include "StreetDriver.cxx"
#include "BaseOutputSaver.cxx"
#include "StreetNetworkChemistry.cxx"
#include "Photochemistry.cxx"

using namespace Polyphemus;

// INCLUDES //
//////////////
  
int main(int argc, char** argv)
{

  TRY;

  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  typedef double real;
  typedef StreetNetworkChemistry<real, Photochemistry<real> > ClassModel;

  StreetDriver<real, ClassModel, BaseOutputSaver<real, ClassModel> >
    Driver(argv[1]);

  Driver.Run();

  END;

  return 0;

}
