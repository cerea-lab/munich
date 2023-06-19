
//////////////
// INCLUDES //

#define SELDONDATA_DEBUG_LEVEL_2

#include "AtmoData.hxx"
#include "StreetDriver.cxx"
#include "BaseOutputSaver.cxx"
#include "StreetNetworkAerosol.cxx"
#include "Aerosol_SSH.cxx"
#include <time.h>


using namespace Polyphemus;

// INCLUDES //
//////////////
  
int main(int argc, char** argv)
{

  TRY;
  clock_t tBegin, tEnd, tPassed_time;
  tBegin = clock();
  
  if (argc != 2)
    {
      string mesg  = "Usage:\n";
      mesg += string("  ") + argv[0] + " [configuration file]";
      cout << mesg << endl;
      return 1;
    }

  typedef double real;
  typedef StreetNetworkAerosol<real, Aerosol_SSH<real> > ClassModel;

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Init(&argc, &argv);    
#endif
  
  StreetDriver<real, ClassModel, BaseOutputSaver<real, ClassModel> >
    Driver(argv[1]);

  Driver.Run();

#ifdef POLYPHEMUS_PARALLEL_WITH_MPI
    MPI_Finalize();    
#endif
  
  tEnd = clock();
  tPassed_time = ((tEnd - tBegin) / (CLOCKS_PER_SEC / 1000));
  // cout<<"Passed time: "<<tPassed_time<<endl;

  END;

  return 0;

}
