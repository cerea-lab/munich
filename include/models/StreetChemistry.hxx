#ifndef POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_HXX

#include "StreetTransport.cxx"

namespace Polyphemus
{

  //////////////
  // INCLUDES //
  //////////////
  
  using namespace std;

  ////////////
  // STREET //
  ////////////


  //! Class that stores all data about a street.
  template<class T>
  class StreetChemistry
  {

  protected:

    //! Emission rate (ug/s)
    Array<T, 1> photolysis_rate_;

  public:

     /*** Constructor and destructor ***/

    StreetChemistry(int street_id, int begin_inter, int end_inter, 
                    T length, T width, T height, int ns_local, int Nr_photolysis);
    virtual ~StreetChemistry();   

  };

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_HXX
#endif
