#ifndef POLYPHEMUS_FILE_MODELS_STREETCHEMITRY_CXX

#include "StreetChemistry.hxx"

namespace Polyphemus
{


  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  //! Main constructor.
  /*!

   */
  template<class T>
  StreetChemistry<T>::StreetChemistry(int street_id, int begin_inter, int end_inter,
                    T length, T width, T height,
                                     int ns_local, int Nr_photolysis):
    Street<T>(street_id, begin_inter, end_inter, length, width, height, ns_local)
  {
    photolysis_rate_.resize(Nr_photolysis);
  }

  //! Destructor
  template<class T>
  StreetChemistry<T>::~StreetChemistry()
  {
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_CXX
#endif

