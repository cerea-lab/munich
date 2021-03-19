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
  class StreetChemistry: public Street<T>
  {

  protected:
    
    //! Photolysis rate (1/s)
    Array<T, 1> photolysis_rate_;
    //! Additional meteo data for chemisty.
    T attenuation_;
    T specific_humidity_;

  public:

     /*** Constructor and destructor ***/

    StreetChemistry(int street_id,
		    int begin_inter,
		    int end_inter, 
		    T length,
		    T width,
		    T height,
		    int typo,
		    int ns_local,
		    int nr_photolysis);
    virtual ~StreetChemistry();   

    /*** Methods ***/
    
    T GetPhotolysisRate(int r) const;
    void SetPhotolysisRate(Array<T, 1> photolysis_rate);
    void SetAttenuation(T attenuation);
    T GetAttenuation() const;
  };

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_HXX
#endif
