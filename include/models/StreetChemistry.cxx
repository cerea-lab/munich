#ifndef POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_CXX

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
  StreetChemistry<T>::StreetChemistry(int street_id,
				      int begin_inter,
				      int end_inter,
				      T length,
				      T width,
				      T height,
				      int typo,
				      int ns_local,
				      int nr_photolysis):
    Street<T>(street_id, begin_inter, end_inter, length, width, height, typo, ns_local)
  {
    photolysis_rate_.resize(nr_photolysis);
    photolysis_rate_ = 0.0;
    attenuation_ = 0.0;
    specific_humidity_ = 0.0;
  }

  //! Destructor
  template<class T>
  StreetChemistry<T>::~StreetChemistry()
  {
  }


  //////////////////////////////////////////
  // ACCESS METHODS FOR STREET ATTRIBUTES //
  //////////////////////////////////////////

  //! Returns the photolysis rate.
  template<class T>
  inline T StreetChemistry<T>::GetPhotolysisRate(int r) const
  {
    return photolysis_rate_(r);
  }

  template<class T>
  inline void StreetChemistry<T>::SetPhotolysisRate(Array<T, 1> photolysis_rate)
  {
    photolysis_rate_ = photolysis_rate;
  }

  //! Sets the meteo data.
  /*!
    \param attenuation the attenuation.
    \param specific_humidity the specific humidity.
  */
  template<class T>
  inline void StreetChemistry<T>::SetAttenuation(T attenuation)
  {
    attenuation_ = attenuation;
  }

  //! Returns the attenuation
  /*!
    \return The attenuation.
  */
  template<class T>
  inline T StreetChemistry<T>::GetAttenuation() const
  {
    return attenuation_;
  }

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODELS_STREETCHEMISTRY_CXX
#endif

