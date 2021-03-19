
#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_API_CXX

#include "API.hxx"

namespace Polyphemus
{
  

  //! Default constructor.
  template<class T>
  API<T>::API()
  {
  }

    //! Destructor.
  template<class T>
  API<T>::~API()
  {
  }

  
  //! Get pointer to a shared library.
  /*!
    \param handle pointer to a shared library.
    \param name name of function symbol in library.
    \param errors_are_fatal abort if true, silently ignore if false
    \return pointer to function in shared library.
  */
  template<class T>
  void *API<T>::get_dl_function_pointer(void *handle,
                                        const char  *name,
                                        bool errors_are_fatal)
  {

    void  *retval = NULL;
    char  *error = NULL;
    char  *name_ = NULL;

    dlerror();    /* Clear any existing error */

    retval = dlsym(handle, name);
    error = dlerror();

    if (error != NULL) { /* Try different symbol names */
      dlerror();    /* Clear any existing error */
      int _size_ = strlen(name) + strlen("_");
      name_ = (char*)malloc((_size_ + 1)*sizeof(char));
      strcpy(name_, name);
      strcat(name_, "_");
      retval = dlsym(handle, name_);
      error = dlerror();
      
      if (error != NULL) cout << "Error on dlsym : " << error << endl;
    
      free(name_);
    }
    
    if (error != NULL && errors_are_fatal) {
      cout << "Error, abort program" << endl;
      exit(1);
    }
    
    return retval;
  }  


  /*----------------------------------------------------------------------------
   * Call a function of the shared library
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  void API<T>::call(void *handle,
                    const char *name)
  {
    
    typedef void* (*_tmp_sshaerosol_t)(void);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct();
  }

  /*----------------------------------------------------------------------------
   * Send a boolean to SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   flag             --> boolean exchanged with the external code
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  void API<T>::send_bool(void *handle,
                         const char *name,
                         bool flag)
  {
  
    typedef void* (*_tmp_sshaerosol_t)(bool*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(&flag);
  }


  /*----------------------------------------------------------------------------
   * Send an integer to SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   flag             --> boolean exchanged with the external code
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  void API<T>::send_int(void *handle,
                        const char *name,
                        int val)
  {
  
    typedef void* (*_tmp_sshaerosol_t)(int*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(&val);
  }

  
  /*----------------------------------------------------------------------------
   * Send a double to SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   val              --> double exchanged with the external code
   *
   *----------------------------------------------------------------------------*/

  template<class T>
  void API<T>::send_double(void *handle,
                           const char *name,
                           T val)
  {

    typedef void* (*_tmp_sshaerosol_t)(T*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    T tmp = val;
    fct(&tmp);
  }


  /*----------------------------------------------------------------------------
   * Receive a boolean from SSH-aerosol, returns it
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *
   * returns the boolean received from the external code
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  bool API<T>::recv_bool(void *handle,
                         const char *name)
  {

    typedef bool (*_tmp_sshaerosol_t)(void);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    bool res = fct();
    
    return res;
  }


  /*----------------------------------------------------------------------------
   * Receive a int from SSH-aerosol, returns it
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *
   * returns the integer received from the external code
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  int API<T>::recv_int(void *handle,
                       const char *name)
  {
    
    typedef int (*_tmp_sshaerosol_t)(void);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    int res = fct();
    
    return res;
  }

  /*----------------------------------------------------------------------------
   * Receive a double from SSH-aerosol, returns it
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *
   * returns the double received from the external code
   *
   *----------------------------------------------------------------------------*/

  template<class T>
  T API<T>::recv_double(void *handle,
                        const char *name)
  {
  
    typedef T (*_tmp_sshaerosol_t)(void);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    T res = fct();
    
    return res;
  }

  /*----------------------------------------------------------------------------
   * Send a char array to SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   array            <-- array exchanged with the external code
   *
   *----------------------------------------------------------------------------*/
    template<class T>
    void API<T>::exchange_char_array(void *handle,
                                     const char *name,
                                     const char *array)
    {

      typedef void* (*_tmp_sshaerosol_t)(const char*);
      _tmp_sshaerosol_t fct =
        (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                    name,
                                                    true);
      fct(array);
    }


  /*----------------------------------------------------------------------------
   * Send a char array to SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   array            <-- array exchanged with the external code
   *
   *----------------------------------------------------------------------------*/
    template<class T>
    void API<T>::exchange_char_array(void *handle,
                                     const char *name,
                                     int id,
                                     char *array)
    {

      typedef void* (*_tmp_sshaerosol_t)(int*, const char*);
      _tmp_sshaerosol_t fct =
        (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                    name,
                                                    true);
      fct(&id, array);
    }
  

  /*----------------------------------------------------------------------------
   * Exchange a double array with SSH-aerosol
   *
   * parameters:
   *   handle           <-- pointer to shared library (result of dlopen)
   *   name             <-- name of function symbol in SSH-aerosol
   *   array            <-> array exchanged with the external code
   *
   *----------------------------------------------------------------------------*/
  template<class T>
  void API<T>::exchange_double_array(void *handle,
                                     const char *name,
                                     T *array)
  {
    typedef void* (*_tmp_sshaerosol_t)(T*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array);
  }

  template<class T>
  void API<T>::exchange_double_array(void *handle,
                                     const char *name,
                                     Array<T, 1> array)
  {
    
    typedef void* (*_tmp_sshaerosol_t)(T*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array.data());
  }

  template<class T>
  void API<T>::exchange_double_array(void *handle,
                                     const char *name,
                                     Array<T, 2> array)
  {

    typedef void* (*_tmp_sshaerosol_t)(T*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array.data());
  }

  template<class T>
  void API<T>::exchange_double_array(void *handle,
                                     const char *name,
                                     Array<T, 3> array)
  {

    typedef void* (*_tmp_sshaerosol_t)(T*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array.data());
  }
  

  template<class T>
  void API<T>::exchange_int_array(void *handle,
                                  const char *name,
                                  Array<int, 1> array)
  {
    
    typedef void* (*_tmp_sshaerosol_t)(int*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array.data());
  }

  template<class T>
  void API<T>::exchange_int_array(void *handle,
                                  const char *name,
                                  Array<int, 2> array)
  {
    typedef void* (*_tmp_sshaerosol_t)(int*);
    _tmp_sshaerosol_t fct =
      (_tmp_sshaerosol_t) get_dl_function_pointer(handle,
                                                  name,
                                                  true);
    fct(array.data());
  }
  
  
}

#define POLYPHEMUS_FILE_MODULES_AEROSOL_API_CXX
#endif
