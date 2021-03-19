#ifndef POLYPHEMUS_FILE_MODULES_AEROSOL_API_HXX

#include <dlfcn.h>

namespace Polyphemus
{

  using namespace std;
  using namespace AtmoData;


  template<class T>
  class API
  {

  public:
    
    /*** Constructor and destructor ***/

    API();
    ~API();
    
    /*** Other methods ***/
    
    void *get_dl_function_pointer(void *handle,
                                  const char *name,
                                  bool errors_are_fatal);

    void call(void *handle,
              const char *name);

    void send_bool(void *handle,
                   const char *name,
                   bool flag);

    void send_int(void *handle,
                  const char *name,
                  int val);
  
    void send_double(void *handle,
                     const char *name,
                     T val);

    bool recv_bool(void *handle,
                   const char *name);

    int recv_int(void *handle,
                 const char *name);

    T recv_double(void *handle,
                  const char *name);

    void exchange_char_array(void *handle,
                             const char *name,
                             const char *array);

    void exchange_char_array(void *handle,
                             const char *name,
                             int id,
                             char *array);
    
    void exchange_double_array(void *handle,
                               const char *name,
                               T *array);

    void exchange_double_array(void *handle,
                               const char *name,
                               Array<T, 1> array);

    void exchange_double_array(void *handle,
                               const char *name,
                               Array<T, 2> array);

    void exchange_double_array(void *handle,
                               const char *name,
                               Array<T, 3> array);

    void exchange_int_array(void *handle,
                            const char *name,
                            Array<int, 1> array);
    
    void exchange_int_array(void *handle,
                            const char *name,
                            Array<int, 2> array);        

  };

} // namespace Polyphemus.


#define POLYPHEMUS_FILE_MODULES_AEROSOL_API_HXX
#endif
