

#ifndef MB_ALLOCATOR_HPP
#define MB_ALLOCATOR_HPP

// Don't #include any MB files here,
// this is just a framework for memory allocations.
#include <memory>

// C malloc/free fucntions
void* mdb_malloc(size_t size, const char* filename=NULL, int linenumber=0);
void* mdb_calloc(size_t size, const char* filename=NULL, int linenumber=0);
void mdb_free(void* ptr);

#define MB_MALLOC(size) mdb_malloc(size, __FILE__, __LINE__)
#define MB_CALLOC(size) mdb_malloc(size, __FILE__, __LINE__)
#define MB_FREE(ptr)   mdb_free  (ptr)

// memory allocation template that works with STL containers
// it keeps track of who is using how much memory

#define DECLARE_MB_ALLOCATOR_TEMPLATE( classname )                             \
template < class T > class classname                                            \
{                                                                               \
  /* the real allocator that we'll use */                                       \
  static std::allocator<T> mAllocator;                                          \
public:                                                                         \
                                                                                \
  static unsigned long mNumBytesAllocated;                                      \
                                                                                \
  /* type definitions */                                                        \
  typedef T              value_type;                                            \
  typedef T*             pointer;                                               \
  typedef const T*       const_pointer;                                         \
  typedef T&             reference;                                             \
  typedef const T&       const_reference;                                       \
  typedef size_t    size_type;                                                  \
  typedef ptrdiff_t difference_type;                                            \
                                                                                \
  /* rebind allocator to type U */                                              \
  template <class U> struct rebind                                              \
  {                                                                             \
    typedef classname <U> other;                                                \
  };                                                                            \
                                                                                \
  /* return addresses of values */                                              \
  pointer address(reference value) const                                        \
  {                                                                             \
    return &value;                                                              \
  }                                                                             \
  const_pointer address(const_reference value) const                            \
  {                                                                             \
    return &value;                                                              \
  }                                                                             \
                                                                                \
  /* constructor and destructor  -- */                                          \
  /* does nothing because allocator has no state */                             \
  classname() {}                                                                \
  classname( const classname & ) {}                                             \
   /* windows doesn't like this constructor */                                  \
  /*template <class U> classname (const classname <U>&) {} */                   \
  ~classname() {}                                                               \
                                                                                \
  /* return the maximum number of elements that can be allocated */             \
  /*size_type max_size() const*/                                                \
  /*{*/                                                                         \
    /*return std::numeric_limits<std::size_t>::max() / sizeof(T);*/             \
  /*}*/                                                                         \
                                                                                \
  /* allocate but don't initialize number of elements of type T */              \
  pointer allocate(size_type num, const void* = 0)                              \
  {                                                                             \
    mNumBytesAllocated+= num*sizeof(T);                                         \
    return mAllocator.allocate(num);                                            \
    /*return (pointer)(::operator new(num*sizeof(T))); */                       \
  }                                                                             \
                                                                                \
  /* initialize elements of allocated storage p with value */                   \
  void construct(pointer p, const T& value)                                     \
  {                                                                             \
    new((void*)p)T(value);                                                      \
  }                                                                             \
                                                                                \
  /* destroy elements of initialized storage p */                               \
  void destroy(pointer p)                                                       \
  {                                                                             \
    p->~T();                                                                    \
  }                                                                             \
                                                                                \
  /* deallocate storage p of deleted elements */                                \
  void deallocate(pointer p, size_type num)                                     \
  {                                                                             \
    mNumBytesAllocated -= num*sizeof(T);                                        \
    mAllocator.deallocate(p,num);                                               \
    /*::operator delete((void*)p);*/                                            \
  }                                                                             \
                                                                                \
};                                                                              \
                                                                                \
/* return that all specializations of this allocator are interchangeable */     \
template <class T1, class T2>                                                   \
bool operator==(const classname <T1>&, const classname <T2>&)                   \
{                                                                               \
  return true;                                                                  \
}                                                                               \
template <class T1, class T2>                                                   \
bool operator!=(const classname <T1>&, const classname <T2>&)                   \
{                                                                               \
  return false;                                                                 \
}                                                                               \



// DEFINE_MB_ALLOCATOR_CLASS macros used to define
// the static members of the allocator class
// takes intermediate class name and template type

// windows and gcc want static template members
// to be instantiated for each template type
#if defined(__GNUC__) || defined(WIN32)
#define DEFINE_MB_ALLOCATOR_CLASS( classname , nametype )                      \
std::allocator<nametype> classname<nametype>::mAllocator;                       \
unsigned long classname<nametype>::mNumBytesAllocated = 0;
#else
#define DEFINE_MB_ALLOCATOR_CLASS( classname , nametype )                      \
template <class T> std::allocator<T> classname<T>::mAllocator;                  \
template <class T> unsigned long classname<T>::mNumBytesAllocated = 0;
#endif


// DECLARE_MB_ALLOCATOR_CLASS macro used to declare
// an allocator class that you'll use
// takes intermediate class name, template type, and name of class you want

#define DECLARE_MB_ALLOCATOR_CLASS( classname, nametype, deftype )             \
DECLARE_MB_ALLOCATOR_TEMPLATE( classname )                                     \
typedef classname < nametype > deftype;


#endif


