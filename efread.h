/* versions of fread and fwrite that swap bytes on little-endian machines */

#ifndef EFREAD_H
#define EFREAD_H

#include <stdio.h>

#ifdef LITTLEENDIAN

size_t efread(void *ptr, size_t size, size_t nitems, FILE *stream);
size_t efwrite(const void *ptr, size_t size, size_t nitems, FILE *stream);

#else  /* ! LITTLEENDIAN */

#define efread fread
#define efwrite fwrite

#endif /* LITTLEENDIAN */

template<class T>
size_t tfread( T &a, FILE *stream)
{
  return efread( &a, sizeof(T), 1, stream);
}

template<class T>
size_t tfwrite( const T &a, FILE *stream)
{
  return efwrite( &a, sizeof(T), 1, stream);
}  

#endif /* MYFREAD_H */
