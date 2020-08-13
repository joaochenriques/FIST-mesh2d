/* versions of fread and fwrite that swap bytes on little-endian machines */

#include <stdio.h>
#include "stdlib.h"
#include <string.h>
#include "efread.h"

#ifdef LITTLEENDIAN

/* copy nitems items of size size from oldbuf to newbuf, reversing the
 * byte order in each item.
 */

static void copyswap(size_t size, size_t nitems, const unsigned char *oldbuf,
	      unsigned char *newbuf)
{
  int i, top, j;
  int total = nitems*size;
  for(i=0; i<total; i+=size) {
    top = i+size-1;
    for(j=0; j<size; j++)
      newbuf[i+j] = oldbuf[top-j];
  }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

size_t efread(void *ptr, size_t size, size_t nitems, FILE *stream) {
  unsigned char *buf;
  size_t status;

  if(size == 1)
    buf = (unsigned char*) ptr;	// don't have to use temp space
  else
    buf = new unsigned char[size*nitems]; // temp space

  status = fread(buf, size, nitems, stream);

  if(size != 1) {
    copyswap(size, nitems, buf, (unsigned char*) ptr);
    delete buf;
  }

  return status;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/

size_t efwrite(const void *ptr, size_t size, size_t nitems, FILE *stream) {
  int status;
  unsigned char *buf;

  if(size == 1)
    buf = (unsigned char *) ptr;
  else {
    buf = new unsigned char[size*nitems]; // temp space
    copyswap(size, nitems, (unsigned char*) ptr, buf);
  }
  status = fwrite(buf, size, nitems, stream);
  if(size != 1)
    delete buf;
  return status;
}

#endif
