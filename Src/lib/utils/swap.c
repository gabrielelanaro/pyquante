#include "swap.h"
#include <string.h>

/* Swap generic function */
void swap(void *vp1, void *vp2, int size)
{
  char buffer[size];
  memcpy(buffer, vp1, size);
  memcpy(vp1,vp2,size);
  memcpy(vp2,buffer,size);
}
