/**
 * swap.h
 * swap utility functions
 */

#ifndef _SWAP_H_
#define _SWAP_H_
/**
 * Generic swap function, can swap everything and should be used in this way:
 *
 * swap(&a,&b,sizeof(a));
 */
void swap(void *vp1, void *vp2, int size);

#endif /* _SWAP_H_ */
