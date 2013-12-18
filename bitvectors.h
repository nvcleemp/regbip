/*
 * Main developer: Nicolas Van Cleemput
 * 
 * Copyright (C) 2013 Nicolas Van Cleemput.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 */

#ifndef BITVECTORS_H
#define	BITVECTORS_H

/** Use of bitvectors */
typedef unsigned int bitv;
typedef unsigned int bitv_size;
#define ZERO 0U
#define ONE 1U
#define EMPTY_SET 0U
#define SINGLETON(el) (ONE << (el))
#define CONTAINS(s, el) ((s) & SINGLETON(el))
#define ADD(s, el) ((s) |= SINGLETON(el))
#define ADD_ALL(s, elements) ((s) |= (elements))
//these will only work of the element is actually in the set
#define REMOVE(s, el) ((s) ^= SINGLETON(el))
#define REMOVE_ALL(s, elements) ((s) ^= (elements))
//the following macros perform an extra step, but will work even if the element is not in the set
#define SAFE_REMOVE(s, el) ADD(s, el); REMOVE(s, el)
#define SAFE_REMOVE_ALL(s, elements) ADD_ALL(s, elements); REMOVE_ALL(s, elements)

#define SETSIZE(s, count) {bitv bitcountTemp = (s); \
bitcountTemp = bitcountTemp - ((bitcountTemp >> 1) & 0x55555555); \
bitcountTemp = (bitcountTemp & 0x33333333) + ((bitcountTemp >> 2) & 0x33333333);\
count = (((bitcountTemp + (bitcountTemp >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;}

#endif	/* BITVECTORS_H */

