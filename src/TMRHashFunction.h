#ifndef TMR_HASH_FUNCTION_H
#define TMR_HASH_FUNCTION_H

/*
  The following file defines has functions used at various
  points in the TMR code

  This code is from the public domain code located here:

  http://burtleburtle.net/bob/c/lookup3.c
 */

#define TMR_ROT_INT(x,k) (((x) << (k)) | ((x) >> (32 - (k))))
#define TMR_MIX_INT(a,b,c) ((a -= c, a ^= TMR_ROT_INT(c, 4), c += b, \
                             b -= a, b ^= TMR_ROT_INT(a, 6), a += c, \
                             c -= b, c ^= TMR_ROT_INT(b, 8), b += a, \
                             a -= c, a ^= TMR_ROT_INT(c,16), c += b, \
                             b -= a, b ^= TMR_ROT_INT(a,19), a += c, \
                             c -= b, c ^= TMR_ROT_INT(b, 4), b += a))
#define TMR_FINAL_HASH(a,b,c) ((c ^= b, c -= TMR_ROT_INT(b,14), \
                                a ^= c, a -= TMR_ROT_INT(c,11), \
                                b ^= a, b -= TMR_ROT_INT(a,25), \
                                c ^= b, c -= TMR_ROT_INT(b,16), \
                                a ^= c, a -= TMR_ROT_INT(c, 4), \
                                b ^= a, b -= TMR_ROT_INT(a,14), \
                                c ^= b, c -= TMR_ROT_INT(b,24)))

/*
  Create a hash value for pairs of usigned integers
*/
inline uint32_t TMRIntegerPairHash( uint32_t u, uint32_t v ){
  uint32_t w = 0;
  TMR_MIX_INT(u, v, w);
  return TMR_FINAL_HASH(u, v, w);
}

/*
  Create a hash value for triplets of unsigned integers
*/
inline uint32_t TMRIntegerTripletHash( uint32_t u, uint32_t v, uint32_t w ){
  TMR_MIX_INT(u, v, w);
  return TMR_FINAL_HASH(u, v, w);
}

/*
  Create a hash value for four-tuples of unsigned integers
*/
inline uint32_t TMRIntegerFourTupleHash( uint32_t u, uint32_t v,
                                         uint32_t w, uint32_t x ){
  TMR_MIX_INT(u, v, w);
  u += x;
  return TMR_FINAL_HASH(u, v, w);
}

/*
  Create a hash value for four-tuples of unsigned integers
*/
inline uint32_t TMRIntegerFiveTupleHash( uint32_t u, uint32_t v,
                                         uint32_t w, uint32_t x,
                                         uint32_t z ){
  TMR_MIX_INT(u, v, w);
  u += x;
  v += z;
  return TMR_FINAL_HASH(u, v, w);
}

#endif // TMR_HASH_FUNCTION_H
