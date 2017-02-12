#ifndef TMR_GEOMETRIC_PREDICATES_H
#define TMR_GEOMETRIC_PREDICATES_H

/*
  Header file for Shewchuk's geometric predicates
*/

TMR_EXTERN_C_BEGIN
double orient2d( double pa[], double pb[], double pc[] );
double incircle( double pa[], double pb[], double pc[], double pd[] );
TMR_EXTERN_C_END

#endif // TMR_GEOMETRIC_PREDICATES_H
