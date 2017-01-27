#ifndef TMR_TRIANGLE_INTERFACE_H
#define TMR_TRIANGLE_INTERFACE_H

/*
  Triangulate the region specified 
*/

void TMR_TriangulateSegments( double *pts, int npts,
                              int **tris, int ntris );

#endif // TMR_TRIANGLE_INTERFACE_H
