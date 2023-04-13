/*
  The following algorithm generates an STL file from an octree
  with a prescribed scalar field.

  The method uses the "marching cubes" algorithm by Paul Bourke that
  was taken from:

  http://paulbourke.net/geometry/polygonise/

  The algorithm is based on look up tables that determine the
  location and orientation of the cutting planes.
*/

#include "TMR_STLTools.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const int edgeTable[256] = {
    0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
    0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a, 0x596, 0x49f,
    0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230,
    0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936,
    0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa,  0x7a6, 0x6af, 0x5a5,
    0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460, 0x569,
    0x663, 0x76a, 0x66,  0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a,
    0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453,
    0x55a, 0x256, 0x35f, 0x55,  0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
    0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc,  0xfcc,
    0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9, 0xac3, 0xbca,
    0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55,
    0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6,
    0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f,
    0x66,  0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af,
    0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0, 0xd30,
    0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636,
    0x13a, 0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895,
    0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190, 0xf00, 0xe09,
    0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a,
    0x203, 0x109, 0x0};

const int triTable[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

const int faceEdgeTable[16] = {15,  158, 61,  172, 107, 250, 89,  200,
                               199, 86,  245, 100, 163, 50,  145, 0};

const int faceTriTable[16][13] = {
    {0, 1, 2, 0, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 4, 2, 4, 1, -1, -1, -1, -1},
    {3, 0, 4, 3, 4, 5, 3, 5, 2, -1, -1, -1, -1},
    {2, 3, 7, 2, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 5, 0, 5, 6, 0, 6, 3, -1, -1, -1, -1},
    {1, 5, 4, 3, 7, 6, 4, 5, 7, 7, 5, 6, -1},
    {0, 4, 3, 3, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {3, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 7, 1, 7, 0, -1, -1, -1, -1},
    {1, 6, 4, 1, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 7, 2, 6, 5, 4, 5, 7, 7, 5, 6, -1},
    {2, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 5, 0, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
};

/*
  The data structures required for the marching cubes algorithm
*/
class Cell {
 public:
  TMRPoint p[8];
  double val[8];
};

/*
  Linearly interpolate the position where an isosurface cuts an edge
  between two vertices, each with their own scalar value
*/
TMRPoint vertex_interp(double isolevel, TMRPoint p1, TMRPoint p2, double valp1,
                       double valp2) {
  double mu;
  TMRPoint p;

  if (fabs(isolevel - valp1) <= 0.0) return (p1);
  if (fabs(isolevel - valp2) <= 0.0) return (p2);
  if (fabs(valp1 - valp2) <= 0.0) return (p1);

  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);

  return (p);
}

/*
  Given a grid cell and an isolevel, calculate the triangular facets
  required to represent the isosurface through the cell.  Return the
  number of triangular facets, the array "triangles" will be loaded
  up with the vertices at most 5 triangular facets.

  0 will be returned if the grid cell is either totally above of
  totally below the isolevel.
*/
int polygonise(Cell grid, double isolevel, TMR_STLTriangle *triangles) {
  int ntriang;
  int cubeindex;
  TMRPoint vertlist[12];

  /*
    Determine the index into the edge table which tells us which
    vertices are inside of the surface
  */
  cubeindex = 0;
  if (grid.val[0] < isolevel) cubeindex |= 1;
  if (grid.val[1] < isolevel) cubeindex |= 2;
  if (grid.val[2] < isolevel) cubeindex |= 4;
  if (grid.val[3] < isolevel) cubeindex |= 8;
  if (grid.val[4] < isolevel) cubeindex |= 16;
  if (grid.val[5] < isolevel) cubeindex |= 32;
  if (grid.val[6] < isolevel) cubeindex |= 64;
  if (grid.val[7] < isolevel) cubeindex |= 128;

  /* Cube is entirely in/out of the surface */
  if (edgeTable[cubeindex] == 0) return (0);

  /* Find the vertices where the surface intersects the cube */
  if (edgeTable[cubeindex] & 1)
    vertlist[0] =
        vertex_interp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
  if (edgeTable[cubeindex] & 2)
    vertlist[1] =
        vertex_interp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
  if (edgeTable[cubeindex] & 4)
    vertlist[2] =
        vertex_interp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3]);
  if (edgeTable[cubeindex] & 8)
    vertlist[3] =
        vertex_interp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0]);
  if (edgeTable[cubeindex] & 16)
    vertlist[4] =
        vertex_interp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
  if (edgeTable[cubeindex] & 32)
    vertlist[5] =
        vertex_interp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
  if (edgeTable[cubeindex] & 64)
    vertlist[6] =
        vertex_interp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7]);
  if (edgeTable[cubeindex] & 128)
    vertlist[7] =
        vertex_interp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4]);
  if (edgeTable[cubeindex] & 256)
    vertlist[8] =
        vertex_interp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
  if (edgeTable[cubeindex] & 512)
    vertlist[9] =
        vertex_interp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
  if (edgeTable[cubeindex] & 1024)
    vertlist[10] =
        vertex_interp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
  if (edgeTable[cubeindex] & 2048)
    vertlist[11] =
        vertex_interp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);

  /* Create the triangles */
  ntriang = 0;
  for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {
    triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i]];
    triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i + 1]];
    triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i + 2]];
    ntriang++;
  }

  return (ntriang);
}

/*
  Given a triangle, compute the face normal
*/
void compute_normal(TMR_STLTriangle tri, double n[]) {
  double a[3], b[3];
  a[0] = tri.p[1].x - tri.p[0].x;
  a[1] = tri.p[1].y - tri.p[0].y;
  a[2] = tri.p[1].z - tri.p[0].z;

  b[0] = tri.p[2].x - tri.p[0].x;
  b[1] = tri.p[2].y - tri.p[0].y;
  b[2] = tri.p[2].z - tri.p[0].z;

  // Compute the cross product
  n[0] = a[1] * b[2] - a[2] * b[1];
  n[1] = a[2] * b[0] - a[0] * b[2];
  n[2] = a[0] * b[1] - a[1] * b[0];

  double norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  n[0] /= norm;
  n[1] /= norm;
  n[2] /= norm;
}

/*
  Store the vertices on a block that correspond to each face
*/
const int face_vertex[6][4] = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
                               {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

/*
  Write out the face
*/
int face_polygonize(TMRPoint *p, double vals[], double cutoff,
                    TMR_STLTriangle *triangles) {
  TMRPoint vertlist[8];

  // Determine the nodes that are below the cutoff
  int faceindex = 0;
  if (vals[0] < cutoff) faceindex |= 1;
  if (vals[1] < cutoff) faceindex |= 2;
  if (vals[2] < cutoff) faceindex |= 4;
  if (vals[3] < cutoff) faceindex |= 8;

  if (faceEdgeTable[faceindex] == 0) return 0;

  // Expand the vertex list
  if (faceEdgeTable[faceindex] & 1) vertlist[0] = p[0];
  if (faceEdgeTable[faceindex] & 2) vertlist[1] = p[1];
  if (faceEdgeTable[faceindex] & 4) vertlist[2] = p[2];
  if (faceEdgeTable[faceindex] & 8) vertlist[3] = p[3];
  if (faceEdgeTable[faceindex] & 16)
    vertlist[4] = vertex_interp(cutoff, p[0], p[1], vals[0], vals[1]);
  if (faceEdgeTable[faceindex] & 32)
    vertlist[5] = vertex_interp(cutoff, p[1], p[2], vals[1], vals[2]);
  if (faceEdgeTable[faceindex] & 64)
    vertlist[6] = vertex_interp(cutoff, p[2], p[3], vals[2], vals[3]);
  if (faceEdgeTable[faceindex] & 128)
    vertlist[7] = vertex_interp(cutoff, p[3], p[0], vals[3], vals[0]);

  // Create the triangles
  int ntri = 0;
  for (int i = 0; faceTriTable[faceindex][i] != -1; i += 3) {
    triangles[ntri].p[0] = vertlist[faceTriTable[faceindex][i]];
    triangles[ntri].p[1] = vertlist[faceTriTable[faceindex][i + 1]];
    triangles[ntri].p[2] = vertlist[faceTriTable[faceindex][i + 2]];
    ntri++;
  }

  return ntri;
}

/*
  Keep a list of all the triangles
*/
class TriangleList {
 public:
  TriangleList(TMR_STLTriangle *tris, int ntris) {
    triangles = tris;
    len = ntris;
    max_len = ntris;
    len_incr = ntris;
  }
  TriangleList(int _max_len) {
    len = 0;
    max_len = _max_len;
    if (max_len < 100) {
      max_len = 100;
    }
    len_incr = max_len;
    triangles = new TMR_STLTriangle[max_len];
  }
  ~TriangleList() { delete[] triangles; }

  // Add a triangle to the list
  void addTriangle(TMR_STLTriangle *tri) {
    if (len >= max_len) {
      max_len += len_incr;
      TMR_STLTriangle *temp = new TMR_STLTriangle[max_len];
      memcpy(temp, triangles, len * sizeof(TMR_STLTriangle));
      delete[] triangles;
      triangles = temp;
    }
    triangles[len] = *tri;
    len++;
  }

  // Get the list of triangles
  void getTriangles(int *ntris, TMR_STLTriangle **tris) {
    *tris = triangles;
    *ntris = len;
  }

  // Write all the triangles to a list
  int writeSTLFile(const char *filename) {
    // Open the file and see if we were successful
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      return 1;
    }

    fprintf(fp, "solid topology\n");
    for (int i = 0; i < len; i++) {
      double n[3];
      compute_normal(triangles[i], n);

      // Write the facet normal
      fprintf(fp, "facet normal %e %e %e\n", n[0], n[1], n[2]);

      // Write the vertex loop
      fprintf(fp, "outer loop\n");
      for (int k = 0; k < 3; k++) {
        fprintf(fp, "vertex %e %e %e\n", triangles[i].p[k].x,
                triangles[i].p[k].y, triangles[i].p[k].z);
      }
      fprintf(fp, "endloop\nendfacet\n");
    }

    fprintf(fp, "endsolid topology\n");
    fclose(fp);

    return 0;
  }

 private:
  int len, max_len, len_incr;
  TMR_STLTriangle *triangles;
};

/*
  Write out the triangular volume elements within the STL file
*/
void add_volume(TriangleList *list, Cell *grid, double cutoff) {
  TMR_STLTriangle triangles[5];
  int ntri = polygonise(*grid, cutoff, triangles);

  // Add the triangles to the list
  for (int k = 0; k < ntri; k++) {
    list->addTriangle(&triangles[k]);
  }
}

/*
  Write out the intersection of the face with the boundaries
*/
void add_faces(TriangleList *list, Cell *grid, double cutoff, int bound[]) {
  // Loop over each of the faces
  for (int face = 0; face < 6; face++) {
    // Extract the points and values from the face
    TMRPoint p[4];
    int b[4];
    double vals[4];

    for (int i = 0; i < 4; i++) {
      p[i] = grid->p[face_vertex[face][i]];
      vals[i] = grid->val[face_vertex[face][i]];
      b[i] = bound[face_vertex[face][i]];
    }

    // Check whether this face is actually on a boundary
    if (b[0] && b[1] && b[2] && b[3]) {
      TMR_STLTriangle triangles[4];
      int ntri = face_polygonize(p, vals, cutoff, triangles);

      // Add the triangles to the list
      for (int k = 0; k < ntri; k++) {
        list->addTriangle(&triangles[k]);
      }
    }
  }
}

/*
  The following code generates an .STL file as output from the
  topology optimization problem.
*/

const int ordering_transform[] = {0, 1, 3, 2, 4, 5, 7, 6};

/**
  Create a list of STL triangles
*/
int TMR_GenerateSTLTriangles(TMROctForest *filter, TACSBVec *x, int x_offset,
                             double cutoff, TriangleList **_list) {
  // Set the return flag
  int fail = 0;

  // Get the MPI communicator
  int mpi_size, mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Retrieve the mesh order
  const int mesh_order = filter->getMeshOrder();

  // Ensure that the values are distributed so that we can access them
  // directly
  x->beginDistributeValues();
  x->endDistributeValues();

  // Get the dependent nodes and weight values
  const int *dep_ptr, *dep_conn;
  const double *dep_weights;
  filter->getDepNodeConn(&dep_ptr, &dep_conn, &dep_weights);

  // Set the maximum length of any of the block sides
  const int32_t hmax = 1 << TMR_MAX_LEVEL;

  // Create the list of Triangles
  TriangleList *list = new TriangleList(4096);

  // Get the block -> face information and the face -> block info.
  // This will be used to determine which faces lie on the boundaries
  // of the domain
  const int *block_face_conn;
  filter->getConnectivity(NULL, NULL, NULL, NULL, NULL, &block_face_conn, NULL,
                          NULL);

  const int *face_block_ptr;
  filter->getInverseConnectivity(NULL, NULL, NULL, NULL, NULL, &face_block_ptr);

  // Get the array of octants
  TMROctantArray *octants;
  int nelems;
  TMROctant *octs;
  filter->getOctants(&octants);
  octants->getArray(&octs, &nelems);

  // Get the connectivity
  const int *conn;
  filter->getNodeConn(&conn);

  // Get the nodal locations from the TMROctree object
  TMRPoint *X;
  filter->getPoints(&X);

  // Allocate space to store the node locations and levelset values
  const int bsize = x->getBlockSize();
  TacsScalar *xvars = new TacsScalar[bsize];
  TacsScalar *levelvals = new TacsScalar[mesh_order * mesh_order * mesh_order];
  TMRPoint *Xe = new TMRPoint[mesh_order * mesh_order * mesh_order];

  // Get the mesh coordinates from the filter
  for (int i = 0; i < nelems; i++) {
    // Compute the side-length of this element
    const int32_t h = 1 << (TMR_MAX_LEVEL - octs[i].level);

    // Get the block id
    int block = octs[i].block;

    // Check if the octree face lies on a boundary or not
    int octree_face_boundary[6];
    for (int k = 0; k < 6; k++) {
      int face = block_face_conn[6 * block + k];
      int nblocks = face_block_ptr[face + 1] - face_block_ptr[face];
      octree_face_boundary[k] = (nblocks == 1);
    }

    // Evaluate whether the element is actually on a face
    octree_face_boundary[0] = octree_face_boundary[0] && (octs[i].x == 0);
    octree_face_boundary[1] =
        octree_face_boundary[1] && (octs[i].x + h == hmax);
    octree_face_boundary[2] = octree_face_boundary[2] && (octs[i].y == 0);
    octree_face_boundary[3] =
        octree_face_boundary[3] && (octs[i].y + h == hmax);
    octree_face_boundary[4] = octree_face_boundary[4] && (octs[i].z == 0);
    octree_face_boundary[5] =
        octree_face_boundary[5] && (octs[i].z + h == hmax);

    // Get the local connectivity
    const int *c = &conn[mesh_order * mesh_order * mesh_order * i];

    // Loop over the nodes within the element
    for (int kk = 0; kk < mesh_order; kk++) {
      for (int jj = 0; jj < mesh_order; jj++) {
        for (int ii = 0; ii < mesh_order; ii++) {
          // Find the local index
          int index = (ii + jj * mesh_order + kk * mesh_order * mesh_order);

          // Find the local node number
          int node = filter->getLocalNodeNumber(c[index]);
          if (node >= 0) {
            // Get the node location
            Xe[index] = X[node];
          } else {
            printf(
                "TMR_GenerateBinFile: Failed at node with block: "
                "%d x %d y: %d z: %d\n",
                octs[i].block, octs[i].x, octs[i].y, octs[i].z);
          }
          if (c[index] >= 0) {
            x->getValues(1, &c[index], xvars);
            levelvals[index] = xvars[x_offset];
          } else {
            int dep = -c[index] - 1;
            levelvals[index] = 0.0;
            for (int jp = dep_ptr[dep]; jp < dep_ptr[dep + 1]; jp++) {
              x->getValues(1, &dep_conn[jp], xvars);
              levelvals[index] += dep_weights[jp] * xvars[x_offset];
            }
          }
        }
      }
    }

    // Loop over each octant in the mesh (if the mesh is higher-order)
    for (int iz = 0; iz < mesh_order - 1; iz++) {
      for (int iy = 0; iy < mesh_order - 1; iy++) {
        for (int ix = 0; ix < mesh_order - 1; ix++) {
          // Set the value of the X/Y/Z locations of the nodes
          Cell cell;

          // Keep track if any of the nodes lies on a boundary surface
          int on_boundary[8];
          int bound = 0;

          // Loop over the nodes within the element
          for (int kk = 0; kk < 2; kk++) {
            for (int jj = 0; jj < 2; jj++) {
              for (int ii = 0; ii < 2; ii++) {
                // Compute the index
                int index = ordering_transform[ii + 2 * jj + 4 * kk];

                // Compute the offset into the local mesh
                int offset = (ix + ii) + (iy + jj) * mesh_order +
                             (iz + kk) * mesh_order * mesh_order;

                if (filter->getInterpType() == TMR_BERNSTEIN_POINTS) {
                  int w = (int)((offset) / (mesh_order * mesh_order));
                  int v = (int)((offset - mesh_order * mesh_order * w) /
                                mesh_order);
                  int u = offset - mesh_order * v - mesh_order * mesh_order * w;
                  double pt[3];
                  pt[0] = -1.0 + 2.0 / (mesh_order - 1.0) * u;
                  pt[1] = -1.0 + 2.0 / (mesh_order - 1.0) * v;
                  pt[2] = -1.0 + 2.0 / (mesh_order - 1.0) * w;

                  int nweights = mesh_order * mesh_order * mesh_order;
                  double N[nweights];
                  filter->evalInterp(pt, N);
                  cell.val[index] = 0.0;
                  for (int nn = 0; nn < nweights; nn++) {
                    cell.val[index] += N[nn] * levelvals[nn];
                  }
                } else {
                  cell.val[index] = levelvals[offset];
                }

                // Set the node location
                cell.p[index].x = Xe[offset].x;
                cell.p[index].y = Xe[offset].y;
                cell.p[index].z = Xe[offset].z;

                int fx0 = (ix == 0 && ii == 0);
                int fx1 = (ix == mesh_order - 2 && ii == 1);
                int fy0 = (iy == 0 && jj == 0);
                int fy1 = (iy == mesh_order - 2 && jj == 1);
                int fz0 = (iz == 0 && kk == 0);
                int fz1 = (iz == mesh_order - 2 && kk == 1);

                // Check if this node lies on a boundary
                on_boundary[index] = ((octree_face_boundary[0] && fx0) ||
                                      (octree_face_boundary[1] && fx1) ||
                                      (octree_face_boundary[2] && fy0) ||
                                      (octree_face_boundary[3] && fy1) ||
                                      (octree_face_boundary[4] && fz0) ||
                                      (octree_face_boundary[5] && fz1));

                // Store whether this node is on the boundary
                bound = bound || on_boundary[index];
              }
            }
          }

          // Write out the surface
          add_volume(list, &cell, cutoff);

          // Write out the faces if they are on the boundary
          if (bound) {
            add_faces(list, &cell, cutoff, on_boundary);
          }
        }
      }
    }
  }

  // Free xvals
  delete[] xvars;
  delete[] levelvals;
  delete[] Xe;

  *_list = list;

  return fail;
}

int TMR_GenerateSTLTriangles(int root, TMROctForest *filter, TACSBVec *x,
                             int x_offset, double cutoff, int *_ntris,
                             TMR_STLTriangle **_tris) {
  // Generate the triangle
  TriangleList *list;
  int fail = TMR_GenerateSTLTriangles(filter, x, x_offset, cutoff, &list);
  if (fail) {
    return fail;
  }

  // Get the local triangles
  int ntris;
  TMR_STLTriangle *tris;
  list->getTriangles(&ntris, &tris);

  // Get the MPI communicator
  int mpi_size, mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  int *counts = new int[mpi_size];
  int *range = new int[mpi_size + 1];
  MPI_Allgather(&ntris, 1, MPI_INT, counts, 1, MPI_INT, comm);

  range[0] = 0;
  for (int i = 0; i < mpi_size; i++) {
    range[i + 1] = counts[i] + range[i];
  }

  // Allocate the receiver only on the root processor
  *_ntris = 0;
  *_tris = NULL;
  if (mpi_rank == root) {
    *_ntris = range[mpi_size];
    *_tris = new TMR_STLTriangle[*_ntris];
  }

  MPI_Gatherv(tris, ntris, TMR_STLTriangle_MPI_type, *_tris, counts, range,
              TMR_STLTriangle_MPI_type, root, comm);
  delete[] range;
  delete[] counts;
  delete list;

  return fail;
}

/*
  Generate a binary file containing the STL triangles
*/
int TMR_GenerateBinFile(const char *filename, TMROctForest *filter, TACSBVec *x,
                        int x_offset, double cutoff) {
  // Generate the triangle
  TriangleList *list;
  int fail = TMR_GenerateSTLTriangles(filter, x, x_offset, cutoff, &list);
  if (fail) {
    return fail;
  }

  // Get the MPI communicator
  int mpi_size, mpi_rank;
  MPI_Comm comm = filter->getMPIComm();
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Get the local triangles
  int ntris;
  TMR_STLTriangle *tris;
  list->getTriangles(&ntris, &tris);

  // Copy the filename to a non-const array
  char *fname = new char[strlen(filename) + 1];
  strcpy(fname, filename);

  int *range = new int[mpi_size + 1];
  range[0] = 0;
  MPI_Allgather(&ntris, 1, MPI_INT, &range[1], 1, MPI_INT, comm);
  for (int i = 0; i < mpi_size; i++) {
    range[i + 1] += range[i];
  }

  // Create the file and write out the information
  MPI_File fp = NULL;
  MPI_File_open(comm, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                &fp);

  if (fp) {
    // Write out the integer indicating how many triangles there are
    if (mpi_rank == 0) {
      MPI_File_write(fp, &range[mpi_size], 1, MPI_INT, MPI_STATUS_IGNORE);
    }

    // Write out all the triangles to the file
    char datarep[] = "native";
    MPI_File_set_view(fp, sizeof(int), TMR_STLTriangle_MPI_type,
                      TMR_STLTriangle_MPI_type, datarep, MPI_INFO_NULL);
    MPI_File_write_at_all(fp, range[mpi_rank], tris, ntris,
                          TMR_STLTriangle_MPI_type, MPI_STATUS_IGNORE);
    MPI_File_close(&fp);
  } else {
    fail = 1;
  }

  delete[] range;
  delete[] fname;
  delete list;

  return fail;
}

/*
  Take the binary file generated from above and convert to the .STL
  data format (in ASCII).
*/
int TMR_ConvertBinToSTL(const char *binfile, const char *stlfile) {
  FILE *fp = fopen(binfile, "rb");
  if (!fp) {
    return 1;
  }

  // Try to read in the number of triangles
  int ntris;
  if (fread(&ntris, sizeof(int), 1, fp) != 1) {
    return 1;
  }

  // Read in the triangles themselves
  TMR_STLTriangle *tris = new TMR_STLTriangle[ntris];
  unsigned int unsigned_ntris = ntris;
  if (fread(tris, sizeof(TMR_STLTriangle), ntris, fp) != unsigned_ntris) {
    delete[] tris;
    return 1;
  }

  TriangleList *list = new TriangleList(tris, ntris);
  list->writeSTLFile(stlfile);
  delete list;

  return 0;
}
