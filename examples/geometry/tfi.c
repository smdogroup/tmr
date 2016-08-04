#include "TMRGeometry.h"
#include "TMROctForest.h"

/*
  The following code is used to generate geometry for a problem
*/

void evalPoint( int32_t u_int, int32_t v_int, int32_t w_int,
                TMRPoint *pt ){
  const double u = 1.0*u_int/(1 << TMR_MAX_LEVEL);
  const double v = 1.0*v_int/(1 << TMR_MAX_LEVEL);
  const double w = 1.0*w_int/(1 << TMR_MAX_LEVEL);

  // Evaluate the x/y/z locations and store them in pt
  pt->x = 1.0*u;
  pt->y = 1.0*v;
  pt->z = 1.0*w*(u + 0.25);
}

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  MPI_Comm comm = MPI_COMM_WORLD;

  // Create the edge look up table
  int conn[] = {0, 1, 2, 3, 4, 5, 6, 7};

  // Create the forest and set the connectivity
  TMROctForest *forest = new TMROctForest(comm);
  forest->setConnectivity(8, conn, 1);

  // Create a random tree and the corresponding nodes/
  forest->createRandomTrees(100, 0, 15);
  forest->balance(1);
  forest->createNodes(2);
  
  // Get the parametric node locations
  int nlevel = 3;
  int nvals = (1 << nlevel) + 1;

  // Set the u-values along the x direction
  int32_t *params = new int32_t[ nvals ];    
  int32_t increment = (1 << (TMR_MAX_LEVEL - nlevel));
  int32_t loc = 0;
  for ( int k = 0; loc <= (1 << TMR_MAX_LEVEL); k++, loc += increment ){
    params[k] = loc;
  }

  // Create the points
  TMRPoint corners[8];
  const int32_t hmax = 1 << TMR_MAX_LEVEL;
  for ( int k = 0; k < 8; k++ ){
    evalPoint(hmax*(k%2), hmax*((k%4)/2), hmax*(k/4), &corners[k]);
  }

  // Set the look up values along the edges of the block
  TMR_EdgeLookup *edges[12];
  TMR_SurfaceLookup *faces[6];
  for ( int k = 0; k < 12; k++ ){
    int32_t t = 0;
    int32_t u = 0, v = 0, w = 0;

    // Compute the fixed u/v/w coordinates
    if (k < 4){
      v = hmax*(k % 2);
      w = hmax*(k/2);
    }
    else if (k < 8){
      u = hmax*(k % 2);
      w = hmax*((k-4)/2);
    }
    else {
      u = hmax*(k % 2);
      v = hmax*((k-8)/2);
    }

    // Allocate the array of points
    TMRPoint *pts = new TMRPoint[ nvals ];
    for ( int i = 0; i < nvals; i++ ){
      if (k < 4){
        evalPoint(params[i], v, w, &pts[i]);
      }
      else if (k < 8){
        evalPoint(u, params[i], w, &pts[i]);
      }
      else {
        evalPoint(u, v, params[i], &pts[i]);
      }
    }
    
    // Set the u/v/w locations
    edges[k] = new TMR_EdgeLookup(params, pts, nvals);
    delete [] pts;    
  }

  // Set the parameter locations along the faces
  int32_t *uparams = new int32_t[ nvals*nvals ];
  int32_t *vparams = new int32_t[ nvals*nvals ];
  for ( int k = 0; k < 6; k++ ){    
    // Set the face points on the surface
    TMRPoint *pts = new TMRPoint[ nvals*nvals ];
    for ( int indx = 0, j = 0; j < nvals; j++ ){
      for ( int i = 0; i < nvals; i++, indx++ ){
        uparams[indx] = params[i];
        vparams[indx] = params[j];

        if (k < 2){
          evalPoint(hmax*(k % 2), uparams[indx], vparams[indx], &pts[indx]);
        }
        else if (k < 4){
          evalPoint(uparams[indx], hmax*(k % 2), vparams[indx], &pts[indx]);
        }
        else {
          evalPoint(uparams[indx], vparams[indx], hmax*(k % 2), &pts[indx]);
        }
      }
    }

    // Set the face points on the surface
    faces[k] = new TMR_SurfaceLookup(uparams, vparams, pts, nvals*nvals);
  }

  // Create the volume parametrization
  TMR_TFIVolume *volume = new TMR_TFIVolume(corners, edges, faces);

  // Write out a file for each processor - bad practice!
  char filename[128];
  sprintf(filename, "volume.dat");
  FILE *fp = fopen(filename, "w");

  // Write the tecplot header
  fprintf(fp, "Variables = X, Y, Z\n");
  fprintf(fp, "ZONE T=TMR N=%d E=%d ", nvals*nvals*nvals, (nvals-1)*(nvals-1)*(nvals-1));
  fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

  for ( int k = 0; k < nvals; k++ ){
    for ( int j = 0; j < nvals; j++ ){
      for ( int i = 0; i < nvals; i++ ){
        // Set the parameter locations
        int32_t u = params[i];
        int32_t v = params[j];
        int32_t w = params[k];

        // Set the point value
        TMRPoint pt;
        int fail = volume->evalPoint(u, v, w, &pt);
        if (fail){
          printf("Missing point %d %d %d\n", u, v, w);
        }

        // Print the result to the file
        fprintf(fp, "%e %e %e\n", pt.x, pt.y, pt.z);
      }
    }
  }

  for ( int k = 0; k < nvals-1; k++ ){
    for ( int j = 0; j < nvals-1; j++ ){
      for ( int i = 0; i < nvals-1; i++ ){
        fprintf(fp, "%d %d %d %d  %d %d %d %d\n",
                i+1 + j*nvals + k*nvals*nvals, i+2 + j*nvals + k*nvals*nvals,
                i+2 + (j+1)*nvals + k*nvals*nvals, i+1 + (j+1)*nvals + k*nvals*nvals,
                i+1 + j*nvals + (k+1)*nvals*nvals, i+2 + j*nvals + (k+1)*nvals*nvals,
                i+2 + (j+1)*nvals + (k+1)*nvals*nvals, i+1 + (j+1)*nvals + (k+1)*nvals*nvals);
      }
    }
  }  

  fclose(fp);

  delete volume;
  delete forest;
  TMRFinalize();
  MPI_Finalize();
  return (0);
}

