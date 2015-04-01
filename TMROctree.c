#include "TMROctree.h"

/*
  Refine the initial octree to the specified depth along all
  coordinate directions 
*/
TMROctree::TMROctree( int refine_level ){
  // Check that the refinement lies within legal bounds
  if (refine_level < 0){
    refine_level = 0;
  }
  else if (refine_level-1 > TMR_MAX_LEVEL){
    refine_level = TMR_MAX_LEVEL-1;
  }

  // Compute the number of octants along each edge in
  // the leaves
  int h = 1 << (TMR_MAX_LEVEL - refine_level);

  // Set the number of octants created
  int nocts = h*h*h;

  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nocts ];

  // Create all the octants required
  int index = 0;
  for ( int z = 0; z < h; z++ ){
    for ( int y = 0; y < h; y++ ){
      for ( int x = 0; x < h; x++ ){
	array[index].x = x;
	array[index].y = y;
	array[index].z = z;
	array[index].level = refine_level;
	index++;
      }
    }
  }

  list = new TMROctantArray(array, nocts);
  list->sort();
  is_sorted = 1;
}

/*
  Generate a random initial octree that can be used for
  testing.
*/
TMROctree::TMROctree( int nrand, int min_level, int max_level ){
  // Create an array of the octants that will be stored
  TMROctant *array = new TMROctant[ nrand ];

  // Generate a random number of octants along random directions
  for ( int i = 0; i < nrand; i++ ){
    int level = min_level + (rand() % (max_level - min_level + 1));

    int h = 1 << (TMR_MAX_LEVEL - level);
    int x = h*(rand() % (1 << level));
    int y = h*(rand() % (1 << level));
    int z = h*(rand() % (1 << level));

    array[i].x = x;
    array[i].y = y;
    array[i].z = z;
    array[i].level = level;
  }

  list = new TMROctantArray(array, nrand);
  list->sortUnique();
  is_sorted = 1;
}

/*
  Create a tree from an array
*/
TMROctree::TMROctree( TMROctantArray *_list, int _is_sorted ){
  list = _list;
  is_sorted = _is_sorted;
}

/*
  Print out the octree to a file for visualization
*/
void TMROctree::printTree( const char * filename ){
  // Iterate through and print out everything  
  FILE * fp = fopen(filename, "w");
  if (fp){
    int size;
    TMROctant *array;
    list->getArray(&array, &size);

    fprintf(fp, "Variables = X, Y, Z\n");
    fprintf(fp, "ZONE T=TMR N=%d E=%d ", 8*size, size);
    fprintf(fp, "DATAPACKING=POINT ZONETYPE=FEBRICK\n");

    // Set the length of the cube
    double dh = 1.0/(1 << TMR_MAX_LEVEL);

    for ( int i = 0; i < size; i++ ){
      int h = 1 << (TMR_MAX_LEVEL - array[i].level);
      int x = array[i].x;
      int y = array[i].y;
      int z = array[i].z;

      fprintf(fp, "%e %e %e\n", x*dh, y*dh, z*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, y*dh, z*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, (y+h)*dh, z*dh);
      fprintf(fp, "%e %e %e\n", x*dh, (y+h)*dh, z*dh);

      fprintf(fp, "%e %e %e\n", x*dh, y*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, y*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", (x+h)*dh, (y+h)*dh, (z+h)*dh);
      fprintf(fp, "%e %e %e\n", x*dh, (y+h)*dh, (z+h)*dh);
    }

    for ( int i = 0; i < size; i++ ){
      for ( int k = 0; k < 8; k++ ){
	fprintf(fp, "%d ", 8*i+k+1);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
}

/*
  Balance the octree in place by scanning across  using a simple algorithm 
*/
void TMROctree::balance( int balance_type ){
  TMROctantHash *hash = new TMROctantHash();

  // Go through the existing list of octants and add up everything for
  // balancing
  TMROctant p, q, neighbor;

  int size;
  TMROctant *array;
  list->getArray(&array, &size);

  TMROctantQueue *queue = new TMROctantQueue();

  for ( int i = 0; i < size; i++ ){
    // Try adding all of the children
    array[i].getSibling(0, &q);
    hash->addOctant(&q);
    
    // Get the parent, and then add the face matching
    // octants around them and their siblings
    if (q.level > 0){
      q.parent(&p);

      for ( int face = 0; face < 6; face++ ){
	p.faceNeighbor(face, &neighbor);
	neighbor.getSibling(0, &q);
	
	int hmax = 1 << TMR_MAX_LEVEL;
	// If we're in bounds, add the neighbor
	if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	    q.x < hmax && q.y < hmax && q.z < hmax){
	  if (hash->addOctant(&q)){
	    queue->push(&q);
	  }
	}
      }
      
      if (balance_type & 2){
	for ( int edge = 0; edge < 12; edge++ ){
	  p.edgeNeighbor(edge, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  int hmax = 1 << TMR_MAX_LEVEL;
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
      if ((balance_type & 2) && (balance_type & 4)){
	for ( int corner = 0; corner < 8; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  int hmax = 1 << TMR_MAX_LEVEL;
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  while (queue->length() > 0){
    q = queue->pop();
   
    if (q.level > 1){
      q.parent(&p);

      for ( int face = 0; face < 6; face++ ){
	p.faceNeighbor(face, &neighbor);
	neighbor.getSibling(0, &q);
	
	int hmax = 1 << TMR_MAX_LEVEL;
	// If we're in bounds, add the neighbor
	if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	    q.x < hmax && q.y < hmax && q.z < hmax){
	  if (hash->addOctant(&q)){
	    queue->push(&q);
	  }
	}
      }

      if (balance_type & 2){
	for ( int edge = 0; edge < 12; edge++ ){
	  p.edgeNeighbor(edge, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  int hmax = 1 << TMR_MAX_LEVEL;
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
      if ((balance_type & 2) && (balance_type & 4)){
	for ( int corner = 0; corner < 8; corner++ ){
	  p.cornerNeighbor(corner, &neighbor);
	  neighbor.getSibling(0, &q);
	  
	  int hmax = 1 << TMR_MAX_LEVEL;
	  // If we're in bounds, add the neighbor
	  if (q.x >= 0 && q.y >= 0 && q.z >= 0 && 
	      q.x < hmax && q.y < hmax && q.z < hmax){
	    if (hash->addOctant(&q)){
	      queue->push(&q);
	    }
	  }
	}
      }
    }
  }

  TMROctantArray *list0 = hash->toArray();
  list0->getArray(&array, &size);

  for ( int i = 0; i < size; i++ ){
    for ( int j = 0; j < 8; j++ ){
      array[i].getSibling(j, &q);
      hash->addOctant(&q);
    }
  }

  list = hash->toArray();
  list->sortUnique();
}


TMROctree *TMROctree::coarsen(){
  list->sortUnique();

  int size;
  TMROctant *array;
  list->getArray(&array, &size);

  // Set the offset to be 2**d-1
  int offset = (1 << 3) - 1;

  TMROctantQueue *queue = new TMROctantQueue();
  int coarsened = 0;

  // Scan through the list, if we have offset octants which
  // all share the same parent, then coarsen the parent
  for ( int i = 0; i < size; i++ ){
    int same_parent = 0;
    if (array[i].level > 0 && 
	array[i].childId() == 0 && 
	i+offset < size &&
	array[i+offset].childId() == offset){

      // Test to see if the children have the same parent
      TMROctant p;
      array[i+offset].getSibling(0, &p);
      if (array[i].compare(&p) == 0){
	same_parent = 1;
	array[i].parent(&p);
	queue->push(&p);

	i += offset;
	coarsened++;
      }
    }

    if (!same_parent){
      queue->push(&array[i]);
    }
  }

  printf("coarsend = %d\n", coarsened);

  TMROctree * tree = new TMROctree(queue->toArray(), 1);
  delete queue;

  return tree;
}
