#include "TMRQuadrant.h"
#include "TMRHashFunction.h"

/*
  Get the child id of the quadrant
*/
int TMRQuadrant::childId(){
  int id = 0;
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);
  
  id = id | ((x & h) ? 1 : 0);
  id = id | ((y & h) ? 2 : 0);

  return id;
}

/*
  Return the sibling of the quadrant
*/
void TMRQuadrant::getSibling( int id, TMRQuadrant *sib ){
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);

  int32_t xr = ((x & h) ? x-h : x);
  int32_t yr = ((y & h) ? y-h : y);

  sib->face = face;
  sib->level = level;
  sib->info = 0;
  sib->x = ((id & 1) ? xr+h : xr);
  sib->y = ((id & 2) ? yr+h : yr);
}

/*
  Get the parent of the quadrant
*/
void TMRQuadrant::parent( TMRQuadrant *p ){
  if (level > 0){
    p->face = face;
    p->level = level-1;
    p->info = 0;
    const int32_t h = 1 << (TMR_MAX_LEVEL - level);

    p->x = x & ~h;
    p->y = y & ~h;
  }
  else {
    p->face = face;
    p->level = 0;
    p->info = 0;
    p->x = x;
    p->y = y;
  }
}

/*
  Get the edge neighbour
*/
void TMRQuadrant::edgeNeighbor( int edge, TMRQuadrant *neighbor ){
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);
  neighbor->face = face;
  neighbor->level = level;
  neighbor->info = 0;

  neighbor->x = x + ((edge == 0) ? -h : (edge == 1) ? h : 0);
  neighbor->y = y + ((edge == 2) ? -h : (edge == 3) ? h : 0);
}

/*
  Get the neighbor along a corner
*/
void TMRQuadrant::cornerNeighbor( int corner, TMRQuadrant *neighbor ){
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);
  neighbor->face = face;
  neighbor->level = level;
  neighbor->info = 0;

  neighbor->x = x + (2*(corner & 1) - 1)*h;
  neighbor->y = y + ((corner & 2) - 1)*h;
}

/*
  Compare two quadrants using the Morton enconding (z-ordering)

  This function returns -1 if self < quadrant, 0 if self == quadrant and 1
  if self > quadrant. Ties are broken by the level of the quadrant such
  that the quadrants will be sorted by location then level.
*/
int TMRQuadrant::compare( const TMRQuadrant *quadrant ) const {
  if (face != quadrant->face){
    return face - quadrant->face;
  }

  uint32_t xxor = x ^ quadrant->x;
  uint32_t yxor = y ^ quadrant->y;
  uint32_t sor = xxor | yxor;

  // If there is no most-significant bit, then we are done
  if (sor == 0){
    return level - quadrant->level;
  }

  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = x - quadrant->x;
  }
  else {
    discrim = y - quadrant->y;
  }

  if (discrim > 0){
    return 1;
  }
  else if (discrim < 0){
    return -1;
  }
  return 0;
}

/*
  Compare two quadrants to determine whether they have the same Morton
  encoding, but may be at different levels.
*/
int TMRQuadrant::comparePosition( const TMRQuadrant *quadrant ) const {
  if (face != quadrant->face){
    return face - quadrant->face;
  }

  uint32_t xxor = x ^ quadrant->x;
  uint32_t yxor = y ^ quadrant->y;
  uint32_t sor = xxor | yxor;

  // Note that here we do not distinguish between levels
  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = x - quadrant->x;
  }
  else {
    discrim = y - quadrant->y;
  }

  if (discrim > 0){
    return 1;
  }
  else if (discrim < 0){
    return -1;
  }
  return 0;
}

/*
  Compare two quadrants to determine whether they have the same Morton
  encoding, but may be at different levels.
*/
int TMRQuadrant::compareNode( const TMRQuadrant *quadrant ) const {
  if (face != quadrant->face){
    return face - quadrant->face;
  }

  uint32_t xxor = x ^ quadrant->x;
  uint32_t yxor = y ^ quadrant->y;
  uint32_t sor = xxor | yxor;

  // Note that here we do not distinguish between levels
  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = x - quadrant->x;
  }
  else {
    discrim = y - quadrant->y;
  }

  if (discrim > 0){
    return 1;
  }
  else if (discrim < 0){
    return -1;
  }
  return info - quadrant->info;
}

/*
  Determine whether the input quadrant is contained within the
  quadrant itself. This can be used to determine whether the given
  quadrant is a descendent of this object.
*/
int TMRQuadrant::contains( TMRQuadrant *quad ){
  const int32_t h = 1 << (TMR_MAX_LEVEL - level);

  // Check whether the quadrant lies within this quadrant
  if ((quad->face == face) && 
      (quad->x >= x && quad->x < x + h) &&
      (quad->y >= y && quad->y < y + h)){
    return 1;
  }
  
  return 0;
}

/*
  Compare two quadrants within the same sub-tree
*/
static int compare_quadrants( const void *a, const void *b ){
  const TMRQuadrant *ao = static_cast<const TMRQuadrant*>(a);
  const TMRQuadrant *bo = static_cast<const TMRQuadrant*>(b);
  
  return ao->compare(bo);
}

/*
  Compare two quadrant positions
*/
static int compare_position( const void *a, const void *b ){
  const TMRQuadrant *ao = static_cast<const TMRQuadrant*>(a);
  const TMRQuadrant *bo = static_cast<const TMRQuadrant*>(b);
  
  return ao->comparePosition(bo);
}

/*
  Compare two quadrant nodes
*/
static int compare_nodes( const void *a, const void *b ){
  const TMRQuadrant *ao = static_cast<const TMRQuadrant*>(a);
  const TMRQuadrant *bo = static_cast<const TMRQuadrant*>(b);
  
  return ao->compareNode(bo);
}

/*
  Store a array of quadrants
*/
TMRQuadrantArray::TMRQuadrantArray( TMRQuadrant *_array, int _size,
                                    int _use_node_index ){
  array = _array;
  size = _size;
  max_size = size;
  is_sorted = 0;
  use_node_index = _use_node_index;
}

/*
  Destroy the quadrant array
*/
TMRQuadrantArray::~TMRQuadrantArray(){
  delete [] array;
}

/*
  Duplicate the array and return the copy
*/
TMRQuadrantArray* TMRQuadrantArray::duplicate(){
  TMRQuadrant *arr = new TMRQuadrant[ size ];
  memcpy(arr, array, size*sizeof(TMRQuadrant));

  TMRQuadrantArray *dup = new TMRQuadrantArray(arr, size);
  dup->is_sorted = is_sorted;
  dup->use_node_index = use_node_index;

  return dup;
}

/*
  Sort the list and remove duplicates from the array of possible
  entries.
*/
void TMRQuadrantArray::sort(){
  if (use_node_index){
    qsort(array, size, sizeof(TMRQuadrant), compare_nodes);

    // Now that the Quadrants are sorted, remove duplicates
    int i = 0; // Location from which to take entries
    int j = 0; // Location to place entries
    
    for ( ; i < size; i++, j++ ){
      while ((i < size-1) && 
             (array[i].compareNode(&array[i+1]) == 0)){
        i++;
      }

      if (i != j){
        array[j] = array[i];
      }
    }

    // The new size of the array
    size = j;
  }
  else {
    qsort(array, size, sizeof(TMRQuadrant), compare_quadrants);

    // Now that the Quadrants are sorted, remove duplicates
    int i = 0; // Location from which to take entries
    int j = 0; // Location to place entries
    
    for ( ; i < size; i++, j++ ){
      while ((i < size-1) && 
             (array[i].comparePosition(&array[i+1]) == 0)){
        i++;
      }

      if (i != j){
        array[j] = array[i];
      }
    }

    // The new size of the array
    size = j;
  }

  is_sorted = 1;
}

/*
  Determine if the array contains the specified quadrant
*/
TMRQuadrant* TMRQuadrantArray::contains( TMRQuadrant *q,
                                         const int use_position ){
  if (!is_sorted){
    is_sorted = 1;
    sort();
  }

  if (use_node_index){
    return (TMRQuadrant*)bsearch(q, array, size, sizeof(TMRQuadrant), 
                                 compare_nodes);
  }
  else {
    // Search for nodes - these will share the same
    if (use_position){
      return (TMRQuadrant*)bsearch(q, array, size, sizeof(TMRQuadrant), 
                                   compare_position);
    }

    // Search the array for an identical element
    return (TMRQuadrant*)bsearch(q, array, size, sizeof(TMRQuadrant), 
                                 compare_quadrants);
  }
}

/*
  Merge the entries of two arrays
*/
void TMRQuadrantArray::merge( TMRQuadrantArray * list ){
  if (!is_sorted){
    sort();
  }
  if (!list->is_sorted){
    list->sort();
  }

  // Keep track of the number of duplicates
  int nduplicates = 0;

  // Scan through the list and determine the number
  // of duplicates 
  int j = 0, i = 0;
  for ( ; i < size; i++ ){
    while ((j < list->size) && 
           (list->array[j].compare(&array[i]) < 0)){
      j++;
    }
    if (j >= list->size){
      break;
    }
    if (array[i].compare(&list->array[j]) == 0){
      nduplicates++;
    }
  }

  // Compute the required length of the new array
  int len = size + list->size - nduplicates; 

  // Allocate a new array if required
  if (len > max_size){
    max_size = len;

    TMRQuadrant *temp = array;
    array = new TMRQuadrant[ max_size ];
    memcpy(array, temp, size*sizeof(TMRQuadrant));

    // Free the old array
    delete [] temp;
  }

  // Set the pointer to the end of the new array
  int end = len-1;

  // Copy the new array back
  i = size-1;
  j = list->size-1;
  while (i >= 0 && j >= 0){
    if (array[i].compare(&list->array[j]) > 0){
      array[end] = array[i];
      end--, i--;
    }
    else if (list->array[j].compare(&array[i]) > 0){
      array[end] = list->array[j];
      end--, j--;
    }
    else { // b[j] == a[i]
      array[end] = array[i];
      end--, j--, i--;
    }      
  }
    
  // Only need to copy over remaining elements from b - if any
  while (j >= 0){
    array[j] = list->array[j];
    j--;
  }

  // Set the new size of the array
  size = len;
}

/*
  Get the array
*/
void TMRQuadrantArray::getArray( TMRQuadrant **_array, int *_size ){
  if (_array){ *_array = array; }
  if (_size){ *_size = size; }
}

/*
  Create an queue of quadrants 
*/
TMRQuadrantQueue::TMRQuadrantQueue(){
  root = tip = NULL;
  num_elems = 0;
}

/*
  Free the queue
*/
TMRQuadrantQueue::~TMRQuadrantQueue(){
  QuadQueueNode *node = root;
  while (node){
    QuadQueueNode *tmp = node;
    node = node->next;
    delete tmp;
  }
}

/*
  Get the length of the quadrant queue
*/
int TMRQuadrantQueue::length(){ 
  return num_elems;
}

/*
  Push a value onto the quadrant queue
*/
void TMRQuadrantQueue::push( TMRQuadrant *quad ){
  if (!tip){
    root = new QuadQueueNode();
    root->quad = *quad;
    tip = root;
  }
  else {
    tip->next = new QuadQueueNode();
    tip->next->quad = *quad;
    tip = tip->next;
  }
  num_elems++;
}

/*
  Pop a value from the quadrant queue
*/
TMRQuadrant TMRQuadrantQueue::pop(){
  if (!root){
    return TMRQuadrant();
  }
  else {
    num_elems--;
    TMRQuadrant temp = root->quad;
    QuadQueueNode *tmp = root;
    root = root->next;
    delete tmp;
    if (num_elems == 0){ tip = NULL; }
    return temp;
  }
}

/*
  Convert the queue to an array 
*/
TMRQuadrantArray* TMRQuadrantQueue::toArray(){
  // Allocate the array
  TMRQuadrant *array = new TMRQuadrant[ num_elems ];
  
  // Scan through the queue and retrieve the quadrants
  QuadQueueNode *node = root;
  int index = 0;
  while (node){
    array[index] = node->quad;
    node = node->next;
    index++;
  }

  // Create the array object
  TMRQuadrantArray *list = new TMRQuadrantArray(array, num_elems);
  return list;
}

/*
  A hash for quadrants

  Note that this isn't a true hash table since it does not associate
  elements with other values. It is used to create unique lists of
  elements and nodes within the quadree mesh.
*/
TMRQuadrantHash::TMRQuadrantHash( int _use_node_index ){
  use_node_index = _use_node_index;
  num_elems = 0;
  num_buckets = min_num_buckets;
  hash_buckets = new QuadHashNode*[ num_buckets ];
  memset(hash_buckets, 0, num_buckets*sizeof(QuadHashNode*));
}

/*
  Free the memory allocated by the quadrant hash
*/
TMRQuadrantHash::~TMRQuadrantHash(){
  // Free all the elements in the hash
  for ( int i = 0; i < num_buckets; i++ ){
    QuadHashNode *node = hash_buckets[i];
    while (node){        
      // Delete the old guy
      QuadHashNode *tmp = node;
      node = node->next;
      delete tmp;
    }
  }

  delete [] hash_buckets;
}

/*
  Covert the hash table to an array
*/
TMRQuadrantArray* TMRQuadrantHash::toArray(){
  // Create an array of quadrants
  TMRQuadrant *array = new TMRQuadrant[ num_elems ];

  // Loop over all the buckets
  for ( int i = 0, index = 0; i < num_buckets; i++ ){
    // Get the hash bucket and extract all the elements from this
    // bucket into the array
    QuadHashNode *node = hash_buckets[i];
    
    while (node){
      array[index] = node->quad;
      index++;
      node = node->next;
    }
  }
  
  // Create an array object and add it to the list
  TMRQuadrantArray *list = new TMRQuadrantArray(array, num_elems, 
                                                use_node_index);
  return list;
}

/*
  Add an quadrant to the hash table. 

  A new quadrant is added only if it is unique within the list of
  objects. The function returns true if the quadrant is added, and false
  if it already exists within the hash table.

  input:
  quad:   the quadrant that may be added to the hash table
  
  returns:
  true if the quadrant is added, false if it is not
*/
int TMRQuadrantHash::addQuadrant( TMRQuadrant *quad ){ 
  if (num_elems > 10*num_buckets){
    // Redistribute the quadrants to new buckets
    int num_old_buckets = num_buckets;
    num_buckets = 2*num_buckets;
    QuadHashNode **new_buckets = new QuadHashNode*[ num_buckets ];
    memset(new_buckets, 0, num_buckets*sizeof(QuadHashNode*));

    // Keep track of the end bucket
    QuadHashNode **end_buckets = new QuadHashNode*[ num_buckets ];
    
    // Redistribute the quadrant nodes based on the new
    // number of buckets within the hash data structure
    for ( int i = 0; i < num_old_buckets; i++ ){
      QuadHashNode *node = hash_buckets[i];
      while (node){
        int bucket = getBucket(&(node->quad));

        // If this is the first new bucket, create
        // the new node
        if (!new_buckets[bucket]){
          new_buckets[bucket] = new QuadHashNode;
          new_buckets[bucket]->quad = node->quad;
          end_buckets[bucket] = new_buckets[bucket];
        }
        else {
          end_buckets[bucket]->next = new QuadHashNode;
          end_buckets[bucket]->next->quad = node->quad;
          end_buckets[bucket] = end_buckets[bucket]->next;
        }
        
        // Delete the old guy
        QuadHashNode *tmp = node;

        // Increment the node to the next hash
        node = node->next;

        // Free the node
        delete tmp;
      }
    }

    delete [] end_buckets;
    delete [] hash_buckets;
    hash_buckets = new_buckets;
  }

  int bucket = getBucket(quad);
  
  // If no quadrant has been added to the bucket, 
  // create a new bucket
  if (!hash_buckets[bucket]){
    hash_buckets[bucket] = new QuadHashNode;
    hash_buckets[bucket]->quad = *quad;
    num_elems++;
    return 1;
  }
  else {
    // Get the head node for the corresponding bucket
    QuadHashNode *node = hash_buckets[bucket];
    if (use_node_index){
      // Search the bucket to see if there's another node in the list
      while (node){
        if (node->quad.compareNode(quad) == 0){
          return 0;
        }
        
        // If the next node does not exist, quit while node 
        // is the last node in the linked list
        if (!node->next){
          break;
        }
        node = node->next;
      }
    }
    else {
      // Perform the same search, except using quads
      while (node){
        if (node->quad.compare(quad) == 0){
          return 0;
        }
        
        if (!node->next){
          break;
        }
        node = node->next;
      }
    }
    
    // Add the quadrant as the last node
    node->next = new QuadHashNode;
    node->next->quad = *quad;
    num_elems++;
  }

  return 1;
}

/*
  Get a bucket to place the quadrant in. 

  This code creates a value based on the quadrant location within the
  mesh and then takes the remainder of the number of buckets.  
*/
int TMRQuadrantHash::getBucket( TMRQuadrant *quad ){
  // The hash value
  uint32_t val = 0;

  if (use_node_index){
    uint32_t u = 0, v = 0, w = 0, x = 0;
    u = quad->face;
    if (quad->x >= 0){
      v = quad->x;
    }
    else {
      v = (1 << (TMR_MAX_LEVEL + 1)) - quad->x;
    }
    if (quad->y >= 0){
      w = quad->y;
    }
    else {
      w = (1 << (TMR_MAX_LEVEL + 1)) - quad->y;
    }
    x = quad->info;

    // Compute the hash value
    val = TMRIntegerFourTupleHash(u, v, w, x);
  }
  else {
    uint32_t u = 0, v = 0, w = 0;
    u = quad->face;
    if (quad->x >= 0){
      v = quad->x;
    }
    else {
      v = (1 << (TMR_MAX_LEVEL + 1)) - quad->x;
    }
    if (quad->y >= 0){
      w = quad->y;
    }
    else {
      w = (1 << (TMR_MAX_LEVEL + 1)) - quad->y;
    }

    // Compute the hash value
    val = TMRIntegerTripletHash(u, v, w);
  }

  // Compute the bucket value
  int bucket = val % num_buckets;
  
  return bucket;
}
