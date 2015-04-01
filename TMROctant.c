#include "TMROctant.h"

/*
  Get the child id of the quadrant
*/
int TMROctant::childId(){
  int id = 0;
  const int h = 1 << (TMR_MAX_LEVEL - level);
  
  id = id | ((x & h) ? 1 : 0);
  id = id | ((y & h) ? 2 : 0);
  id = id | ((z & h) ? 4 : 0);

  return id;
}

/*
  Return the sibling of the octant
*/
void TMROctant::getSibling( int id, TMROctant *sib ){
  const int h = 1 << (TMR_MAX_LEVEL - level);

  int xr = ((x & h) ? x-h : x);
  int yr = ((y & h) ? y-h : y);
  int zr = ((z & h) ? z-h : z);

  sib->level = level;
  sib->x = ((id & 1) ? xr+h : xr);
  sib->y = ((id & 2) ? yr+h : yr);
  sib->z = ((id & 4) ? zr+h : zr);
}

/*
  Get the parent of the octant
*/
void TMROctant::parent( TMROctant *p ){
  if (level > 0){
    p->level = level-1;
    const int h = 1 << (TMR_MAX_LEVEL - level);

    p->x = x & ~h;
    p->y = y & ~h;
    p->z = z & ~h;
  }
}

/*
  Get the face neighbour
*/
void TMROctant::faceNeighbor( int face, TMROctant *neighbor ){
  neighbor->level = level;
  const int h = 1 << (TMR_MAX_LEVEL - level);

  neighbor->x = x + ((face == 0) ? -h : (face == 1) ? h : 0);
  neighbor->y = y + ((face == 2) ? -h : (face == 3) ? h : 0);
  neighbor->z = z + ((face == 4) ? -h : (face == 5) ? h : 0);
}

/*
  Get the neighbor along an edge
*/
void TMROctant::edgeNeighbor( int edge, TMROctant *neighbor ){
  const int h = 1 << (TMR_MAX_LEVEL - level);
  neighbor->level = level;
  
  if (edge < 4){
    // Edges parallel to the x-direction
    neighbor->x = x;
    if (edge % 2 == 0){
      neighbor->y = y-h;
    }
    else {
      neighbor->y = y+h;
    }

    if (edge < 2){
      neighbor->z = z-h;
    }
    else {
      neighbor->z = z+h;
    }    
  }
  else if (edge < 8){
    // Edges parallel to the y-direction
    neighbor->y = y;

    if (edge % 2 == 0){
      neighbor->x = x-h;
    }
    else {
      neighbor->x = x+h;
    }

    if (edge < 6){
      neighbor->z = z-h;
    }
    else {
      neighbor->z = z+h;
    }    
  }
  else {
    // Edges parallel to the z-direction
    neighbor->z = z;
  
    if (edge % 2 == 0){
      neighbor->x = x-h;
    }
    else {
      neighbor->x = x+h;
    }

    if (edge < 10){
      neighbor->y = y-h;
    }
    else {
      neighbor->y = y+h;
    }
  }
}

/*
  Get the neighbor along a corner
*/
void TMROctant::cornerNeighbor( int corner, TMROctant *neighbor ){
  const int h = 1 << (TMR_MAX_LEVEL - level);
  neighbor->level = level;

  neighbor->x = x + (2*(corner & 1) - 1)*h;
  neighbor->y = y + ((corner & 2) - 1)*h;
  neighbor->z = z + ((corner & 4)/2 - 1)*h;
}

/*
  Compare two octants using the Morton enconding (z-ordering)
*/
int TMROctant::compare( const TMROctant *octant ) const {
  unsigned int xxor = x ^ octant->x;
  unsigned int yxor = y ^ octant->y;
  unsigned int zxor = z ^ octant->z;
  unsigned int sor = xxor | yxor | zxor;

  // If there is no most-significant bit, then we are done
  if (sor == 0){
    return level - octant->level;
  }

  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = x - octant->x;
  }
  else if (yxor > (sor ^ yxor)){
    discrim = y - octant->y;
  }
  else {
    discrim = z - octant->z;
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
  Compare two octants to determine whether they have the same Morton
  encoding, but may be at different levels.
*/
int TMROctant::compareEncoding( const TMROctant *octant ) const {
  unsigned int xxor = x ^ octant->x;
  unsigned int yxor = y ^ octant->y;
  unsigned int zxor = z ^ octant->z;
  unsigned int sor = xxor | yxor | zxor;

  // Note that here we do not distinguish between levels
  // Check for the most-significant bit
  int discrim = 0;
  if (xxor > (sor ^ xxor)){
    discrim = x - octant->x;
  }
  else if (yxor > (sor ^ yxor)){
    discrim = y - octant->y;
  }
  else {
    discrim = z - octant->z;
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
  Compare two octants within the same sub-tree
*/
static int compare_octants( const void *a, const void *b ){
  const TMROctant *ao = static_cast<const TMROctant*>(a);
  const TMROctant *bo = static_cast<const TMROctant*>(b);
  
  return ao->compare(bo);
}

/*
  Store a array of octants
*/
TMROctantArray::TMROctantArray( TMROctant *_array, int _size ){
  array = _array;
  size = _size;
  max_size = size;
  is_sorted = is_unique = 0;
}

/*
  Destroy the octant array
*/
TMROctantArray::~TMROctantArray(){
  delete [] array;
}

/*
  Just sort the array of entries without uniquifying them
*/
void TMROctantArray::sort(){
  qsort(array, size, sizeof(TMROctant), compare_octants);
  is_sorted = 1;
}

/*
  Sort the list and remove duplicates from the array of possible
  entries.
*/
void TMROctantArray::sortUnique(){
  qsort(array, size, sizeof(TMROctant), compare_octants);

  // Now that the Octants are sorted, remove duplicates
  int i = 0; // Location from which to take entries
  int j = 0; // Location to place entries
  
  for ( ; i < size; i++, j++ ){
    while ((i < size-1) && 
	   (array[i].compareEncoding(&array[i+1]) == 0)){
      i++;
    }

    if (i != j){
      array[j] = array[i];
    }
  }

  // The new size of the array
  size = j;

  is_sorted = is_unique = 1;
}

/*
  Merge the entries of two arrays
*/
void TMROctantArray::merge( TMROctantArray * list ){
  if (!is_unique){
    sortUnique();
  }
  if (!list->is_unique){
    list->sortUnique();
  }

  // Keep track of the number of duplicates
  int ndup = 0;

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
      ndup++;
    }
  }

  // Compute the required length of the new array
  int len = size + list->size - ndup; 

  // Allocate a new array if required
  if (len > max_size){
    max_size = len;

    TMROctant *temp = array;
    array = new TMROctant[ max_size ];
    memcpy(array, temp, size*sizeof(TMROctant));

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
void TMROctantArray::getArray( TMROctant **_array, int *_size ){
  *_array = array;
  *_size = size;
}


TMROctantQueue::TMROctantQueue(){
  root = tip = NULL;
  num_elems = 0;
}

int TMROctantQueue::length(){ 
  return num_elems;
}

void TMROctantQueue::push( TMROctant *oct ){
  if (!tip){
    root = new OcQueueNode();
    root->oct = *oct;
    tip = root;
  }
  else {
    tip->next = new OcQueueNode();
    tip->next->oct = *oct;
    tip = tip->next;
  }
  num_elems++;
}

TMROctant TMROctantQueue::pop(){
  if (!root){
    return TMROctant();
  }
  else {
    num_elems--;
    TMROctant temp = root->oct;
    root = root->next;
    return temp;
  }
}

TMROctantArray* TMROctantQueue::toArray(){
  TMROctant *array = new TMROctant[ num_elems ];
  
  OcQueueNode *node = root;
  int index = 0;
  while (node){
    array[index] = node->oct;
    node = node->next;
    index++;
  }

  TMROctantArray *list = new TMROctantArray(array, num_elems);
  return list;
}


TMROctantHash::TMROctantHash(){
  num_elems = 0;
  num_buckets = min_num_buckets;
  hash_buckets = new OcHashNode*[ num_buckets ];
  memset(hash_buckets, 0, num_buckets*sizeof(OcHashNode*));
}

TMROctantArray * TMROctantHash::toArray(){
  TMROctant *array = new TMROctant[ num_elems ];

  int max_len = 0;
  int min_len = 0;
  int av_len = 0;
  int min_index = -1;

  for ( int i = 0, index = 0; i < num_buckets; i++ ){
    OcHashNode *node = hash_buckets[i];
    
    int len = 0;
    while (node){
      array[index] = node->oct;
      index++;
      len++;
      node = node->next;
    }

    av_len += len;
    if (i == 0){
      min_len = len;
    }

    min_len = (len < min_len ? len : min_len);
    max_len = (len > max_len ? len : max_len);
    if (min_len == len){
      min_index = i;
    }
  }
  
  printf("maximum hash array length: %d\n", max_len);
  printf("minimum hash array length: %d  index: %d\n", min_len, min_index);
  printf("average hash array length: %d\n", av_len/num_buckets);

  TMROctantArray *list = new TMROctantArray(array, num_elems);
  return list;
}

int TMROctantHash::addOctant( TMROctant *oct ){
  int bucket = getBucket(oct);
  
  if (!hash_buckets[bucket]){
    hash_buckets[bucket] = new OcHashNode;
    hash_buckets[bucket]->oct = *oct;
    num_elems++;
    return 1;
  }
  else {
    OcHashNode *node = hash_buckets[bucket];
    while (node){
      if (node->oct.compare(oct) == 0){
	return 0;
      }
      
      if (!node->next){
	break;
      }
      node = node->next;
    }
    
    node->next = new OcHashNode;
    node->next->oct = *oct;
    node->next->next = NULL;
    num_elems++;
  }

  return 1;
}

// Get the buckets to place the octant in
int TMROctantHash::getBucket( TMROctant *oct ){
  int rx = oct->x >> (TMR_MAX_LEVEL - oct->level); 
  int ry = oct->y >> (TMR_MAX_LEVEL - oct->level); 
  int rz = oct->z >> (TMR_MAX_LEVEL - oct->level); 

  return rx*ry*rz % num_buckets;
}
