#include "TMRBase.h"
#include "TMR_STLTools.h"

/*
  This is a conversion tool that converts the output .bstl file (which
  is the results of a parallel I/O and is in bindary) to a regular
  .stl file.
*/
int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);
  TMRInitialize();

  // Convert the extension (if any) on the input file to a .stl
  // extension since this will convert successful
  if (argc == 1){
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }

  char *infile = new char[ strlen(argv[1])+1 ];
  strcpy(infile, argv[1]);

  // Set the output file
  char *outfile = new char[ strlen(infile)+5 ];
  int len = strlen(infile);
  int i = len-1;
  for ( ; i >= 0; i-- ){
    if (infile[i] == '.'){ break; }     
  }
  if (i == 0){ 
    i = len-1; 
  }
  strcpy(outfile, infile);
  strcpy(&outfile[i], ".stl");

  if (strcmp(infile, outfile) != 0){
    TMR_ConvertBinToSTL(infile, outfile);
  }

  delete [] infile;
  delete [] outfile;

  TMRFinalize();
  MPI_Finalize();
  return (0);
}
