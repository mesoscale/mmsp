// File:    mmsp2snow.cpp
// Purpose: reads MMSP grid containing vector of doubles, uses first field (phase) and second field (undercooling)
//          to create a graphic representation of the solidified structure
// Output:  grayscale portable network graphics
// Depends: MMSP, DevIL image library, zlib

// Questions/Comments to kellet@rpi.edu (Trevor Keller)

// DevIL usage after http://bobobobo.wordpress.com/2009/03/02/how-to-use-openil-to-generate-and-save-an-image/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <zlib.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "MMSP.hpp"
#include "devil_cpp_wrapper.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  if ( argc != 3 )
  {
    cout << "Usage: " << argv[0] << " data.dat output.csv\n";
    return ( 1 );
  }

  // file open error check
  std::ifstream input(argv[1]);
  if (!input) {
    std::cerr<<"File input error: could not open "<<argv[1]<<".\n\n";
    exit(-1);
  }

  // read data type
  std::string type;
  getline(input,type,'\n');

  // grid type error check
  if (type.substr(0,4)!="grid") {
    std::cerr<<"File input error: file does not contain grid data."<<std::endl;
    exit(-1);
  }
  else if (type.substr(5,6)!="vector")
  {
    std::cerr<<"File input error: grid does not contain vector data."<<std::endl;
    exit(-1);
  }
  else if (type.substr(12,6)!="double")
  {
    std::cerr<<"File input error: vector data does not contain double-precision floats."<<std::endl;
    exit(-1);
  }

  // read grid dimension
  int dim;
  input>>dim;

  if (dim!=2)
  {
    std::cerr<<"ERROR: Snowflakes must be 2D.\n"<<std::endl;
    exit(-1);
  }

  // Initialize image
  ilInit();

  // construct grid object
  MMSP::grid<2,MMSP::vector<double> > grid(argv[1]);
  int bytesPerPx = 1; // set to 1 for 8-bit grayscale
  int width = MMSP::g1(grid,0)-MMSP::g0(grid,0);
  int height = MMSP::g1(grid,1)-MMSP::g0(grid,1);
  unsigned int theSize = height*width*bytesPerPx; // number of char's
  unsigned char* imData = new unsigned char[theSize];

  #ifdef DEBUG
  cout<<"Grid contains "<<MMSP::nodes(grid)<<" nodes. Allocated array of "<<theSize<<" bytes."<<endl;
  #endif

  // Populate the image from MMSP data.
  for (int i=0; i<MMSP::nodes(grid); ++i)
  {
    MMSP::vector<double> data = grid(i);
    unsigned char phi(255*data[0]);
    unsigned char T=255*(0.5-data[1]);
    for (int j=0; j<bytesPerPx; ++j) imData[ i+j ] = (0.375*255<phi&&phi<0.625*255)?0:((phi>0.625*255)?255:255-T);
  }

  ILenum Error;
  ILuint imageID = ilGenImage() ;
  ilBindImage(imageID);
  ilTexImage(width, height, 1, bytesPerPx, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData);
  Error = ilGetError();
  if (Error!=IL_NO_ERROR) cout<<"Error making image: "<<iluErrorString(Error)<<endl;
  ilEnable(IL_FILE_OVERWRITE);
  ilSave( IL_PNG, argv[2] ) ;
  Error = ilGetError();
  if (Error!=IL_NO_ERROR) cout<<"Error saving image: "<<iluErrorString(Error)<<endl;

  delete [] imData;

  return 0;
}
