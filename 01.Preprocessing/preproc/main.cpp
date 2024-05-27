//
// Reachtools2 for Paired-Tag V2 preprocessing
// Chenxu Zhu (czhu@nygenome.org)
// 06/22/2023
//

#include <iostream>
#include <string>
#include "preprocessing.h"
#include "cxstring.h"

using namespace std;

string version = "2024.04.20";

void help();

int main(int argc, char** argv){
  if(argc < 2){
    help();
    return 0;
  }

  string mod(argv[1]);


  // preprocessing DNA, extract reads
  if(mod == "combine_384plex"){
    string par1;
    if(argc<3){
      par1 = "help";
    }
    else{
      par1 = argv[2];
    }
    preprocessing::combine_384plex(par1);
    return 0;
  }

  if(mod == "convert"){
    string par1;
    if(argc<3){
      par1 = "help";
    }
    else{
      par1 = argv[2];
    }
    preprocessing::convert(par1);
    return 0;
  }




  return 0;
}


void help(){
  cout << "Preprocessing of Paired-seq/Tag and Paired-Damage-seq datasets." << "\nVersion: " << version << endl;
  cout << endl;
  cout << "Usage:\n\nStep 1: combine function:" << endl;
  cout << "[reachtools2 combine_384plex infile_prefix]\tProcess Paired-Tag/Paired-seq/Paired-damage-seq infile_prefix_R1.fq.gz & infile_prefix_R2.fq.gz.\n\t\t\t\t\t\tFor 384 plex format (8bp + 8bp + 4bp barcodes). \n\t\t\t\t\t\tOutput DNA, RNA, undermined fastq files. " << endl;
  cout << endl;

  cout << "\nStep 2: convert files:" << endl;
  cout << "[reachtools2 convert prefix_BC.sam]\t\tConvert the barcode mapped SAM file to fastq files.\n\t\t\t\t\t\tFor 384 plex format (8bp + 8bp + 4bp barcodes)." << endl;


  return;
}
