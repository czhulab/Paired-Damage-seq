#!/usr/bin/env
if g++ -std=c++0x -pthread main.cpp cxstring.cpp preprocessing.cpp -o ../preproc
then 
  echo make_finished
fi
