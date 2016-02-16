#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

// int N=100;

double drand(){
  return double(rand() / (RAND_MAX + 1.0));
}


int intrand(int N){
  return int ( drand() * N);
}
