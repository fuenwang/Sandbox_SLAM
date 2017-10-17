#include "Depth.h"

DepthMap::DepthMap(int h, int w){
    //this->depth = (MAT*)malloc(sizeof(MAT));
    this->depth = new MAT(h, w);
    //this->depth->MAT(h, w);
}

DepthMap::~DepthMap(){
    //this->depth->~MAT();
    delete this->depth;
}

void DepthMap::Print(){
    this->depth->print();
}
