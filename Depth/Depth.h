#ifndef DEPTH_H
#define DEPTH_H
#include "../Matrix/MAT.h"

class DepthMap{
    private:
        MAT *depth;
    public:
        DepthMap(int h, int w);
        ~DepthMap();
        void Print();
};



#endif
