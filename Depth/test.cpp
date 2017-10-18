#include "Depth.h"


int main(){
    DepthMap t(5, 9);
    //t.Print();
    //MAT a(3, 6);

    t.Print();
    t(4, 3) = 10000;
    t.Print();
    return 0;
}
