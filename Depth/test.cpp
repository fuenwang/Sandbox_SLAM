#include "Depth.h"


int main(){
    DepthMap t(720, 1280);
    MAT a(720, 1280);
    //t.Print();
    //MAT a(3, 6);

    //t.Print();
    t(4, 3) = 10000;
    a.T();
    //t.Print();
    return 0;
}
