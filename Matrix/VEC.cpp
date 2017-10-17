// Vector CPP
// ID: 102061149
// Name: Fu-En Wang

#include "VEC.h"

VEC::VEC(int n){
    this->dimension = n;
    this->val = new double[n]();
}

VEC::VEC(const VEC &data){
    this->dimension = data.dimension;
    this->val = new double[data.dimension]();
    for(int i=0; i<this->dimension; i++)
        this->val[i] = data.val[i];
}

VEC::~VEC(){
    delete this->val;
}

int VEC::dim() const{
    return this->dimension; 
}

void VEC::sort(){
    double tmp;
    for(int i=0; i <= this->dimension; i++){
        for(int j=i+1; j<this->dimension; j++){
            if(this->val[i] > this->val[j]) {
                tmp = this->val[i];
                this->val[i] = this->val[j];
                this->val[j] = tmp;
            }
        }
    }
}

double VEC::sum(){
    double sum=0;
    for(int i=0; i<this->dimension; i++)
        sum += this->val[i];
    return sum;
}
void VEC::print(){
    for(int i=0; i<this->dimension; i++)
        printf("%g ",this->val[i]);
    printf("\n");
    return;
}

void VEC::checkDim(const VEC &data){
    if(this->dimension != data.dimension){
        printf("Dimension dismatch!\n");
        exit(0);
    }
}

bool VEC::operator==(VEC &data){
    if(this->dimension == data.dimension){
        for(int i=0; i<this->dimension; i++)
            if(this->val[i] != data.val[i])
                return false;
        return true;
    }
    else
        return false;
}

void VEC::operator=(const VEC &data){
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        this->val[i] = data.val[i];
    return;
}

void VEC::operator=(double num){
    for(int i=0; i<this->dimension; i++)
        this->val[i] = num;
    return;
}
void VEC::operator=(int num){
    for(int i=0; i<this->dimension; i++)
        this->val[i] = num;
    return;
}

void VEC::operator=(const double *a){
    for(int i=0; i<this->dimension; i++){
        this->val[i] = a[i];
    }
    return;
}

void VEC::operator+=(const VEC &data){
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        this->val[i] += data.val[i];
}

void VEC::operator-=(const VEC &data){
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        this->val[i] -= data.val[i];
}

void VEC::operator*=(const VEC &data){
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        this->val[i] *= data.val[i];
}

void VEC::operator/=(const VEC &data){
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        this->val[i] /= data.val[i];
}
void VEC::operator/=(double num){
    for(int i=0; i<this->dimension; i++)
        this->val[i] /= num;
}


double& VEC::operator[](int index) const{
    if(index >= this->dimension || index < 0){
        printf("array[index] out of bound\n");
        exit(0);
    }
    return this->val[index];
}

VEC VEC::operator+(const VEC &data){
    VEC out(this->dimension);
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        out.val[i] = this->val[i] + data.val[i];
    return out;
}

VEC VEC::operator-(const VEC &data){
    VEC out(this->dimension);
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        out.val[i] = this->val[i] - data.val[i];
    return out;

}

double VEC::operator*(const VEC &data){
    double sum = 0;
    //VEC out(this->dimension);
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        sum  += this->val[i] * data.val[i];
    return sum;
}

VEC VEC::operator/(const VEC &data){
    VEC out(this->dimension);
    //this->checkDim(data);
    for(int i=0; i<this->dimension; i++)
        out.val[i] = this->val[i] / data.val[i];
    return out;
}

VEC VEC::operator+(double num){
    VEC out(this->dimension);
    for(int i=0; i<this->dimension; i++)
        out[i] = this->val[i] + num;
    return out;
}

VEC VEC::operator-(double num){
    VEC out(this->dimension);
    for(int i=0; i<this->dimension; i++)
        out[i] = this->val[i] - num;
    return out;
}

VEC VEC::operator*(double num){
    VEC out(this->dimension);
    for(int i=0; i<this->dimension; i++)
        out[i] = this->val[i] * num;
    return out;
}

VEC VEC::operator/(double num){
    VEC out(this->dimension);
    for(int i=0; i<this->dimension; i++)
        out[i] = this->val[i] / num;
    return out;
}

VEC operator+(double num, const VEC &vec){
    VEC out(vec.dimension);
    for(int i=0; i<vec.dimension; i++)
        out.val[i] = num + vec.val[i];
    return out;
}
VEC operator-(double num, const VEC &vec){
    VEC out(vec.dimension);
    for(int i=0; i<vec.dimension; i++)
        out.val[i] = num - vec.val[i];
    return out;
}

VEC operator*(double num, const VEC &vec){
    VEC out(vec.dimension);
    for(int i=0; i<vec.dimension; i++)
        out.val[i] = num * vec.val[i];
    return out;
}

VEC operator/(double num, const VEC &vec){
    VEC out(vec.dimension);
    for(int i=0; i<vec.dimension; i++)
        out.val[i] = num / vec.val[i];
    return out;
}

double Lagrange(double x, VEC &XDATA, VEC &YDATA){
    int n = YDATA.dim();
    VEC NS(YDATA);
    for(int k=1; k<n; k++){
        for(int j=0; j<n-k; j++){
            NS[j] = ((x - XDATA[j]) * NS[j+1] - (x - XDATA[k+j]) * NS[j]) / (XDATA[j+k] - XDATA[j]);
        }
    }
    return NS[0];
}

double Integrate(int order, VEC &X, VEC &Y){
    const double *weight;
    switch(order){
        case 1:
            weight = W1;
            break;
        case 2:
            weight = W2;
            break;
        case 3:
            weight = W3;
            break;
        case 4:
            weight = W4;
            break;
        case 5:
            weight = W5;
            break;
        case 6:
            weight = W6;
            break;
        default:
            printf("The order is incorrect\n");
            exit(0);
    }
    int shape = Y.dim();
    double sum = 0;
    for(int index=0; index<shape-1; index += order){
        for(int k=0; k<=order; k++){
            sum += weight[k] * Y[index + k];
        }
    }
    return (X[1] - X[0]) * sum;
}

double Integrate(double (*func)(double), int order, int nblocks, double start, double end){
    const double *weight;
    switch(order){
        case 1:
            weight = W1;
            break;
        case 2:
            weight = W2;
            break;
        case 3:
            weight = W3;
            break;
        case 4:
            weight = W4;
            break;
        case 5:
            weight = W5;
            break;
        case 6:
            weight = W6;
            break;
        default:
            printf("The order is incorrect\n");
            exit(0);
    }
    double step = (end - start) / nblocks;
    //cout << step << endl;
    double total_sum = 0;
    double block_sum;
    double fx;
    double h = step / order;
    double xi = start;
    for(int block = 0; block < nblocks; block++){
        block_sum = 0;
        //cout << xi << endl;
        for(int i=0; i<=order; i++){
            fx = func(xi + i*h);
            block_sum += weight[i] * fx;
        }
        xi += step;
        total_sum += h * block_sum;
    }
    return total_sum;
}

int Newton(double (*func)(double), double (*func_der)(double), double &x, int maxIter, double tol){
    double tmp;
    for(int it = 1; it <= maxIter; it++){
        tmp = func(x);
        double frac = tmp / func_der(x);
        x -= frac;
        if(fabs(tmp) < tol){
            return it;
        }
    }
    return maxIter+1;
}

int Bisection(double (*func)(double), double &a, double &b, double &x, int maxIter, double tol){
    x = (a + b) / 2;
    for(int it = 1; it <= maxIter; it++){
        if(func(x) * func(a) <= 0){
            b = x;
        }
        else{
            a = x;
        }
        x = (a + b) / 2;
        if(fabs(x - a) <= tol){
            return it;
        }
    }
    return maxIter + 1;
}

int NewtonPoly(VEC &a, VEC &x, int maxIter, double tol){
    int order = a.dim() - 1;
    VEC b(order);
    VEC c(order-1);
    double x_i = x[order-1];
    double err;
    double b_minus_1;
    double c_minus_1;
    for(int n=order; n>=1; n--){
        err = tol + 1;
        for(int it=1; it <=maxIter && err >= tol; it++){
            b[n - 1] = a[n];
            if(n >= 2){
                c[n - 2] = b[n - 1];
            }
            for(int j = n - 2; j >= 0; j--){
                b[j] = a[j+1] + x_i * b[j+1];
            }
            for(int j = n - 3; j >= 0; j--){
                c[j] = b[j+1] + x_i * c[j+1];
            }
            b_minus_1 = a[0] + x_i * b[0];
            if(n >= 2){
                c_minus_1 = b[0] + x_i * c[0];
            }
            else{
                c_minus_1 = b[0];
            }
            x_i -= b_minus_1 / c_minus_1;
            err = fabs(b_minus_1);
        }
        x[n-1] = x_i;
        a = b;
    }
    return 0;
}

void FindMax(VEC &i, VEC &j){
    for(int idx = 0; idx < i.dim(); idx++){
        if(i[idx] > j[idx])
            j[idx] = i[idx];
    }
};
void FindMin(VEC &i, VEC &j){
    for(int idx = 0; idx < i.dim(); idx++){
        if(i[idx] < j[idx])
            j[idx] = i[idx];
    }
};






