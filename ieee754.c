#include <stdio.h>
#include <math.h>

#define BIAS 1023

    double double2ieee(double value, int *signal, int *expoent, double *mantissa )
    {
        if( value == 0.0)
        {
            *signal = 0;
            *expoent = 0;
            *mantissa = 0;
        }
        else
        {
            *signal = value < 0? value = fabs(value), 1 : 0; 
            *expoent = floor(log2(value));
            *mantissa = value/pow(2,*expoent) - 1.0;
            *expoent += BIAS; 
        }
    }

    int main()
    {
        int s, e;
        double m;
        double2ieee(-9.5, &s, &e, &m);
        printf("\n%d, %d, %lf\n", s, e, m);
        return 0;
    }
