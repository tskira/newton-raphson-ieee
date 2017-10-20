/* 
 * Trampo MC:
 *
 * Aplicar newton raphson para calculo da raiz de um numero
 * no padrao ieee 754 de 64 bits
 * 
 * Para compilar:
 * gcc ieee754.c -lm
 */
 
#include <stdio.h>
#include <math.h>

#define BIAS 1023

typedef struct Ieee
{
    int sign;
    int expoent;
    double mantissa;
} IeeeStandart;

    void Double2Ieee(double value, IeeeStandart *convert_value )
    {
        if( value == 0.0)
        {
            convert_value->sign = 0;
            convert_value->expoent = 0;
            convert_value->mantissa = 0.0;
        }
        else
        {
            convert_value->sign = value < 0? value = fabs(value), 1 : 0; 
            convert_value->expoent = floor(log2(value));
            convert_value->mantissa = value/pow(2,convert_value->expoent) - 1.0;
            convert_value->expoent += BIAS; 
        }
    }

    int CheckZero(IeeeStandart x)
    {
        return (x.sign == 0 & x.expoent == 0 & x.mantissa == 0.0);
    }

    void Multipication(IeeeStandart x, IeeeStandart y, IeeeStandart *resp)
    {
        int normalize_exp;
        double normalize_mantissa;

        if(CheckZero(x) || CheckZero(y))
        {
            resp->sign = 0;
            resp->expoent = 0;
            resp->mantissa = 0.0;
        }
        else
        {
            resp->sign = x.sign ^ y.sign;
            resp->mantissa = (x.mantissa + 1.0) * (y.mantissa + 1.0) - 1.0;
            resp->expoent = (x.expoent - BIAS) + (y.expoent) ;
        }
    }

    int main()
    {
        IeeeStandart x1, x2, x3;
        Double2Ieee(0.085, &x1);
        Double2Ieee(0.085, &x2);
        Multipication(x1, x2, &x3);
        printf("\n%d, %d, %lf\n", x3.sign, x3.expoent, x3.mantissa);
        return 0;
    }
