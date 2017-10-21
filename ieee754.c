/* 
 * Trampo MC:
 *
 * Aplicar newton raphson para calculo da raiz de um numero
 * no padrao ieee 754 de 64 bits
 * 
 * O padrao ieee, para 64 bits sera representado por uma struct
 * As funcoes de manipulacao implementadas sao:
 *  Multiplication
 *  Addition
 *  Division
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

    void Addition(IeeeStandart x, IeeeStandart y, IeeeStandart *resp)
    {
        /* TODO:
         *
         * Resolver caso em que x ou y == 0;
         * 
         */
         
        int exp_difference;
        int normalize_mantissa;

        x.mantissa += 1.0;
        y.mantissa += 1.0;

        if ( x.expoent >= y.expoent)
        {
            resp->sign = x.sign;
            exp_difference = x.expoent - y.expoent;
            resp->expoent = x.expoent;
            y.expoent += exp_difference;
            y.mantissa /= pow(2, exp_difference);
        }
        else
        {
            resp->sign = y.sign;
            exp_difference = y.expoent - x.expoent;
            resp->expoent = y.expoent;
            x.expoent += exp_difference;
            x.mantissa /= pow(2, exp_difference);
        }

        resp->mantissa = x.sign == y.sign? x.mantissa + y.mantissa : x.mantissa - y.mantissa;
        normalize_mantissa = floor(log2(resp->mantissa));
        resp->mantissa = resp->mantissa/pow(2, normalize_mantissa);
        resp->expoent += normalize_mantissa;
        resp->mantissa -= 1.0;
    }

    void Division(IeeeStandart x, IeeeStandart y, IeeeStandart *resp)
    {
        int normalize_mantissa;

        resp->sign = x.sign ^ y.sign;
        x.mantissa += 1.0;
        y.mantissa += 1.0;
        resp->mantissa = x.mantissa/y.mantissa;
        normalize_mantissa = floor(log2(resp->mantissa));
        resp->mantissa /= pow(2, normalize_mantissa);
        resp->mantissa -= 1.0;
        resp->expoent = x.expoent - y.expoent + normalize_mantissa + BIAS;
    }

    int main()
    {
        IeeeStandart x1, x2, x3;
        Double2Ieee(127.03125, &x1);
        Double2Ieee(16.9375, &x2);
        Division(x1, x2, &x3);
        printf("\n%d, %d, %lf\n", x3.sign, x3.expoent, x3.mantissa);
        return 0;
    }
