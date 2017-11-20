/* 
 * Trampo MC:
 * 
 * Thiago Kira
 * ra78750
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
 * 
 * bibliografia:
 * https://acassis.wordpress.com/2012/02/18/descobrindo-a-raiz-quadrada-com-metodo-newton-raphson/
 * http://www.rfwireless-world.com/Tutorials/floating-point-tutorial.html
 * http://class.ece.iastate.edu/arun/Cpre305/ieee754/ie4.html
 * http://babbage.cs.qc.cuny.edu/IEEE-754.old/Decimal.html
 */
 
#include <stdio.h>
#include <math.h>

#define BIAS 1023

typedef struct Ieee
{
    int sign;
    int expoent;
    double mantissa;
} IeeeStandard;

    void Double2Ieee(double value, IeeeStandard *convert_value )
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

    int CheckZero(IeeeStandard x)
    {
        return (x.sign == 0 & x.expoent == 0 & x.mantissa == 0.0);
    }

    void Multiplication(IeeeStandard x, IeeeStandard y, IeeeStandard *resp)
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

    void Addition(IeeeStandard x, IeeeStandard y, IeeeStandard *resp)
    {
        /* TODO:
         *
         * Tratar caso em que x ou y == 0;
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

    void Division(IeeeStandard x, IeeeStandard y, IeeeStandard *resp)
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
    
    void Copy(IeeeStandard x, IeeeStandard *y)
    {
        y->sign = x.sign;
        y->expoent = x.expoent;
        y->mantissa = x.mantissa;
    }

    /*
     * Funcao Newton-Rapshon
     * 
     * Aplicacao do metodo newton rapshon para calculo da raiz quadrada
     * utilizando o padrao ieee para 64 bits. 
     * 
     * x(k+1) = (xk^2 + X) / 2*xk
     * resp  =  factor1   / factor2
     * 
     * @param:
     * x = valor a ser calculado
     * resp = resposta(referencia)
     * k = numero de iteracoes desejadas
     * 
     */
    void NewtonRaphson(IeeeStandard x, IeeeStandard *resp, int k)
    {
        IeeeStandard xk;
        IeeeStandard dois;
        IeeeStandard factor1, factor2;

        Double2Ieee(2.0, &dois);
        
        /* 
         * Determinar o xk inicial:
         * Quando o chute inicial eh muito ruim a funcao nao converge
         * Por exemplo, a funcao padrao que utilizei foi x0 = x/2
         * Para valores mais precisos (ex 25, 36 ...) a funcao conver-
         * ge normalmente
         * Para o caso do numero de avogrado a funcao apresentou pro - 
         * blemas, portanto foi nescessario roubar o chute inicial.
         */
        
        /* Chute padrao: */
        Division(x, dois, &xk);
    

        /* Chute roubado: */         
        Double2Ieee(700000000000.0, &xk);        


        while( k-- > 0)
        {
            Multiplication(xk, xk, &factor1);
            Addition(factor1, x, &factor1);
            Multiplication(dois, xk, &factor2);
            Division(factor1, factor2, resp);
            Copy(*resp, &xk);
        }
    }

    /* 
     * expoente - 1023
     */
    double SquareRoot(double mantissa, int exp_exc, double erro)
    {
        int impar = 0;

        if (exp_exc & 1)
        {
            impar = 1;
            exp_exc--;
        }

        mantissa += 1.0;    
       double x1 = mantissa;
       double x0 = 0.0;
       exp_exc /= 2;
       
       do
       {
           x0 = x1;
           x1 = x0 + (mantissa - (x0 * x0))/ (2*x0); 
           printf("\n[x0] : %.20lf \t [x1] : %.20lf ", x0, x1);
       } while (fabs(x1 - x0) > erro);
       printf("\nexp : %d", exp_exc);
       if(impar) return(x1 * sqrt(2));
       return x1;
    }


    int main()
    {
        IeeeStandard x1, x2;

        // 6 * 10²³, numero de avogrado no ieee standard:
        /*
        x1.sign = 0;
        x1.expoent = 1028;
        x1.mantissa = 0.125;
        */
        Double2Ieee(36.0, &x1);
        //NewtonRaphson(x1, &x2, 4);
        //printf("\n[S]:%d - [E]:%d - [M]%.12lf\n", x2.sign, x2.expoent, x2.mantissa);

        printf("\nResp : %.20lf\n", SquareRoot(x1.mantissa, (x1.expoent - BIAS), 1e-15));
        return 0;
    }
