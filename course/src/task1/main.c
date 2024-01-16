//          Computing  formal  solution  to  interval  linear  systems 
//          of the fixed-point form form x = Cx+d  in Kaucher complete 
//          interval  arithmetic by subddifferential Newton method 
//          with a special starting approximation.
//          (c) Sergey P. Shary, 1996                                    

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#pragma warning(disable : 4996)

int   DIM, * ip;
double* C, * d, * F, * x, * xx;


char solver()

//   solving a point linear system of the double dimension by Gaussian elimination 
{
    int   ier, i, j, k, m;
    double   p, q;

    ier = 0;     *(ip + DIM - 1) = 1;

    for (k = 0; k < DIM - 1; k++)
    {
        m = k; p = fabs(*(F + k * DIM + k));
        for (i = k + 1; i < DIM; i++)
        {
            q = fabs(*(F + i * DIM + k));
            if (p < q) { m = i; p = q; }
        }
        *(ip + k) = m; p = *(F + m * DIM + k);
        if (m != k)
        {
            *(ip + DIM - 1) = -*(ip + DIM - 1);
            *(F + m * DIM + k) = *(F + k * DIM + k);
            *(F + k * DIM + k) = p;
        }
        if (p == 0)
        {
            ier = k; *(ip + DIM - 1) = 0; break;
        }
        p = 1. / p;
        for (i = k + 1; i < DIM; i++)
            *(F + i * DIM + k) *= -p;
        for (j = k + 1; j < DIM; j++)
        {
            p = *(F + m * DIM + j);
            *(F + m * DIM + j) = *(F + k * DIM + j);
            *(F + k * DIM + j) = p;
            if (p != 0)
                for (i = k + 1; i < DIM; i++)
                    *(F + i * DIM + j) += *(F + i * DIM + k) * p;
        }
    }

    if (ier != 0 || *(F + DIM * DIM - 1) == 0)
    {
        printf("\n\r   The matrix of the system is not absolutely regular.");
        printf("\n\r   Try to change it a little. \n");
        return 1;
    }

    for (k = 0; k < DIM - 1; k++)
    {
        m = *(ip + k); q = *(xx + m); *(xx + m) = *(xx + k); *(xx + k) = q;
        for (i = k + 1; i < DIM; i++) *(xx + i) += *(F + i * DIM + k) * q;
    }

    for (j = 0; j < DIM - 1; j++)
    {
        k = DIM - j - 1;  *(xx + k) /= *(F + k * DIM + k);  q = -*(xx + k);
        for (i = 0; i < k; i++) *(xx + i) += *(F + i * DIM + k) * q;
    }

    *xx /= *F;  return 0;
}


int main(int argc, char* argv[])

{
    char   DataFile[64];
    int   i, j, l, m, ni, Dim;
    double   g0, g1, h0, h1, p, q, r, s0, s1, t0, t1, u0, u1, v0, v1, Eps, Tau;
    FILE* fi, * fopen();
    long   IterLim;

    printf("\r\n  ***   ");
    printf("Computing  algebraic  solutions  of  interval  linear  systems  ");
    printf("  *** \r\n   ** ");
    printf("   of the form  x = Cx + d  by  subdifferential  Newton  method ");
    printf("    ** \r\n    * ");
    printf("                    (c) Sergey P. Shary, 1995                ");
    printf("       * \n\n");

    
    strcpy(DataFile, "data.txt");

    //  reading input data of the problem,
    //              forming the "middle" point matrix F 
    //                          of the interval system under solution  

    fi = fopen(DataFile, "r");
    if (fi == NULL)
    {
        printf("\n\r    Cannot open the data file `%1s' \n", DataFile);
        goto final;
    }

    if (fscanf(fi, "%d %lf %lf %lu", &Dim, &Eps, &Tau, &IterLim) < 4)
    {
        printf("\n\r    Error reading the parameters of the algorithm! \n");
        goto final;
    }

    if (Dim < 1)
    {
        printf("\n\r    Incorrect dimension of the equation! \n");
        goto final;
    }

    if (Tau <= 0 || Tau > 1)
    {
        printf("\n\r    Incorrect damping factor! \n");
        goto final;
    }

    if (Eps < 0)
    {
        printf("\n\r    Incorrect relative defect! \n");
        goto final;
    }

    DIM = 2 * Dim;
    ip = calloc(DIM, sizeof(int));
    F = calloc(DIM * DIM, sizeof(double));
    x = calloc(DIM, sizeof(double));
    xx = calloc(DIM, sizeof(double));
    C = calloc(Dim * DIM, sizeof(double));
    d = calloc(DIM, sizeof(double));

    if (ip == NULL || F == NULL || x == NULL || xx == NULL || C == NULL || d == NULL)
    {
        printf("\n\r    Not enough room to allocate the data! \n");
        goto final;
    }

    for (i = 0; i < Dim; i++)
        for (j = 0; j < Dim; j++)
        {
            if (fscanf(fi, "%lf %lf", &u0, &u1) < 2)
            {
                printf("\n\r    Error reading the interval matrix! \n");
                goto final;
            }
            *(C + i * DIM + 2 * j) = u0; *(C + i * DIM + 2 * j + 1) = u1;
            p = -0.5 * (u0 + u1); if (i == j) p += 0.5;
            if (p >= 0.)
            {
                *(F + i * DIM + j) = p; *(F + (i + Dim) * DIM + j + Dim) = p;
                *(F + (i + Dim) * DIM + j) = 0; *(F + i * DIM + j + Dim) = 0;
            }
            else
            {
                *(F + i * DIM + j) = 0; *(F + (i + Dim) * DIM + j + Dim) = 0;
                *(F + (i + Dim) * DIM + j) = p; *(F + i * DIM + j + Dim) = p;
            }
        }

    for (i = 0; i < Dim; i++)
    {
        if (fscanf(fi, "%lf %lf", &v0, &v1) < 2)
        {
            printf("\n\r    Error reading the right-hand side vector!\n");
            goto final;
        }
        *(xx + i) = v0;  *(xx + i + Dim) = v1;
        *(d + i) = v0;   *(d + i + Dim) = v1;
    }

    fclose(fi);

    //      displaying the input data of the problem 

    printf("\r\n Dimension of the equation = %u", Dim);
    printf("\r\n The interval matrix C of the system: \n\n\r");

    for (i = 0; i < Dim; i++)
    {
        for (j = 0; j < Dim; j++)
            printf(" [ %- .3f, %- .3f ]", *(C + i * DIM + 2 * j), *(C + i * DIM + 2 * j + 1));
        printf("\n\r");
    }

    printf("\n The interval vector b: \n\n\r");

    for (i = 0; i < Dim; i++)
        printf(" [ %- .3f, %- .3f ]", *(d + i), *(d + i + Dim));

    printf("\n\n");

    //  solving the embedded "middle" point system,
    //          we find the starting iteration of the algorithm  

    if (solver()) goto final;

    //      launching the iteration 
    //                  of the subdifferential Newton method 
    ni = 0;

    do {
        printf("\n\n\r Iteration number: %d\n", ni);

        for (i = 0; i < Dim; i++)
        {
            u0 = *(xx + i);  u1 = *(xx + i + Dim);
            printf("\r\n%15d%7s", i + 1, ".     [");
            printf(" %16.11f", u0); printf(" ,");
            printf("%16.11f", u1);  printf(" ]");
            if (u0 <= u1)  printf("     ->\n");
            else  printf("     <-\n");
        }
        ni++; r = 0;

        for (i = 0; i < DIM; i++)
        {
            *(x + i) = *(xx + i);
            for (j = 0; j < DIM; j++) *(F + i * DIM + j) = (i != j) ? 0 : 1;
        }

        for (i = 0; i < Dim; i++)
        {
            s0 = *(x + i);    s1 = *(x + i + Dim);
            for (j = 0; j < Dim; j++)
            {
                g0 = *(C + i * DIM + 2 * j); g1 = *(C + i * DIM + 2 * j + 1);
                h0 = *(x + j); h1 = *(x + j + Dim);

                if (g0 * g1 > 0) l = (g0 > 0) ? 0 : 2; else l = (g0 <= g1) ? 1 : 3;
                if (h0 * h1 > 0) m = (h0 > 0) ? 1 : 3; else m = (h0 <= h1) ? 2 : 4;

                switch (4 * l + m)
                {
                case 1: t0 = g0 * h0; t1 = g1 * h1;
                    *(F + i * DIM + j) -= g0; *(F + (i + Dim) * DIM + j + Dim) -= g1;
                    break;
                case 2: t0 = g1 * h0; t1 = g1 * h1;
                    *(F + i * DIM + j) -= g1; *(F + (i + Dim) * DIM + j + Dim) -= g1;
                    break;
                case 3: t0 = g1 * h0; t1 = g0 * h1;
                    *(F + i * DIM + j) -= g1; *(F + (i + Dim) * DIM + j + Dim) -= g0;
                    break;
                case 4: t0 = g0 * h0; t1 = g0 * h1;
                    *(F + i * DIM + j) -= g0; *(F + (i + Dim) * DIM + j + Dim) -= g0;
                    break;
                case 5: t0 = g0 * h1; t1 = g1 * h1;
                    *(F + i * DIM + j + Dim) -= g0; *(F + (i + Dim) * DIM + j + Dim) -= g1;
                    break;
                case 6: u0 = g0 * h1; v0 = g1 * h0;  u1 = g0 * h0; v1 = g1 * h1;
                    if (u0 < v0) { t0 = u0; *(F + i * DIM + j + Dim) -= g0; }
                    else { t0 = v0; *(F + i * DIM + j) -= g1; }
                    if (u1 > v1) { t1 = u1; *(F + (i + Dim) * DIM + j) -= g0; }
                    else { t1 = v1; *(F + (i + Dim) * DIM + j + Dim) -= g1; }
                    break;
                case 7: t0 = g1 * h0; t1 = g0 * h0;
                    *(F + i * DIM + j) -= g1; *(F + (i + Dim) * DIM + j) -= g0;
                    break;
                case 8: t0 = 0; t1 = 0;
                    break;
                case 9: t0 = g0 * h1; t1 = g1 * h0;
                    *(F + i * DIM + j + Dim) -= g0; *(F + (i + Dim) * DIM + j) -= g1;
                    break;
                case 10: t0 = g0 * h1; t1 = g0 * h0;
                    *(F + i * DIM + j + Dim) -= g0; *(F + (i + Dim) * DIM + j) -= g0;
                    break;
                case 11: t0 = g1 * h1; t1 = g0 * h0;
                    *(F + i * DIM + j + Dim) -= g1; *(F + (i + Dim) * DIM + j) -= g0;
                    break;
                case 12: t0 = g1 * h1; t1 = g1 * h0;
                    *(F + i * DIM + j + Dim) -= g1; *(F + (i + Dim) * DIM + j) -= g1;
                    break;
                case 13: t0 = g0 * h0; t1 = g1 * h0;
                    *(F + i * DIM + j) -= g0; *(F + (i + Dim) * DIM + j) -= g1;
                    break;
                case 14: t0 = 0; t1 = 0;
                    break;
                case 15: t0 = g1 * h1; t1 = g0 * h1;
                    *(F + i * DIM + j + Dim) -= g1; *(F + (i + Dim) * DIM + j + Dim) -= g0;
                    break;
                case 16: u0 = g0 * h0; v0 = g1 * h1;  u1 = g0 * h1; v1 = g1 * h0;
                    if (u0 > v0) { t0 = u0; *(F + i * DIM + j) -= g0; }
                    else { t0 = v0; *(F + i * DIM + j + Dim) -= g1; }
                    if (u1 < v1) { t1 = u1; *(F + (i + Dim) * DIM + j + Dim) -= g0; }
                    else { t1 = v1; *(F + (i + Dim) * DIM + j) -= g1; }
                    break;
                }
                s0 -= t0;     s1 -= t1;
            }
            *(xx + i) = t0 = s0 - *(d + i);   *(xx + i + Dim) = t1 = s1 - *(d + i + Dim);
            t0 = fabs(t0); t1 = fabs(t1);     r += (t0 > t1) ? t0 : t1;
        }

        //      solving the point linear system 
        //                      with the subgradient matrix F    

        if (solver()) goto final;

        q = 0;
        for (i = 0; i < DIM; i++)
        {
            *(xx + i) = *(x + i) - *(xx + i) * Tau;
            q += fabs(*(xx + i));
        }

        if (q == 0) q = 1;

    } while (r / q > Eps && ni < IterLim);


    if (ni >= IterLim)
    {
        printf("\r    Doesn't the method diverge?! \n");
        printf("\r    The number of iterations alloted is exhausted!\n");
    }

    printf("\r Formal solution of the equation:\n");

    for (i = 0; i < Dim; i++)
    {
        u0 = *(xx + i);  u1 = *(xx + i + Dim);
        printf("\r\n%15d%7s", i + 1, ".     [");
        printf(" %16.11f", u0); printf(" ,");
        printf("%16.11f", u1);  printf(" ]");
        if (u0 <= u1)  printf("     ->");
        else  printf("     <-");
    }

    printf("\n\n\r Number of iterations = %d", ni);
    printf("\n\r 1-norm of the defect of the approximate solution = %.3e\n", r);

    final:;

}
