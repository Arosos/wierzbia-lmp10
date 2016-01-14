#include "makespl.h"
#include <math.h>

void make_spl (points_t *pts, spline_t *spl)
{
	int i, j, m;
	int n = pts->n;
	double *x = pts->x, *y = pts->y;
	if (n % 2 == 0)
		m = n/2-1;
	else
		m = (n-1)/2-1;
	double p = x[0], k = x[n-1], a0, a[m], b[m], sumay, sumaa, sumab, c, d = 2.0/n, e;
	sumay = sumaa = sumab = 0.0;
	for (i = 0 ; i < n ; i++)
		sumay += y[i];
	a0 = sumay/n;
	if (alloc_spl(spl, n) == 0)
	{
		spl->n = m;
		for (i = 0 ; i < m ; i++)
		{
			for (j = 0 ; j < n ; j++)
			{
				c = 2.0*M_PI*j*i/n;
				sumaa += y[j] * cos(c);
				sumab += y[j] * sin(c);
			}
			a[i] = d*sumaa;
			b[i] = d*sumab;
			sumaa = sumab = 0.0;
		}
		for(i = 0 ; i < n ; i++)
		{
			double xx = spl->x[i] = p + i*(k-p)/(m-1);
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (j = 0 ; j < m ; j++)
			{
				e = 2.0*M_PI*j/n;
				spl->f[i] += a[j]*cos(e*xx) + b[j]*sin(e*xx);
				spl->f1[i] += (-1)*e*a[j]*sin(e*xx) + e*b[j]*cos(e*xx);
				spl->f2[i] += (-1)*e*e*a[j]*cos(e*xx) + (-1)*e*e*b[j]*sin(e*xx);
				spl->f3[i] += e*e*e*a[j]*sin(e*xx) + (-1)*e*e*e*b[j]*cos(e*xx);
			}
		}
	}
}
