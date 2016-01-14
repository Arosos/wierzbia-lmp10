#include "makespl.h"
#include <math.h>

void make_spl (points_t *pts, spline_t *spl)
{
	int i, j, m;
	int n = pts->n;
	double *x = pts->x;
	double *y = pts->y;
	m = (n-1)/2;
	if (n % 2 == 0)
		m = (n/2)-1;
	else
		m = (n-1)/2-1;
	double a0, a[m], b[m], sumay, sumaa, sumab, c, d = 2.0/n;
	sumay = sumaa = sumab = 0.0;
	for (i = 0 ; i < n ; i++)
		sumay += y[i];
	a0 = sumay/n;
	if (alloc_spl(spl, n) == 0)
	{
		spl->x = x;
		spl->y = y;
	}
	for (i = 0 ; i < m ; i++)
	{
		spl->a[i] = spl->b[i] = 0.0;
		for (j = 0 ; j < n ; j++)
		{
			c = 2.0*M_PI*(j+1)*(i+1)/n;
			sumaa += y[j] * cos(c);
			sumab += y[j] * sin(c);
		}
		spl->a[i] = d*sumaa;
		spl->b[i] = d*sumab;
	}
}
