#include "stdafx.h"
#include <vector>
#include <iostream>  
#include <cmath>

using namespace std;

// ====================================================================

const double g_e(exp(1.0));
const double g_pi(acos(-1.0));

extern double *g_oShift, *g_m;
extern vector<double> g_y, g_z;

// ====================================================================

void TwoPeaksFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Two-Peak Trap */
void FiveUnevenFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Five-Uneven-Peak Trap */
void EqualMinFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Equal Minima */
void DecreaseMinFunc(vector<double> &, double &, int, double, int, int);							/* Expanded Decreasing Minima */
void UnevenMinFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Uneven Minima */
void HimmelblauFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Himmelblau¡¯s Function */
void CamelbackFunc(vector<double> &, double &, int, double, int, int);								/* Expanded Six-Hump Camel Back */
void VincentFunc(vector<double> &, double &, int, double, int, int);								/* Modified Vincent Function */

void SphereFunc(vector<double> &, double &, int, double *, double *, double, int, int);				/* Sphere */
void EllipsFunc(vector<double> &, double &, int, double *, double *, double, int, int);				/* Ellipsoidal */
void BentCigarFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Bent_Cigar */
void DiscusFunc(vector<double> &, double &, int, double *, double *, double, int, int);				/* Discus */
void DifPowersFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Different Powers */
void RosenbrockFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Rosenbrock's */
void AckleyFunc(vector<double> &, double &, int, double *, double *, double, int, int);				/* Ackley's */
void RastriginFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Rastrigin's */
void WeierstrassFunc(vector<double> &, double &, int, double *, double *, double, int, int);		/* Weierstrass's */
void GriewankFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Griewank's */
void SchwefelFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Schwefel's */
void KatsuuraFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Katsuura */
void GrieRosenFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* Griewank-Rosenbrock */
void Escaffer6Func(vector<double> &, double &, int, double *, double *, double, int, int);			/* Expanded Scaffer's F6 */
void HappyCatFunc(vector<double> &, double &, int, double *, double *, double, int, int);			/* HappyCat */
void HgBatFunc(vector<double> &, double &, int, double *, double *, double, int, int);				/* HGBat */

void Cf01(vector<double> &, double &, int, int);													/* Composition Function 1 */
void Cf02(vector<double> &, double &, int, int);													/* Composition Function 2 */
void Cf03(vector<double> &, double &, int, int);													/* Composition Function 3 */
void Cf04(vector<double> &, double &, int, int);													/* Composition Function 4 */
void Cf05(vector<double> &, double &, int, int);													/* Composition Function 5 */
void Cf06(vector<double> &, double &, int, int);													/* Composition Function 6 */
void Cf07(vector<double> &, double &, int, int);													/* Composition Function 7 */

void ShiftFunc(vector<double> &, vector<double> &, int, double *);
void RotateFunc(vector<double> &, vector<double> &, int, double *);
void SRFunc(vector<double> &, vector<double> &, int, double *, double *, double, int, int);			/* shift and rotate */
void CfCal(vector<double> &, double &, int, double *, double *, double *, int);

void Func(vector<double> &, double &, int, int);

// ====================================================================

void Func(vector<double> &x, double &f, int nx, int funcNum)
{
	switch (funcNum)
	{
	case 1:
		TwoPeaksFunc(x, f, nx, 1.0, 1, 1);
		f += 100.0;
		break;
	case 2:
		FiveUnevenFunc(x, f, nx, 1.0, 1, 1);
		f += 200.0;
		break;
	case 3:
		EqualMinFunc(x, f, nx, 0.05, 1, 1);
		f += 300.0;
		break;
	case 4:
		DecreaseMinFunc(x, f, nx, 0.05, 1, 1);
		f += 400.0;
		break;
	case 5:
		UnevenMinFunc(x, f, nx, 0.05, 1, 1);
		f += 500.0;
		break;
	case 6:
		HimmelblauFunc(x, f, nx, 0.2, 1, 1);
		f += 600.0;
		break;
	case 7:
		CamelbackFunc(x, f, nx, 0.05, 1, 1);
		f += 700.0;
		break;
	case 8:
		VincentFunc(x, f, nx, 0.2, 1, 1);
		f += 800.0;
		break;
	case 9:
		Cf01(x, f, nx, 1);
		f += 900.0;
		break;
	case 10:
		Cf02(x, f, nx, 1);
		f += 1000.0;
		break;
	case 11:
		Cf03(x, f, nx, 1);
		f += 1100.0;
		break;
	case 12:
		Cf04(x, f, nx, 1);
		f += 1200.0;
		break;
	case 13:
		Cf05(x, f, nx, 1);
		f += 1300.0;
		break;
	case 14:
		Cf06(x, f, nx, 1);
		f += 1400.0;
		break;
	case 15:
		Cf07(x, f, nx, 1);
		f += 1500.0;
		break;
	case 16:
		GriewankFunc(x, f, nx, g_oShift, g_m, 1.0, 1, 0);
		f += 1600.0;
		break;
	case 17:
		AckleyFunc(x, f, nx, g_oShift, g_m, 1.0, 1, 0);
		f += 1700.0;
		break;
	case 18:
		RosenbrockFunc(x, f, nx, g_oShift, g_m, 1.0, 1, 0);
		f += 1800.0;
		break;
	case 19:
		RastriginFunc(x, f, nx, g_oShift, g_m, 1.0, 1, 0);
		f += 1900.0;
		break;
	case 20:
		Escaffer6Func(x, f, nx, g_oShift, g_m, 1.0, 1, 0);
		f += 2000.0;
		break;
	}
}

void TwoPeaksFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Two-Peak Trap */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 20.0;  // shift to orgin
		if ((g_z.at(i) < 15.0)&(g_z.at(i) >= 0.0))
			f += -(160.0 / 15.0)*(15.0 - g_z.at(i));
		else if ((g_z.at(i) <= 20.0)&(g_z.at(i) >= 15.0))
			f += -40.0*(g_z.at(i) - 15.0);
		else if (g_z.at(i) < 0.0)
			f += -160.0 + pow(g_z.at(i), 2.0);
		else
			f += -200.0 + pow(g_z.at(i) - 20.0, 2.0);
	}

	f += 200.0*nx;
}

void FiveUnevenFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Five-Uneven-Peak Trap */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		if (g_z.at(i) < 0)
			f += -200.0 + pow(g_z.at(i), 2.0);
		else if (g_z.at(i) < 2.5)
			f += -80.0*(2.5 - g_z.at(i));
		else if (g_z.at(i) < 5.0)
			f += -64.0*(g_z.at(i) - 2.5);
		else if (g_z.at(i) < 7.5)
			f += -160.0 + pow(g_z.at(i), 2.0);
		else if (g_z.at(i) < 12.5)
			f += -28.0*(g_z.at(i) - 7.5);
		else if (g_z.at(i) < 17.5)
			f += -28.0*(17.5 - g_z.at(i));
		else if (g_z.at(i) < 22.5)
			f += -32.0*(g_z.at(i) - 17.5);
		else if (g_z.at(i) < 27.5)
			f += -32.0*(27.5 - g_z.at(i));
		else if (g_z.at(i) <= 30.0)
			f += -80.0*(g_z.at(i) - 27.5);
		else
			f += -200.0 + pow(g_z.at(i) - 30.0, 2.0);
	}

	f += 200.0*nx;
}

void EqualMinFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Equal Minima */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 0.1;  // shift to orgin
		if ((g_z.at(i) <= 1.0)&(g_z.at(i) >= 0.0))
			f += -pow((sin(5 * g_pi*g_z.at(i))), 6.0);
		else
			f += pow(g_z.at(i), 2.0);
	}

	f += 1.0*nx;
}

void DecreaseMinFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Decreasing Minima */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 0.1;  // shift to orgin
		if ((g_z.at(i) <= 1.0)&(g_z.at(i) >= 0.0))
			f += -exp(-2.0*log(2.0)*pow((g_z.at(i) - 0.1) / 0.8, 2.0))*pow(sin(5.0*g_pi*g_z.at(i)), 6.0);
		else
			f += pow(g_z.at(i), 2.0);
	}

	f += 1.0*nx;
}

void UnevenMinFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Uneven Minima */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 0.079699392688696;		  // pow(0.15, 4.0 / 3.0);  // shift to orgin
		if ((g_z.at(i) <= 1.0)&(g_z.at(i) >= 0.0))
			f -= pow(sin(5.0*g_pi*(pow(g_z.at(i), 0.75) - 0.05)), 6.0);
		else
			f += pow(g_z.at(i), 2.0);
	}

	f += 1.0*nx;
}

void HimmelblauFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Himmelblau¡¯s Function */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i < nx - 1; i += 2)
	{
		g_z.at(i) += 3.0;
		g_z.at(i + 1) += 2.0;  // shift to orgin

		f += pow((g_z.at(i) * g_z.at(i) + g_z.at(i + 1) - 11.0), 2.0) + pow((g_z.at(i) + g_z.at(i + 1) * g_z.at(i + 1) - 7.0), 2.0);
	}
}

void CamelbackFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Expanded Six-Hump Camel Back */
{
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i < nx - 1; i += 2)
	{
		g_z.at(i) += 0.089842;
		g_z.at(i + 1) += -0.712656;  // shift to orgin
		f += ((4.0 - 2.1*pow(g_z.at(i), 2.0) + pow(g_z.at(i), 4.0) / 3.0)*pow(g_z.at(i), 2.0) + g_z.at(i) * g_z.at(i + 1) + ((-4.0 + 4.0*pow(g_z.at(i + 1), 2.0))*pow(g_z.at(i + 1), 2.0)))*4.0;
	}

	f += 4.126514*nx / 2.0;
}

void VincentFunc(vector<double> &x, double &f, int nx, double shRate, int sFlag, int rFlag) /* Modified Vincent Function */
{
	// orginal bound [0.25, 10], optima=[0.333; 0.6242; 1.1701; 2.1933; 4.1112; 7.7063]
	f = 0.0;

	SRFunc(x, g_z, nx, g_oShift, g_m, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 4.1112;  // shift to orgin
		if ((g_z.at(i) >= 0.25)&(g_z.at(i) <= 10.0))
			f += -sin(10.0*log(g_z.at(i)));
		else if (g_z.at(i) < 0.25)
			f += pow(0.25 - g_z.at(i), 2.0) - sin(10.0*log(2.5));
		else
			f += pow(g_z.at(i) - 10, 2.0) - sin(10.0*log(10.0));
	}

	f = f / nx;
	f += 1.0;
}

void SphereFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Sphere */
{
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
		f += g_z.at(i) * g_z.at(i);
}

void EllipsFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Ellipsoidal */
{
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
		f += pow(10.0, 6.0*i / (nx - 1))*g_z.at(i) * g_z.at(i);
}

void BentCigarFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Bent_Cigar */
{
	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	f = g_z.at(0) * g_z.at(0);

	for (auto i(1); i != nx; ++i)
		f += pow(10.0, 6.0)*g_z.at(i) * g_z.at(i);
}

void DiscusFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Discus */
{
	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	f = pow(10.0, 6.0)*g_z.at(0) * g_z.at(0);

	for (auto i(1); i != nx; ++i)
		f += g_z.at(i) * g_z.at(i);
}

void DifPowersFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Different Powers */
{
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
		f += pow(fabs(g_z.at(i)), 2.0 + 4.0*i / (nx - 1.0));

	f = pow(f, 0.5);
}

void RosenbrockFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Rosenbrock's */
{
	double tmp1, tmp2;
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */
	g_z.at(0) += 1.0;  // shift to orgin

	for (auto i(0); i != nx - 1; ++i)
	{
		g_z.at(i + 1) += 1.0;  // shift to orgin
		tmp1 = g_z.at(i) * g_z.at(i) - g_z.at(i + 1);
		tmp2 = g_z.at(i) - 1.0;

		f += 100.0*tmp1*tmp1 + tmp2*tmp2;
	}
}

void AckleyFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Ackley's  */
{
	auto sum1(0.0), sum2(0.0);

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		sum1 += g_z.at(i) * g_z.at(i);
		sum2 += cos(2.0*g_pi*g_z.at(i));
	}

	sum1 = -0.2*sqrt(sum1 / nx);
	sum2 /= nx;
	f = g_e - 20.0*exp(sum1) - exp(sum2) + 20.0;
}

void WeierstrassFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Weierstrass's  */
{
	auto kMax(20);
	double sum, sum2, a(0.5), b(3.0);

	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (auto j(0); j <= kMax; ++j)
		{
			sum += pow(a, j)*cos(2.0*g_pi*pow(b, j)*(g_z.at(i) + 0.5));
			sum2 += pow(a, j)*cos(2.0*g_pi*pow(b, j)*0.5);
		}
		f += sum;
	}

	f -= nx*sum2;
}

void GriewankFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Griewank's  */
{
	auto s(0.0), p(1.0);

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		s += g_z.at(i) * g_z.at(i);
		p *= cos(g_z.at(i) / sqrt(1.0 + i));
	}

	f = 1.0 + s / 4000.0 - p;
}

void RastriginFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Rastrigin's  */
{
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
		f += (g_z.at(i) * g_z.at(i) - 10.0*cos(2.0*g_pi*g_z.at(i)) + 10.0);
}

void SchwefelFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Schwefel's  */
{
	double tmp;
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) += 4.209687462275036E+002;
		if (g_z.at(i) > 500)
		{
			f -= (500.0 - fmod(g_z.at(i), 500.0))*sin(pow(500.0 - fmod(g_z.at(i), 500.0), 0.5));
			tmp = (g_z.at(i) - 500.0) / 100.0;
			f += tmp*tmp / nx;
		}
		else if (g_z.at(i) < -500)
		{
			f -= (-500.0 + fmod(fabs(g_z.at(i)), 500.0))*sin(pow(500.0 - fmod(fabs(g_z.at(i)), 500.0), 0.5));
			tmp = (g_z.at(i) + 500.0) / 100.0;
			f += tmp*tmp / nx;
		}
		else
			f -= g_z.at(i) * sin(pow(fabs(g_z.at(i)), 0.5));
	}

	f += 4.189828872724338E+002*nx;
}

void KatsuuraFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Katsuura  */
{
	double temp, tmp1, tmp2, tmp3;
	f = 1.0;

	tmp3 = pow(1.0*nx, 1.2);

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		temp = 0.0;
		for (auto j(1); j <= 32; ++j)
		{
			tmp1 = pow(2.0, j);
			tmp2 = tmp1*g_z.at(i);
			temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
		}

		f *= pow(1.0 + (i + 1.0)*temp, 10.0 / tmp3);
	}

	tmp1 = 10.0 / nx / nx;

	f = f*tmp1 - tmp1;
}

void GrieRosenFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Griewank-Rosenbrock  */
{
	double temp, tmp1, tmp2;
	f = 0.0;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */
	g_z.at(0) += 1.0;  // shift to orgin

	for (auto i(0); i != nx - 1; ++i)
	{
		g_z.at(i + 1) += 1.0;  // shift to orgin
		tmp1 = g_z.at(i) * g_z.at(i) - g_z.at(i + 1);
		tmp2 = g_z.at(i) - 1.0;
		temp = 100.0*tmp1*tmp1 + tmp2*tmp2;

		f += (temp*temp) / 4000.0 - cos(temp) + 1.0;
	}

	tmp1 = g_z.at(nx - 1) * g_z.at(nx - 1) - g_z.at(0);
	tmp2 = g_z.at(nx - 1) - 1.0;
	temp = 100.0*tmp1*tmp1 + tmp2*tmp2;

	f += (temp*temp) / 4000.0 - cos(temp) + 1.0;
}

void Escaffer6Func(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* Expanded Scaffer¡¯s F6  */
{
	double temp1, temp2;

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	f = 0.0;

	for (auto i(0); i != nx - 1; ++i)
	{
		temp1 = sin(sqrt(g_z.at(i) * g_z.at(i) + g_z.at(i + 1) * g_z.at(i + 1)));
		temp1 = temp1*temp1;
		temp2 = 1.0 + 0.001*(g_z.at(i) * g_z.at(i) + g_z.at(i + 1) * g_z.at(i + 1));

		f += 0.5 + (temp1 - 0.5) / (temp2*temp2);
	}

	temp1 = sin(sqrt(g_z.at(nx - 1) * g_z.at(nx - 1) + g_z.at(0) * g_z.at(0)));
	temp1 = temp1*temp1;
	temp2 = 1.0 + 0.001*(g_z.at(nx - 1) * g_z.at(nx - 1) + g_z.at(0) * g_z.at(0));

	f += 0.5 + (temp1 - 0.5) / (temp2*temp2);
}

void HappyCatFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
																													 /* original global optimum: [-1, -1, ..., -1] */
{
	auto alpha(1.0 / 8.0), r2(0.0), sumZ(0.0);

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) = g_z.at(i) - 1.0;		  // shift to orgin
		r2 += g_z.at(i) * g_z.at(i);
		sumZ += g_z.at(i);
	}

	f = pow(fabs(r2 - nx), 2.0*alpha) + (0.5*r2 + sumZ) / nx + 0.5;
}

void HgBatFunc(vector<double> &x, double &f, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
																												  /* original global optimum: [-1, -1, ..., -1] */
{
	double alpha(1.0 / 4.0), r2(0.0), sumZ(0.0);

	SRFunc(x, g_z, nx, os, mr, shRate, sFlag, rFlag); /* shift and rotate */

	for (auto i(0); i != nx; ++i)
	{
		g_z.at(i) = g_z.at(i) - 1.0;  // shift to orgin
		r2 += g_z.at(i) * g_z.at(i);
		sumZ += g_z.at(i);
	}

	f = pow(fabs(pow(r2, 2.0) - pow(sumZ, 2.0)), 2 * alpha) + (0.5*r2 + sumZ) / nx + 0.5;
}

void Cf01(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 1 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 20, 10, 20, 10, 20, 10, 20, 10, 20 };
	double bias[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	i = 0;
	SphereFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	i = 1;
	SphereFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	i = 2;
	EllipsFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 3;
	EllipsFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 4;
	BentCigarFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 5;
	BentCigarFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 6;
	DiscusFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+4;
	i = 7;
	DiscusFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+4;
	i = 8;
	DifPowersFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+5;
	i = 9;
	DifPowersFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+5;

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf02(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 3 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	double bias[10] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90 };

	i = 0;
	EllipsFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+5;
	i = 1;
	EllipsFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+5;
	i = 2;
	DifPowersFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 3;
	DifPowersFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 4;
	BentCigarFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 5;
	BentCigarFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+6;
	i = 6;
	DiscusFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+4;
	i = 7;
	DiscusFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 1E+4;
	i = 8;
	SphereFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	i = 9;
	SphereFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf03(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 4 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 };
	double bias[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	i = 0;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10;
	i = 1;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10;
	i = 2;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10;
	i = 3;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10;
	i = 4;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10;
	i = 5;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10;
	i = 6;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 7;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 8;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 9;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf04(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 5 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 10, 20, 20, 30, 30, 40, 40, 50, 50 };
	double bias[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	i = 0;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10;
	i = 1;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10;
	i = 2;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10;
	i = 3;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10;
	i = 4;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10;
	i = 5;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10;
	i = 6;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 7;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 8;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 9;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf05(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 6 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	double bias[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	i = 0;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10.0;
	i = 1;
	HgBatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 2;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 3;
	AckleyFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 10.0;
	i = 4;
	WeierstrassFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-3, 1, rFlag);
	fit[i] *= 2.5;
	i = 5;
	KatsuuraFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] /= 1E+3;
	i = 6;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 7;
	GrieRosenFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 2.5;
	i = 8;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 9;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf06(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 7 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 10, 20, 20, 30, 30, 40, 40, 50, 50 };
	double bias[10] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180 };

	i = 0;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 1;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 2;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 3;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 4;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 5;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 6;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 7;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);
	i = 8;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 9;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void Cf07(vector<double> &x, double &f, int nx, int rFlag) /* Composition Function 8 */
{
	int i, cfNum(10);
	double fit[10];
	double delta[10] = { 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
	double bias[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	i = 0;
	RosenbrockFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 2.048E-2, 1, rFlag);
	fit[i] /= 10.0;
	i = 1;
	HgBatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 2;
	RastriginFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.12E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 3;
	AckleyFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] /= 10.0;
	i = 4;
	WeierstrassFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-3, 1, rFlag);
	fit[i] *= 2.5;
	i = 5;
	KatsuuraFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] /= 1E+3;
	i = 6;
	Escaffer6Func(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 1.0, 1, rFlag);
	fit[i] *= 100;
	i = 7;
	GrieRosenFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 2.5;
	i = 8;
	HappyCatFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 5.0E-2, 1, rFlag);
	fit[i] *= 10.0;
	i = 9;
	SchwefelFunc(x, fit[i], nx, &g_oShift[i*nx], &g_m[i*nx*nx], 10.0, 1, rFlag);

	CfCal(x, f, nx, delta, bias, fit, cfNum);
}

void ShiftFunc(vector<double> &x, vector<double> &xShift, int nx, double *os)
{
	for (auto i(0); i != nx; ++i)
		xShift.at(i) = x.at(i) - os[i];
}

void RotateFunc(vector<double> &x, vector<double> &xRot, int nx, double *mr)
{
	for (auto i(0); i != nx; ++i)
	{
		xRot.at(i) = 0;
		for (auto j(0); j != nx; ++j)
			xRot.at(i) += x.at(j) * mr[i*nx + j];
	}
}

void SRFunc(vector<double> &x, vector<double> &srX, int nx, double *os, double *mr, double shRate, int sFlag, int rFlag) /* shift and rotate */
{
	if (sFlag == 1)
	{
		if (rFlag == 1)
		{
			ShiftFunc(x, g_y, nx, os);

			for (auto i(0); i != nx; ++i)  // shrink to the orginal search range
				g_y.at(i) *= shRate;

			RotateFunc(g_y, srX, nx, mr);
		}
		else
		{
			ShiftFunc(x, srX, nx, os);

			for (auto i(0); i != nx; ++i)  // shrink to the orginal search range
				srX[i] *= shRate;
		}
	}
	else
	{
		if (rFlag == 1)
		{
			for (auto i(0); i != nx; ++i)  // shrink to the orginal search range
				g_y.at(i) = x.at(i) * shRate;

			RotateFunc(g_y, srX, nx, mr);
		}
		else
		{
			for (auto i(0); i != nx; ++i)  // shrink to the orginal search range
				srX[i] = x.at(i) * shRate;
		}
	}
}

void CfCal(vector<double> &x, double &f, int nx, double *delta, double *bias, double *fit, int cfNum)
{
	vector<double> w(cfNum);
	auto wMax(0.0), wSum(0.0);

	for (auto i(0); i != cfNum; ++i)
	{
		fit[i] += bias[i];
		w.at(i) = 0;

		for (auto j(0); j != nx; ++j)
			w.at(i) += pow(x.at(j) - g_oShift[i*nx + j], 2.0);

		if (w.at(i) != 0)
			w.at(i) = pow(1.0 / w.at(i), 0.5)*exp(-w.at(i) / 2.0 / nx / pow(delta[i], 2.0));
		else
			w.at(i) = DBL_MAX;
		if (w.at(i) > wMax)
			wMax = w.at(i);
	}

	for (auto i(0); i != cfNum; ++i)
		wSum += w.at(i);

	if (wMax == 0)
	{
		for (auto i(0); i != cfNum; ++i)
			w.at(i) = 1;

		wSum = cfNum;
	}

	f = 0.0;

	for (auto i(0); i != cfNum; ++i)
		f = f + w.at(i) / wSum*fit[i];
}