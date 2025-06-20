/* Spherical harmonics functions (unnormalized) for Laplace series */

/*
Copyright (C) 2006  Mikko Kaasalainen, Josef Durech
Modified by Shuai Feng during 2021 to 2023

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "declarations.h"

void sphfunc(int ndir, double *at, double *af,
	int Mmax, int Lmax,
	double** Ts, double** Tc, double** Fs, double** Fc, double*** Pleg, double** Dsph)
{
	int i, j, m, l, n, k, ibot;

	double aleg[MAX_LM + 1][MAX_LM + 1][MAX_LM + 1];

	for (i = 1; i <= ndir; i++)
	{
		Ts[i][0] = 1;
		Tc[i][0] = 1;
		Ts[i][1] = sin(at[i]);
		Tc[i][1] = cos(at[i]);
		for (j = 2; j <= Lmax; j++)
		{
			Ts[i][j] = Ts[i][1] * Ts[i][j - 1];
			Tc[i][j] = Tc[i][1] * Tc[i][j - 1];
		}
		Fs[i][0] = 0;
		Fc[i][0] = 1;
		for (j = 1; j <= Mmax; j++)
		{
			Fs[i][j] = sin(j * af[i]);
			Fc[i][j] = cos(j * af[i]);
		}
	}

	for (m = 0; m <= Lmax; m++)
		for (l = 0; l <= Lmax; l++)
			for (n = 0; n <= Lmax; n++)
				aleg[n][l][m] = 0;

	aleg[0][0][0] = 1;
	aleg[1][1][0] = 1;

	for (l = 1; l <= Lmax; l++)
		aleg[0][l][l] = aleg[0][l - 1][l - 1] * (2 * l - 1);

	for (m = 0; m <= Mmax; m++)
		for (l = m + 1; l <= Lmax; l++)
		{
			if ((2 * ((l - m) / 2)) == (l - m))
			{
				aleg[0][l][m] = -(l + m - 1) * aleg[0][l - 2][m] / (1 * (l - m));
				ibot = 2;
			}
			else
				ibot = 1;

			if (l != 1)
				for (n = ibot; n <= l - m; n = n + 2)
					aleg[n][l][m] = ((2 * l - 1) * aleg[n - 1][l - 1][m] - (l + m - 1) * aleg[n][l - 2][m]) / (1 * (l - m));
		}

	for (i = 1; i <= ndir; i++)
	{
		k = 0;
		for (m = 0; m <= Mmax; m++)
			for (l = m; l <= Lmax; l++)
			{
				Pleg[i][l][m] = 0;
				if ((2 * ((l - m) / 2)) == (l - m))
					ibot = 0;
				else
					ibot = 1;

				for (n = ibot; n <= l - m; n = n + 2)
					Pleg[i][l][m] = Pleg[i][l][m] + aleg[n][l][m] * Tc[i][n] * Ts[i][m];
				k++;
				Dsph[i][k] = Fc[i][m] * Pleg[i][l][m];
				if (m != 0)
				{
					k++;
					Dsph[i][k] = Fs[i][m] * Pleg[i][l][m];
				}
			}
	}
}

