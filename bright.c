/* this computes brightness and its derivatives w.r.t. parameters */

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
#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include "constants.h"
#include "declarations.h"

double bright(double *ee, double *ee0, double t, double *cg,
	double *dyda, int ncoef,int Nphpar,int Numfac,
	double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg,double Phi_0,
	double* fTime, int dSegmentOfExposuret, double** fCoorAS0, double** fCoorAE0)
{
	int ncoef0, i, j, k,
		*incl;

	double cos_alpha, br, cl, cls, alpha, sum, sum0, dnom,
		e[4], e0[4],
		* php, * dphp,
		* mu, * mu0, * s, * dbr,
		* dsmu, * dsmu0,
		de[4][4], de0[4][4], tmat[4][4],
		dtm[4][4][4];

	/*fast rotating*/
	int k3;
	double* fdyda = vector_double(MAX_N_PAR + 1);
	double fbr = 0;
	/*fast rotating*/

	incl = vector_int(MAX_N_FAC + 1);

	php = vector_double(N_PHOT_PAR + 1);
	dphp = vector_double(N_PHOT_PAR + 1);
	mu = vector_double(MAX_N_FAC + 1);
	mu0 = vector_double(MAX_N_FAC + 1 + 1);
	s = vector_double(MAX_N_FAC + 1);
	dbr = vector_double(MAX_N_FAC + 1);
	dsmu = vector_double(MAX_N_FAC + 1);
	dsmu0 = vector_double(MAX_N_FAC + 1);

	ncoef0 = ncoef - 2 - Nphpar;
	cl = exp(cg[ncoef - 1]); /* Lambert */
	cls = cg[ncoef];       /* Lommel-Seeliger */

	/*fast rotating*/
	for (i = 1; i <= ncoef0 - 3; i++)
	{
		dyda[i] = 0;
	}
	for (k = 1; k <= 3; k++)
	{
		dyda[ncoef0 - 3 + k] = 0;
	}
	br = 0;
	/*fast rotating*/

	for (k3 = 0; k3 < dSegmentOfExposuret; k3++)
	{
		cos_alpha = dot_product(fCoorAE0[k3], fCoorAS0[k3]);
		//printf("%f\n", cos_alpha);
		alpha = acos(cos_alpha);
		for (i = 1; i <= Nphpar; i++)
			php[i] = cg[ncoef0 + i];

		phasec(dphp, alpha, php, Scale); /* computes also Scale */

		matrix(cg[ncoef0], fTime[k3], tmat, dtm, Blmat, Dblm, Phi_0);	//¼ÆËãtmat¡¢dtm

		fbr = 0;
		/* Directions (and ders.) in the rotating system */
		for (i = 1; i <= 3; i++)
		{
			e[i] = 0;
			e0[i] = 0;
			for (j = 1; j <= 3; j++)
			{
				e[i] = e[i] + tmat[i][j] * fCoorAE0[k3][j];
				e0[i] = e0[i] + tmat[i][j] * fCoorAS0[k3][j];
				de[i][j] = 0;
				de0[i][j] = 0;
				for (k = 1; k <= 3; k++)
				{
					de[i][j] = de[i][j] + dtm[j][i][k] * fCoorAE0[k3][k];
					de0[i][j] = de0[i][j] + dtm[j][i][k] * fCoorAS0[k3][k];
				}
			}
		}

		/*Integrated brightness (phase coeff. used later) */
		for (i = 1; i <= Numfac; i++)
		{
			incl[i] = 0;
			mu[i] = e[1] * Nor[i][1] + e[2] * Nor[i][2] + e[3] * Nor[i][3];
			mu0[i] = e0[1] * Nor[i][1] + e0[2] * Nor[i][2] + e0[3] * Nor[i][3];
			if ((mu[i] > TINY) && (mu0[i] > TINY))
			{
				incl[i] = 1;
				dnom = mu[i] + mu0[i];
				s[i] = mu[i] * mu0[i] * (cl + cls / dnom);
				fbr = fbr + Area[i] * s[i];
				dsmu[i] = cls * pow(mu0[i] / dnom, 2) + cl * mu0[i];
				dsmu0[i] = cls * pow(mu[i] / dnom, 2) + cl * mu[i];
				dbr[i] = Darea[i] * s[i];
			}
		}

		/* Derivatives of brightness w.r.t. g-coeffs */
		for (i = 1; i <= ncoef0 - 3; i++)
		{
			fdyda[i] = 0;
			for (j = 1; j <= Numfac; j++)
				if (incl[j] == 1)
					fdyda[i] = fdyda[i] + dbr[j] * Dg[j][i];
			fdyda[i] = *Scale * fdyda[i];
		}
		/*   printf("%f \n", dyda[1]);   */
		   /* Ders. of brightness w.r.t. rotation parameters */
		for (k = 1; k <= 3; k++)
		{
			fdyda[ncoef0 - 3 + k] = 0;
			for (i = 1; i <= Numfac; i++)
				if (incl[i] == 1)
				{
					sum = 0;
					sum0 = 0;
					for (j = 1; j <= 3; j++)
					{
						sum = sum + Nor[i][j] * de[j][k];
						sum0 = sum0 + Nor[i][j] * de0[j][k];
					}
					fdyda[ncoef0 - 3 + k] = fdyda[ncoef0 - 3 + k] + Area[i] * (dsmu[i] * sum + dsmu0[i] * sum0);
				}
			fdyda[ncoef0 - 3 + k] = *Scale * fdyda[ncoef0 - 3 + k];
		}

		/* Ders. of br. w.r.t. phase function params. */
		for (i = 1; i <= Nphpar; i++)
			fdyda[ncoef0 + i] = fbr * dphp[i];

		/* Ders. of br. w.r.t. cl, cls */
		fdyda[ncoef - 1] = 0;
		fdyda[ncoef] = 0;
		for (i = 1; i <= Numfac; i++)
			if (incl[i] == 1)
			{
				fdyda[ncoef - 1] = fdyda[ncoef - 1] + mu[i] * mu0[i] * Area[i];
				fdyda[ncoef] = fdyda[ncoef] + Area[i] * mu[i] * mu0[i] / (mu[i] + mu0[i]);
			}
		fdyda[ncoef - 1] = *Scale * fdyda[ncoef - 1] * cl;
		fdyda[ncoef] = *Scale * fdyda[ncoef];

		/*fast rotating*/
		for (i = 1; i <= ncoef0 - 3; i++)
		{
			dyda[i] += fdyda[i];
		}
		for (k = 1; k <= 3; k++)
		{
			dyda[ncoef0 - 3 + k] += fdyda[ncoef0 - 3 + k];
		}
		/*fast rotating*/

		/* Scaled brightness */
		fbr = fbr * *Scale;

		/*fast rotating*/
		br += fbr;
		/*fast rotating*/

	}
	//printf("%f\n", br);

	/*fast rotating*/
	br = br / dSegmentOfExposuret;
	for (i = 1; i <= ncoef0 - 3; i++)
	{
		dyda[i] /= dSegmentOfExposuret;
		//printf("%f\n", dyda[i]);
	}
	for (k = 1; k <= 3; k++)
	{
		dyda[ncoef0 - 3 + k] /= dSegmentOfExposuret;
		//printf("%f\n", dyda[ncoef0 - 3 + k]);
	}
	/*fast rotating*/

	deallocate_vector((void*)incl);
	deallocate_vector((void*)php);
	deallocate_vector((void*)dphp);
	deallocate_vector((void*)mu);
	deallocate_vector((void*)mu0);
	deallocate_vector((void*)s);
	deallocate_vector((void*)dbr);
	deallocate_vector((void*)dsmu);
	deallocate_vector((void*)dsmu0);
	deallocate_vector((void*)fdyda);

	return(br);
}
