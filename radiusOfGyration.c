#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include "postprocess.h"

// Compute and plot radius of gyration (2D) of individual chains and the distribution in melt.

typedef struct bounds_float
{
	float min, max;
} BOUNDS_FLOAT;

typedef struct coordinates
{
	float x, y, z;
} COORDINATES;

typedef struct RgStats
{
	float average, standardDeviation;
} RG_STATS;

void *reinitializeCOORDINATES (COORDINATES **array, int nAtomsPerChain)
{
	for (int i = 0; i < nAtomsPerChain; ++i)
	{
		(*array)[i].x = 0;
		(*array)[i].y = 0;
		(*array)[i].z = 0;
	}
}

void *computeRG2D (float **radiusOfGyration2d_all, struct lammpsdump *dumpinfo, int nAtoms, struct isotactic_polypropylene *monomers, int nChains, int nMonomers, int nTimeframes)
{
	float *centerOfMass_X, *centerOfMass_Y, *centerOfMass_Z, sum_coords_X = 0, sum_coords_Y = 0, sum_coords_Z = 0, sumRadius = 0, distance;

	centerOfMass_X = (float *) malloc (nChains * sizeof (float));
	centerOfMass_Y = (float *) malloc (nChains * sizeof (float));
	centerOfMass_Z = (float *) malloc (nChains * sizeof (float));
	int nAtomsPerChain = nAtoms / nChains;

	COORDINATES *currentChain;
	currentChain = (COORDINATES *) malloc (nAtomsPerChain * sizeof (COORDINATES));
	int currentUA = 0, CH_pos, CH2_pos, CH3_pos, currentChain_id = 0;

	reinitializeCOORDINATES (&currentChain, nAtomsPerChain);

	// Storing coordinates of all UAs in current chain and calculating center of mass
	for (int i = 0; i < nMonomers; ++i)
	{
		currentChain[currentUA].x = dumpinfo[monomers[i].CH].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH].z;
		currentUA++;

		currentChain[currentUA].x = dumpinfo[monomers[i].CH3].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH3].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH3].z;
		currentUA++;

		currentChain[currentUA].x = dumpinfo[monomers[i].CH2].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH2].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH2].z;
		currentUA++;

		if ((monomers[i].chain_id_number > currentChain_id) || ((i+1) == nMonomers))
		{
			sum_coords_X = 0;
			sum_coords_Y = 0;
			sum_coords_Z = 0;

			for (int j = 0; j < nAtomsPerChain; ++j)
			{
				sum_coords_X += currentChain[j].x;
				sum_coords_Y += currentChain[j].y;
				sum_coords_Z += currentChain[j].z;
			}

			centerOfMass_X[currentChain_id] = sum_coords_X / nAtomsPerChain;
			centerOfMass_Y[currentChain_id] = sum_coords_Y / nAtomsPerChain;
			centerOfMass_Z[currentChain_id] = sum_coords_Z / nAtomsPerChain;

			currentChain_id++;
			currentUA = 0;
			reinitializeCOORDINATES (&currentChain, nAtomsPerChain);
		}
	}

	// Calculating radius of gyration along XZ plane (ignoring Y axis)
	currentChain_id = 0;
	currentUA = 0;
	reinitializeCOORDINATES (&currentChain, nAtomsPerChain);

	for (int i = 0; i < nMonomers; ++i)
	{
		currentChain[currentUA].x = dumpinfo[monomers[i].CH].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH].z;
		currentUA++;

		currentChain[currentUA].x = dumpinfo[monomers[i].CH3].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH3].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH3].z;
		currentUA++;

		currentChain[currentUA].x = dumpinfo[monomers[i].CH2].x;
		currentChain[currentUA].y = dumpinfo[monomers[i].CH2].y;
		currentChain[currentUA].z = dumpinfo[monomers[i].CH2].z;
		currentUA++;

		if ((monomers[i].chain_id_number > currentChain_id) || ((i+1) == nMonomers))
		{
			for (int j = 0; j < nAtomsPerChain; ++j)
			{
				distance = sqrt (
					pow ((currentChain[j].x - centerOfMass_X[currentChain_id]), 2) + 
					pow ((currentChain[j].z - centerOfMass_Z[currentChain_id]), 2));
				sumRadius += distance;
			}
			(*radiusOfGyration2d_all)[currentChain_id + (nChains * (nTimeframes - 1))] = sumRadius / nAtomsPerChain;

			currentChain_id++;
			currentUA = 0;
			sumRadius = 0;
			reinitializeCOORDINATES (&currentChain, nAtomsPerChain);
		}
	}

	free (centerOfMass_X);
	free (centerOfMass_Y);
	free (centerOfMass_Z);
	free (currentChain);
}

RG_STATS *computeAverageRG2D (float *radiusOfGyration2d_all, int nTimeframes, int nChains)
{
	FILE *output;
	output = fopen ("distribution.rg2d", "w");
	RG_STATS *radiusOfGyration;
	radiusOfGyration = (RG_STATS *) malloc (nChains * sizeof (RG_STATS));

	float *sumRG;
	sumRG = (float *) calloc (nChains, sizeof (float));

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			sumRG[j] += radiusOfGyration2d_all[(j + (i * nChains))];
		}
	}

	for (int i = 0; i < nChains; ++i)
	{
		radiusOfGyration[i].average = sumRG[i] / nTimeframes;
	}

	float sd_numerator, *sd_summation, sd_denominator;
	sd_summation = (float *) calloc (nChains, sizeof (float));

	for (int i = 0; i < nTimeframes; ++i)
	{
		for (int j = 0; j < nChains; ++j)
		{
			sd_numerator = pow ((radiusOfGyration2d_all[(j + (i * nChains))] - radiusOfGyration[j].average), 2);
			sd_summation[j] += sd_numerator;
		}
	}

	for (int i = 0; i < nChains; ++i)
	{
		radiusOfGyration[i].standardDeviation = sqrt (sd_summation[i] / (nTimeframes - 1));
		// To do:
		// Randomize the %d value around 1.
		fprintf(output, "%d %f %f\n", (1), radiusOfGyration[i].average, radiusOfGyration[i].standardDeviation);
	}

	fclose (output);
	return radiusOfGyration;
}

BOUNDS_FLOAT findBoundsRG2D (RG_STATS *radiusOfGyration, int arraySize)
{
	BOUNDS_FLOAT boundary;

	for (int i = 0; i < arraySize; ++i)
	{
		if (i == 0)
		{
			boundary.max = radiusOfGyration[i].average;
			boundary.min = radiusOfGyration[i].average;
		}
		else
		{
			if (radiusOfGyration[i].average > boundary.max)
				boundary.max = radiusOfGyration[i].average;
			if (radiusOfGyration[i].average < boundary.min)
				boundary.min = radiusOfGyration[i].average;
		}
	}

	return boundary;
}

void printFrequency (RG_STATS *radiusOfGyration, float binWidth, int nChains)
{
	FILE *outputFrequency;
	outputFrequency = fopen ("frequency.rg2d", "w");

	BOUNDS_FLOAT boundary;
	boundary = findBoundsRG2D (radiusOfGyration, nChains);

	int binMin = (int) boundary.min, binMax = (int) boundary.max + 1, nBins = (int) ((binMax - binMin) / binWidth);
	int *binFrequency;
	binFrequency = (int *) calloc (nBins, sizeof (int));

	for (int bin = binMin; bin < binMax; )
	{
		bin += binWidth;

		for (int i = 0; i < nChains; ++i)
		{
			if ((radiusOfGyration[i].average > bin) && (radiusOfGyration[i].average <= (int) (bin + binWidth)))
			{
				binFrequency[bin - binMin]++;
			}
		}
	}

	for (int i = 0; i < nBins; ++i)
	{
		fprintf(outputFrequency, "%d %d\n", (i + binMin), binFrequency[i]);
	}

	fclose (outputFrequency);
}

int main(int argc, char const *argv[])
{
	FILE *input;
	char *inputfilename, lineString[5000];
	int nAtoms, nTimeframes = 0, nChains, nMonomers;
	float *radiusOfGyration2d_all;

	inputfilename = (char *) malloc (2000 * sizeof (char));

	if (argc == 2)
		input = fopen (argv[1], "r");
	else if (argc == 1)
	{
		inputfilename = getInputFileName();
		input = fopen (inputfilename, "r");
	}

	nAtoms = getNatoms (inputfilename);
	struct lammpsdump *dumpinfo;
	struct isotactic_polypropylene *monomers;

	do
	{
		nTimeframes++;
		fprintf(stdout, "%s: %d  \n", "nTimeframes", nTimeframes);
		dumpinfo = lscanf (input, nAtoms);
		monomers = classifyiPP (nAtoms, dumpinfo, &nMonomers, &nChains, nTimeframes);
		radiusOfGyration2d_all = (float *) realloc (radiusOfGyration2d_all, (nChains * nTimeframes) * sizeof (radiusOfGyration2d_all));
		computeRG2D (&radiusOfGyration2d_all, dumpinfo, nAtoms, monomers, nChains, nMonomers, nTimeframes);

		if (nTimeframes > 100)
			break;

	} while (dumpinfo[0].status);

	// Computing average and standard deviation
	RG_STATS *radiusOfGyration;
	radiusOfGyration = (RG_STATS *) malloc (nChains * sizeof (RG_STATS));

	radiusOfGyration = computeAverageRG2D (radiusOfGyration2d_all, nTimeframes, nChains);

	float binWidth = 1;
	printFrequency (radiusOfGyration, binWidth, nChains);

	return 0;
}