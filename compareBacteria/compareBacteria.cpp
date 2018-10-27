#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
	long* second;
	long one_l[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;


	void InitVectors()
	{
		vector = new long[M];
		second = new long[M1];
		stochastic = new long double[M];
		memset(vector, 0, M * sizeof(long));
		memset(second, 0, M1 * sizeof(long));
		memset(one_l, 0, AA_NUMBER * sizeof(long));
		memset(stochastic, 0, M * sizeof(long double));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i = 0; i < LEN - 1; i++)
		{
			short enc = encode(buffer[i]);
			one_l[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		second[indexs]++;
	}

	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		one_l[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		vector[index]++;
		total++;
		indexs = (indexs % M2) * AA_NUMBER + enc;
		second[indexs]++;
	}

public:
	long* vector;
	long double* stochastic;

	Bacteria(char* filename)
	{
		FILE * bacteria_file = fopen(filename, "r");
		InitVectors();

		char ch;
		while ((ch = fgetc(bacteria_file)) != EOF)
		{
			if (ch == '>')
			{
				while (fgetc(bacteria_file) != '\n'); // skip rest of line

				char buffer[LEN - 1];
				fread(buffer, sizeof(char), LEN - 1, bacteria_file);
				init_buffer(buffer);
			}
			else if (ch != '\n')
				cont_buffer(ch);
		}
		fclose(bacteria_file);

#pragma omp parallel for
	for (long i = 0; i < M; i++) {
		double p1 = (double)second[i / AA_NUMBER] / (total + complement);
		double p2 = (double)one_l[i % AA_NUMBER] / total_l;
		double p3 = (double)second[i % M1] / (total + complement);
		double p4 = (double)one_l[i / M1] / total_l;
		stochastic[i] = total * (p1 * p2 + p3 * p4) / 2;
	}
}

	
	~Bacteria()
	{
		delete vector;
		delete second;
		delete stochastic;
	}
	
	double stochastic_compute(long i)
	{
		double p1 = (double)second[i / AA_NUMBER] / (total + complement);
		double p2 = (double)one_l[i % AA_NUMBER] / total_l;
		double p3 = (double)second[i % M1] / (total + complement);
		double p4 = (double)one_l[i / M1] / total_l;
		return total * (p1 * p2 + p3 * p4) / 2;
	}
};

void ReadInputFile(char* input_name)
{
	FILE* input_file = fopen(input_name, "r");
	fscanf(input_file, "%d", &number_bacteria);
	bacteria_name = new char*[number_bacteria];

	for (long i = 0; i < number_bacteria; i++)
	{
		bacteria_name[i] = new char[20];
		fscanf(input_file, "%s", bacteria_name[i]);
		strcat(bacteria_name[i], ".faa");
	}
	fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1 = 0;
	double vector_len2 = 0;

	for (long i = 0; i < M; i++)
	{

		double t1;
		if (b1->stochastic[i] > EPSILON) {
			t1 = (b1->vector[i] - b1->stochastic[i]) / b1->stochastic[i];
			vector_len1 += (t1 * t1);
		}
		else
			t1 = 0;

		double t2;
		if (b2->stochastic[i] > EPSILON) {
			t2 = (b2->vector[i] - b2->stochastic[i]) / b2->stochastic[i];
			vector_len2 += (t2 * t2);
		}
		else
			t2 = 0;

		correlation += t1 * t2;
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void CompareAllBacteria()
{
	for (int i = 0; i < number_bacteria - 1; i++)
	{
		Bacteria* b1 = new Bacteria(bacteria_name[i]);

		for (int j = i + 1; j < number_bacteria; j++)
		{
			Bacteria* b2 = new Bacteria(bacteria_name[j]);
			double correlation = CompareBacteria(b1, b2);
			printf("%03d %03d -> %.10lf\n", i, j, correlation);
			delete b2;
		}
		delete b1;
	}
}

void main(int argc, char * argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile(argv[1]);
	CompareAllBacteria();

	time_t t2 = time(NULL);
	printf("time elapsed: %d seconds\n", t2 - t1);
}
