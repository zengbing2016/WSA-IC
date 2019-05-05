// WSA-IC.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Algorithm.h"
#include <direct.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

// ====================================================================

int g_n;								// function dimension
int g_p;								// population size
unsigned long g_evalNum;				// evaluation number
int g_runNum;							// number of runs

double g_fErr;							// fitness error to distinguish global optima and local optima
double g_dThr;							// distance threshold to distinguish different optima

vector<pair<double, double>> g_vBd;		// variables bounds

auto g_totFuncNum(0);					// total number of functions

extern int g_sThr;						// stability threshold
extern double g_fThr;

// ====================================================================

vector<vector<double>> g_gloOpt;		// the positions of global optima in each run over g_runNum independent runs
vector<double> g_gloFit;				// the fitness of global optima in each run over g_runNum independent runs

vector<vector<double>> g_gloFit_;

double *g_oShift, *g_m;
vector<double> g_y, g_z;

int g_cfNums[] = { 0, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 1, 1, 1, 1, 1 };

clock_t g_start, g_finish;

// ====================================================================

bool CalcFuncNum()
{
	ifstream paraIn("input/parameters/parameters.txt");

	if (!paraIn.is_open())
	{
		cout << "The parameters file is nonexistent, or cannot open the parameters file for reading!" << endl;
		return false;
	}
	else
	{
		char str[256];
		memset(str, '0', 256);

		while (!paraIn.eof())
		{
			paraIn.getline(str, 100);
			if (str + 255 != find(str, str + 255, 'F'))
				g_totFuncNum++;
		}

		paraIn.close();
	}

	return true;
}

bool Input(int funcNum)
{
	// Load parameters from files	
	ifstream paraIn("input/parameters/parameters.txt");

	if (!paraIn.is_open())
	{
		cout << "The parameters file is nonexistent, or cannot open the parameters file for reading!" << endl;
		return false;
	}
	else
	{
		char str[256], str_[256];

		sprintf_s(str_, "F%d", funcNum);

		while (!paraIn.eof())
		{
			paraIn.getline(str, 100);
			if (!strcmp(str, str_))
				break;
		}

		if (paraIn.eof())
		{
			cout << "There are not matching parameters!" << endl;
			return false;
		}
		else
		{
			paraIn.seekg(19L, ios_base::cur);
			paraIn.getline(str, 100);
			g_n = atoi(str);
			g_sThr = g_n * 100;

			paraIn.seekg(16L, ios_base::cur);
			paraIn.getline(str, 100);
			g_p = atoi(str);

			paraIn.seekg(18L, ios_base::cur);
			paraIn.getline(str, 100);
			g_evalNum = atoi(str);

			paraIn.seekg(15L, ios_base::cur);
			paraIn.getline(str, 100);
			g_runNum = atoi(str);

			paraIn.seekg(14L, ios_base::cur);
			paraIn.getline(str, 100);
			g_fThr = g_fErr = atof(str);

			paraIn.seekg(19L, ios_base::cur);
			paraIn.getline(str, 100);
			g_dThr = atof(str);

			g_vBd.clear();
			g_vBd.resize(g_n);
			paraIn.seekg(12L, ios_base::cur);
			paraIn.getline(str, 100);
			for (auto i(g_vBd.begin()); i != g_vBd.end(); ++i)
				(*i).first = atof(str);
			paraIn.seekg(12L, ios_base::cur);
			paraIn.getline(str, 100);
			for (auto i(g_vBd.begin()); i != g_vBd.end(); ++i)
				(*i).second = atof(str);

			paraIn.close();
		}
	}

	g_y.resize(g_n);
	g_z.resize(g_n);

	string fNum = to_string(funcNum);
	string _str;

	if (funcNum <= 20)
	{
		fNum = to_string(funcNum);
		_str = "input/data/shift_data_" + fNum + ".txt";

		ifstream shiftDataIn(_str);

		if (!shiftDataIn.is_open())
		{
			cout << "The shift_data file is nonexistent, or cannot open the shift_data file for reading!" << endl;
			return false;
		}
		else
		{
			g_oShift = new double[g_n*g_cfNums[funcNum]];
			if (g_oShift == NULL)
			{
				cout << "Error: there is insufficient memory available!" << endl;
				return false;
			}

			for (auto i(0); i != g_cfNums[funcNum]; ++i)
			{
				for (auto j(0); j != g_n; ++j)
					shiftDataIn >> g_oShift[i*g_n + j];
			}

			shiftDataIn.close();
		}
	}

	if (funcNum <= 15)
	{
		/* Load Matrix M*/
		_str = "input/data/M_";
		_str.append(fNum + "_D");
		fNum = to_string(g_n);
		_str.append(fNum + ".txt");

		ifstream matrixIn(_str);

		if (!matrixIn.is_open())
		{
			cout << "The Matrix file is nonexistent, or cannot open the Matrix file for reading!" << endl;
			return false;
		}
		else
		{
			g_m = new double[g_cfNums[funcNum] * g_n*g_n];
			if (g_m == NULL)
			{
				cout << "Error: there is insufficient memory available!" << endl;
				return false;
			}

			for (auto i(0); i != g_cfNums[funcNum] * g_n*g_n; ++i)
				matrixIn >> g_m[i];

			matrixIn.close();
		}
	}

	return true;
}

void Output(int funcNum, int runNum)
{
	_mkdir("output");
	string path("output\\function_");
	string fNum = to_string(funcNum);
	path.append(fNum);
	_mkdir(path.c_str());

	path = "output/function_" + fNum + "/global optima.txt";

	if (0 == runNum)
	{
		ofstream gloOptPosAndFitOut(path);

		if (g_gloOpt.size() > 1)
			gloOptPosAndFitOut << "Find " << g_gloOpt.size() << " global optima" << endl;
		else
			gloOptPosAndFitOut << "Find " << g_gloOpt.size() << " global optimum" << endl;
		gloOptPosAndFitOut << setiosflags(ios::scientific) << setprecision(38);

		for (unsigned i(0); i != g_gloOpt.size(); ++i)
		{
			gloOptPosAndFitOut << "optimum_" << i + 1 << ":" << ends;
			for (auto j(0); j != g_n; ++j)
				gloOptPosAndFitOut << g_gloOpt.at(i).at(j) << ends;
			gloOptPosAndFitOut << endl << "fitness_" << i + 1 << ":" << ends << g_gloFit.at(i) << endl;
		}
		gloOptPosAndFitOut << endl;

		gloOptPosAndFitOut.close();
	}
	else
	{
		ofstream gloOptPosAndFitOut(path, ios_base::out | ios_base::app);

		if (g_gloOpt.size() > 1)
			gloOptPosAndFitOut << "Find " << g_gloOpt.size() << " global optima" << endl;
		else
			gloOptPosAndFitOut << "Find " << g_gloOpt.size() << " global optimum" << endl;
		gloOptPosAndFitOut << setiosflags(ios::scientific) << setprecision(38);

		for (unsigned i(0); i != g_gloOpt.size(); ++i)
		{
			gloOptPosAndFitOut << "optimum_" << i + 1 << ":" << ends;
			for (auto j(0); j != g_n; ++j)
				gloOptPosAndFitOut << g_gloOpt.at(i).at(j) << ends;
			gloOptPosAndFitOut << endl << "fitness_" << i + 1 << ":" << ends << g_gloFit.at(i) << endl;
		}
		gloOptPosAndFitOut << endl;

		gloOptPosAndFitOut.close();
	}

	g_gloFit_.push_back(g_gloFit);

	g_gloOpt.clear();
	g_gloFit.clear();
}

void Output_(int funcNum)
{
	auto sr(0.0);
	vector<double> nof(g_gloFit_.size(), 0.0);
	auto anof(0.0);
	auto _nof(0);

	for (unsigned i(0); i != g_gloFit_.size(); ++i)
	{
		for (unsigned j(0); j != g_gloFit_.at(i).size(); ++j)
			g_gloFit_.at(i).at(j) -= funcNum*100.0;
	}

	for (unsigned i(0); i != g_gloFit_.size(); ++i)
	{
		for (unsigned j(0); j != g_gloFit_.at(i).size(); ++j)
		{
			if (g_gloFit_.at(i).at(j) <= g_fErr)
			{
				nof.at(i) += 1.0;
				_nof++;
			}
		}

		if (1 == funcNum || 4 == funcNum || 10 == funcNum || 14 == funcNum || 15 == funcNum || 16 == funcNum || 17 == funcNum || 18 == funcNum || 19 == funcNum || 20 == funcNum)
		{
			if (1 == _nof)
				sr += 1.0;
		}
		else if (9 == funcNum || 11 == funcNum || 12 == funcNum || 13 == funcNum)
		{
			if (10 == _nof)
				sr += 1.0;
		}
		else if (2 == funcNum)
		{
			if (32 == _nof)
				sr += 1.0;
		}
		else if (3 == funcNum)
		{
			if (625 == _nof)
				sr += 1.0;
		}
		else if (5 == funcNum)
		{
			if (125 == _nof)
				sr += 1.0;
		}
		else if (6 == funcNum)
		{
			if (16 == _nof)
				sr += 1.0;
		}
		else if (7 == funcNum)
		{
			if (8 == _nof)
				sr += 1.0;
		}
		else if (8 == funcNum)
		{
			if (216 == _nof)
				sr += 1.0;
		}

		_nof = 0;
	}
	sr /= g_gloFit_.size();

	for (auto i(nof.begin()); i != nof.end(); ++i)
		anof += *i;
	anof /= g_gloFit_.size();

	auto stdNof(0.0);
	for (auto i(nof.begin()); i != nof.end(); ++i)
		stdNof += pow(*i - anof, 2);
	stdNof /= g_gloFit_.size();
	stdNof = sqrt(stdNof);

	vector<double> gloOptFit(g_gloFit_.size());
	double aver(0.0), stdDev(0.0), temp;

	for (unsigned i(0); i != g_gloFit_.size(); ++i)
	{
		temp = 0.0;

		for (unsigned j(0); j != g_gloFit_.at(i).size(); ++j)
			temp += g_gloFit_.at(i).at(j);

		temp /= g_gloFit_.at(i).size();
		gloOptFit.at(i) = temp;
		aver += temp;
	}
	aver /= g_gloFit_.size();

	for (unsigned i(0); i != gloOptFit.size(); ++i)
		stdDev += pow(gloOptFit.at(i) - aver, 2);
	stdDev = sqrt(stdDev / gloOptFit.size());

	string path;
	string fNum = to_string(funcNum);
	path = "output/function_" + fNum + "/result analysis.txt";
	ofstream resultAnalysisOut(path);
	resultAnalysisOut << setiosflags(ios::scientific) << setprecision(38);

	resultAnalysisOut << "SR: " << sr << endl << endl;
	resultAnalysisOut << "ANOF: " << anof << endl;
	resultAnalysisOut << "StdNof: " << stdNof << endl << endl;
	resultAnalysisOut << "Aver: " << aver << endl;
	resultAnalysisOut << "Std.: " << stdDev << endl;
	resultAnalysisOut.close();

	path = "output/function_" + fNum + "/runtime.txt";
	ofstream runTimeOut(path);
	runTimeOut << double(g_finish - g_start) / g_runNum << "ms" << endl;
	runTimeOut.close();
}

int _tmain(int argc, _TCHAR* argv[])
{
	srand((unsigned)(time(NULL) + rand()));

	if (!CalcFuncNum())
		return -1;

	for (auto i(0); i != g_totFuncNum; )
	{
		++i;

		if (!Input(i))
			return -1;

		g_start = clock();

		for (auto j(0); j != g_runNum; ++j)
		{
			Algorithm agrth(g_p);
			agrth.search(i, j);
		}

		g_finish = clock();

		Output_(i);

		g_gloFit_.clear();
		delete[] g_oShift;
		g_oShift = nullptr;
		delete[] g_m;
		g_m = nullptr;
		g_y.clear();
		g_z.clear();
	}

	system("pause");

	return 0;
}