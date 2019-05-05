#include <vector>

using namespace std;

class Whale
{
public:
	Whale();

	void updateFit(int);
	void setPos(int, double);
	void _setPos(int, double);

	double getPos(int id)
	{
		return m_pos.at(id);
	}

	double getFit()
	{
		return m_fit;
	}

	vector<double> & getWhl()
	{
		return m_pos;
	}

	int m_count;

private:
	double m_fit;
	vector<double> m_pos;
};

class Algorithm
{
public:
	Algorithm(int);

	void search(int, int);
	void initialize(int);
	bool move(int, int);
	void calcBetterNearest(int);
	void judWthGloOpt(int);
	void rcdGloInSwarm();

	void outputOptima(int, int);

private:
	vector<Whale> m_whale;
	vector<vector<double>> m_gloOpt;
	vector<double> m_gloFit;
	double m_curGloFit;
	int m_trgWhl;
	unsigned long m_evalFlag;
};

double Distance(const vector<double> &, const vector<double> &);

template<typename T> inline double Rnd(T low, T uper)
{
	return ((double)rand() / RAND_MAX)*(uper - low) + low;
}