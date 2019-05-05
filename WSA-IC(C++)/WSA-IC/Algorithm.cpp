#include "stdafx.h"
#include "Algorithm.h"
#include <map>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>

// =========================================================

extern int g_n;
extern int g_p;
extern unsigned long g_evalNum;

extern vector<pair<double, double>> g_vBd;
extern double g_dThr;							// distance threshold

extern vector<vector<double>> g_gloOpt;
extern vector<double> g_gloFit;

// =========================================================

// parameters
int g_sThr;										// stability threshold
double g_fThr;									// fitness threshold

// =========================================================

extern void Func(vector<double> &, double &, int, int);
extern void Output(int, int);

//=========================================================

Whale::Whale() :m_pos(g_n)
{
}

void Whale::updateFit(int funFlag)
{
	Func(m_pos, m_fit, g_n, funFlag);
}

void Whale::setPos(int id, double value)
{
	if (value < g_vBd.at(id).first)
		m_pos.at(id) = g_vBd.at(id).first;
	else if (value > g_vBd.at(id).second)
		m_pos.at(id) = g_vBd.at(id).second;
	else
		m_pos.at(id) = value;
}

void Whale::_setPos(int id, double value)
{
	m_pos.at(id) = value;
}

// =========================================================

Algorithm::Algorithm(int whlNum) :m_whale(whlNum), m_curGloFit(DBL_MAX), m_evalFlag(0)
{
}

void Algorithm::search(int funFlag, int runFlag)
{
	for (auto i(0); i != g_p; ++i)
	{
		initialize(i);
		m_whale.at(i).updateFit(funFlag);
	}

	while (m_evalFlag != g_evalNum)
	{
		for (auto i(0); i != g_p; ++i)
		{
			calcBetterNearest(i);

			if (g_p != m_trgWhl)
			{
				if (!move(i, funFlag))
					break;
			}
			else
			{
				if (m_whale.at(i).m_count != g_sThr)
					m_whale.at(i).m_count++;
				else
				{
					judWthGloOpt(i);
					
					initialize(i);
					m_whale.at(i).updateFit(funFlag);

					if (++m_evalFlag == g_evalNum)
						break;
				}
			}
		}
	}

	rcdGloInSwarm();
	outputOptima(funFlag, runFlag);
}

void Algorithm::initialize(int whlId)
{
	static double random(0.0);

	for (auto i(0); i != g_n; ++i)
	{
		random = Rnd(g_vBd.at(i).first, g_vBd.at(i).second);
		m_whale.at(whlId)._setPos(i, random);
	}

	m_whale.at(whlId).m_count = 0;
}

bool Algorithm::move(int whlId, int funFlag)
{
	static auto temp(0.0);

	Whale whlTemp(m_whale.at(whlId));

	for (auto i(0); i != g_n; ++i)
	{
		temp = whlTemp.getPos(i) + Rnd(0.0, 2.0) * (m_whale.at(m_trgWhl).getPos(i) - whlTemp.getPos(i));
		whlTemp.setPos(i, temp);
	}

	whlTemp.updateFit(funFlag);

	if (whlTemp.getFit() < m_whale.at(whlId).getFit())
	{
		m_whale.at(whlId) = whlTemp;
		m_whale.at(whlId).m_count = 0;
	}
	else
	{
		if (m_whale.at(whlId).m_count != g_sThr)
			m_whale.at(whlId).m_count++;
		else
		{
			judWthGloOpt(whlId);

			if (++m_evalFlag != g_evalNum)
			{
				initialize(whlId);
				m_whale.at(whlId).updateFit(funFlag);
			}
			else
				return false;
		}
	}

	if (++m_evalFlag == g_evalNum)
		return false;
	else
		return true;
}

void Algorithm::calcBetterNearest(int whlId)
{
	static double distTemp, betterNearest;
	betterNearest = DBL_MAX;

	m_trgWhl = g_p;

	for (auto i(0); i != g_p; ++i)
	{
		if (m_whale.at(i).getFit() < m_whale.at(whlId).getFit())
		{
			distTemp = Distance(m_whale.at(i).getWhl(), m_whale.at(whlId).getWhl());

			if (distTemp < betterNearest)
			{
				m_trgWhl = i;
				betterNearest = distTemp;
			}
		}
	}
}

void Algorithm::judWthGloOpt(int whlId)
{
	if (m_whale.at(whlId).getFit() < m_curGloFit)
	{
		if (m_curGloFit - m_whale.at(whlId).getFit() > g_fThr)
		{
			m_gloOpt.clear();
			m_gloFit.clear();
		}

		m_curGloFit = m_whale.at(whlId).getFit();

		m_gloOpt.push_back(m_whale.at(whlId).getWhl());
		m_gloFit.push_back(m_whale.at(whlId).getFit());
	}
	else
	{
		if (m_whale.at(whlId).getFit() - m_curGloFit <= g_fThr)
		{
			m_gloOpt.push_back(m_whale.at(whlId).getWhl());
			m_gloFit.push_back(m_whale.at(whlId).getFit());
		}
	}
}

void Algorithm::rcdGloInSwarm()
{
	for (auto i(0); i != g_p; ++i)
	{
		if (m_whale.at(i).getFit() < m_curGloFit)
		{
			if (m_curGloFit - m_whale.at(i).getFit() > g_fThr)
			{
				m_gloOpt.clear();
				m_gloFit.clear();
			}

			m_curGloFit = m_whale.at(i).getFit();

			m_gloOpt.push_back(m_whale.at(i).getWhl());
			m_gloFit.push_back(m_whale.at(i).getFit());
		}
		else
		{
			if (m_whale.at(i).getFit() - m_curGloFit <= g_fThr)
			{
				m_gloOpt.push_back(m_whale.at(i).getWhl());
				m_gloFit.push_back(m_whale.at(i).getFit());
			}
		}
	}
}

void Algorithm::outputOptima(int funFlag, int runFlag)
{
	bool flag;

	for (auto i(0); i != m_gloOpt.size(); ++i)
	{
		if (m_gloFit.at(i) - m_curGloFit <= g_fThr)
		{
			flag = true;

			for (auto j(0); j != g_gloOpt.size(); ++j)
			{
				if (Distance(m_gloOpt.at(i), g_gloOpt.at(j)) <= g_dThr)
				{
					if (m_gloFit.at(i) < g_gloFit.at(j))
					{
						g_gloOpt.at(j) = m_gloOpt.at(i);
						g_gloFit.at(j) = m_gloFit.at(i);
					}

					flag = false;
					break;
				}
			}

			if (flag)
			{
				g_gloOpt.push_back(m_gloOpt.at(i));
				g_gloFit.push_back(m_gloFit.at(i));
			}
		}
	}

	Output(funFlag, runFlag);
}

double Distance(const vector<double> &pos0, const vector<double> &pos1)
{
	auto dist(0.0);

	for (auto i(0); i != g_n; ++i)
		dist += pow((pos0.at(i) - pos1.at(i)), 2);
	dist = sqrt(dist);

	return dist;
}