/* Copyright (c) 2020, Dyssol Development Team. All rights reserved. This file is part of Dyssol. See LICENSE file for license information. */

#pragma once

#include <fstream>

template<typename T>
bool writerCSV1D(const std::vector<T>& _vec, const std::string _name, std::ofstream& _csvfile)
{
	if (!_csvfile.is_open()) return false;

	_csvfile << _name << "\n";
	for (size_t i = 0; i < _vec.size(); ++i)
	{
		_csvfile << std::to_string(_vec[i]) << "\n";
	}
	return true;
}

inline bool writerCSV2D(const std::vector<std::vector<std::string>>& _matr, const std::string& _name, std::ofstream& _csvfile)
{
	if (!_csvfile.is_open()) return false;

	_csvfile << _name << "\n";
	for (size_t i = 0; i < _matr.size(); ++i)
	{
		for (size_t j = 0; j < _matr[i].size(); ++j)
		{
			_csvfile << _matr[i][j];
			if (j < _matr[i].size() - 1)
			{
				_csvfile << ";";
			}
		}
		_csvfile << "\n";
	}
	return true;
}

inline bool writerCSV2D(const std::vector<std::vector<double>>& _matr, const std::string& _name, std::ofstream& _csvfile)
{
	if (!_csvfile.is_open()) return false;

	_csvfile << _name << "\n";
	for (size_t i = 0; i < _matr.size(); ++i)
	{
		for (size_t j = 0; j < _matr[i].size(); ++j)
		{
			_csvfile << std::to_string(_matr[i][j]);
			if (j < _matr[i].size() - 1)
			{
				_csvfile << ";";
			}
		}
		_csvfile << "\n";
	}
	return true;
}

inline void transpose(std::vector<std::vector<std::string>>& _matr)
{
	if (_matr.empty())
		return;

	std::vector<std::vector<std::string>> trans_vec(_matr[0].size(), std::vector<std::string>());

	for (size_t i = 0; i < _matr.size(); i++)
	{
		for (size_t j = 0; j < _matr[i].size(); j++)
		{
			trans_vec[j].push_back(_matr[i][j]);
		}
	}

	_matr = trans_vec;    // <--- reassign here
}

inline std::vector<std::string> vecDouble2Str(const std::vector<double>& _vec)
{
	std::vector<std::string> vecStr;
	vecStr.reserve(_vec.size());
	for (const double d : _vec)
		vecStr.push_back(std::to_string(d));
	return vecStr;
}

inline std::vector<std::string> push_back_vec(const std::vector<std::string>& _vec1, const std::vector<std::string>& _vec2)
{
	std::vector<std::string> vVec = _vec1;
	for (const auto& v : _vec2)
		vVec.push_back(v);
	return vVec;
}

inline std::vector<std::string> strVector(const std::string& _name, const std::vector<double>& _vec)
{
	std::vector<std::string> res;
	res.push_back(_name);
	for (const double d : _vec)
		res.push_back(std::to_string(d));
	return res;
}

inline std::string BoolToOnOff(bool _var)
{
	return _var ? "on" : "off";
}

template<typename T>
void Clear1D(std::vector<T>& _vec)
{
	_vec.clear();
}

template<typename T>
void Clear2D(std::vector<std::vector<T>>& _vec)
{
	for (auto& v1 : _vec)
	{
		v1.clear();
	}
	_vec.clear();
}

template<typename T>
void Clear3D(std::vector<std::vector<std::vector<T>>>& _vec)
{
	if (!_vec.empty())
	{
		for (auto& v1 : _vec)
		{
			for (auto& v2 : v1)
			{
				v2.clear();
			}
			v1.clear();
		}
		_vec.clear();
	}
}

template<typename T>
void Clear4D(std::vector<std::vector<std::vector<std::vector<T>>>>& _vec)
{
	for (auto& v1 : _vec)
	{
		for (auto& v2 : v1)
		{
			for (auto& v3 : v2)
			{
				v3.clear();
			}
			v2.clear();
		}
		v1.clear();
	}
	_vec.clear();
}

/// Vector handling ///
inline void Resize3D(std::vector<std::vector<std::vector<double>>>& _vec, size_t _sz1, size_t _sz2, size_t _sz3)
{
	Clear3D(_vec);
	_vec.resize(_sz1, std::vector<std::vector<double>>(_sz2, std::vector<double>(_sz3, 0.0)));
}

inline void Resize4D(std::vector<std::vector<std::vector<std::vector<double>>>>& _vec, size_t _sz1, size_t _sz2, size_t _sz3, size_t _sz4)
{
	Clear4D(_vec);
	_vec.resize(_sz1, std::vector<std::vector<std::vector<double>>>(_sz2, std::vector<std::vector<double>>(_sz3, std::vector<double>(_sz4, 0.0))));
}

/// Vector functions ///
inline void GetValues1D(const double* _vals, std::vector<double>& _vec, const std::vector<size_t>& _ind, const size_t _subGridPoints = 0)
{
	if (_ind.empty())
		return;

	if (!_vec.empty())
		_vec.clear();

	const size_t gridPoints = _ind.size();

	for (size_t i = 0; i < gridPoints - 1; ++i)
	{
		const double val_i0 = _vals[_ind[i]];
		const double val_i1 = _vals[_ind[i + 1]];
		const double C = (val_i1 - val_i0) / static_cast<double>(_subGridPoints + 1);
		for (size_t j = 0; j < _subGridPoints + 1; ++j)
		{
			double val = val_i0 + static_cast<double>(j) * C;
			_vec.push_back(val);
		}
	}
	// last element
	_vec.push_back(_vals[_ind[gridPoints - 1]]);
}

inline void SetValues1D(double* _vals, const std::vector<double>& _vec, const std::vector<size_t>& _ind, const size_t& _subGridPoints = 0)
{
	if (_ind.empty() || _vec.empty())
		return;

	const size_t gridPointsFine = _vec.size();
	const size_t gridPoints = _ind.size();

	for (size_t i = 0; i < gridPoints - 1; i++)
	{
		const size_t iFine = i * (_subGridPoints + 1);
		_vals[_ind[i]] = _vec[iFine];
	}
	// last element
	_vals[_ind[gridPoints - 1]] = _vec[gridPointsFine - 1];
}

inline void GetValues3D(const double* _vals, std::vector<std::vector<std::vector<double>>>& _vec, const std::vector<std::vector<std::vector<size_t>>>& _ind, const size_t _subGridPoints = 0)
{
	if (_ind.empty())
		return;

	if (!_vec.empty())
		Clear3D(_vec);

	for (size_t i = 0; i < _ind.size(); ++i)
	{
		_vec.emplace_back();
		for (size_t j = 0; j < _ind[i].size(); ++j)
		{
			_vec[i].emplace_back();
			GetValues1D(_vals, _vec[i][j], _ind[i][j], _subGridPoints);
		}
	}
}

inline void SetValues3D(double* _vals, const std::vector<std::vector<std::vector<double>>>& _vec, const std::vector<std::vector<std::vector<size_t>>>& _ind, const size_t _subGridPoints = 0)
{
	if (_ind.empty() || _vec.empty())
		return;

	for (size_t i = 0; i < _ind.size(); ++i)
	{
		for (size_t j = 0; j < _ind[i].size(); ++j)
		{
			SetValues1D(_vals, _vec[i][j], _ind[i][j], _subGridPoints);
		}
	}
}

inline std::vector<double> Linspace(double _min, double _max, size_t _num)
{
	std::vector<double> vec(_num);
	for (size_t i = 0; i < _num; ++i)
	{
		vec[i] = static_cast<double>(i) / static_cast<double>(_num - 1) * (_max - _min) + _min;
	}
	return vec;
}

inline std::vector<double> Logspace(int _min, int _max, size_t _num)
{
	std::vector<double> vec = Linspace(_min, _max, _num);
	for (size_t i = 0; i < _num; ++i)
	{
		vec[i] = pow(10, vec[i]);
	}
	return vec;
}

inline std::vector<double> VecGamma(const std::vector<double>& _t, double _n)
{
	std::vector<double> gamma(_t.size());
	for (size_t i = 0; i < _t.size(); ++i)
	{
		gamma[i] = pow(_t[i], _n) * exp(-_t[i]);
	}
	return gamma;
}

inline double IntegrateTrapez(const std::vector<double>& _grid, const std::vector<double>& _x)
{
	double Int = 0;
	for (size_t i = 0; i < _grid.size() - 1; ++i)
	{
		Int += 0.5 * (_x[i] + _x[i + 1]) * (_grid[i + 1] - _grid[i]);
	}
	return Int;
}

inline double Gamma(double _k)
{
	const std::vector<double> vec = Logspace(-12, 3, 10000);
	const std::vector<double> gamma = VecGamma(vec, _k - 1);
	return IntegrateTrapez(vec, gamma);
}

inline std::vector<double> TISModel(std::vector<double> _tau, double _tau_m, double _k, double _gamma)
{
	std::vector<double> vec(_tau.size());
	vec[0] = 1e-6;
	for (size_t i = 1; i < _tau.size(); ++i)
	{
		vec[i] = 1.0 / _tau_m * pow(_tau[i] / _tau_m, _k - 1) * pow(_k, _k) / _gamma * exp(-_tau[i] * _k / _tau_m);
	}
	return vec;
}

inline double IntegrateTIS(double _tau_max, double _tau_m, double _k)
{
	const double gamma = Gamma(_k);
	const std::vector<double> tau = Linspace(0, _tau_max, 100);
	const std::vector<double> TIS = TISModel(tau, _tau_m, _k, gamma);
	return IntegrateTrapez(tau, TIS);
}