#pragma once

#include "utils.h"

// Loading vector from file
void load_data_to(size_t n, std::vector<double>& input_v, std::ifstream& from_file, size_t event_number)
{
	if (event_number != 13)
	{
		for (size_t i = 0; i <= n; i++)
		{
			double value = 0.0;
			from_file >> value;
			input_v.push_back(value);
		}
	}
	else
	{
		if (event_number == 13)
		{
			for (size_t i = 1; i <= 1000; i++)
			{
				input_v.push_back(static_cast<double>(i));
			}
		}
		else
		{
			std::cout << "No such test\n";
		}
	}
}

void assembly(
	size_t n, size_t m,
	double gamma,
	std::vector<size_t>& ig,
	std::vector<double>& di,
	std::vector<double>& gg,
	std::vector<std::vector<double>>& B,
	std::vector<std::vector<double>>& C,
	std::vector<double>& F,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{  // Заполняем профиль
	ig[0] = 1;
	ig[1] = 1;
	for (size_t i = 2; i < m+1; i++)
	{
		ig[i] = ig[i - 1] + i - 1;
	}
	gg.resize(ig[m]-1);
	for (size_t i = 0; i < m; i++)
	{
		G[i] = 0.0;
		di[i] = 0.0;
	}

	for (size_t i = 1; i <= n; i++)
	{
		local_build(n, i - 1, gamma, B, C, F, r, event_number);
		for (size_t j = 0; j < 3; j++)
		{
			di[2 * (i - 1) + j] += B[j][j] + gamma * C[j][j];
			G[2 * (i - 1) + j] += F[j];
		}
		gg[ig[2 * i - 1] - 1] = B[1][0] + gamma * C[1][0];
		gg[ig[2 * i] - 1] = B[2][0] + gamma * C[2][0];
		gg[ig[2 * i]] = B[2][1] + gamma * C[2][1];
	}
}

void local_build(
	size_t n,
	size_t k,
	double gamma,
	std::vector<std::vector<double>>& B,
	std::vector<std::vector<double>>& C,
	std::vector<double>& F,
	std::vector<double>& r,
	size_t event_number
)
{
	std::vector<double> tmp(3);

	double h = r[k + 1] - r[k];
	double rk = r[k];

	tmp[0] = lambda(r[k], event_number);
	tmp[1] = lambda((r[k] + r[k + 1]) / 2.0, event_number);
	tmp[2] = lambda(r[k + 1], event_number);

	B[0][0] = tmp[0] * (2. / 15. * h + 37. / 30. * rk)
		+ tmp[1] * (1. / 3. * h + 6. / 5. * rk)
		+ tmp[2] * (1. / 30. * h + (-1. / 10.) * rk);

	B[0][1] = tmp[0] * ((-2. / 15.) * h + (-22. / 15.) * rk)
		+ tmp[1] * ((-4. / 15.) * h + (-16. / 15.) * rk)
		+ tmp[2] * ((-4. / 15.) * h + (-2. / 15.) * rk);

	B[0][2] = tmp[0] * (7. / 30. * rk)
		+ tmp[1] * ((-1. / 15.) * h + (-2. / 15.) * rk)
		+ tmp[2] * (7. / 30. * h + 7. / 30. * rk);

	B[1][1] = tmp[0] * (8. / 5. * rk)
		+ tmp[1] * (16. / 15. * h + 32. / 15. * rk)
		+ tmp[2] * (8. / 5. * h + 8. / 5. * rk);

	B[1][2] = tmp[0] * (2. / 15. * h + (-2. / 15.) * rk)
		+ tmp[1] * ((-4. / 5.) * h + (-16. / 15.) * rk)
		+ tmp[2] * ((-4. / 3.) * h + (-22. / 15.) * rk);

	B[2][2] = tmp[0] * ((-2. / 15.) * h + (-1. / 10.) * rk)
		+ tmp[1] * (13. / 15. * h + 6. / 5. * rk)
		+ tmp[2] * (11. / 10. * h + 37. / 30. * rk);

	B[1][0] = B[0][1];
	B[2][0] = B[0][2];
	B[2][1] = B[1][2];

	C[0][0] = 2. / 15. * rk + h / 60.;
	C[0][1] = rk / 15.;
	C[0][2] = -1. / 30. * rk + -1. / 60. * h;

	C[1][0] = C[0][1];
	C[1][1] = 8. / 15. * rk + 4. / 15. * h;
	C[1][2] = rk / 15. + h / 15.;

	C[2][0] = C[0][2];
	C[2][1] = C[1][2];
	C[2][2] = 2. / 15. * rk + 7. / 60. * h;

	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			B[i][j] /= h;
			C[i][j] *= h;
		}
	}

	tmp[0] = func(r[k], gamma, event_number);
	tmp[1] = func(((r[k] + r[k + 1]) / 2.), gamma, event_number);
	tmp[2] = func(r[k + 1], gamma, event_number);

	for (size_t i = 0; i < 3; i++)
	{
		F[i] = 0.;
		for (size_t j = 0; j < 3; j++)
		{
			F[i] += C[i][j] * tmp[j];
		}
	}
}

double func(double r, double gamma, size_t event_number)
{
	if (event_number == 1)
	{
		return -4.;
	}
	if (event_number == 2)
	{
		return gamma * r * r;
	}
	if (event_number == 3 || event_number == 4)
	{
		return -6. * r + gamma * r * r;
	}
	if (event_number >= 5 && event_number <= 8 || event_number == 13)
	{
		return -2. + gamma * r;
	}
	if (event_number == 9 || event_number == 10)
	{
		return -12. * r * r + gamma * r * r * r;
	}
	if (event_number == 11)
	{
		return -8. * r * r + gamma * r * r;
	}
	if (event_number == 12)
	{
		return -10. * r * r * r + gamma * r * r;
	}
	if (event_number > 13 || event_number == 0)
	{
		std::cout << "No such test";
	}
}

double lambda(double r, size_t event_number)
{
	if (event_number == 1)
	{
		return 1.0;
	}
	if (event_number == 2)
	{
		return 0.0;
	}
	if (event_number >= 3 && event_number <= 10 || event_number == 13)
	{
		return r;
	}
	if (event_number == 11)
	{
		return r * r;
	}
	if (event_number == 12)
	{
		return r * r * r;
	}
	if (event_number > 13 || event_number == 0)
	{
		std::cout << "No such test";
	}
}

double exact(double r, size_t event_number)
{
	if (event_number >= 1 && event_number <= 4 || event_number == 11 || event_number == 12)
	{
		return r * r;
	}
	if (event_number >= 5 && event_number <= 8 || event_number == 13)
	{
		return r;
	}
	if (event_number >= 9 && event_number <= 10)
	{
		return r * r * r;
	}
	if (event_number > 13 || event_number == 0)
	{
		std::cout << "No such test";
	}
}

void region(
	size_t n, size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	if (event_number >= 1 && event_number <= 5 || event_number == 9 || event_number == 10)
	{
		first_boundary_left(di, G, r, event_number);
		first_boundary_right(n, m, di, G, r, event_number);
	}
	if (event_number == 6)
	{
		first_boundary_left(di, G, r, event_number);
		second_boundary_right(m, di, G, r, event_number);
	}
	if (event_number == 7)
	{
		second_boundary_left(di, G, r, event_number);
		first_boundary_right(n, m, di, G, r, event_number);
	}
	if (event_number == 8)
	{
		third_boundary_left( di, G, r, event_number);
		second_boundary_right(m, di, G, r, event_number);
	}
	if (event_number == 11)
	{
		second_boundary_left(di, G, r, event_number);
		first_boundary_right(n, m, di, G, r, event_number);
	}
	if (event_number == 12)
	{
		second_boundary_left(di, G, r, event_number);
		first_boundary_right(n, m, di, G, r, event_number);
	}
	if (event_number == 13)
	{
		first_boundary_left(di, G, r, event_number);
		second_boundary_right(m, di, G, r, event_number);
	}
}

void SLAU_LLT(
	size_t n, size_t m,
	std::vector<size_t>& ig,
	std::vector<double>& di,
	std::vector<double>& gg,
	std::vector<double>& G
)
{
	LLT_decompose(n, G, gg, ig, di);

	//прямой ход
	forward_prop(n, G, gg, ig, di);

	//обратный ход
	backward_prop(n, m, G, gg, ig, di);
}

void LLT_decompose(
	size_t n,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di
)
{
	di[0] = sqrt(di[0]);
	for (size_t i = 1; i <= n; i++)
	{
		gg[ig[2 * i - 1] - 1] /= di[2 * i - 2];
		di[2 * i - 1] = sqrt(di[2 * i - 1] - (gg[ig[2 * i - 1] - 1] * gg[ig[2 * i - 1] - 1]));
		gg[ig[2 * i] - 1] /= di[2 * i - 2];
		gg[ig[2 * i]] = (gg[ig[2 * i]] - gg[ig[2 * i - 1] - 1] * gg[ig[2 * i] - 1]) / di[2 * i - 1];
		di[2 * i] = sqrt(di[2 * i] - gg[ig[2 * i] - 1] * gg[ig[2 * i] - 1] - gg[ig[2 * i]] * gg[ig[2 * i]]);
	}
}

void forward_prop(
	size_t n,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di
)
{
	G[0] /= di[0];
	for (size_t i = 1; i <= n; i++)
	{
		G[2 * i - 1] = (G[2 * i - 1] - gg[ig[2 * i - 1] - 1] * G[2 * i - 2]) / di[2 * i - 1];
		G[2 * i] = (G[2 * i] - gg[ig[2 * i] - 1] * G[2 * i - 2] - gg[ig[2 * i]] * G[2 * i - 1]) / di[2 * i];
	}
}

void backward_prop(
	size_t n,
	size_t m,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di)
{
	G[m - 1] /= di[m - 1];
	for (size_t i = n; i > 0; i--)
	{
		G[2 * i - 1] = (G[2 * i - 1] - gg[ig[2 * i]] * G[2 * i]) / di[2 * i - 1];
		G[2 * i - 2] = (G[2 * i - 2] - gg[ig[2 * i] - 1] * G[2 * i] - gg[ig[2 * i - 1] - 1] * G[2 * i - 1]) / di[2 * i - 2];
	}
}

void export_data_to(
	size_t n,
	size_t m,
	std::vector<double>& u,
	std::vector<double>& r,
	std::ofstream& output,
	size_t event_number
)
{
	output << r[0] << " " << exact(r[0], event_number) << " " << u[0] << "\n";
	for (size_t i = 1; i <= n; i++)
	{
		double rr = (r[i - 1] + r[i]) / 2.0;
		output << rr << " " << exact(rr, event_number) << " " << u[2 * i - 1] << "\n";
		output << r[i] << " " << exact(r[i], event_number) << " " << u[2 * i] << "\n";
	}
}