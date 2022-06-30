#include "boundary_cond.h"

void first_boundary_left(
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	di[0] = 1.e+30;
	G[0] = exact(r[0], event_number) * 1.e+30;
}

void first_boundary_right(
	size_t n,
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	di[m - 1] = 1.e+30;
	G[m - 1] = exact(r[n], event_number) * 1.e+30;
}

void second_boundary_left(
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	double theta = 2.0;
	di[0] += theta;
}

void second_boundary_right(
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	double theta = -10.0;
	di[m - 1] += theta;
}

void third_boundary_right(
	size_t n,
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	double ub = 5.0;
	double beta = 2.0;
	di[m - 1] += ub * beta;
	G[m - 1] += exact(r[n], event_number) * beta;
}

void third_boundary_left(
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
)
{
	double ub = 1.5;
	double beta = 2.0;
	di[0] += ub * beta;
	G[0] += exact(r[0], event_number) * beta;
}