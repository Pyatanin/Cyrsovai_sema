#pragma once

#include <vector>
#include <fstream>
#include <numeric>
#include <cmath>
#include <iostream>

#include "boundary_cond.h"

// Loading vector from file
void load_data_to(size_t n, std::vector<double>& input_v, std::ifstream& from_file,size_t event_number);

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
);

void local_build(
	size_t n,
	size_t k,
	double gamma,
	std::vector<std::vector<double>>& B,
	std::vector<std::vector<double>>& C,
	std::vector<double>& F,
	std::vector<double>& r,
	size_t event_number
);

double func(double r, double gamma,size_t event_number);
double lambda(double r,size_t event_number);
double exact(double r,size_t event_number);

void region(
	size_t n, size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);

void SLAU_LLT(
	size_t n, size_t m,
	std::vector<size_t>& ig,
	std::vector<double>& di,
	std::vector<double>& gg,
	std::vector<double>& G
);

void LLT_decompose(
	size_t n,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di
);

void forward_prop(
	size_t n,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di
);

void backward_prop(
	size_t n,
	size_t m,
	std::vector<double>& G,
	std::vector<double>& gg,
	std::vector<size_t>& ig,
	std::vector<double>& di);

void export_data_to(
	size_t n, size_t m,
	std::vector<double>& u,
	std::vector<double>& r,
	std::ofstream& output,
	size_t event_number
);
