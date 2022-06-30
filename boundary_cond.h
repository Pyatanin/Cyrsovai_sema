#pragma once

#include <vector>

#include "utils.h"

void first_boundary_left(std::vector<double>& di,
	std::vector<double>& G, 
	std::vector<double>& r,
	size_t event_number
);

void first_boundary_right(
	size_t n,
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);

void second_boundary_left(
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);

void second_boundary_right(
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);

void third_boundary_left(
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);

void third_boundary_right(
	size_t n,
	size_t m,
	std::vector<double>& di,
	std::vector<double>& G,
	std::vector<double>& r,
	size_t event_number
);
