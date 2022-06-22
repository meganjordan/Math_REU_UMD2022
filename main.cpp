#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cstdio>
#include <complex>

#include <thread>
#include <mutex>

#include <Eigen/Dense>

using Eigen::MatrixX;
using item = long double;
using dimension = int;
using mat = Eigen::MatrixX<std::complex<item>>;
using namespace std::complex_literals;

void vMatrix(mat& v, const dimension j, const dimension k, const std::complex<item> t) {
	// Ensures values don't go past the max index of the matrix
	if (j != k) {
		// u sub j =
		// 2(2k-2j)(2j)
		// * 2(2k-2j+1)(2j-1)
		v(j - 1, j) = -t * (item)(
			(2 * ((2 * k) - (2 * j)) * (2 * j))
			* (2 * ((2 * k) - (2 * j) + 1) * ((2 * j) - 1))
			);
		v(j, j - 1) = -std::conj(t);
	}
	// Is in the center, index goes further with this function
	// d sub j = 
	// 2(2k-2j+1)(2j-1)
	// + |t|^2
	// * 2(2j-2)(2k-2j+2)
	v(j - 1, j - 1) = (item)(
		(2 * ((2 * k) - (2 * j) + 1) * ((2 * j) - 1))
		+ std::pow(std::abs(t), 2)
		* (2 * ((2 * j) - 2) * ((2 * k) - (2 * j) + 2))
		);
}

void wMatrix(mat& w, const dimension j, const dimension k, const std::complex<item> t) {
	// Ensures values don't go past the max index of the matrix
	if (j != k) {
		// u sub j =
		// 2(2k-2j)(2j)
		// * 2(2k-2j-1)(2j+1)
		w(j - 1, j) = -t * (item)(
			(2 * ((2 * k) - (2 * j)) * (2 * j))
			* (2 * ((2 * k) - (2 * j) - 1) * ((2 * j) + 1))
			);
		// l sub j =
		// - t bar
		w(j, j - 1) = -std::conj(t);
	}
	// Is in the center, index goes further with this function
	// d sub j = 
	// 2(2k-2j+1)(2j-1)
	// * |t|^2
	// + 2(2k-2j)(2j)
	w(j - 1, j - 1) = (item)(
		((2 * ((2 * k) - (2 * j) + 1) * ((2 * j) - 1))
			* std::pow(std::abs(t), 2))
		+ (2 * ((2 * k) - (2 * j)) * (2 * j))
		);
}

/*
* @param k the dimension of the matrix
* @param t the time (complex)
*/
std::pair<mat,mat> tridiag(dimension k, std::complex<item> t, bool useV) {
	if (std::abs(t) >= 1) {
		std::cerr << "t must be a complex number with an absolute value less than 1."
			<< std::endl << "It is currently " << std::abs(t) << std::endl;
		return std::make_pair(mat::Zero(k, k), mat::Zero(k, k));
	}

	mat v = mat::Zero(k, k);
	mat w = mat::Zero(k, k);

	std::vector<std::thread> vThreads, wThreads;

	for (dimension j = 1; j <= k; j++) {
		if (useV) {
			vThreads.push_back(std::thread(vMatrix, std::ref(v), j, k, t));
		}
		else {
			wThreads.push_back(std::thread(wMatrix, std::ref(w), j, k, t));
		}
	}

	for (auto& t : vThreads) {
		t.join();
	}
	for (auto& t : wThreads) {
		t.join();
	}

	return std::make_pair(v,w);
}

mat toSymmetrical(mat m) {
	std::vector<std::thread> threads;

	auto f = [&m](const dimension j) {
		m(j, j - 1) = std::pow((m(j - 1, j) * m(j, j - 1)), 0.5);
		m(j - 1, j) = m(j, j - 1);
	};

	for (dimension j = 1; j < m.cols(); j++) {
		threads.push_back(std::thread(f, j));
	}

	for (auto& t : threads) {
		t.join();
	}

	return m;
}

int main() {
	Eigen::setNbThreads(48);

	dimension k = 28;
	std::complex<double> t = 0.25i;
	std::pair<mat, mat> m = tridiag(k, t, true);
	std::cout << m.first.eigenvalues() << std::endl;
}