#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cstdio>
#include <complex>

#include <thread>
#include <mutex>

#include <Eigen/Eigenvalues>

using Eigen::MatrixX;
using item = long double;
using dimension = int;
using mat = Eigen::MatrixX<std::complex<item>>;
using vec = Eigen::VectorX<std::complex<item>>;
using namespace std::complex_literals;

/*
* Changes cells in the W matrix according to the 2017 research
* http://dx.doi.org/10.2140/involve.2019.12.125
* 
* Only a separate method for the sake of multithreading, do not call on its own.
* 
* @param v V Matrix
* @param j Iterator which indicates which cells to modify
* @param k Scale of the matrix
* @param t Time
*/
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

/*
* Changes cells in the W matrix according to the 2017 research
* http://dx.doi.org/10.2140/involve.2019.12.125
* 
* Only a separate method for the sake of multithreading, do not call on its own.
* 
* @param w W Matrix
* @param j Iterator which indicates which cells to modify
* @param k Scale of the matrix
* @param t Time
*/
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
* Forms a tridiagonal matrix according to the 2017 research
* http://dx.doi.org/10.2140/involve.2019.12.125
* 
* @param k Scale of the matrix
* @param t Time (complex)
* @param useV Decides whether to create the V Matrix or the W Matrix. If true, V. If false, W.
* @return The computed matrix
*/
mat tridiag(const dimension k, const std::complex<item> t, const bool useV) {

	if (std::abs(t) >= 1) {
		std::cerr << "t must be a complex number with an absolute value less than 1."
			<< std::endl << "It is currently " << std::abs(t) << std::endl;
		return mat::Zero(k, k);
	}

	mat m = mat::Zero(k, k);

	std::vector<std::thread> threads;

	for (dimension j = 1; j <= k; j++) {
		if (useV) {
			threads.push_back(std::thread(vMatrix, std::ref(m), j, k, t));
		}
		else {
			threads.push_back(std::thread(wMatrix, std::ref(m), j, k, t));
		}
	}

	for (auto& t : threads) {
		t.join();
	}

	// Multiply the matrix by the h constant.
	m *= (1 + std::pow(std::real(std::abs(t)), 2)) 
		/ std::pow((1 - std::pow(std::real(std::abs(t)), 2)), 2);

	return m;
}

/*
* Exports the eigenvalues of a matrix to an output file.
* 
* @param m Matrix of which to export the eigenvalues
* @param filename Name of file to export eigenvalues to
*/
void exportEigens(const mat m, const mat n, item t, const std::string filename) {
	std::ofstream file;
	file.open(filename);

	vec e = m.eigenvalues();
	vec f = n.eigenvalues();

	for (int i = 0; i < e.size(); i++) {
		file << "  " << std::to_string(std::real(e(i))) << "  " << 0.0 << "\n"
		<< "  "  << std::to_string(std::real(f(i))) << "  " << 0.0 << "\n";
	}

	file.close();
}

/*
* Makes an asymmetrical matrix symmetrical.
* 
* @param m Asymmetictrical tridiagonal matrix 
* @return Symmetrical tridiagonal matrix
*/
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

int main(int argc, char** argv) {
	//Eigen::setNbThreads(48);
	
	std::vector<std::thread> threads;
	
	dimension k = 4;
	if (argc > 2) {
		k = std::atoi(argv[2]);
	}

	std::string directory;
	if (argc > 1) {
		directory = argv[1];
	}

	// How many t values to view eigenvalues of
	item tIter = 0.001, tInit = 0.7;
	if (argc > 3) {
		tIter = std::atof(argv[3]);
	}
	if (argc > 4)
	{
		tInit = std::atof(argv[4]);
	}
	for (item t = tInit; t >= 0.0; t -= tIter) {
		threads.push_back(std::thread(exportEigens, tridiag(k, t, true), tridiag(k, t, false), t,
			directory + std::to_string(k) + "k - " + std::to_string(t) + "t.txt"));
	}

	for (auto& t : threads) {
		t.join();
	}

	exportEigens(tridiag(k, 0.0, true), tridiag(k, 0.0, false), t, directory + std::to_string(k) + "k - 0.000000t.txt");

	return 0;
}
