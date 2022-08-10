#include <iostream>
#include <vector>
#include <fstream>
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
 * @return The computed matrices, V as 'first' and W as 'second'
 */
std::pair<mat,mat> tridiag(const dimension k, const std::complex<item> t) {

	if (std::abs(t) >= 1) {
		std::cerr << "t must be a complex number with an absolute value less than 1."
			<< std::endl << "It is currently " << std::abs(t) << std::endl;
		return std::make_pair(mat::Zero(k, k),mat::Zero(k, k));
	}

	mat v = mat::Zero(k, k);
	mat w = mat::Zero(k, k);

	std::vector<std::thread> threads;

	for (dimension j = 1; j <= k; j++) {
		threads.push_back(std::thread(vMatrix, std::ref(v), j, k, t));
		threads.push_back(std::thread(wMatrix, std::ref(w), j, k, t));
	}

	for (auto& t : threads) {
		t.join();
	}

	// Multiply the matrices by the h constant.
	double h = (1 + std::pow(std::real(std::abs(t)), 2)) 
		/ std::pow((1 - std::pow(std::real(std::abs(t)), 2)), 2);
	v *= h;
	w *= h;

	return std::make_pair(v, w);
}

/*
 * Exports the eigenvalues of a matrix to an output file.
 * 
 * @param p The computed matrices of which to export the eigenvalues, V as 'first' and W as 'second'
 * @param k Max dimension of eigenvalues
 * @param filename Name of file to export eigenvalues to
 */
void exportEigens(const dimension kInit, item t, const std::string filename) {
	std::ofstream file;
	file.open(filename);

	for (dimension kIter = kInit; kIter > 0; kIter--) {
		std::pair<mat,mat> p = tridiag(kIter, t);

		vec vEigens = p.first.eigenvalues();
		vec wEigens = p.second.eigenvalues();

		// We are assuming that the two inputted vectors have the same dimensions
		for (int i = 0; i < vEigens.size(); i++) {
			file << "  " << std::to_string(std::real(vEigens(i))) << "  " << 0.0 << "\n";
		} for (int i = 0; i < wEigens.size(); i++) {
			file << "  "  << std::to_string(std::real(wEigens(i))) << "  " << 0.0 << "\n";
		}
	}
	file.close();
}

/*
 * Makes an asymmetrical matrix symmetrical.
 * 
 * param m Asymmetictrical tridiagonal matrix 
 * return Symmetrical tridiagonal matrix
 *
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
 */

int main(int argc, char** argv) {
	//Eigen::setNbThreads(48);

	std::vector<std::thread> threads;

	dimension kInit = 4;
	if (argc > 2) {
		kInit = std::atoi(argv[2]);
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
	for (item t = tInit; t >= -tIter; t -= tIter) {
		for (dimension kIter = kInit; kIter > 0; kIter--) {
			std::string strT = std::to_string(t);
			if (t < 0.0) {
				strT = std::to_string(0.0);
			}
			std::cout << strT.substr(0,8);
			std::cout << "\b\b\b\b\b\b\b\b\b";
			threads.push_back(std::thread(exportEigens, kIter, t,
						directory + std::to_string(kIter) + "k_" + strT + "t.txt"));
		}
	}	
	std::cout << std::endl;
	for (auto& t : threads) {
		t.join();
	}

	return 0;
}
