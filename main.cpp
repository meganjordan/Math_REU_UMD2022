#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cstdio>
#include <thread>
#include <mutex>

using namespace std;

namespace ev {
    // variable that controls the name of the output files.
    std::string filename = "evallist";
    std::mutex m;
}

/*
* Method that runs a "choose" function or a combination.
* Credit to rathbhupendra on geeksforgeeks.org for providing a low time-complexity method. O(n).
* https://www.geeksforgeeks.org/binomial-coefficient-dp-9/
* @param n Number of things being chosen from
* @param k Number of things selected from the set of size n
* @return Number of possible combinations of k things from a set of size n
*/
int nChoosek(int n, int k) {
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;

    // Calculate value of
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}

/*
* Finds an eigenvalue
* @param p Index of the matrix columns
* @param q Index of the matrix rows
* @param n Number of dimensions
* @return The eigenvalue
*/
int eval(int p, int q, int n) {
	return (2 * p * q) + ((n - 1)*(p + q));
}

/*
* Finds the multiplicity of an eigenvalue
* @param p Index of the matrix columns
* @param q Index of the matrix rows
* @param n Number of dimensions
* @return The multiplicity
*/
long double mult(int p, int q, int n) {
	return (double)(((p+q) / (n - 1)) + 1) 
        * (double)nChoosek(p+n-2, n-2) 
        * (double)nChoosek(q+n-2, n-2);
}

/*
* Finds the eigenvalue and the multiplicity
* @param p Index of the matrix columns
* @param q Index of the matrix rows
* @param n Number of dimensions
* @return A vector of the eigenvalue listed with multiplicity.
*/
std::vector<int> evalwmult(int p, int q, int n) {
    std::vector<int> v;
    long double multiplicity = mult(p, q, n);
    for (int i = 0; i < multiplicity; i++) {
        int temp = eval(p, q, n);
        v.push_back(temp);
        // Display values in command line
        //std::cout << temp << ", ";
    }
    //std::cout << std::endl;
    return v;
}

/*
* Loops through many p and q values in eigenvalue equation.
* @param n Number of dimensions
* @param loopVal Bounds of a loop through possible values of p and q. Loop goes from 0 to loopVal
* @return Vector of all discovered eigenvalues with mutliplicity
*/
std::unordered_map<int, size_t> evallist(int n, int loopVal) {
    std::vector<int> v;
    std::vector<std::thread> threads;
    std::unordered_map<int, size_t> umap;
    
    // Unused code for eigenvalue storage
    std::ofstream tmp;
    tmp.open(ev::filename + "_tmp.txt");

    // Section of code to be threaded
    auto f = [n, &v, &tmp, &umap](int i, int j) {
        std::vector<int> temp = evalwmult(i, j, n);

        ev::m.lock();

        //std::cout << i << "i " << j << "j - " << temp.at(0) << "," << temp.size() << " - ";

        for (int k = 0; k < temp.size(); k++) {
            tmp << temp.at(0) << ", ";
        }

        // Increase unordered map that associates a key eigenvalue with the value of the multiplicity
        umap[temp.at(0)] += temp.size();
        //std::cout << umap.at(temp.at(0)) << std::endl;

        ev::m.unlock();
    };

    for (int p = 0; p <= loopVal; p++) {
        for (int q = 0; 
            q <= (loopVal-(n-1)*p)/((2*p)+(n-1)); 
        q++) {
                threads.push_back(std::thread(f, p, q));
        }
    }

    for (auto& t : threads) {
        t.join();
    }

    tmp.close();
    return umap;
}

/*
* Unused method that gives the number of values in a vector within a certain range.
* Could be done more efficiently by sorting v and subtracting the first index of min from the last index of max.
* @param v Vector of values
* @param min Minimum value of range (inclusive)
* @param max Maximum value of range (inclusive)
* @return Number of values between min and max
*/
size_t nCount(std::vector<int> v, int min, int max) {
    int counter = 0;
    for (int i = 0; i < v.size(); i++) {
        if (v[i] <= max && v[i] >= min)
            counter++;
    }
    return counter;
}

int main() {
    std::cout << "What should the range of values (from 0 to N) be?" << std::endl;
    int n2 = 100;
    std::cin >> n2;
    std::cout << std::endl;

    std::cout << "What should the dimension (n) be?" << std::endl; 
    int dimension = 2;
    std::cin >> dimension;
    std::cout << std::endl;

    std::cout << "What should I name the CSV file?" << std::endl << "(Don't include the file type)" << std::endl;
    std::cin >> ev::filename;
    std::cout << std::endl;
    std::ofstream output, rawOutput;
    output.open("./output/" + ev::filename + ".csv");
    rawOutput.open("./output/" + ev::filename + "_raw.csv");

    // a map that counts the number of occurences of an eigenvalue.
    // in other words, key = eigenvalue, value = multiplicity
    std::unordered_map<int, size_t> originalUmap = evallist(dimension, n2);
    std::unordered_map<int, size_t> umap = originalUmap;
    
    // changes map to have the value represent the number of occurences of eigenvalues at or below the key eigenvalue.
    for (int i = n2; i >= 0; i--) {
        for (int j = 0; j < i; j++) {
            if (umap.find(j) != umap.end()) {
                umap[i] += umap.at(j);
                //std::cout << i << "," << j << " - " << umap.at(i) << "," << umap.at(j) << std::endl;
            }
        }
    }
    
    // prints umap and sends data to csv file.
    for (int i = 0; i <= n2; i++) {
        std::cout << i << "   -   ";
        output << i << ",";
        rawOutput << i << ",";

        // For step function
        if (umap.find(i) != umap.end()) {
            std::cout << umap.at(i);
            output << umap.at(i) << std::endl;
        } else {
            std::cout << 0;
            output << 0 << std::endl;
        }
        
        // For multiplicity
        if (originalUmap.find(i) != originalUmap.end()) {
            rawOutput << originalUmap.at(i) << std::endl;
        }
        else {
            rawOutput << 0 << std::endl;
        }
        std::cout << std::endl;
    }
    output.close();
    rawOutput.close();
    return 0;
}