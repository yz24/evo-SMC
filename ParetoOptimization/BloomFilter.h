#include <iostream>
#include <bitset>
#include <random>
#include <set>

using namespace std;

const int FILTER_SIZE = 10000;
const int NUM_HASHES = 3;

class BloomFilter {

public:
    BloomFilter() {
        filter.reset();
    }

    void insert(set<int> s) {
        for (auto i : s) {
            for (int j = 0; j < NUM_HASHES; ++j) {
                filter[hash(i, j)] = true;
            }
        }
    }

    bool contains(set<int> s) const {
        bitset<FILTER_SIZE> temp_filter;
        temp_filter.reset();

        for (auto i : s) {
            for (int j = 0; j < NUM_HASHES; ++j) {
                temp_filter[hash(i, j)] = true;
            }
        }

        // Check if any of the bits in the temporary filter match the Bloom filter
        return (temp_filter & filter) == temp_filter;
    }

    void clear() {
        filter.reset();
    }

private:
    bitset<FILTER_SIZE> filter;
    hash<int> h;

    int hash(int i, int seed) const {
        // Use a random seed for each hash function to minimize collisions
        std::mt19937 mt(seed);
        return h(i ^ mt());
    }
};

//int main() {
//    BloomFilter bf;
//
//    // Insert some sets of integers into the Bloom filter
//    set<int> myset1 = {42, 1337, 9001};
//    bf.insert(myset1);
//
//    set<int> myset2 = {1234, 5678, 9012};
//    bf.insert(myset2);
//
//    // Check if some sets of integers match any of the sets already inserted
//    set<int> myset3 = {42, 1337, 9001};
//    cout << bf.contains(myset3) << endl;  // Output: 1
//
//    set<int> myset4 = {42, 1337};
//    cout << bf.contains(myset4) << endl;  // Output: 0
//
//    set<int> myset5 = {1234, 5678, 9012};
//    cout << bf.contains(myset5) << endl;  // Output: 1
//
//    set<int> myset6 = {1234, 5678};
//    cout << bf.contains(myset6) << endl;  // Output: 0
//
//    // Clear the Bloom filter
//    bf.clear();
//}
