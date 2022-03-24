#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>
using namespace std;
ifstream f ("teste.txt");
long long max(vector<long long> v, long long n)
{
    long long i, maxim = 0;
    for(i = 0; i < n; i++)
        if(v[i] > maxim)
        maxim = v[i];
    return maxim;
}
void countingSort(vector<long long> &v, long long n) {
    long long maxim = max(v,n), i;
    vector<long long> output(n, 0), freq(maxim + n, 0);
    for (i = 0; i < n; i++) {
        freq[v[i]]++;
    }
    for (i = 1; i <= maxim; i++) {
        freq[i] += freq[i - 1];
    }
    for (i = n - 1; i >= 0; i--) {
        output[freq[v[i]] - 1] = v[i];
        freq[v[i]]--;
    }
    v.clear();
    for(auto it: output)
        v.push_back(it);
    output.clear();
}
int e_sortat(vector <long long> v, vector <long long> v1, int dimensiune)
{

    for(long long i = 0; i < dimensiune; i++)
        if(v[i] != v1[i])
            return 0;
    return 1;
}
void merge(vector<long long> &vec, vector<long long> &v1, vector<long long> &v2) {
    long long p1 = 0, p2 = 0;
    while (p1 < v1.size() && p2 < v2.size()) {
        if (v1.at(p1) < v2.at(p2))
            vec.push_back(v1.at(p1++));
        else
            vec.push_back(v2.at(p2++));
    }
    while (p1 < v1.size()) {
        vec.push_back(v1.at(p1++));
    }
    while (p2 < v2.size()) {
        vec.push_back(v2.at(p2++));
    }
}

void mergeSort(vector<long long> &v) {
    if (v.size() <= 1)
        return;

    auto mijloc = v.begin() + v.size() / 2;
    vector<long long> v1(v.begin(), mijloc);
    vector<long long> v2(mijloc, v.end());
    mergeSort(v1);
    mergeSort(v2);
    v.resize(0);
    merge(v, v1, v2);

}

void shellSort(vector<long long> &v, long long dimensiune)
{
    for(long long gap = dimensiune/2; gap > 0; gap /= 2)
    {
        for(long long i = gap; i < dimensiune; i++)
        {
            long long temp = v[i], j;
            for(j = i; j >= gap && v[j-gap] > temp; j -= gap)
                v[j] = v[j-gap];
            v[j] = temp;
        }
    }
}


void countSort(vector <long long> &v, long long dimensiune, long long exp) {
    long long i;
    vector<long long> output(dimensiune, 0), freq(10, 0);
    for (i = 0; i < dimensiune; i++) {
        freq[(v[i] / exp) % 10]++;
    }
    for (i = 1; i < 10; i++) {
        freq[i] += freq[i - 1];
    }
    for (i = dimensiune - 1; i >= 0; i--) {
        output[freq[(v[i] / exp) % 10] - 1] = v[i];
        freq[(v[i] / exp) % 10]--;
    }
    v.resize(0);
    for (i = 0; i < dimensiune; i++)
        v.push_back(output[i]);
}
void radixsort(vector <long long> &v, long long dimensiune)
{
    long long maxim = v[0], i, exp;
    for (i = 1; i < dimensiune; i++)
        if(v[i] > maxim) maxim = v[i];
    for (exp = 1; maxim / exp > 0; exp *= 10)
        countSort(v, dimensiune, exp);
}

void countSort2(vector <long long> &v, long long dimensiune, long long exp) {
    long long i;
    vector<long long> output(dimensiune, 0), freq(10, 0);
    for (i = 0; i < dimensiune; i++) {
        freq[(v[i] / exp) % 2]++;
    }
    for (i = 1; i < 2; i++) {
        freq[i] += freq[i - 1];
    }
    for (i = dimensiune - 1; i >= 0; i--) {
        output[freq[(v[i] / exp) % 2] - 1] = v[i];
        freq[(v[i] / exp) % 2]--;
    }
    v.resize(0);
    for (i = 0; i < dimensiune; i++)
        v.push_back(output[i]);
}
void radixsort2(vector <long long> &v, long long dimensiune)
{
    long long maxim = v[0], i, exp;
    for (i = 1; i < dimensiune; i++)
        if(v[i] > maxim) maxim = v[i];
    for (exp = 1; maxim / exp > 0; exp *= 2)
        countSort2(v, dimensiune, exp);
}
void countSort3(vector <long long> &v, long long dimensiune, long long exp) {
    long long i, p = pow(2,16);
    vector<long long> output(dimensiune, 0), freq(p, 0);
    for (i = 0; i < dimensiune; i++) {
        freq[(v[i] / exp) % p]++;
    }
    for (i = 1; i < p; i++) {
        freq[i] += freq[i - 1];
    }
    for (i = dimensiune - 1; i >= 0; i--) {
        output[freq[(v[i] / exp) % p] - 1] = v[i];
        freq[(v[i] / exp) % p]--;
    }
    v.resize(0);
    for (i = 0; i < dimensiune; i++)
        v.push_back(output[i]);
}
void radixsort3(vector <long long> &v, long long dimensiune)
{
    long long maxim = v[0], i, exp;
    for (i = 1; i < dimensiune; i++)
        if(v[i] > maxim) maxim = v[i];
    for (exp = 1; maxim / exp > 0; exp *= pow(2, 16))
        countSort3(v, dimensiune, exp);
}

void heapify(vector<long long> &v, long long dimensiune, long long i) {
    long long largest = i, left = 2 * i + 1, right = 2 * i + 2;
    if (left < dimensiune && v[left] > v[largest]) {
        largest = left;
    }
    if (right < dimensiune && v[right] > v[largest]) {
        largest = right;
    }
    if (largest != i) {
        swap(v[i], v[largest]);
        heapify(v, dimensiune, largest);
    }
}
void heapSort(vector<long long> &v, long long dimensiune) {
    long long i;
    for (i = dimensiune / 2 - 1; i >= 0; i--)
        heapify(v, dimensiune, i);
    for (i = dimensiune - 1; i >= 0; i--) {
        swap(v[0], v[i]);
        heapify(v, i, 0);
    }
}
int main() {
    srand(time(0));
    int T;
    f >> T;
    long long n, m;
    vector<long long> v, v1, v2, v3, v4, v5, v6, v7;
    for (int j = 0; j < T; j++) {
        f >> n >> m;
        for (long long i = 0; i < n; i++) {
            long long x = rand();
            x = x*x*x;
            x = x%m;
            v.push_back(x);
            v1.push_back(x);
            v2.push_back(x);
            v3.push_back(x);
            v4.push_back(x);
            v5.push_back(x);
            v6.push_back(x);
            v7.push_back(x);
        }
        cout << "N = " << n << " M = " << m << endl;
        cout << "Stl Sort: ";
        clock_t tStart = clock();
        sort(v.begin(), v.end());
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << endl;

        cout << "Shellsort: ";
        tStart = clock();
        shellSort(v1, v.size());
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v1, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << "MergeSort: ";
        tStart = clock();
        mergeSort(v2);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v2, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << "Radix Sort(baza 2): ";
        tStart = clock();
        radixsort2(v3, n);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v3, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;


        cout << "Radix Sort(baza 10): ";
        tStart = clock();
        radixsort(v4, n);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v4, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << "Heap Sort: ";
        tStart = clock();
        heapSort(v6, n);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v6, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << "Radix Sort(baza 2^16): ";
        tStart = clock();
        radixsort3(v7, n);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v7, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << "Counting sort: ";
        tStart = clock();
        countingSort(v5, n);
        cout.precision(4);
        cout << fixed << float(clock() - tStart) / CLOCKS_PER_SEC << " | ";
        if (e_sortat(v, v5, n)) cout << "Sortat cu succes!" << endl;
        else cout << "Sortarea a esuat!" << endl;

        cout << endl;
        v.clear();
        v1.clear();
        v2.clear();
        v3.clear();
        v4.clear();
        v5.clear();
        v6.clear();
        v7.clear();
        }

}

