#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <QDebug>
#include <QImage>
#include <vector>
#include <math.h>
#include <random>
#include <algorithm>
#include <future>
#include <tuple>
#include <set>

using std::vector;

const int WH = 512;
const int SIZE = 64;
const int CNTSQ = WH / SIZE;
const int CORES = 8;
const int WAVES = 3;
const int GENS = 50;
const int POP = 60;
const int LIVE = 20;
const int ODDS = 5;
const int MRATE = 5;

class Solver
{
private:
    struct square {
        vector<vector<vector<int>>> rgb;
        square() {
            rgb.resize(3, vector<vector<int>>(SIZE, vector<int>(SIZE)));
            for (int k = 0; k < 3; ++k)
            for (int i = 0; i < SIZE; ++i)
            for (int j = 0; j < SIZE; ++j)
                rgb[k][i][j] = 0;
        }
    };
    QImage image;
    vector<vector<square>> PICTURE;
    vector<vector<vector<int>>> DISS;
public:
    Solver(QImage image);
    vector<int> solve();
    void selection(vector<vector<int>>& population);
    void breeding(vector<vector<int>>& population);
    double get_diss(const vector<int>& v);
    std::pair<QImage, double> get_image(bool flag = 0);
};

#endif // SOLVER_H
