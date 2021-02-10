#include "solver.h"
#define cbr(a) (a) * (a) * (a)

Solver::Solver(QImage image) : image(image)
{
    std::srand(std::time(0));
    DISS.resize(4, vector<vector<int>>(CNTSQ * CNTSQ, vector<int>(CNTSQ * CNTSQ)));
    PICTURE.resize(CNTSQ, vector<square>(CNTSQ));
    for (int i = 0; i < CNTSQ; ++i)
    for (int j = 0; j < CNTSQ; ++j)
    for (int x = 0; x < SIZE; ++x)
    for (int y = 0; y < SIZE; ++y) {
        QColor color = image.pixelColor(i * SIZE + x, j * SIZE + y);
        PICTURE[i][j].rgb[0][x][y] = color.red();
        PICTURE[i][j].rgb[1][x][y] = color.green();
        PICTURE[i][j].rgb[2][x][y] = color.blue();
    }
    vector<int> help = {};
    for (int k1 = 0; k1 < CNTSQ * CNTSQ; ++k1) {
        for (int k2 = 0; k2 < CNTSQ * CNTSQ; ++k2) {
            if (k1 == k2) continue;
            int i1 = k1 / CNTSQ; int j1 = k1 % CNTSQ;
            int i2 = k2 / CNTSQ; int j2 = k2 % CNTSQ;
            for (int dir = 0 ; dir < 4; ++dir) {
                double mb = 0;
                for (int x = 0; x < SIZE; ++x) {
                    for (int c = 0; c < 3; ++c) {
                        if (dir == 0) {
                               mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][x][0] - PICTURE[i2][j2].rgb[c][x][SIZE - 1])));
                        } else if (dir == 1) {
                                mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][0][x] - PICTURE[i2][j2].rgb[c][SIZE - 1][x])));
                        } else if (dir == 2) {
                                mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][x][SIZE - 1] - PICTURE[i2][j2].rgb[c][x][0])));
                        } else if (dir == 3) {
                                mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][SIZE - 1][x] - PICTURE[i2][j2].rgb[c][0][x])));
                        }
                    }
                }
                DISS[dir][k1][k2] = std::sqrt(mb);
            }
        }
    }
}

double Solver::get_diss(const vector<int>& v) {
    double res = 0;
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
        if (i % CNTSQ != CNTSQ - 1) res += DISS[2][v[i]][v[i + 1]];
        if (i / CNTSQ != CNTSQ - 1) res += DISS[3][v[i]][v[i + CNTSQ]];
    }
    return res;
}

void vertical(vector<int>& v, bool flag) {
    vector<int> temp(CNTSQ * CNTSQ);
    for (int i = 0; i < CNTSQ * CNTSQ; i += CNTSQ) {
        for (int j = 0; j < CNTSQ - 1; ++j) {
            temp[i + j + flag] = v[i + j + (flag ^ 1)];
        }
        temp[i + ((flag ^ 1) ? CNTSQ - 1 : 0)] = v[i + (flag ? CNTSQ - 1 : 0)];
    }
    v = temp;
}

void horizontal(vector<int>& v, bool flag) {
    vector<int> temp(CNTSQ * CNTSQ);
    for (int i = 0; i < CNTSQ; ++i) {
        for (int j = 0; j < CNTSQ - 1; ++j) {
            temp[i + (j + (flag ^ 1)) * CNTSQ] = v[i + (j + flag) * CNTSQ];
        }
        temp[i + (flag ? (CNTSQ - 1) * CNTSQ : 0)] = v[i + ((flag ^ 1) ? (CNTSQ - 1) * CNTSQ : 0)];
    }
    v = temp;
}

void mutate(vector<int>& v) {
    std::srand(std::time(0));
    int r = rand() % 100;
    if (r < 50) {
        if (r % 2 == 0) vertical(v, 1);
        else vertical(v, 0);
    } else {
        if (r % 2 == 0) horizontal(v, 1);
        else horizontal(v, 0);
    }
}

void Solver::selection(vector<vector<int>>& population) {
    vector<std::pair<double, int>> best;
    for (int i = 0; i < POP; ++i) {
        best.push_back({get_diss(population[i]), i});
    }
    std::sort(best.begin(), best.end());
    for (int i = 0; i < LIVE; ++i) {
        std::swap(population[i], population[best[i].second]);
    }
    for (int i = LIVE; i < LIVE + ODDS; ++i) {
        int idx = rand() % (POP - i) + (i - 1);
        std::swap(population[i], population[idx]);
    }
    for (int i = 0; i < LIVE + ODDS; ++i) {
        if (rand() % 100 < MRATE) {
            mutate(population[i]);
        }
    }
}

void greedy(vector<int>& can, vector<int>& ready, const vector<vector<vector<int>>>& DISS) {
    if (can.size() == 0) return;
    std::tuple<double, int, int> ans = {1e15, can[0], -1};
    vector<int> help = {-1, -CNTSQ, 1, CNTSQ};
    for (int what : can) {
        for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
            if (ready[i] != -1) continue;
            if (std::get<2>(ans) == -1) ans = {1e15, what, i};
            double mb = 1e15;
            for (int dir = 0; dir < 4; ++dir) {
                if ((i % CNTSQ == 0 && dir == 0) || (i % CNTSQ == CNTSQ - 1 && dir == 2) ||
                    (i / CNTSQ == 0 && dir == 1) || (i / CNTSQ == CNTSQ - 1 && dir == 3)) continue;
                if (ready[i + help[dir]] != -1) {
                    if (mb == 1e15) mb = 0;
                    mb += DISS[dir][what][ready[i + help[dir]]];
                }
            }
            if (mb < std::get<0>(ans)) {
                ans = {mb, what, i};
            }
        }
    }
    can.erase(std::find(can.begin(), can.end(), std::get<1>(ans)));
    ready[std::get<2>(ans)] = std::get<1>(ans);
}

void divider(vector<vector<int>>& population, int l, int r, const vector<vector<vector<int>>>& DISS) {
    std::srand(std::time(0));
    for (int i = l; i <= r; ++i) {
        int idx_a = rand() % (LIVE + ODDS);
        int idx_b = rand() % (LIVE + ODDS - 1);
        if (idx_b >= idx_a) ++idx_b;
        vector<int> can;
        for (int i = 0; i < CNTSQ * CNTSQ; ++i) can.push_back(i);
        std::random_shuffle(can.begin(), can.end());
        vector<int> ready(CNTSQ * CNTSQ, -1);
        while (can.size()) {
            greedy(can, ready, DISS);
        }
        population[i] = ready;
    }
}

void Solver::breeding(vector<vector<int>>& population) {
    vector<std::future<void>> futures;
    int block = (POP - (LIVE + ODDS)) / CORES;
    int ost = (POP - (LIVE + ODDS)) % CORES;
    vector<int> blocks(CORES, block);
    for (int i = 0; i < CORES; ++i) {
        if (ost == 0) break;
        ++blocks[i]; --ost;
    }
    int l = LIVE + ODDS;
    for (int i = 0; i < CORES; ++i) {
        if (blocks[i] == 0) break;
        int r = l + blocks[i] - 1;
        futures.push_back(std::async(std::launch::async, divider, ref(population), l, r, ref(DISS)));
        l += blocks[i];
    }
    for (auto& x : futures) x.get();
}

vector<int> Solver::solve() {
    vector<int> answer;
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) answer.push_back(i);
    for (int i = 0; i < WAVES; ++i) {
        vector<vector<int>> population;
        for (int i = 0; i < POP; ++i) {
            vector<int> temp = answer;
            std::random_shuffle(temp.begin(), temp.end());
            population.push_back(answer);
        }
        for (int j = 0; j < GENS; ++j) {
            selection(population);
            breeding(population);
        }
        if (get_diss(population[0]) < get_diss(answer)) {
            answer = population[0];
        }
    }
    return answer;
}

std::pair<QImage, double> Solver::get_image(bool flag) {
    vector<int> answer;
    if (flag == 0) answer = solve();
    else for (int i = 0; i < CNTSQ * CNTSQ; ++i) answer.push_back(i);
    double score = get_diss(answer);
    qDebug() << "answer : "<< score;
    QImage result(WH, WH, QImage::Format_RGB32);
    for (int k = 0; k < CNTSQ * CNTSQ; ++k) {
        int i1 = k / CNTSQ; int j1 = k % CNTSQ;
        int i = answer[k] / CNTSQ; int j = answer[k] % CNTSQ;
        for (int x = 0; x < SIZE; ++x) {
            for (int y = 0; y < SIZE; ++y) {
                QColor color(PICTURE[i][j].rgb[0][x][y], PICTURE[i][j].rgb[1][x][y], PICTURE[i][j].rgb[2][x][y]);
                result.setPixelColor(i1 * SIZE + x, j1 * SIZE + y, color);
            }
        }
    }
    return {result, score};
}
