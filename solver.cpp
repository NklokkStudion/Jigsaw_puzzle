#include "solver.h"
#define cbr(a) (a) * (a) * (a)

Solver::Solver(QImage image) : image(image)
{
    std::srand(std::time(0));
    DISS.resize(4, vector<vector<int>>(CNTSQ * CNTSQ, vector<int>(CNTSQ * CNTSQ)));
    BDISS.resize(4, vector<double>(CNTSQ * CNTSQ, -1));
    BBODY.resize(4, vector<int>(CNTSQ * CNTSQ, -1));
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
    for (int k1 = 0; k1 < CNTSQ * CNTSQ; ++k1) {
        for (int k2 = 0; k2 < CNTSQ * CNTSQ; ++k2) {
            if (k1 == k2) continue;
            int i1 = k1 / CNTSQ;
            int j1 = k1 % CNTSQ;
            int i2 = k2 / CNTSQ;
            int j2 = k2 % CNTSQ;
            for (int dir = 0 ; dir < 4; ++dir) {
                double mb = 0;
                for (int x = 0; x < SIZE; ++x) {
                    if (dir == 0) {
                       for (int c = 0; c < 3; ++c) {
                           mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][x][0] - PICTURE[i2][j2].rgb[c][x][SIZE - 1])));
                       }
                    } else if (dir == 1) {
                        for (int c = 0; c < 3; ++c) {
                            mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][0][x] - PICTURE[i2][j2].rgb[c][SIZE - 1][x])));
                        }
                    } else if (dir == 2) {
                        for (int c = 0; c < 3; ++c) {
                            mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][x][SIZE - 1] - PICTURE[i2][j2].rgb[c][x][0])));
                        }
                    } else if (dir == 3) {
                        for (int c = 0; c < 3; ++c) {
                            mb += (double)(cbr(abs(PICTURE[i1][j1].rgb[c][SIZE - 1][x] - PICTURE[i2][j2].rgb[c][0][x])));
                        }
                    }
                }
                mb = std::sqrt(mb);
                DISS[dir][k1][k2] = mb;
                if (BDISS[dir][k1] == -1 || mb < BDISS[dir][k1]) {
                    BDISS[dir][k1] = mb;
                    BBODY[dir][k1] = k2;
                }
            }
        }
    }
}

double Solver::get_diss(const vector<int>& v) {
    double res = 0;
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
        if (i % CNTSQ != CNTSQ - 1) {
            res += DISS[2][v[i]][v[i + 1]];
        }
        if (i / CNTSQ != CNTSQ - 1) {
            res += DISS[3][v[i]][v[i + CNTSQ]];
        }
    }
    return res;
}

void mutate(vector<int>& v) {
    std::srand(std::time(0));
    int r = rand() % 100;
    if (r < 50) {
        if (r % 2 == 0) { // DOWN
            vector<int> temp(CNTSQ * CNTSQ);
            for (int i = 0; i < CNTSQ * CNTSQ; i += CNTSQ) {
                for (int j = 0; j < CNTSQ - 1; ++j) {
                    temp[i + 1 + j] = v[i + j];
                }
                temp[i] = v[i + CNTSQ - 1];
            }
            v = temp;
        } else { // UP
            vector<int> temp(CNTSQ * CNTSQ);
            for (int i = 0; i < CNTSQ * CNTSQ; i += CNTSQ) {
                for (int j = 0; j < CNTSQ - 1; ++j) {
                    temp[i + j] = v[i + j + 1];
                }
                temp[i + CNTSQ - 1] = v[i];
            }
            v = temp;
        }
    } else {
        if (r % 2 == 0) { // LEFT
            vector<int> temp(CNTSQ * CNTSQ);
            for (int i = 0; i < CNTSQ; ++i) {
                for (int j = 0; j < CNTSQ - 1; ++j) {
                    temp[i + j * CNTSQ] = v[i + (j + 1) * CNTSQ];
                }
                temp[i + (CNTSQ - 1) * CNTSQ] = v[i];
            }
            v = temp;
        } else { // RIGHT
            vector<int> temp(CNTSQ * CNTSQ);
            for (int i = 0; i < CNTSQ; ++i) {
                for (int j = 0; j < CNTSQ - 1; ++j) {
                    temp[i + (j + 1) * CNTSQ] = v[i + j * CNTSQ];
                }
                temp[i] = v[i + (CNTSQ - 1) * CNTSQ];
            }
            v = temp;
        }
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

void phase_1(const vector<int>& a, const vector<int>& b, vector<int>& can, vector<int>& ready) {
    std::srand(std::time(0));
    std::set<int> have;
    for (int x : can) have.insert(x);
    vector<std::pair<int, int>> best;
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
        if (a[i] == b[i] && have.count(a[i]) == 1 && ready[i] == -1) {
            best.push_back({a[i], i});
        }
    }
    if (best.size() == 0) return;
    int idx = rand() % (int)(best.size());
    can.erase(std::find(can.begin(), can.end(), best[idx].first));
    ready[best[idx].second] = best[idx].first;
}

void phase_2(const vector<int>& a, const vector<int>& b, vector<int>& can, vector<int>& ready, const vector<vector<int>>& BBODY) {
    std::srand(std::time(0));

    std::set<int> have;
    for (int x : can) have.insert(x);

    vector<std::tuple<int, int, int, int>> best;
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
        for (int dir = 0; dir < 2; ++dir) {
            if (ready[i] != -1) continue;
            if (i % CNTSQ == 0 && dir == 0) continue;
            if (i / CNTSQ == 0 && dir == 1) continue;
            if (dir == 0) {
                if (ready[i - 1] != -1) continue;
                if (have.count(a[i]) == 0) continue;
                if (have.count(a[i - 1]) == 0) continue;
                if (BBODY[0][a[i]] == a[i - 1] && BBODY[2][a[i - 1]] == a[i]) best.push_back({a[i - 1], a[i], i - 1, i});
            } else if (dir == 1) {
                if (ready[i - CNTSQ] != -1) continue;
                if (have.count(a[i]) == 0) continue;
                if (have.count(a[i - CNTSQ]) == 0) continue;
                if (BBODY[1][a[i]] == a[i - CNTSQ] && BBODY[3][a[i - CNTSQ]] == a[i]) best.push_back({a[i - CNTSQ], a[i], i - CNTSQ, i});
            }
        }
    }
    for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
        for (int dir = 0; dir < 2; ++dir) {
            if (ready[i] != -1) continue;
            if (i % CNTSQ == 0 && dir == 0) continue;
            if (i / CNTSQ == 0 && dir == 1) continue;
            if (dir == 0) {
                if (ready[i - 1] != -1) continue;
                if (have.count(b[i]) == 0) continue;
                if (have.count(b[i - 1]) == 0) continue;
                if (BBODY[0][b[i]] == b[i - 1] && BBODY[2][b[i - 1]] == b[i]) best.push_back({b[i - 1], b[i], i - 1, i});
            } else if (dir == 1) {
                if (ready[i - CNTSQ] != -1) continue;
                if (have.count(b[i]) == 0) continue;
                if (have.count(b[i - CNTSQ]) == 0) continue;
                if (BBODY[1][b[i]] == b[i - CNTSQ] && BBODY[3][b[i - CNTSQ]] == b[i]) best.push_back({b[i - CNTSQ], b[i], i - CNTSQ, i});
            }
        }
    }
    if (best.size() == 0) return;
    int idx = rand() % (int)(best.size());
    can.erase(std::find(can.begin(), can.end(), std::get<0>(best[idx])));
    can.erase(std::find(can.begin(), can.end(), std::get<1>(best[idx])));
    ready[std::get<2>(best[idx])] = std::get<0>(best[idx]);
    ready[std::get<3>(best[idx])] = std::get<1>(best[idx]);
}

void phase_3(vector<int>& can, vector<int>& ready, const vector<vector<vector<int>>>& DISS) {
    if (can.size() == 0) return;
    int bwhat = can[0];
    int bidx = -1;
    double bvalue = 1e15;
    for (int what : can) {
        for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
            if (ready[i] != -1) continue;
            if (bidx == -1) {
                bidx = i;
                bwhat = what;
                bvalue = 1e15;
            }
            double mb = 1e15;
            for (int dir = 0; dir < 4; ++dir) {
                if (i % CNTSQ == 0 && dir == 0) continue;
                if (i / CNTSQ == 0 && dir == 1) continue;
                if (i % CNTSQ == CNTSQ - 1 && dir == 2) continue;
                if (i / CNTSQ == CNTSQ - 1 && dir == 3) continue;
                if (dir == 0 && ready[i - 1] != -1) {
                    if (mb == 1e15) mb = 0;
                    mb += DISS[dir][what][ready[i - 1]];
                }
                if (dir == 1 && ready[i - CNTSQ] != -1) {
                    if (mb == 1e15) mb = 0;
                    mb += DISS[dir][what][ready[i - CNTSQ]];
                }
                if (dir == 2 && ready[i + 1] != -1) {
                    if (mb == 1e15) mb = 0;
                    mb += DISS[dir][what][ready[i + 1]];
                }
                if (dir == 3 && ready[i + CNTSQ] != -1) {
                    if (mb == 1e15) mb = 0;
                    mb += DISS[dir][what][ready[i + CNTSQ]];
                }
            }
            if (mb < bvalue) {
                bvalue = mb;
                bidx = i;
                bwhat = what;
            }
        }
    }
    can.erase(std::find(can.begin(), can.end(), bwhat));
    ready[bidx] = bwhat;
}

void divider(vector<vector<int>>& population, int l, int r, int nogen, const vector<vector<vector<int>>>& DISS, const vector<vector<int>>& BBODY) {
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
            //if (nogen >= FPHASE) phase_1(population[idx_a], population[idx_b], can, ready);
            //phase_2(population[idx_a], population[idx_b], can, ready, BBODY);
            phase_3(can, ready, DISS);
        }
        population[i] = ready;
    }
}

void Solver::breeding(vector<vector<int>>& population, int nogen) {
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
        futures.push_back(std::async(std::launch::async, divider, ref(population), l, r, nogen, ref(DISS), ref(BBODY)));
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
            breeding(population, j);
        }
        if (get_diss(population[0]) < get_diss(answer)) {
            answer = population[0];
        }
    }
    return answer;
}

std::pair<QImage, double> Solver::get_image(bool flag) {
    vector<int> answer;
    if (flag == 0) {
        answer = solve();
    } else {
        for (int i = 0; i < CNTSQ * CNTSQ; ++i) {
            answer.push_back(i);
        }
    }
    double score = get_diss(answer);
    qDebug() << "answer : "<< score;
    QImage result(WH, WH, QImage::Format_RGB32);
    for (int k = 0; k < CNTSQ * CNTSQ; ++k) {
        int i1 = k / CNTSQ;
        int j1 = k % CNTSQ;
        int i = answer[k] / CNTSQ;
        int j = answer[k] % CNTSQ;
        for (int x = 0; x < SIZE; ++x) {
            for (int y = 0; y < SIZE; ++y) {
                QColor color(PICTURE[i][j].rgb[0][x][y], PICTURE[i][j].rgb[1][x][y], PICTURE[i][j].rgb[2][x][y]);
                result.setPixelColor(i1 * SIZE + x, j1 * SIZE + y, color);
            }
        }
    }
    return {result, score};
}
