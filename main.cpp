#include "mainwindow.h"

#include <QApplication>
#include <QLabel>
#include <QImage>
#include <QDebug>

struct Timer
{
    QString name; clock_t tm;
    Timer(QString name) : name(name), tm(clock()){}
    ~Timer(){qDebug() << name << " : " <<(double)(clock()-tm)/CLOCKS_PER_SEC << "ms";}
};

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
//    w.set_image("img_64//1200.png");
//    w.solve(0);
    int cnt = 10;
    {
        Timer t("Time");
        for (int i = 1200; i < 1200 + cnt; ++i) {
            qDebug() << "------" << i << "-----";
            w.set_image("img_64//" + QString::number(i) + ".png");
            w.solve(1);
            w.set_image("img_64//" + QString::number(i) + ".png");
            w.solve();
            w.set_image("img_64_ready//" + QString::number(i) + ".png");
            w.solve(1);
            qDebug() << "-----------------";
        }
    }
    double avg_percent = 0;
    for (int i = 0; i < (int)w.scores.size(); i += 3) {
        double percent = (w.scores[i] - w.scores[i + 1]) / (w.scores[i] - w.scores[i + 2]);
        avg_percent += percent;
    }
    avg_percent /= cnt;
    avg_percent *= 100;
    qDebug() << "AVG : " << avg_percent << "%";
    w.show();
    return a.exec();
}
