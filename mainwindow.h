#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "solver.h"

#include <QMainWindow>
#include <QLabel>
#include <QImage>
#include <QDebug>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void set_image(const QString& name);
    void solve(bool flag = 0);
    std::vector<double> scores;
private:
    Ui::MainWindow *ui;
    QImage image_start;
    QLabel *label_start;
    QImage image_finish;
    QLabel *label_finish;
};
#endif // MAINWINDOW_H
