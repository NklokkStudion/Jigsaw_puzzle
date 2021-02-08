#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setFixedSize(1044, 512);
    label_start = new QLabel(this);
    label_finish = new QLabel(this);
    label_start->setGeometry(0, 0, 512, 512);
    label_finish->setGeometry(532, 0, 512, 512);
}

void MainWindow::set_image(const QString& name) {
    image_start =  QImage(name);
    image_start.scaled(512, 512);
    label_start->setPixmap(QPixmap::fromImage(image_start, Qt::AutoColor));
    label_start->show();
}

void MainWindow::solve(bool flag) {
    Solver algorithm(image_start);
    std::pair<QImage, double> getted = algorithm.get_image(flag);
    image_finish = getted.first;
    scores.push_back(getted.second);
    label_finish->setPixmap(QPixmap::fromImage(image_finish, Qt::AutoColor));
    label_finish->show();
}

MainWindow::~MainWindow()
{
    delete ui;
}

