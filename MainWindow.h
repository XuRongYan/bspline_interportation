//
// Created by 徐溶延 on 2020/2/3.
//

#ifndef BSPLINE_INTERPORTATION_MAINWINDOW_H
#define BSPLINE_INTERPORTATION_MAINWINDOW_H
#include <QMainWindow>
#include "BSpline.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = 0);

    ~MainWindow() override;

    bool eventFilter(QObject *watched, QEvent *event) override;

private:
    Ui::MainWindow *ui;
    void init();
    void setBackGroundColors();
    void handleMousePressEvent(QMouseEvent *e);
    void handleMouseMoveEvent(QMouseEvent *e);
    void handleMouseReleaseEvent(QMouseEvent *e);
    void draw(QPaintEvent *e);
    void drawValuePoints(QPaintEvent *e);
    void drawValueLines(QPaintEvent *e);
    void drawControlPoints(QPaintEvent *e);
    void drawControlLines(QPaintEvent *e);
    void drawBSpline(QPaintEvent *e);
    void drawKnotPoints(QPaintEvent *e);
    void checkMeetPoint(QMouseEvent *e);
protected:
    std::vector<QPointF> QPoints;
    std::vector<double> vecV;
    BSpline bSpline;
    void keyReleaseEvent(QKeyEvent *event) override;
};


#endif //BSPLINE_INTERPORTATION_MAINWINDOW_H
