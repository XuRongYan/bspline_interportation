//
// Created by 徐溶延 on 2020/2/3.
//

#include "MainWindow.h"
#include <QEvent>
#include <QPaintEvent>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QDebug>
#include <QPainter>
#include <QPen>
#include <Eigen/Dense>
#include <iostream>
#include "ui_MainWindow.h"
#include "bspline_interportation.h"
#include "BSpline.h"

MainWindow::MainWindow(QWidget *parent) :
        QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    setBackGroundColors();
    ui->widget_panel->installEventFilter(this);
}


MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::init() {

}

void MainWindow::setBackGroundColors() {
    QPalette panelPalette = ui->widget_panel->palette();
    panelPalette.setColor(QPalette::Active, QPalette::Window, Qt::white);
    ui->widget_panel->setPalette(panelPalette);
}

void MainWindow::handleMousePressEvent(QMouseEvent *e) {

}

void MainWindow::handleMouseMoveEvent(QMouseEvent *e) {

}

void MainWindow::handleMouseReleaseEvent(QMouseEvent *e) {
    QPointF pointF(e->x(), e->y());
    qDebug() << "add point" << pointF;
    QPoints.push_back(pointF);
    vecV.push_back(pointF.x());
    vecV.push_back(pointF.y());
    ui->widget_panel->update();
}

void MainWindow::draw(QPaintEvent *e) {
    drawValuePoints(e);
    //drawValueLines(e);
    //drawControlLines(e);
    //drawControlPoints(e);
    drawBSpline(e);
    //drawKnotPoints(e);
}

void MainWindow::drawValuePoints(QPaintEvent *e) {
    QPainter controlPointPainter(ui->widget_panel);
    QPen controlPointsPen(Qt::red);
    controlPointsPen.setCapStyle(Qt::RoundCap);
    controlPointsPen.setWidth(5);
    controlPointPainter.setPen(controlPointsPen);
    controlPointPainter.drawPoints(QPoints.data(), QPoints.size());
}

void MainWindow::drawValueLines(QPaintEvent *e) {
    QPainter controlLinePainter(ui->widget_panel);
    QPen controlLinePen(Qt::red);
    controlLinePen.setStyle(Qt::DotLine);
    controlLinePainter.setPen(controlLinePen);
    for (int i = 0; i < (int)QPoints.size() - 1; i++) {
        controlLinePainter.drawLine(QPoints[i], QPoints[i + 1]);
    }
}

void MainWindow::drawControlPoints(QPaintEvent *e) {
    qDebug() << "drawControlPoints";
    QPainter controlPointPainter(ui->widget_panel);
    QPen controlPointsPen(Qt::blue);
    controlPointsPen.setWidth(10);
    controlPointsPen.setCapStyle(Qt::SquareCap);
    controlPointPainter.setPen(controlPointsPen);
    std::vector<QPointF> controlPoints = bSpline.getControlPoints();
    controlPointPainter.drawPoints(controlPoints.data(), controlPoints.size());
}

void MainWindow::drawControlLines(QPaintEvent *e) {
    QPainter controlLinePainter(ui->widget_panel);
    QPen controlLinePen(Qt::red);
    controlLinePen.setStyle(Qt::DotLine);
    controlLinePainter.setPen(controlLinePen);
    std::vector<QPointF> controlPoints = bSpline.getControlPoints();
    for (int i = 0; i < (int)controlPoints.size() - 1; i++) {
        controlLinePainter.drawLine(controlPoints[i], controlPoints[i + 1]);
    }
}

void MainWindow::drawBSpline(QPaintEvent *e) {
    QPainter curvePainter(ui->widget_panel);
    curvePainter.setPen(QPen(Qt::black));
    std::vector<QPointF> points = bSpline.getPoints();
    int first = 0, last = (int)points.size() - 1;
    for (int i = 0; i < (int)points.size() - 1; i++) {
        QPointF zero(0, 0);
        if (points[i] == zero) {
            first = i + 1;
            continue;
        } else if (points[i + 1] == zero) {
            last = i;
            break;
        }
        curvePainter.drawLine(points[i], points[i + 1]);
    }
    if (!points.empty()) curvePainter.drawLine(points[last], points[first]);
}

void MainWindow::checkMeetPoint(QMouseEvent *e) {

}

void MainWindow::keyReleaseEvent(QKeyEvent *event) {
    if (event->key() == Qt::Key_Space ) {
        qDebug() << "enter";
        QPoints.push_back(QPoints[0]);
        vecV.push_back(vecV[0]);
        vecV.push_back(vecV[1]);
        Eigen::Matrix2Xd V = Eigen::Map<Eigen::Matrix2Xd>(vecV.data(), 2, vecV.size() / 2);
        Eigen::Matrix2Xd Vout;
        std::vector<double> knots;
        xry_mesh::interporlate_bspline1(V, Vout, knots);
        std::vector<QPointF> controlPoints;
        std::cout << "Vout:\n" << Vout << std::endl;
        controlPoints.reserve(Vout.cols());
        for (int i = 0; i < Vout.cols(); i++) {
            controlPoints.emplace_back(Vout(0, i), Vout(1, i));
        }
        qDebug() << controlPoints;
        bSpline = BSpline(3, controlPoints, knots);
        bSpline.calc();
        ui->widget_panel->update();
    } else {
        QWidget::keyReleaseEvent(event);
    }
}

bool MainWindow::eventFilter(QObject *watched, QEvent *event) {
    //监听
    if (watched == ui->widget_panel) {
        if (event->type() == QEvent::Paint) {
            qDebug() << "paint";
            draw(dynamic_cast<QPaintEvent *>(event));
            return true;
        }
        if (event->type() == QEvent::MouseButtonPress) {
            handleMousePressEvent(dynamic_cast<QMouseEvent *>(event));
            return true;
        }

        if (event->type() == QEvent::MouseMove) {
            handleMouseMoveEvent(dynamic_cast<QMouseEvent *>(event));
            return true;
        }

        if (event->type() == QEvent::MouseButtonRelease) {
            handleMouseReleaseEvent(dynamic_cast<QMouseEvent *>(event));
            return true;
        }

    }
    return QObject::eventFilter(watched, event);
}

void MainWindow::drawKnotPoints(QPaintEvent *e) {
    qDebug() << "drawControlPoints";
    QPainter controlPointPainter(ui->widget_panel);
    QPen controlPointsPen(Qt::green);
    controlPointsPen.setWidth(10);
    controlPointPainter.setPen(controlPointsPen);
    std::vector<QPointF> controlPoints = bSpline.getKnotPoints();
    controlPointPainter.drawPoints(controlPoints.data(), controlPoints.size());
}


