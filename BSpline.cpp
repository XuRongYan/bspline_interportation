//
// Created by 徐溶延 on 2020/1/31.
//

#include "BSpline.h"
#include <cstdlib>
#include <QDebug>
#include <utility>
#include <cmath>

int BSpline::getOrder() const {
    return order;
}

void BSpline::setOrder(int order) {
    BSpline::order = order;
}

int BSpline::getControlNum() const {
    return controlNum;
}

void BSpline::setControlNum(int controlNum) {
    BSpline::controlNum = controlNum;
}

int BSpline::getKnotNum() const {
    return knotNum;
}

void BSpline::setKnotNum(int knotNum) {
    BSpline::knotNum = knotNum;
}

BSpline::KnotType BSpline::getType() const {
    return type;
}

void BSpline::setType(BSpline::KnotType type) {
    BSpline::type = type;
}

const std::vector<QPointF> &BSpline::getControlPoints() const {
    return controlPoints;
}

void BSpline::setControlPoints(const std::vector<QPointF> &controlPoints) {
    BSpline::controlPoints = controlPoints;
}

const std::vector<double> &BSpline::getKnots() const {
    return knots;
}

void BSpline::setKnots(const std::vector<double> &knots) {
    BSpline::knots = knots;
}

BSpline::BSpline(int order, BSpline::KnotType type, std::vector<QPointF> controlPoints) : order(order),
                                                                                          type(type),
                                                                                          controlPoints(std::move(
                                                                                                  controlPoints)),
                                                                                          pointNum(100) {
    init();
}

BSpline::BSpline(int order, const std::vector<QPointF> &controlPoints, const std::vector<double> &knots) : order(order),
                                                                                                           controlPoints(
                                                                                                                   controlPoints),
                                                                                                           knots(knots) {
    P.clear();
    ts.clear();
    pointNum = 500;
    controlNum = controlPoints.size();
    knotNum = knots.size() + 2;
    this->knots.insert(this->knots.begin(), 0);
    this->knots.push_back(1);
    double step = 1.0 / pointNum;
    for (int i = 0; i < pointNum; i++) {
        ts.push_back(step * i);
    }
    ts.push_back(1);
    int idx = 1;
    for (auto it = ts.begin(); it != ts.end(); it++) {
        if (fabs(*it - knots[idx]) < step && fabs(*it - knots[idx]) > 1e-12 && idx < knots.size() - 1) {
            it = ts.insert(it, knots[idx++]);
        } else if (knots[idx] < 1e-12) {
            idx++;
        } else if (knots[idx] == 1) {
            break;
        }
    }
    points = std::vector<QPointF>(ts.size(), QPointF(0, 0));

}

double BSpline::cox_de_Boor(int i, int p, double t) {
    if (p == 0) {
        if (t < knots[i] || t >= knots[i + 1]) return 0;
        return 1;
    } else {
        if (knots[i + p] - knots[i] <= 0 && knots[i + p + 1] - knots[i + 1] > 0) {
            return ((knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1])) * cox_de_Boor(i + 1, p - 1, t);
        } else if (knots[i + p] - knots[i] > 0 && knots[i + p + 1] - knots[i + 1] <= 0) {
            return ((t - knots[i]) / (knots[i + p] - knots[i])) * cox_de_Boor(i, p - 1, t);
        } else if (knots[i + p] - knots[i] <= 0 && knots[i + p + 1] - knots[i + 1] <= 0) {
            return 0;
        }
        return ((t - knots[i]) / (knots[i + p] - knots[i])) * cox_de_Boor(i, p - 1, t) +
               ((knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1])) * cox_de_Boor(i + 1, p - 1, t);
    }
}

int BSpline::calc() {
    points.clear();
    knotPoints.clear();
    for (double t : ts) {
        QPointF qPointF(0, 0);
        //if (t >= knots[0] && t <= knots[knots.size() - 1] && knotNum == knots.size()) {
        if (!isClosed) {
            if (t >= knots[0] && t <= knots[knots.size() - 1] && knotNum == knots.size()) {
                for (int i = 0; i < controlPoints.size(); i++) {
                    qPointF += cox_de_Boor(i, order, t) * controlPoints[i];
                }
            }
        } else {
            if (t >= knots[0] && t <= knots[knots.size() - 1]) {
                for (int i = 0; i < controlPoints.size(); i++) {
                    qPointF += cox_de_Boor(i, order, t) * controlPoints[i];
                }
            }
        }

        qDebug() << "t=" << t << qPointF << (t >= knots[order] && t <= knots[knots.size() - order - 1] ? "show" : "hide");
        points.push_back(qPointF);
    }
    //for (double t : knots) qDebug() << t;
    if (type == SEQ || isClosed) {
        for (int i = order; i < knots.size() - order; i++) {
            QPointF qPointF(0, 0);
            for (int j = 0; j < controlPoints.size(); j++) {
                qPointF += cox_de_Boor(j, order, knots[i]) * controlPoints[j];
            }
            knotPoints.push_back(qPointF);
        }
    } else {

        for (int i = 0; i < knots.size(); i++) {
            QPointF qPointF(0, 0);
            if (knots[i] == 1) knots[i] -= 0.000001;
            for (int j = 0; j < controlPoints.size(); j++) {
                qPointF += cox_de_Boor(j, order, knots[i]) * controlPoints[j];
            }
            knotPoints.push_back(qPointF);
        }
    }
    printInfo();
    return 0;
}

BSpline::BSpline() {}

void BSpline::deleteAll() {
    controlPoints.clear();
    notifyUpdateAll();
}

void BSpline::deleteLastPoint() {
    if (controlPoints.empty())
        return;
    controlPoints.pop_back();
    notifyUpdateAll();
}

void BSpline::addPoint(const QPointF &qPointF) {
    controlPoints.push_back(qPointF);
    notifyAddPoint();
}

int BSpline::getPointNum() const {
    return pointNum;
}

void BSpline::setPointNum(int pointNum) {
    BSpline::pointNum = pointNum;
}

const std::vector<QPointF> &BSpline::getPoints() const {
    return points;
}

void BSpline::setPoints(const std::vector<QPointF> &points) {
    BSpline::points = points;
}

const std::vector<int> &BSpline::getP() const {
    return P;
}

void BSpline::setP(const std::vector<int> &p) {
    P = p;
}

const std::vector<double> &BSpline::getTs() const {
    return ts;
}

void BSpline::setTs(const std::vector<double> &ts) {
    BSpline::ts = ts;
}

void BSpline::setControlPoint(int i, const QPointF &pointF) {
    controlPoints[i] = pointF;
    notifyUpdateAll();
}

void BSpline::notifyUpdateAll() {
    init();
    calc();
}

void BSpline::notifyPointChange() {

}

void BSpline::notifyAddPoint() {
    init();
    calc();
}

void BSpline::notifyDeletePoint() {

}

void BSpline::init() {
    knots.clear();
    P.clear();
    ts.clear();
    pointNum = 100;
    controlNum = controlPoints.size();
    knotNum = controlNum + order + 1;
    if (type == KnotType::SEQ || isClosed) {
        //if (isClosed) knotNum -= order;
        double step = 1.0 / (knotNum - 1);
        for (int i = 0; i < knotNum - 1; i++) {
            knots.push_back(i * step);
        }
        knots.push_back(1);
    } else {
        for (int i = 0; i < order + 1; i++) knots.push_back(0);
        double step = 1.0 / (knotNum - 2 * (order + 1) + 1);
        for (int i = 1; i <= (knotNum - 2 * (order + 1)); i++) {
            knots.push_back(i * step);
        }
        for (int i = 0; i < order + 1; i++) knots.push_back(1.0);
    }
    for (int i = 0; i < knotNum - 1; i++) {
        P.push_back(1);
    }
    double step = 1.0 / pointNum;
    for (int i = 0; i < pointNum; i++) {
        ts.push_back(step * i);
    }
    ts.push_back(1);
    points = std::vector<QPointF>(ts.size(), QPointF(0, 0));
}

const std::vector<QPointF> &BSpline::getKnotPoints() const {
    return knotPoints;
}

bool BSpline::isClosed1() const {
    return isClosed;
}

void BSpline::setIsClosed(bool isClosed) {
    BSpline::isClosed = isClosed;
}

void BSpline::toOpen() {
    if (!isClosed) return;
    isClosed = false;
    for (int i = 0; i < order; i++) {
        controlPoints.pop_back();
    }
    notifyUpdateAll();
}

void BSpline::toClose() {
    if (isClosed) return;
    isClosed = true;
    //增加p个点
    for (int i = 0; i < order; i++) {
        controlPoints.push_back(controlPoints[i]);
    }
    notifyUpdateAll();
}

void BSpline::printInfo() {
    qDebug() << "n = " << controlPoints.size();
    qDebug() << "m = " << knots.size();
    qDebug() << knots;
    //qDebug() << points;
}
