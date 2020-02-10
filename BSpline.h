//
// Created by 徐溶延 on 2020/1/31.
//

#ifndef B_SPLINE_QT_BSPLINE_H
#define B_SPLINE_QT_BSPLINE_H

#include <QPointF>
#include <vector>



class BSpline {
public:
    enum KnotType{
        SEQ,
        CLAMPED,
    };

    BSpline(int order, KnotType type, std::vector<QPointF> controlPoints);

    BSpline();

    BSpline(int order, const std::vector<QPointF> &controlPoints, const std::vector<double> &knots);

    int getOrder() const;

    void setOrder(int order);

    int getControlNum() const;

    void setControlNum(int controlNum);

    int getKnotNum() const;

    void setKnotNum(int knotNum);

    KnotType getType() const;

    void setType(KnotType type);

    const std::vector<QPointF> &getControlPoints() const;

    void setControlPoints(const std::vector<QPointF> &controlPoints);

    void setControlPoint(int i, const QPointF &pointF);

    const std::vector<double > &getKnots() const;

    void setKnots(const std::vector<double> &knots);

    int getPointNum() const;

    void setPointNum(int pointNum);

    const std::vector<QPointF> &getPoints() const;

    void setPoints(const std::vector<QPointF> &points);

    bool isClosed1() const;

    void setIsClosed(bool isClosed);

    const std::vector<int> &getP() const;

    void setP(const std::vector<int> &p);

    const std::vector<double> &getTs() const;

    void setTs(const std::vector<double> &ts);

    const std::vector<QPointF> &getKnotPoints() const;

    int calc();

    void addPoint(const QPointF &qPointF);

    void deleteLastPoint();

    void deleteAll();

    void notifyUpdateAll();

    void notifyPointChange();

    void notifyAddPoint();

    void notifyDeletePoint();

    void toClose();

    void toOpen();

private:
    int order;
    int controlNum;
    int knotNum;
    int pointNum;
    bool isClosed;
    KnotType type;
    std::vector<QPointF> controlPoints;
    std::vector<QPointF> points;
    std::vector<QPointF> knotPoints;
    std::vector<double> knots;
    std::vector<int> P;
    std::vector<double > ts;
    double cox_de_Boor(int i, int p, double t);
    void init();
    void printInfo();
};


#endif //B_SPLINE_QT_BSPLINE_H
