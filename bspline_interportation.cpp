//
// Created by 徐溶延 on 2020/2/3.
//

#include "bspline_interportation.h"
#include <Eigen/Sparse>
#include <vector>
#include <QDebug>
#include <iostream>

namespace xry_mesh {
    std::vector<double> knots;
    std::vector<double> t;

    int solve_tri_diagnal(const size_t n, const double *l,
                          const double *d, const double *u, double *x) {
        const double eps = 1e-8;
        Eigen::VectorXd dd = Eigen::Map<const Eigen::VectorXd>(d, n);
        for (size_t i = 1; i < n; i++) {
            if (fabs(dd[i - 1]) < eps) return 0;

            double fac = l[i] / dd[i - 1];
            dd[i] -= fac * u[i - 1];
            x[i] -= fac * x[i - 1];
        }
        if (fabs(dd[n - 1]) < eps) {
            return 0;
        }
        x[n - 1] /= dd[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            if (fabs(dd[i]) < eps) {
                return 0;
            }
            x[i] = (x[i] - u[i] * x[i + 1]) / dd[i];
        }
        return 1;
    }

    double cox_de_Boor(int i, int p, double t) {
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

    int prepareValuePointKnot(const Eigen::Matrix2Xd &Q) {
        t.clear();
        double length = 0;
        for (int i = 0; i < (int) Q.cols() - 1; i++) {
            Eigen::Vector2d p1, p2;
            p1 = Q.col(i);
            p2 = Q.col(i + 1);
            double subLen = (p1 - p2).norm();
            t.push_back(length);
            length += subLen;
        }
        t.push_back(length);
        for (auto &v : t) {
            v /= length;
        }
        return 0;
    }

    double prepareKnots(const Eigen::Matrix2Xd &Q) {
        int knotNum = 5 + Q.cols();
        double step = 1.0 / (knotNum - 6 - 1);
        assert(step > 0);
        for (int i = 0; i < 3; i++) knots.push_back(0);
        for (int i = 0; i < knotNum - 6; i++) {
            knots.push_back(i * step);
        }
        for (int i = 0; i < 3; i++) knots.push_back(1);
        return step;
    }

    int prepareMatrix(int size, double step,
                      std::vector<double> &l,
                      std::vector<double> &d,
                      std::vector<double> &u) {
        double h1, h2, hn_1, hn, f1, f2, f3, g1, g2, g3;
        h1 = t[1] - t[0];
        h2 = t[2] - t[1];
        hn_1 = t[t.size() - 2] - t[t.size() - 3];
        hn = t.back() - t[t.size() - 2];
        f1 = h1 + h2;
        f2 = -(2 * h1 + h2);
        f3 = h1;
        g1 = hn_1 + hn;
        g2 = -(2 * hn + hn_1);
        g3 = hn;
        d.push_back(1);
        d.push_back(f2);
        u.push_back(0);
        l.push_back(0);
        std::vector<double> tt;
        tt.push_back(0);
        tt.push_back(0);
        tt.insert(tt.end(), t.begin(), t.end());
        tt.push_back(1);
        tt.push_back(1);
        for (int i = 2; i < size; i++) {
            double N1, N2, N3;
            N1 = (t[i + 1] - t[i]) * (t[i + 1] - t[i]) / ((t[i + 1] - t[i - 2]) * (t[i + 1] - t[i - 1]));
            N2 = (t[i] - t[i - 2]) * (t[i + 1] - t[i]) / ((t[i + 1] - t[i - 2]) * (t[i + 1] - t[i - 1])) +
                 (t[i + 2] - t[i]) * (t[i] - t[i - 1]) / ((t[i + 2] - t[i - 1]) * (t[i + 1] - t[i - 1]));
            N3 = (t[i] - t[i - 1]) * (t[i] - t[i - 1]) / ((t[i + 2] - t[i - 1]) * (t[i + 1] - t[i - 1]));
            d.push_back(N2);
            l.push_back(N1);
            u.push_back(N3);
        }
        d.push_back(g2);
        d.push_back(1);
        l.push_back(1);

        return 0;
    }

    int prepareSparse(Eigen::SparseMatrix<double> &G, const Eigen::Matrix2Xd &Q, std::vector<double > &vecKnots) {
        std::vector<Eigen::Triplet<double> > triplets;
        double h1, h2, hn_1, hn, f1, f2, f3, g1, g2, g3;
        int n = Q.cols();
        h1 = t[1] - t[0];
        h2 = t[2] - t[1];
        hn_1 = t[t.size() - 2] - t[t.size() - 3];
        hn = t.back() - t[t.size() - 2];
        f1 = -3.0 / h1;
        f2 = 3.0 / h1;
        f3 = 0;
        g1 = 3.0 / hn;
        g2 = -3.0 / hn;
        g3 = 0;
        triplets.emplace_back(0, 0, f1);
        triplets.emplace_back(0, 1, f2);
        triplets.emplace_back(0, 2, f3);
        triplets.emplace_back(n + 1, n + 1, g1);
        triplets.emplace_back(n + 1, n, g2);
        triplets.emplace_back(n + 1, n - 1, g3);
        triplets.emplace_back(1, 0, 1);
        triplets.emplace_back(n, n + 1, 1);
        std::vector<double> tt;
        tt.push_back(0);
        tt.push_back(0);
        tt.insert(tt.end(), t.begin(), t.end());
        tt.push_back(1);
        tt.push_back(1);
        knots = tt;
        knots.push_back(1);
        knots.insert(knots.begin(), 0);
        vecKnots = tt;
        qDebug() << "knots:" << knots;
        for (int i = 2; i < n; i++) {
            double N1, N2, N3, N1_, N2_, N3_;
            N1 = (tt[i + 2] - tt[i + 1]) * (tt[i + 2] - tt[i + 1]) / ((tt[i + 2] - tt[i - 1]) * (tt[i + 2] - tt[i]));
            N2 = (tt[i + 1] - tt[i - 1]) * (tt[i + 2] - tt[i + 1]) / ((tt[i + 2] - tt[i - 1]) * (tt[i + 2] - tt[i])) +
                 (tt[i + 3] - tt[i + 1]) * (tt[i + 1] - tt[i]) / ((tt[i + 3] - tt[i]) * (tt[i + 2] - tt[i]));
            N3 = (tt[i + 1] - tt[i]) * (tt[i + 1] - tt[i]) / ((tt[i + 3] - tt[i]) * (tt[i + 2] - tt[i]));
            N1_ = cox_de_Boor(i - 1, 3, knots[i + 2]);
            N2_ = cox_de_Boor(i, 3, knots[i + 2]);
            N3_ = cox_de_Boor(i + 1, 3, knots[i + 2]);
            triplets.emplace_back(i, i - 1, N1);
            triplets.emplace_back(i, i, N2);
            triplets.emplace_back(i, i + 1, N3);
        }
        G.setFromTriplets(triplets.begin(), triplets.end());
        return 0;
    }


    int interporlate_bspline(const Eigen::Matrix2Xd &Q,
                             Eigen::Matrix2Xd &Vout) {
        Eigen::SparseMatrix<double, Eigen::ColMajor> G(Q.cols(), Q.cols());
        std::vector<double> l, d, u, x, y, vecV;
        prepareValuePointKnot(Q);
        double step = prepareKnots(Q);
        prepareMatrix(Q.cols(), step, l, d, u);
        x.push_back(0);
        for (int i = 0; i < Q.cols(); i++) {
            x.push_back(Q(0, i));
        }
        x.push_back(0);
        y.push_back(0);
        for (int i = 0; i < Q.cols(); i++) {
            y.push_back(Q(1, i));
        }
        y.push_back(0);

        std::swap(x[0], x[1]);
        std::swap(x[x.size() - 2], x[x.size() - 1]);
        std::swap(y[0], y[1]);
        std::swap(y[y.size() - 2], y[y.size() - 1]);
        int res1, res2;
        res1 = solve_tri_diagnal(d.size(), l.data(), d.data(), u.data(), x.data());
        res2 = solve_tri_diagnal(d.size(), l.data(), d.data(), u.data(), y.data());
        if (res1 == 0) qDebug() << "solve x failed";
        if (res2 == 0) qDebug() << "solve y failed";
        qDebug() << x;
        qDebug() << y;
        for (int i = 0; i < x.size(); i++) {
            vecV.push_back(x[i]);
            vecV.push_back(y[i]);
        }
        Vout = Eigen::Map<Eigen::Matrix2Xd>(vecV.data(), 2, vecV.size() / 2);
        return 0;
    }

/**
 *
 * @param Q 型值点
 * @param startSlop
 * @param endSlop
 * @param Vout 3次B样条插值点
 * @return
 */
    int interporlate_bspline1(const Eigen::Matrix2Xd &Q,
                              Eigen::Matrix2Xd &Vout,
                              std::vector<double > &vecKnots) {
        Eigen::SparseMatrix<double, Eigen::ColMajor> G(Q.cols() + 2, Q.cols() + 2);
        prepareValuePointKnot(Q);
        prepareKnots(Q);
        Eigen::Vector2d slope1 = (t[1] + t[2] - 2 * t[0]) * (Q.col(1) - Q.col(0)) - (t[1] - t[0]) / (t[2] - t[1]) * (Q.col(2) - Q.col(1));
        int n = Q.cols() - 1;
        int n_1 = Q.cols() - 2;
        int n_2 = Q.cols() - 3;
        Eigen::Vector2d slope2 = (t[n_1] + t[n_2] - 2 * t[n]) * (Q.col(n_1) - Q.col(n)) - (t[n_1] - t[n]) / (t[n_2] - t[n_1]) * (Q.col(n_2) - Q.col(n_1));
        prepareSparse(G, Q, vecKnots);
//        Eigen::Vector2d slope = Q.col(1) - Q.col(Q.cols() - 2);
//        slope.normalize();

        std::cout << G << std::endl;
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > lu(G);
        Eigen::Matrix2Xd _Q(2, Q.cols() + 2);
        _Q.col(0) =  slope1;
        for (int i = 0; i < Q.cols(); i++) {
            _Q.col(i + 1) = Q.col(i);
        }
        _Q.col(Q.cols() + 1) = slope2;
        Eigen::MatrixXd Vt = lu.solve(_Q.transpose());
        Vout = Vt.transpose();
        double error = (G * Vout.transpose() - _Q.transpose()).norm();
        qDebug() << "error =" << error;
        return 0;
    }
}

