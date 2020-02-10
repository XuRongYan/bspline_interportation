//
// Created by 徐溶延 on 2020/2/3.
//

#ifndef BSPLINE_INTERPORTATION_BSPLINE_INTERPORTATION_H
#define BSPLINE_INTERPORTATION_BSPLINE_INTERPORTATION_H

#include <Eigen/Dense>

namespace xry_mesh {
    int interporlate_bspline(const Eigen::Matrix2Xd &Q,
                             Eigen::Matrix2Xd &Vout);

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
                              std::vector<double> &vecKnots);
}


#endif //BSPLINE_INTERPORTATION_BSPLINE_INTERPORTATION_H
