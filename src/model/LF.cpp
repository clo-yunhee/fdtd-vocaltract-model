#include "LF.h"

#include <array>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/tools/roots.hpp>
#include <iostream>

namespace math = boost::math;
using namespace boost::math::constants;
using namespace boost::math::quadrature;
using namespace models;
using namespace models::precomp;
using boost::math::cos_pi;
using boost::math::sin_pi;

double LF::evaluate(double t) const {
    static constexpr double T0 = 1;

    double dg;

    if (t < 0 || t > 1) {
        t = -std::fmod(std::abs(t), 1.0);
    }

    if (t <= m_Te) {
        dg = -m_Ee * std::exp(m_alpha * (t - m_Te)) * sin_pi(t / m_Tp) /
             sin_pi(m_Te / m_Tp);
    } else if (!std::isinf(m_epsilon) && !std::isnan(m_epsilon)) {
        dg = -m_Ee / (m_epsilon * m_Ta) *
             (std::exp(-m_epsilon * (t - m_Te)) -
              std::exp(-m_epsilon * (T0 - m_Te)));
    } else {
        dg = 0;
    }
    return dg;
}

double LF::evaluateAntiderivative(double t) const {
    static constexpr double T0 = 1;

    double g;

    if (t < 0 || t > 1) {
        t = -std::fmod(std::abs(t), 1.0);
    }

    if (t <= m_Te) {
        g = -(m_Ee * std::exp(-m_alpha * m_Te)) / sin_pi(m_Te / m_Tp) * 1.0 /
            (m_alpha * m_alpha + pi_sqr<double>() / (m_Tp * m_Tp)) *
            (pi<double>() / m_Tp +
             m_alpha * std::exp(m_alpha * t) * sin_pi(t / m_Tp) -
             pi<double>() / m_Tp * std::exp(m_alpha * t) * cos_pi(t / m_Tp));
    } else {
        // Start integrating from Te.
        // Doing it this way means we are integrating a smooth function.
        const auto& fn = [this](double t) -> double { return evaluate(t); };

        g = m_gTe + gauss_kronrod<double, 31>::integrate(fn, m_Te, t, 3, 1e-3);
    }
    return g;
}

void LF::setRd(const double Rd) {
    static constexpr double Ee = 1;
    static constexpr double T0 = 1;

    // Just interpolate the pre-computed values.
    const int index = std::floor((Rd - Rd_min) / Rd_step);

    m_Ee = Ee;
    std::tie(m_Te, m_Tp, m_Ta, m_alpha, m_epsilon) = Rd_table[index];

    // Store Ug(Te)
    m_gTe = evaluateAntiderivative(m_Te);
}

double LF::Te() const { return m_Te; }
