#ifndef SOURCEMODEL__MODELS_LF_H
#define SOURCEMODEL__MODELS_LF_H

#include <array>
#include <tuple>

namespace models {
class LF {
   public:
    double evaluate(double t) const;
    double evaluateAntiderivative(double t) const;

    void setRd(double Rd);

    // Specific to LF model, get Te directly, used for when Rd is used.
    double Te() const;

   private:
    void fitParameters(double Ee, double T0, double Te, double Tp, double Ta);

    double m_Ee;  // calculated for E0 = 1
    double m_Te;  // = Oq * T0
    double m_Tp;  // = am * Oq * T0
    double m_Ta;  // = Qa * (1 - Oq) * T0

    double m_epsilon;
    double m_alpha;

    double m_gTe;
};

namespace precomp {
#include "LF_precomputed_Rd_double.inc.h"
}  // namespace precomp
}  // namespace models

#endif  //  SOURCEMODEL__MODELS_LF_H