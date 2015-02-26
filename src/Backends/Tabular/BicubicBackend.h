#ifndef BICUBICBACKEND_H
#define BICUBICBACKEND_H

#include "TabularBackends.h"

namespace CoolProp
{

/** \brief This class implements bicubic interpolation, as very clearly laid out by
 * the page on wikipedia: http://en.wikipedia.org/wiki/Bicubic_interpolation
 * 
 * Essentially you have an already-inverted matrix that you need to multiply 

\f[
A^{-1} =  \left[ \begin{array}{*{16}c} 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 \\
 -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 \\
 9 & -9 & -9 & 9 & 6 & 3 & -6 & -3 & 6 & -6 & 3 & -3 & 4 & 2 & 2 & 1 \\
 -6 & 6 & 6 & -6 & -3 & -3 & 3 & 3 & -4 & 4 & -2 & 2 & -2 & -2 & -1 & -1 \\
 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 \\
 -6 & 6 & 6 & -6 & -4 & -2 & 4 & 2 & -3 & 3 & -3 & 3 & -2 & -1 & -2 & -1 \\
 4 & -4 & -4 & 4 & 2 & 2 & -2 & -2 & 2 & -2 & 2 & -2 & 1 & 1 & 1 & 1 \end{array} \right]
 \f]

 \f[
 x = \frac{h-h_i}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial x}{\partial h} = \frac{1}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial h}{\partial x} = h_{i+1}-h_i
 \f]

 \f[
 y = \frac{p-p_j}{p_{j+1}-p_j}
 \f]
 \f[
 \frac{\partial y}{\partial p} = \frac{1}{p_{j+1}-p_j}
 \f]
  \f[
 \frac{\partial p}{\partial y} = p_{j+1}-p_j
 \f]

 \f[
 \frac{\partial f}{\partial x} = \frac{\partial f}{\partial h}\cdot\frac{\partial h}{\partial x}
 \f]
 \
*/
class BicubicBackend : public TabularBackend
{
    public:
        std::string backend_name(void){return "BicubicBackend";}
        BicubicBackend(shared_ptr<CoolProp::AbstractState> AS) : TabularBackend (AS) {}
};

double do_one();
}

#endif // BICUBICBACKEND_H
