#include "ReducingFunctions.h"

namespace CoolProp{

long double ReducingFunction::d_ndTrdni_dxj__constxi(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        long double s = 0;
        for (std::size_t k = 0; k < N; k++)
        {
            s += x[k]*d2Trdxidxj(x,j,k, xN_flag);
        }
        return d2Trdxidxj(x,i,j, xN_flag)-dTrdxi__constxj(x,j, xN_flag)-s;
    }
    else if (xN_flag == XN_DEPENDENT){
        if (j == N-1){ return 0;}
        long double s = 0;
        for (std::size_t k = 0; k < N-1; k++)
        {
            s += x[k]*d2Trdxidxj(x,k,j, xN_flag);
        }
        long double val = d2Trdxidxj(x,j,i,xN_flag)-dTrdxi__constxj(x,j, xN_flag)-s;
        return val;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
    
}
long double ReducingFunction::d_ndrhorbardni_dxj__constxi(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    long double s = 0;
    for (std::size_t k = 0; k < N; k++)
    {
        s += x[k]*d2rhormolardxidxj(x,j,k, xN_flag);
    }
    return d2rhormolardxidxj(x,j,i, xN_flag)-drhormolardxi__constxj(x,j, xN_flag)-s;
}
long double ReducingFunction::ndrhorbardni__constnj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        long double summer_term1 = 0;
        for (std::size_t j = 0; j < N; j++)
        {
            summer_term1 += x[j]*drhormolardxi__constxj(x,j, xN_flag);
        }
        return drhormolardxi__constxj(x,i, xN_flag)-summer_term1;
    }
    else if (xN_flag == XN_DEPENDENT){
        long double summer_term1 = 0;
        for (std::size_t k = 0; k < N-1; ++k)
        {
            summer_term1 += x[k]*drhormolardxi__constxj(x, k, xN_flag);
        }
        return drhormolardxi__constxj(x, i, xN_flag)-summer_term1;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
}
long double ReducingFunction::ndTrdni__constnj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        // GERG Equation 7.54
        long double summer_term1 = 0;
        for (std::size_t j = 0; j < N; j++)
        {
            summer_term1 += x[j]*dTrdxi__constxj(x,j, xN_flag);
        }
        return dTrdxi__constxj(x,i, xN_flag)-summer_term1;
    }
    else if (xN_flag == XN_DEPENDENT){
        long double summer_term1 = 0;
        for (std::size_t k = 0; k < N-1; ++k)
        {
            summer_term1 += x[k]*dTrdxi__constxj(x, k, xN_flag);
        }
        return dTrdxi__constxj(x, i, xN_flag)-summer_term1;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
}

long double GERG2008ReducingFunction::Tr(const std::vector<long double> &x)
{
    return Yr(x, beta_T, gamma_T, T_c, Yc_T);
}
long double GERG2008ReducingFunction::dTrdxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    return dYrdxi__constxj(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
long double GERG2008ReducingFunction::d2Trdxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    return d2Yrdxi2__constxj(x, i, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
long double GERG2008ReducingFunction::d2Trdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    return d2Yrdxidxj(x, i, j, beta_T, gamma_T, T_c, Yc_T, xN_flag);
}
long double GERG2008ReducingFunction::rhormolar(const std::vector<long double> &x)
{
    return 1/Yr(x, beta_v, gamma_v, v_c, Yc_v);
}
long double GERG2008ReducingFunction::drhormolardxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -pow(rhormolar(x),2)*dvrmolardxi__constxj(x, i, xN_flag);
}
long double GERG2008ReducingFunction::dvrmolardxi__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    return dYrdxi__constxj(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
long double GERG2008ReducingFunction::d2vrmolardxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    return d2Yrdxi2__constxj(x, i, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
long double GERG2008ReducingFunction::d2vrmolardxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    return d2Yrdxidxj(x, i, j, beta_v, gamma_v, v_c, Yc_v, xN_flag);
}
long double GERG2008ReducingFunction::d2rhormolardxi2__constxj(const std::vector<long double> &x, std::size_t i, x_N_dependency_flag xN_flag)
{
    long double rhor = this->rhormolar(x);
    long double dvrbardxi = this->dvrmolardxi__constxj(x,i, xN_flag);
    return 2*pow(rhor,3)*pow(dvrbardxi,2)-pow(rhor,2)*this->d2vrmolardxi2__constxj(x,i, xN_flag);
}
long double GERG2008ReducingFunction::d2rhormolardxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double rhor = this->rhormolar(x);
    double dvrbardxi = this->dvrmolardxi__constxj(x,i, xN_flag);
    double dvrbardxj = this->dvrmolardxi__constxj(x,j, xN_flag);
    return 2*pow(rhor,3)*dvrbardxi*dvrbardxj-pow(rhor,2)*this->d2vrmolardxidxj(x,i,j, xN_flag);
}

long double GERG2008ReducingFunction::Yr(const std::vector<long double> &x, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc)
{
    long double Yr = 0;
    for (std::size_t i = 0; i < N; i++)
    {
        double xi = x[i];
        Yr += xi*xi*Yc[i];
        
        // The last term is only used for the pure component, as it is sum_{i=1}^{N-1}sum_{j=1}^{N}
        if (i==N-1){ break; }

        for (std::size_t j = i+1; j < N; j++)
        {
            Yr += c_Y_ij(i, j, beta, gamma, Y_c_ij)*f_Y_ij(x, i, j, beta);
        }
    }
    return Yr;
}
long double GERG2008ReducingFunction::dYrdxi__constxj(const std::vector<long double> &x, std::size_t i,  const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        long double xi = x[i];
        long double dYr_dxi = 2*xi*Yc[i];
        for (std::size_t k = 0; k < i; k++)
        {
            dYr_dxi += c_Y_ij(k,i,beta,gamma,Y_c_ij)*dfYkidxi__constxk(x,k,i,beta);
        }
        for (std::size_t k = i+1; k < N; k++)
        {
            dYr_dxi += c_Y_ij(i,k,beta,gamma,Y_c_ij)*dfYikdxi__constxk(x,i,k,beta);
        }
        return dYr_dxi;
    }
    else if (xN_flag == XN_DEPENDENT){
        // Table S1 from Gernert, 2014, supplemental information
        if (i == N-1){return 0.0;}
        long double dYr_dxi = 2*x[i]*Yc[i] - 2*x[N-1]*Yc[N-1];
        for (std::size_t k = 0; k < i; k++)
        {
            dYr_dxi += c_Y_ij(k, i, beta, gamma, Y_c_ij)*dfYkidxi__constxk(x,k,i,beta);
        }
        for (std::size_t k = i+1; k < N-1; k++)
        {
            dYr_dxi += c_Y_ij(i, k, beta, gamma, Y_c_ij)*dfYikdxi__constxk(x,i,k,beta);
        }
        double beta_Y_iN = beta[i][N-1], xN = x[N-1];
        dYr_dxi += c_Y_ij(i, N-1, beta, gamma, Y_c_ij)*(xN*(x[i]+xN)/(pow(beta_Y_iN,2)*x[i]+xN)+(1-beta_Y_iN*beta_Y_iN)*x[i]*xN*xN/pow(beta_Y_iN*beta_Y_iN*x[i]+xN, 2));
        for (std::size_t k = 0; k < N-1; ++k)
        {
            double beta_Y_kN = beta[k][N-1];
            dYr_dxi += c_Y_ij(k, N-1, beta, gamma, Y_c_ij)*(-x[k]*(x[k]+xN)/(pow(beta_Y_kN,2)*x[k]+xN)+(1-beta_Y_kN*beta_Y_kN)*xN*x[k]*x[k]/pow(beta_Y_kN*beta_Y_kN*x[k]+xN, 2));        
        }
        return dYr_dxi;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
}
long double GERG2008ReducingFunction::d2Yrdxi2__constxj(const std::vector<long double> &x, std::size_t i,  const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        // See Table B9 from Kunz Wagner 2012 (GERG 2008)
        long double d2Yr_dxi2 = 2*Yc[i];
        for (std::size_t k = 0; k < i; k++)
        {
            d2Yr_dxi2 += c_Y_ij(k,i,beta,gamma,Y_c_ij)*d2fYkidxi2__constxk(x,k,i,beta);
        }
        for (std::size_t k = i+1; k < N; k++)
        {
            d2Yr_dxi2 += c_Y_ij(i,k,beta,gamma,Y_c_ij)*d2fYikdxi2__constxk(x,i,k,beta);
        }
        return d2Yr_dxi2;
    }
    else if (xN_flag == XN_DEPENDENT){
        // Table S1 from Gernert, 2014, supplemental information
        if (i == N-1){return 0.0;}
        long double d2Yr_dxi2 = 2*Yc[i] + 2*Yc[N-1];
        for (std::size_t k = 0; k < i; k++)
        {
            d2Yr_dxi2 += c_Y_ij(k, i, beta, gamma, Y_c_ij)*d2fYkidxi2__constxk(x,k,i,beta);
        }
        for (std::size_t k = i+1; k < N-1; k++)
        {
            d2Yr_dxi2 += c_Y_ij(i, k, beta, gamma, Y_c_ij)*d2fYikdxi2__constxk(x,i,k,beta);
        }
        double beta_Y_iN = beta[i][N-1], xN = x[N-1];
        d2Yr_dxi2 += 2*c_Y_ij(i, N-1, beta, gamma, Y_c_ij)*(-(x[i]+xN)/(pow(beta_Y_iN,2)*x[i]+xN)+(1-beta_Y_iN*beta_Y_iN)*(xN*xN/pow(beta_Y_iN*beta_Y_iN*x[i]+xN, 2)+((1-beta_Y_iN*beta_Y_iN)*x[i]*xN*xN-beta_Y_iN*beta_Y_iN*x[i]*x[i]*xN)/pow(beta_Y_iN*beta_Y_iN*x[i]+xN, 3)));
        for (std::size_t k = 0; k < N-1; ++k)
        {
            double beta_Y_kN = beta[k][N-1];
            d2Yr_dxi2 += 2*c_Y_ij(k, N-1, beta, gamma, Y_c_ij)*x[k]*x[k]*(1-beta_Y_kN*beta_Y_kN)/pow(pow(beta_Y_kN,2)*x[k]+xN, 2)*(xN/(beta_Y_kN*beta_Y_kN*x[k]+xN)-1);
        }
        return d2Yr_dxi2;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
}

long double GERG2008ReducingFunction::d2Yrdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c_ij, const std::vector<long double> &Yc, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        if (i == j)
        {
            return d2Yrdxi2__constxj(x, i, beta, gamma, Y_c_ij, Yc, xN_flag);
        }
        else
        {
            // See Table B9 from Kunz Wagner 2012 (GERG 2008)
            return c_Y_ij(i, j, beta, gamma, Y_c_ij)*d2fYijdxidxj(x, i, j, beta);
        }
    }
    else if (xN_flag == XN_DEPENDENT){
        // Table S1 from Gernert, 2014, supplemental information
        if (j == N-1 || i == N-1){ return 0.0; }
        if (i == j){ return d2Yrdxi2__constxj(x, i, beta, gamma, Y_c_ij, Yc, xN_flag); }
        long double xN = x[N-1];
        long double d2Yr_dxidxj = 2*Yc[N-1];
        d2Yr_dxidxj += c_Y_ij(i, j, beta, gamma, Y_c_ij)*d2fYijdxidxj(x,i,j,beta);
        
        for (std::size_t k = 0; k < N-1; k++)
        {
            long double beta_Y_kN = beta[k][N-1];
            d2Yr_dxidxj += 2*c_Y_ij(k, N-1, beta, gamma, Y_c_ij)*x[k]*x[k]*(1-beta_Y_kN*beta_Y_kN)/pow(beta_Y_kN*beta_Y_kN*x[k]+xN,2)*(xN/(beta_Y_kN*beta_Y_kN*x[k]+xN)-1);
        }
        {
            long double beta_Y_iN = beta[i][N-1];
            d2Yr_dxidxj +=  c_Y_ij(i, N-1, beta, gamma, Y_c_ij)*((1-beta_Y_iN*beta_Y_iN)*(2*x[i]*xN*xN/pow(beta_Y_iN*beta_Y_iN*x[i]+xN, 3)-x[i]*xN/pow(beta_Y_iN*beta_Y_iN*x[i]+xN, 2)) - (x[i]+xN)/(beta_Y_iN*beta_Y_iN*x[i]+xN));
        }
        {
            long double beta_Y_jN = beta[j][N-1];
            d2Yr_dxidxj += -c_Y_ij(j, N-1, beta, gamma, Y_c_ij)*((1-beta_Y_jN*beta_Y_jN)*(2*x[j]*x[j]*xN*beta_Y_jN*beta_Y_jN/pow(beta_Y_jN*beta_Y_jN*x[j]+xN, 3)-x[j]*xN/pow(beta_Y_jN*beta_Y_jN*x[j]+xN, 2)) + (x[j]+xN)/(beta_Y_jN*beta_Y_jN*x[j]+xN));
        }
        return d2Yr_dxidxj;
    }
    else{
        throw ValueError(format("xN dependency flag invalid"));
    }
}

long double GERG2008ReducingFunction::dfYkidxi__constxk(const std::vector<long double> &x, std::size_t k, std::size_t i, const STLMatrix &beta)
{
    double xk = x[k], xi = x[i], beta_Y = beta[k][i];
    return xk*(xk+xi)/(beta_Y*beta_Y*xk+xi)+xk*xi/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi));
}
long double GERG2008ReducingFunction::dfYikdxi__constxk(const std::vector<long double> &x, std::size_t i, std::size_t k, const STLMatrix &beta)
{
    double xk = x[k], xi = x[i], beta_Y = beta[i][k];
    return xk*(xi+xk)/(beta_Y*beta_Y*xi+xk)+xi*xk/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk));
}
long double GERG2008ReducingFunction::c_Y_ij(std::size_t i, std::size_t j, const STLMatrix &beta, const STLMatrix &gamma, const STLMatrix &Y_c)
{
    return 2*beta[i][j]*gamma[i][j]*Y_c[i][j];
}
long double GERG2008ReducingFunction::f_Y_ij(const std::vector<long double> &x, std::size_t i, std::size_t j, const STLMatrix &beta)
{
    double xi = x[i], xj = x[j], beta_Y = beta[i][j];
    return xi*xj*(xi+xj)/(beta_Y*beta_Y*xi+xj);
}
long double GERG2008ReducingFunction::d2fYikdxi2__constxk(const std::vector<long double> &x, std::size_t i, std::size_t k, const STLMatrix &beta)
{
    double xi = x[i], xk = x[k], beta_Y = beta[i][k];
    return 1/(beta_Y*beta_Y*xi+xk)*(1-beta_Y*beta_Y*(xi+xk)/(beta_Y*beta_Y*xi+xk))*(2*xk-xi*xk*2*beta_Y*beta_Y/(beta_Y*beta_Y*xi+xk));
}
long double GERG2008ReducingFunction::d2fYkidxi2__constxk(const std::vector<long double> &x, std::size_t k, std::size_t i, const STLMatrix &beta)
{
    double xi = x[i], xk = x[k], beta_Y = beta[k][i];
    return 1/(beta_Y*beta_Y*xk+xi)*(1-(xk+xi)/(beta_Y*beta_Y*xk+xi))*(2*xk-xk*xi*2/(beta_Y*beta_Y*xk+xi));
}
long double GERG2008ReducingFunction::d2fYijdxidxj(const std::vector<long double> &x, std::size_t i, std::size_t j, const STLMatrix &beta)
{
    double xi = x[i], xj = x[j], beta_Y = beta[i][j], beta_Y2 = beta_Y*beta_Y;
    return (xi+xj)/(beta_Y2*xi+xj) + xj/(beta_Y2*xi+xj)*(1-(xi+xj)/(beta_Y2*xi+xj))
        +xi/(beta_Y2*xi+xj)*(1-beta_Y2*(xi+xj)/(beta_Y2*xi+xj))
        -xi*xj/pow(beta_Y2*xi+xj,2)*(1+beta_Y2-2*beta_Y2*(xi+xj)/(beta_Y2*xi+xj));
}

} /* namespace CoolProp */
