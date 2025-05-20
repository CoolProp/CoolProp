#include <Eigen/Dense>

template<typename Callable>
auto romberg_diff(Callable& func, double x, std::size_t order=2, double h=0.1){
   
   // Initialize the table
   auto r = Eigen::ArrayXd(order, order);
   
   // Compute the first column using the central difference formula
   for (auto i = 0; i < order; ++i){
       r(i, 0) = (func(x + h) - func(x - h)) / (2 * h);
       h /= 2.0;
   }
   
   // Apply Richardson extrapolation
   for (auto i = 1; i < order; ++i){
       for (auto j = i; j < order; ++j){
           double fouri = pow(4, i);
           r(j, i) = (fouri * r(j, i-1) - r(j-1, i-1)) / (fouri - 1);
       }
   }
   
   return r(order - 1, order - 1);
}