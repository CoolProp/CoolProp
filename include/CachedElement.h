/*
 * CachedElement.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef CACHEDELEMENT_H_
#define CACHEDELEMENT_H_

#include <algorithm>
#include <array>
#include "CoolPropTools.h"
#include "DataStructures.h"
#include "CPnumerics.h"

namespace CoolProp {

/*!
 A class that contains the magic to cache a value.
 
 Includes an "=" assignment operator and casting to boolean
 so you can do something like::
 
 double CoolPropStateClassSI::d3phir_dTau3(double tau, double delta){
 if (cache.d3phir_dTau3)    {
 return cache.d3phir_dTau3;
 } else {
 cache.d3phir_dTau3 = pFluid->d3phir_dTau3(tau,delta);
 return cache.d3phir_dTau3;
 }
 };
 */


class CachedElement
{
    
private:
    bool is_cached = false;
    double value;
    
public:
    /// Default constructor
    CachedElement() {
        this->clear();
    };
    
    /// Function to carry out the caching
    void _do_cache(double value) {
        this->value = value;
        this->is_cached = true;
    }
    
    /// Assignment operator - sets the value and sets the flag
    void operator=(const double& value) {
        _do_cache(value);
    };
    
    /// Cast to boolean, for checking if cached
    operator bool() {
        return is_cached;
    };
    
    /// Cast to double, for returning value
    operator double() {
        if (is_cached) {
            return static_cast<double>(value);
        } else {
            throw std::exception();
        }
    }
#ifndef COOLPROPDBL_MAPS_TO_DOUBLE
    operator CoolPropDbl() {
        if (is_cached) {
            return value;
        } else {
            throw std::exception();
        }
    }
#endif
    /// Clear the flag and the value
    void clear() {
        is_cached = false;
        this->value = _HUGE;
    };
    double& pt() {
        return this->value;
    }
};

template<typename NumType>
class CacheArrayElement
{
    
private:
    NumType& value;
    bool& is_cached;
    
public:
    
    // Constructor with value
    CacheArrayElement(NumType& val, bool& is_cached) : value(val), is_cached(is_cached) {};
    
    /// Function to carry out the caching
    void _do_cache(double value) {
        this->value = value;
        this->is_cached = true;
    }
    
    /// Assignment operator - sets the value and sets the flag
    void operator=(const double& value) {
        _do_cache(value);
    };
    
    /// Cast to boolean, for checking if cached
    operator bool() {
        return is_cached;
    };
    
    /// Cast to double, for returning value
    operator double() {
        if (is_cached) {
            return static_cast<double>(value);
        } else {
            throw std::exception();
        }
    }
#ifndef COOLPROPDBL_MAPS_TO_DOUBLE
    operator CoolPropDbl() {
        if (is_cached) {
            return value;
        } else {
            throw std::exception();
        }
    }
#endif
    /// Clear the flag and the value
    void clear() {
        is_cached = false;
        this->value = _HUGE;
    };
    NumType& pt() {
        return this->value;
    }
};


template<int N>
class CacheArray{

private:

    std::size_t inext = 0;
    std::array<double, N> m_values = create_filled_array<double, N>(_HUGE);
    std::array<bool, N> m_cached = create_filled_array<bool, N>(false);

   public:
    void clear(){
        memset(m_values.data(), 0, sizeof(m_values));
        memset(m_cached.data(), false, sizeof(m_cached));
    }
    auto factory(std::size_t i){
        return CacheArrayElement<double>(m_values[i], m_cached[i]);
    }
    auto next(){
        if (inext > N){
            throw ValueError("No more cache elements available");
        }
        auto el = factory(inext);
        inext++;
        return el;
    }
};

} /* namespace CoolProp */
#endif /* CACHEDELEMENT_H_ */
