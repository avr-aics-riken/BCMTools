/*
 * BCMTools
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

///
/// @file Vec3.h
/// @brief 3次元ベクトル テンプレートクラス.
/// 

#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <cmath>
#include <stdexcept>

#ifdef BCMT_NAMESPACE
namespace BCMT_NAMESPACE {
#endif


/// 3次元ベクトル テンプレートクラス.
template <typename T>
class Vec3 {
public:

    T x, y, z;

    Vec3(T s = 0) : x(s), y(s), z(s) {}

    Vec3(T x, T y, T z) : x(x), y(y), z(z) {}

    Vec3(const T v[3]) : x(v[0]), y(v[1]), z(v[2]) {}

    Vec3(const Vec3& v) : x(v.x), y(v.y), z(v.z) {}

    Vec3& operator=(const Vec3& v) {
        x = v.x; y = v.y; z = v.z;
        return *this;
    }

    Vec3& operator=(T s) {
        x = y = z = s;
        return *this;
    }

    const T& operator[](int i) const {
        switch (i) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw std::out_of_range("Vec3 index must be 0, 1, or 2");
        }
    }

    Vec3& operator+=(const Vec3& v) {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& v) {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    Vec3& operator*=(T s) {
        x *= s; y *= s; z *= s;
        return *this;
    }

    friend const Vec3 operator+(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
    }

    friend const Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
        return Vec3(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
    }

    const Vec3 operator*(T s) const {
        return Vec3(x*s, y*s, z*s);
    }

    const Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    bool operator==(const Vec3& v) const {
        return x == v.x && y == v.y && z == v.z;
    }

    bool operator!=(const Vec3& v) const {
        return !(*this == v);
    }

    double length() const {
        return sqrt(x*x + y*y + z*z);
    }

    Vec3& normalize() {
        double len = length();
        if (len != 0.0) return *this *= (1.0/len);
        else return *this;
    }

    Vec3& normalize(double& len) {
        len = length();
        if (len != 0.0) return *this *= (1.0/len);
        else return *this;
    }

    Vec3& projectX() {
        x = 0;
        return *this;
    }

    Vec3& projectY() {
        y = 0;
        return *this;
    }

    Vec3& projectZ() {
        z = 0;
        return *this;
    }


};

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const ::Vec3<T>& v) {
    return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

template <typename T>
inline std::istream& operator>>(std::istream& is, ::Vec3<T>& v) {
    return is >> v.x >> v.y >> v.z;
}


/// 3次元整数ベクトルクラス.
typedef ::Vec3<int> Vec3i;

/// 3次元実数ベクトルクラス.
typedef ::Vec3<double> Vec3r;


inline double dotProduct(const Vec3r& u, const Vec3r& v) {
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline const Vec3r crossProduct(const Vec3r& u, const Vec3r& v) {
    return Vec3r(u.y*v.z - u.z*v.y,
                 u.z*v.x - u.x*v.z,
                 u.x*v.y - u.y*v.x);
}

inline double distance(const Vec3r& u, const Vec3r& v) {
    return (u - v).length();
}


#ifdef BCMT_NAMESPACE
} // namespace BCMT_NAMESPACE
#endif

#endif // VEC3_H
