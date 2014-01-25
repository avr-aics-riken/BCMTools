/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __MATH_VEC_H__
#define __MATH_VEC_H__

#include <cmath>

template <typename T>
class vec2
{
public:
	vec2() { m_vec[0] = m_vec[1] = 0; }

	vec2(const T x, const T y) {
		m_vec[0] = x;
		m_vec[1] = y;
	}

	vec2(const vec2 &v ){
		m_vec[0] = v[0];
		m_vec[1] = v[1];
	}

	void operator =(const vec2 &v){
		m_vec[0] = v[0];
		m_vec[1] = v[1];
	}

	T operator [](const int i) const {
		return m_vec[i];	
	}

	vec2 operator +(const vec2 &v) const {
		return vec2(m_vec[0] + v[0], m_vec[1] + v[1]);
	}

	vec2 operator -(const vec2 &v) const {
		return vec2(m_vec[0] - v[0], m_vec[1] - v[1]);
	}
private:
	T m_vec[2];
};

template <typename T>
class vec3
{
public:
	vec3() {
		m_vec[0] = m_vec[1] = m_vec[2] = 0;		
	}
	
	vec3(const T x, const T y, const T z) {
		m_vec[0] = x;
		m_vec[1] = y;
		m_vec[2] = z;
	}

	vec3(const vec3 &v) {
		m_vec[0] = v[0];
		m_vec[1] = v[1];
		m_vec[2] = v[2];
	}


	void operator =(const vec3 &v) {
		m_vec[0] = v[0];
		m_vec[1] = v[1];
		m_vec[2] = v[2];
	}

	T operator [](const int i) const {
		return m_vec[i];
	}

	vec3 operator +(const vec3 &v) const {
		return vec3(m_vec[0] + v[0], m_vec[1] + v[1], m_vec[2] + v[2]);
	}

	vec3 operator -(const vec3 &v) const {
		return vec3(m_vec[0] - v[0], m_vec[1] - v[1], m_vec[2] - v[2]);
	}

	T dot(const vec3 &v) const {
		return (m_vec[0] * v[0] + m_vec[1] * v[1] + m_vec[2] * v[2]);
	}

	vec3 cross(const vec3 &v) const {
		return vec3(m_vec[1] * v[2] - m_vec[2] * v[1],
		            m_vec[2] * v[0] - m_vec[0] * v[2],
					m_vec[0] * v[1] - m_vec[1] * v[0]);
	}

	T norm() const
	{
		return sqrt(m_vec[0] * m_vec[0] + m_vec[1] * m_vec[1] + m_vec[2] * m_vec[2]);
	}

	T dot(const vec3 &v1, const vec3 &v2) const {
		return v1.dot(v2);
	}

	vec3 cross(const vec3 &v1, const vec3 &v2) const {
		return v1.cross(v2);
	}

	void normalize(){
		T _norm = this->norm();

		if(_norm != 0){
			m_vec[0] /= _norm;
			m_vec[1] /= _norm;
			m_vec[2] /= _norm;
		}
	}

	T x(){ return m_vec[0]; }
	T y(){ return m_vec[1]; }
	T z(){ return m_vec[2]; }

	T* ptr(){ return &(m_vec[0]); }
private:
	T m_vec[3];
};


template <typename T>
class vec4
{
public:
	vec4() {
		m_vec[0] = m_vec[1] = m_vec[2] = m_vec[3] = 0;		
	}
	
	vec4(const T x, const T y, const T z, const T w) {
		m_vec[0] = x;
		m_vec[1] = y;
		m_vec[2] = z;
		m_vec[3] = w;
	}

	vec4(const vec4 &v) {
		m_vec[0] = v[0];
		m_vec[1] = v[1];
		m_vec[2] = v[2];
		m_vec[3] = v[3];
	}

	vec4(const vec3<T> &v, const T w){
		m_vec[0] = v[0];
		m_vec[1] = v[1];
		m_vec[2] = v[2];
		m_vec[3] = w;
	}


	void operator =(const vec4 &v) {
		m_vec[0] = v[0];
		m_vec[1] = v[1];
		m_vec[2] = v[2];
		m_vec[3] = v[3];
	}

	T operator [](const int i) const {
		return m_vec[i];
	}

	vec4 operator +(const vec4 &v) const {
		return vec4(m_vec[0] + v[0], m_vec[1] + v[1], m_vec[2] + v[2], m_vec[3] + v[3]);
	}

	vec4 operator -(const vec4 &v) const {
		return vec4(m_vec[0] - v[0], m_vec[1] - v[1], m_vec[2] - v[2], m_vec[3] - v[3]);
	}

	T norm() const
	{
		return sqrt(m_vec[0] * m_vec[0] + m_vec[1] * m_vec[1] + m_vec[2] * m_vec[2] + m_vec[3] * m_vec[3]);
	}

	void normalize(){
		T _norm = this->norm();

		if(_norm != 0){
			m_vec[0] /= _norm;
			m_vec[1] /= _norm;
			m_vec[2] /= _norm;
			m_vec[3] /= _norm;
		}
	}

	T x(){ return m_vec[0]; }
	T y(){ return m_vec[1]; }
	T z(){ return m_vec[2]; }
	T w(){ return m_vec[3]; }

	T* ptr(){ return &(m_vec[0]); }
private:
	T m_vec[4];
};


#endif // __MATH_VEC_H__
