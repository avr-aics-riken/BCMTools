/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __MATH_MAT_H__
#define __MATH_MAT_H__

#include "vec.h"
#include <cmath>
#include <cstring>

#include <iostream>

template <typename T>
class mat3
{
public:
	mat3(){
		identity();
	}

	mat3(const mat3 &m){
		memcpy(&(m_mat[0]), m.ptr(), 9 * sizeof(T));
	}

	mat3(const T* p){
		memcpy(&(m_mat[0]), p, 9 * sizeof(T));
	}

	mat3(const vec3<T>& v1, const vec3<T>& v2, const vec3<T>& v3){
		m_mat[0] = v1[0]; m_mat[1] = v2[0]; m_mat[2] = v3[0];
		m_mat[3] = v1[1]; m_mat[4] = v2[1]; m_mat[5] = v3[1];
		m_mat[6] = v1[2]; m_mat[7] = v2[2]; m_mat[8] = v3[2];
	}

	void operator =(const mat3 m){
		memcpy(&(m_mat[0]), m.ptr(), 9 * sizeof(T));
	}

	T operator [](const int i) const {
		return m_mat[i];
	}

	T at(int u, int v) const {
		return m_mat[u + 3 * v];
	}

	T& operator [](const int i){
		return m_mat[i];
	}

	T& at(int u, int v){
		return m_mat[u + 3 * v];
	}

	mat3 operator +(const mat3 &m) const {
		T ret_mat[9] = {
			m_mat[0] + m[0], m_mat[1] + m[1], m_mat[2] + m[2],
			m_mat[3] + m[3], m_mat[4] + m[4], m_mat[5] + m[5],
			m_mat[6] + m[6], m_mat[7] + m[7], m_mat[8] + m[8]
		};

		return mat3(&ret_mat[0]);
	}

	mat3 operator -(const mat3 &m) const {
		T ret_mat[9] = {
			m_mat[0] - m[0], m_mat[1] - m[1], m_mat[2] - m[2],
			m_mat[3] - m[3], m_mat[4] - m[4], m_mat[5] - m[5],
			m_mat[6] - m[6], m_mat[7] - m[7], m_mat[8] - m[8]
		};

		return mat3(&ret_mat[0]);
	}

	mat3 operator *(const mat3 &m) const {
		T ret_mat[9] = {0};
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				for(int k = 0; k < 3; k++){
					ret_mat[j + i * 3] += m_mat[k + i * 3] * m[j + k * 3];
				}
			}
		}
		
		return mat3(&ret_mat[0]);
	}

	vec3<T> operator *(const vec3<T> &v) const {
		T ret_vec[3] = {0};
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				ret_vec[i] += m_mat[j + i * 3] * v[j];
			}
		}
		return vec3<T>(ret_vec[0], ret_vec[1], ret_vec[2]);
	}

	mat3 operator *(const T s) const {
		T ret_mat[9] = {
			m_mat[0] * s, m_mat[1] * s, m_mat[2] * s,
			m_mat[3] * s, m_mat[4] * s, m_mat[5] * s,
			m_mat[6] * s, m_mat[7] * s, m_mat[8] * s
		};

		return mat3(ret_mat);
	}

	mat3 operator /(const T s) const {
		T ret_mat[9] = {
			m_mat[0] / s, m_mat[1] / s, m_mat[2] / s,
			m_mat[3] / s, m_mat[4] / s, m_mat[5] / s,
			m_mat[6] / s, m_mat[7] / s, m_mat[8] / s
		};

		return mat3(ret_mat);
	}

	void identity(){
		m_mat[0] = 1; m_mat[1] = 0; m_mat[2] = 0;
		m_mat[3] = 0; m_mat[4] = 1; m_mat[5] = 0;
		m_mat[6] = 0; m_mat[7] = 0; m_mat[8] = 1;
	}

	void transpose()
	{
		T back[9];
		memcpy(&back[0], &m_mat[0], 9 * sizeof(T));

		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				m_mat[j + i * 3] = back[i + j * 3];
			}
		}
	}

	//T det();
	//void inv();

	const T* ptr() const {
		return &(m_mat[0]);
	}

	T* ptr() {
		return &(m_mat[0]);
	}

	void print()
	{
		for(int c = 0; c < 3; c++){
			for(int r = 0; r < 3; r++){
				std::cout.width(10);
				std::cout << m_mat[r + c * 3] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

private:
	T m_mat[9];
};


template <typename T>
class mat4
{
public:
	mat4(){
		identity();
	}

	mat4(const mat4 &m){
		memcpy(&(m_mat[0]), m.ptr(), 16 * sizeof(T));
	}

	mat4(const T* p){
		memcpy(&(m_mat[0]), p, 16 * sizeof(T));
	}

	mat4(const vec4<T>& v1, const vec4<T>& v2, const vec4<T>& v3, const vec4<T> v4){
			m_mat[ 0] = v1[0]; m_mat[ 1] = v2[0]; m_mat[ 2] = v3[0]; m_mat[ 3] = v4[0];
			m_mat[ 4] = v1[1]; m_mat[ 5] = v2[1]; m_mat[ 6] = v3[1]; m_mat[ 7] = v4[1];
			m_mat[ 8] = v1[2]; m_mat[ 9] = v2[2]; m_mat[10] = v3[2]; m_mat[11] = v4[2];
			m_mat[12] = v1[3]; m_mat[13] = v2[3]; m_mat[14] = v3[3]; m_mat[15] = v4[3];
	}

	void operator =(const mat4 m){
		memcpy(&(m_mat[0]), m.ptr(), 16 * sizeof(T));
	}

	T operator [](const int i) const {
		return m_mat[i];
	}

	T at(int u, int v) const {
		return m_mat[u + 4 * v];
	}

	T& operator [](const int i) {
		return m_mat[i];
	}

	T& at(int u, int v) {
		return m_mat[u + 4 * v];
	}

	mat4 operator +(const mat4 &m) const {
		T ret_mat[16] = { 0 };
		for(int i = 0; i < 16; i++){
			ret_mat[i] = m_mat[i] + m[i];
		}

		return mat4(&ret_mat[0]);
	}

	mat4 operator -(const mat4 &m) const {
		T ret_mat[16] = { 0 };
		for(int i = 0; i < 16; i++){
			ret_mat[i] = m_mat[i] - m[i];
		}

		return mat4(&ret_mat[0]);
	}

	mat4 operator *(const mat4 &m) const {
		T ret_mat[16] = {0};
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					ret_mat[j + i * 4] += m_mat[k + i * 4] * m[j + k * 4];
				}
			}
		}
		
		return mat4(&ret_mat[0]);
	}

	vec3<T> operator *(const vec4<T> &v) const {
		T ret_vec[4] = {0};
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				ret_vec[i] += m_mat[j + i * 4] * v[j];
			}
		}
		return vec3<T>(ret_vec[0], ret_vec[1], ret_vec[2], ret_vec[3]);
	}

	mat4 operator *(const T s) const {
		T ret_mat[16] = { 0 };
		for(int i = 0; i < 16; i++){
			ret_mat[i] = m_mat[i] * s;
		}

		return mat4(ret_mat);
	}

	mat4 operator /(const T s) const {
		T ret_mat[16] = { 0 };
		for(int i = 0; i < 16; i++){
			ret_mat[i] = m_mat[i] / s;
		}

		return mat4(ret_mat);
	}

	void identity(){
		m_mat[ 0] = 1; m_mat[ 1] = 0; m_mat[ 2] = 0; m_mat[ 3] = 0;
		m_mat[ 4] = 0; m_mat[ 5] = 1; m_mat[ 6] = 0; m_mat[ 7] = 0;
		m_mat[ 8] = 0; m_mat[ 9] = 0; m_mat[10] = 1; m_mat[11] = 0;
		m_mat[12] = 0; m_mat[13] = 0; m_mat[14] = 0; m_mat[15] = 1;
	}

	void transpose()
	{
		T back[16];
		memcpy(&back[0], &m_mat[0], 9 * sizeof(T));

		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				m_mat[j + i * 4] = back[i + j * 4];
			}
		}
	}

	//T det();
	//void inv();

	const T* ptr() const {
		return &(m_mat[0]);
	}

	T* ptr() {
		return &(m_mat[0]);
	}


	void print()
	{
		for(int c = 0; c < 4; c++){
			for(int r = 0; r < 4; r++){
				std::cout.width(10);
				std::cout << m_mat[r + c * 4] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

private:
	T m_mat[16];
};


#endif // __MATH_MAT_H__
