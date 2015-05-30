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

#ifndef __UGL_H__
#define __UGL_H__

#include "type.h"
#include "mat.h"

#ifdef __APPLE__ 
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif



class BBox {
public:
	BBox(vec3<float> org, vec3<float> size)
	{
		this->org  = org;
		this->size = size;
	}
	
	const vec3<float> GetMin() const { return org; }
	const vec3<float> GetMax() const { return org + size; }
	const vec3<float> GetCenter() const { return org + vec3<float>(size[0] / 2.0, size[1] / 2.0, size[2] / 2.0); }
	const float GetRadius() const { return (sqrt(size[0] * size[0] + size[1] * size[1] + size[2] * size[2])) / 2.0; }

private:
	vec3<float> org;
	vec3<float> size;
};

class ViewController
{
public:
	ViewController( vec2<float>& windowSize, const BBox& bbox )
	 : m_windowSize(windowSize), 
	   m_bbox(bbox),
	   m_aspect(windowSize[0] / windowSize[1])
	{
		SetBoundingBox(m_bbox);
		SetPerspective();
	}

	~ViewController(){ }

	void SetBoundingBox( const BBox& bbox ){
		m_bbox = bbox;
		m_scale = bbox.GetRadius();
		const vec3<float> center = bbox.GetCenter();
		const vec3<float> lookat(center[0], center[1], center[2] + bbox.GetRadius() * 2.0);
		SetLookAt( lookat, center );
		SetRotateCenter( center );

		m_near = bbox.GetRadius() * 0.005f;
		m_far  = bbox.GetRadius() * 4.0f;
	}

	void ProjectionUpdate()
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		switch(m_mode){
			case PERSPECTIVE : {
				gluPerspective( m_fov, m_aspect, m_near, m_far );
				break;
			}
			case ORTHO : {
				// TODO
				//glOrtho(m_ortho[0], m_ortho[1], m_ortho[2], m_ortho[3], m_near, m_far);
				break;
			}
			default : {
				break;
			}
		}
		glMatrixMode(GL_MODELVIEW);
	}

	
	void SetPerspective(const float fov = 40.0)
	{
		m_fov = fov;
		m_mode = PERSPECTIVE;
		ProjectionUpdate();
	}

	void SetOrtho()
	{
		m_mode = ORTHO;
		ProjectionUpdate();
	}
	
	void ResizeWindow(const vec2<float>& windowSize )
	{
		m_windowSize = windowSize;
		m_aspect = m_windowSize[0] / m_windowSize[1];
		ProjectionUpdate();
	}

	void ViewUpdate()
	{
		if( m_mode == ORTHO ){
			ProjectionUpdate();
		}

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt( m_position[0], m_position[1], m_position[2], m_lookAt[0], m_lookAt[1], m_lookAt[2], 0.0, 1.0, 0.0);

		glMultMatrixf( m_translate.ptr() );

		mat4<float> trans;
		trans.at(0, 3) = m_rotateCenter[0];
		trans.at(1, 3) = m_rotateCenter[1];
		trans.at(2, 3) = m_rotateCenter[2];
		glMultMatrixf( trans.ptr() );

		glMultMatrixf( m_rotate.ptr() );

		trans.at(0, 3) = -m_rotateCenter[0];
		trans.at(1, 3) = -m_rotateCenter[1];
		trans.at(2, 3) = -m_rotateCenter[2];
		glMultMatrixf( trans.ptr() );
	}

	
	void SetLookAt(  const vec3<float>& position, const vec3<float>& lookAt = vec3<float>(0.0, 0.0, 0.0) )
	{
		m_position = position;
		m_lookAt   = lookAt;
	}

	void SetRotateCenter( const vec3<float>& center ){
		m_rotateCenter = center;
	}

	void SetTranslate( const vec3<float>& translate )
	{
		m_translate.at(0, 3) += translate[0] * m_scale;
		m_translate.at(1, 3) += translate[1] * m_scale;
		m_translate.at(2, 3) += translate[2] * m_scale;

		// TODO m_ortho
	}

	void SetRotateX( const float rad )
	{
		mat4<float> rot;
		rot.at(1, 1) =  cos(rad);
		rot.at(1, 2) = -sin(rad);
		rot.at(2, 1) =  sin(rad);
		rot.at(2, 2) =  cos(rad);

		m_rotate = m_rotate * rot;
	}

	void SetRotateY( const float rad )
	{
		mat4<float> rot;
		rot.at(0, 0) =  cos(rad);
		rot.at(0, 2) =  sin(rad);
		rot.at(2, 0) = -sin(rad);
		rot.at(2, 2) =  cos(rad);

		m_rotate = m_rotate * rot;	
	}

	void SetRotateZ( const float rad )
	{
		mat4<float> rot;
		rot.at(0, 0) =  cos(rad);
		rot.at(0, 1) = -sin(rad);
		rot.at(1, 0) =  sin(rad);
		rot.at(1, 1) =  cos(rad);

		m_rotate = m_rotate * rot;
	}

private:

	enum {
		PERSPECTIVE,
		ORTHO,
	} m_mode;

	float m_near;
	float m_far;

	vec2<float> m_windowSize;
	BBox        m_bbox;
	
	vec3<float> m_position;
	vec3<float> m_lookAt;
	
	mat4<float> m_translate;
	mat4<float> m_rotate;
	
	vec3<float> m_rotateCenter;
	float m_scale;

	// View Parameter
	float m_aspect;
	// Perspective Params
	float m_fov;
	// OrthoParams
	float m_ortho[4]; // left, right, bottom, top

};


#endif // __UGL_H__

