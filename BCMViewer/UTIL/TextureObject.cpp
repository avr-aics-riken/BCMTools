/*
###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

#include "TextureObject.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include <stdio.h>
#include <string.h>


namespace {

	void getFormat(int color_bit, int color_component,
				   unsigned int& gltexformat, unsigned int& internalformat, unsigned int& internalsize)
	{
		gltexformat    = GL_RGBA32F_ARB;
		internalformat = GL_RGB;
		internalsize   = GL_UNSIGNED_BYTE;
		if (color_bit == 128)
		{
			gltexformat    = GL_RGBA32F_ARB;
			internalformat = GL_RGBA;
			internalsize   = GL_FLOAT;
		}
		else if (color_bit == 32)
		{
			if (color_component == 1)
			{
				gltexformat = GL_LUMINANCE32F_ARB;
				internalformat = GL_LUMINANCE;
				internalsize   = GL_FLOAT;
			}
			else
			{
				gltexformat = GL_RGBA;
				internalformat = GL_BGRA;              //BGRA
				internalsize   = GL_UNSIGNED_BYTE;
			}
		}
		else if (color_bit == 24)
		{
			gltexformat = GL_RGB;
			internalformat = GL_RGB;
			internalsize   = GL_UNSIGNED_BYTE;
		}
		else if(color_bit == 8)
		{
			gltexformat    = GL_LUMINANCE;
			internalformat = GL_LUMINANCE;
			internalsize    = GL_UNSIGNED_BYTE;
		}
	}

} // namespace


TextureObject::TextureObject() :
	m_width(0), m_height(0), m_color_component(3), m_color_bit(32), m_depth_bit(-1), m_texname(0)
{
}

TextureObject::TextureObject(int width, int height, int colorComponent, int colorBit, int depthBit) :
	m_width(width), m_height(height), m_color_component(colorComponent), m_color_bit(colorBit), m_depth_bit(depthBit), m_texname(0)
{
	create();
}

TextureObject::~TextureObject()
{
	release();
}

void TextureObject::create()
{
	if (!m_texname)
	{
		GLuint save_texname;
		glGetIntegerv(GL_TEXTURE_BINDING_2D, &(GLint &)save_texname);
		glGenTextures(1, &m_texname);
		glBindTexture(GL_TEXTURE_2D, m_texname);

		unsigned int gltexformat, internalformat, internalsize;
		getFormat(m_color_bit, m_color_component,
				  gltexformat, internalformat, internalsize);

		glTexImage2D(
			GL_TEXTURE_2D, 0,
			gltexformat,
			m_width, m_height, 0,
			internalformat,
			internalsize, 0
		);
		GLenum err = glGetError();
		if (err != GL_NO_ERROR)
			printf("[%s:%d]CreateTexture error=%d\n",__FILE__, __LINE__, err);

		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, save_texname);
	}
}

void TextureObject::BindTexture()
{
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &(GLint &)m_savedTex);
	glBindTexture(GL_TEXTURE_2D, m_texname);
}

void TextureObject::UnbindTexture()
{
	glBindTexture(GL_TEXTURE_2D, m_savedTex);
}

unsigned int  TextureObject::GetTexture() const
{
	return m_texname;
}

void TextureObject::release()
{
	if (m_texname) {
		//printf("TextureObject::release texture %d\n", m_texname);
		glDeleteTextures(1, &m_texname);
		m_texname = 0;
	}
}

int TextureObject::Width() const
{
	return m_width;
}

int TextureObject::Height() const
{
	return m_height;
}

void TextureObject::Resize(int w, int h, int colorComponent, int colorBit, int depthBit)
{
	if (w != m_width || h != m_height)
	{
		release();
		m_width      = w;
		m_height     = h;
		m_color_component = colorComponent;
		m_color_bit  = colorBit;
		m_depth_bit  = depthBit;
		create();
	}
}
void TextureObject::Resize(int w, int h)
{
	Resize(w, h, m_color_component, m_color_bit, m_depth_bit);
}

bool TextureObject::WriteImage(const float* pImage)
{
	if (m_texname)
	{
		GLuint save_texname;
		glGetIntegerv(GL_TEXTURE_BINDING_2D, &(GLint &)save_texname);
		//glGenTextures(1, &m_texname);
		glBindTexture(GL_TEXTURE_2D, m_texname);
		if (m_color_bit == 128
		|| (m_color_bit == 32 && m_color_component == 1))
		{
			unsigned int gltexformat, internalformat, internalsize;
			getFormat(m_color_bit, m_color_component,
					  gltexformat, internalformat, internalsize);

			glTexImage2D(
				GL_TEXTURE_2D, 0,
				gltexformat,
				m_width, m_height, 0,
				internalformat,
				internalsize, pImage
			);
		}
		else
		{
			return false;
		}
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, save_texname);


		return true;
	}
	return false;
}


bool TextureObject::WriteImage(const int* pImage)
{
	if (!m_texname)
		return false;

	GLuint save_texname;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &(GLint &)save_texname);
	//glGenTextures(1, &m_texname);
	glBindTexture(GL_TEXTURE_2D, m_texname);
	if (m_color_bit == 32 && m_color_component == 4)
	{
		unsigned int gltexformat, internalformat, internalsize;
		getFormat(m_color_bit, m_color_component,
				  gltexformat, internalformat, internalsize);

		glTexImage2D(GL_TEXTURE_2D, 0,
					 gltexformat,
					 m_width, m_height, 0,
					 internalformat,
					 internalsize, pImage
					 );
		GLenum err = glGetError();
		err = err;
	}
	else
	{
		return false;
	}
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, save_texname);

	return true;
}


bool TextureObject::WriteImage(const int* pImage, const int width, const int height)
{
	if (!m_texname)
		return false;

	GLuint save_texname;
	glGetIntegerv(GL_TEXTURE_BINDING_2D, &(GLint &)save_texname);
	//glGenTextures(1, &m_texname);
	glBindTexture(GL_TEXTURE_2D, m_texname);
	if (m_color_bit == 32 && m_color_component == 4)
	{
		unsigned int gltexformat, internalformat, internalsize;
		getFormat(m_color_bit, m_color_component,
				  gltexformat, internalformat, internalsize);

		/*
		glTexImage2D(GL_TEXTURE_2D, 0,
					 gltexformat,
					 m_width, m_height, 0,
					 internalformat,
					 internalsize, pImage
					 );
		*/
		glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, internalformat, internalsize, pImage);
		GLenum err = glGetError();
		err = err;
	}
	else
	{
		return false;
	}
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, save_texname);

	return true;
}
