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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include "SceneGraph.h"
#include "Render.h"
#include "ShaderProgramObject.h"
#include "TextureObject.h"

namespace {
	static const char* fsVShader = 
		"void main(void)"
		"{"
		"	gl_Position    = ftransform();"
		"	gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;"
		"}"
		;

	static const char* fsFShader =
		"uniform sampler2D GridTex;"
		""
		"void main(void)"
		"{"
		"	gl_FragColor = texture2D(GridTex, gl_TexCoord[0].xy);"
		"}"
		;

	static const char* lsVShader =
		"void main(void)"
		"{"
		"	gl_FrontColor = vec4(0.0, 0.0, 0.0, 1.0);"
		"	gl_Position = ftransform();"
		"}"
		;
	
	static const char* lsFShader =
		"void main(void)"
		"{"
		"	gl_FragColor = gl_Color;"
		"	gl_FragDepth = gl_FragCoord.z * 0.99999;"
		"}"
		;

};

namespace Render {

	SGRender::SGRender()
	  : m_faceShader(NULL), m_lineShader(NULL)
	{
		m_faceShader = new ProgramObject();

		ShaderObject fsvs, fsfs;
		if( !fsvs.LoadFromMemory(fsVShader, ShaderObject::VERTEX_SHADER) ){
			printf("VS ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}
		if( !fsfs.LoadFromMemory(fsFShader, ShaderObject::FRAGMENT_SHADER) ){
			printf("FS ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}

		if( !m_faceShader->Link(fsvs, fsfs) ){
			printf("LINK ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}

		m_lineShader = new ProgramObject();

		ShaderObject lsvs, lsfs;
		if( !lsvs.LoadFromMemory(lsVShader, ShaderObject::VERTEX_SHADER) ){
			printf("VS ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}
		if( !lsfs.LoadFromMemory(lsFShader, ShaderObject::FRAGMENT_SHADER) ){
			printf("FS ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}
		if( !m_lineShader->Link(lsvs, lsfs) ){
			printf("LINK ERR [%s %s:%d]\n", __func__, __FILE__, __LINE__);
		}

	}

	SGRender::~SGRender()
	{
		if(m_faceShader){
			m_faceShader->Release();
			delete m_faceShader;
		}
	}
	
	void SGRender::RenderLine( const SG::Geometry<SG::VertexLineFormat>* geometry )
	{
		m_lineShader->Bind();

		glBindBuffer(GL_ARRAY_BUFFER, geometry->GetVertexVBO());
		glVertexPointer(3, GL_FLOAT, geometry->GetStride(), 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->GetIndexVBO());
		glDrawElements(GL_LINES, geometry->GetNumIndices(), GL_UNSIGNED_INT, 0);

		m_lineShader->Unbind();
	}

	void SGRender::RenderFace( const SG::Geometry<SG::VertexFaceFormat>* geometry )
	{
		m_faceShader->Bind();
		m_faceShader->SetUniform("GridTex", 0);

		glBindBuffer(GL_ARRAY_BUFFER, geometry->GetVertexVBO());
		glVertexPointer(3, GL_FLOAT, geometry->GetStride(), 0);
		glNormalPointer(GL_FLOAT, geometry->GetStride(), SG::VertexFaceFormat::GetNormalOffset());
		glTexCoordPointer(2, GL_FLOAT, geometry->GetStride(), SG::VertexFaceFormat::GetTexCoordOffset());

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, geometry->GetIndexVBO());
		glDrawElements(GL_TRIANGLES, geometry->GetNumIndices(), GL_UNSIGNED_INT, 0);

		m_faceShader->Unbind();
	}

} // namespace Render
