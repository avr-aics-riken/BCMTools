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

#ifndef __BCM_VIEWER_SCENE_GRAPH_H__
#define __BCM_VIEWER_SCENE_GRAPH_H__

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


class NonCopyable
{
protected:
	NonCopyable() { }
	~NonCopyable() { }
private:
	NonCopyable( const NonCopyable&);
	NonCopyable& operator = (const NonCopyable&);
};

namespace SG {

	struct IndexFormat
	{
		GLuint i;
		IndexFormat() : i(0) { }
		IndexFormat(GLuint _i) : i(_i) { }
	};

	struct VertexLineFormat
	{
		GLfloat x, y, z; 
		VertexLineFormat() : x(.0f), y(.0f), z(.0f) {}
		VertexLineFormat(GLfloat _x, GLfloat _y, GLfloat _z) : x(_x), y(_y), z(_z){}
	};

	struct VertexFaceFormat
	{
		GLfloat vx, vy, vz;
		GLfloat nx, ny, nz;
		GLfloat tx, ty;
		VertexFaceFormat() : vx(.0f), vy(.0f), vz(.0f), nx(.0f), ny(.0f), nz(.0f), tx(.0f), ty(.0f) {}
		VertexFaceFormat(GLfloat _vx, GLfloat _vy, GLfloat _vz, GLfloat _nx, GLfloat _ny, GLfloat _nz, GLfloat _tx, GLfloat _ty) 
		  : vx(_vx), vy(_vy), vz(_vz), nx(_nx), ny(_ny), nz(_nz), tx(_tx), ty(_ty) {}

		static const GLubyte* GetNormalOffset()   { return (GLubyte*)NULL + sizeof(GLfloat) * 3; }
		static const GLubyte* GetTexCoordOffset() { return (GLubyte*)NULL + sizeof(GLfloat) * 6; }
	};

	template<typename VertexFormat> class Geometry : private NonCopyable
	{
	public:

		Geometry( size_t maxVertices, size_t maxIndices ) : m_maxVertices(maxVertices), m_maxIndices(maxIndices)
		{

			m_v = new VertexFormat[m_maxVertices];
			m_i = new IndexFormat[m_maxIndices];

			size_t numBuf = sizeof(VertexFormat) * m_maxVertices;

			glGenBuffers(2, m_vbo);
			glBindBuffer(GL_ARRAY_BUFFER, GetVertexVBO());
			glBufferData(GL_ARRAY_BUFFER, numBuf, NULL, GL_STATIC_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, GetIndexVBO());
			glBufferData(GL_ARRAY_BUFFER, sizeof(IndexFormat) * m_maxIndices, NULL, GL_STATIC_DRAW);

			m_numIndices   = 0;
			m_numVertices  = 0;
			SetVisible(true);
		}

		virtual ~Geometry()
		{ 
			glDeleteBuffers(2, m_vbo);
			delete [] m_v;
			delete [] m_i;
		}

		const GLsizei GetStride() const { return sizeof(VertexFormat); }

		const GLuint GetVertexVBO() const { return m_vbo[0]; }
		const GLuint GetIndexVBO()  const { return m_vbo[1]; }

		const size_t GetNumIndices()   const { return m_numIndices;   }
		const size_t GetNumVertices()  const { return m_numVertices;  }

		void SetVisible(const bool visible) { m_visible = visible; }
		bool GetVisible() const { return m_visible; }

		bool SetVertices( const VertexFormat* vert, const size_t num ){
			if( num > m_maxVertices ){
				printf("num : %ld m_maxVertices :%ld\n", num, m_maxVertices);
				fprintf(stderr, "[ERR] %s [%s:%d]\n", __func__, __FILE__, __LINE__);
				return false;
			}
			memcpy(m_v, vert, sizeof(VertexFormat) * num);
			UploadVertices(num);
		}

		bool SetIndices( const IndexFormat *idx, const size_t num ){
			if( num > m_maxIndices ){
				printf("num : %ld m_maxIndices :%ld\n", num, m_maxIndices);
				fprintf(stderr, "[ERR] %s [%s:%d]\n", __func__, __FILE__, __LINE__);
				return false;
			}
			memcpy(m_i, idx, sizeof(IndexFormat) * num);
			UploadIndices(num);
		}
		
		VertexFormat* GetVertexBuffer(){ return m_v; }
		IndexFormat*  GetIndexBuffer() { return m_i; }

		void UploadVertices(size_t numVertices) {
			m_numVertices = numVertices;
			glBindBuffer(GL_ARRAY_BUFFER, GetVertexVBO());
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(VertexFormat) * m_numVertices, m_v);
		}

		void UploadIndices(size_t numIndices) {
			m_numIndices = numIndices;
			glBindBuffer(GL_ARRAY_BUFFER, GetIndexVBO());
			glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(IndexFormat) * m_numIndices, m_i);

			//for(int i = 0; i < m_numIndices; i++){
			//	printf("%d\n", m_i[i]);
			//}
		}

	private:
		const size_t m_maxVertices;
		const size_t m_maxIndices;
		size_t m_numVertices;
		size_t m_numIndices;

		IndexFormat  *m_i;
		VertexFormat *m_v;

		bool m_visible;
		GLuint m_vbo[2];
	};

} // namespace SG


#endif // __BCM_VIEWER_SCENE_GRAPH_H__ 
