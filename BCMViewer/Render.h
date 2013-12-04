/*
 *  BCMViewer
 *
 *  Copyright 2012 SGI Japan, Ltd.
 *
 */

#ifndef __BCM_VIEWER_RENDER_H__
#define __BCM_VIEWER_RENDER_H__

class ProgramObject;

namespace SG {
	struct IndexFormat;
	struct VertexLineFormat;
	struct VertexFaceFormat;

	template<typename T> class Geometry;
} // namespace SG

namespace Render {

	class SGRender
	{
	public:
		SGRender();
		~SGRender();
		
		void RenderLine( const SG::Geometry<SG::VertexLineFormat>* geometry );
		void RenderFace( const SG::Geometry<SG::VertexFaceFormat>* geometry );

	private:
		ProgramObject *m_faceShader;
		ProgramObject *m_lineShader;

	};

} // namespace Render


#endif // _BCM_VIEWER_RENDER_H__

