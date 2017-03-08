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
