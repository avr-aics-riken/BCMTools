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

#ifndef GL_SHADERPROGRAMOBJECT_H
#define GL_SHADERPROGRAMOBJECT_H

#include <string>

/*! @brief ShaderObject�N���X
 *
 * @note GLSL�pShaderObject�N���X.VertexShader, FragmentShader����
 * GLSL�v���O�����̃��[�h�ƃR���p�C�����s��
 */
class ShaderObject
{
public:
	enum SHADERTYPE
	{
		FRAGMENT_SHADER = 0x8B30, // GL_FRAGMENT_SHADER
		VERTEX_SHADER   = 0x8B31  // GL_VERTEX_SHADER
	};

	/*! @brief �R���X�g���N�^
	 */
	ShaderObject();

	/*! @brief �f�X�g���N�^
	 */
	~ShaderObject();

	/*! @brief �\�[�X�t�@�C�����[�h����
	 *  @param[in] filename �t�@�C����
	 *  @param[in] shaderType �V�F�[�_�^�C�v
	 *  @return ���������ꍇtrue, ���s�����ꍇfalse
	 */
	bool LoadFromFile(const std::string& filename, SHADERTYPE shaderType);

	/*! @brief �\�[�X�t�@�C�����[�h����
	 *  @param[in] programSource �\�[�X�v���O����
	 *  @param[in] shaderType �V�F�[�_�^�C�v
	 *  @return ���������ꍇtrue, ���s�����ꍇfalse
	 */
	bool LoadFromMemory(const std::string& programSource, SHADERTYPE shaderType);

	/*! @brief ��������
	 *  @return void
	 */
	void Release();

	/*! @brief �V�F�[�_id�擾
	 *  @return �V�F�[�_id
	 */
	unsigned int GetShader() const
	{
		return m_shader;
	}

private:
	unsigned int m_shader;
	std::string	m_source;
};

/*! @brief ProgramObject�N���X
 *
 * @note GLSL�pProgramObject�N���X.
 * ShaderObject�N���X�𗘗p��GLSL�v���O�����������N����
 */

class ProgramObject
{
public:

	/*! @brief �R���X�g���N�^
	 */
	ProgramObject();

	/*! @brief �R���X�g���N�^�i�����N�@�\���j
	 *  @param[in] vertexShader vertex�v���O����
	 *  @param[in] fragmentShader fragment�v���O����
	 */
	ProgramObject(const ShaderObject& vertexShader, const ShaderObject& fragmentShader);

	/*! @brief ShaderObject�̃����N
	 *  @param[in] vertexShader vertex�v���O����
	 *  @param[in] fragmentShader fragment�v���O����
	 *  @return ���������ꍇtrue, ���s�����ꍇfalse
	 */
	bool Link(const ShaderObject& vertexShader, const ShaderObject& fragmentShader);

	/*! @brief ProgramObject�̃o�C���h
	 *  @return void
	 */
	void Bind();

	/*! @brief ProgramObject�̃A���o�C���h
	 *  @return void
	 */
	void Unbind();

	/*! @brief ProgramObject�̉���
	 *  @return void
	 */
	void Release();

	/*! @brief Uniform�ϐ��̐ݒ�
	 *  @return void
	 */
	// 1i - 4i
	void SetUniform(const char* name, const int i0);
	void SetUniform(const char* name, const int i0, const int i1);
	void SetUniform(const char* name, const int i0, const int i1, const int i2);
	void SetUniform(const char* name, const int i0, const int i1, const int i2, const int i3);
	void SetUniform(const char* name, const int num, const int* i_array);
	// 1f - 4f
	void SetUniform(const char* name, const float f0);
	void SetUniform(const char* name, const float f0, const float f1);
	void SetUniform(const char* name, const float f0, const float f1, const float f2);
	void SetUniform(const char* name, const float f0, const float f1, const float f2, const float f3);
	void SetUniform(const char* name, const int num, const float* f_array);
	// matrix
	void SetUniformMatrix2x2(const char* name, const int count, const bool transpose, const float* val);
	void SetUniformMatrix3x3(const char* name, const int count, const bool transpose, const float* val);
	void SetUniformMatrix4x4(const char* name, const int count, const bool transpose, const float* val);


private:
	unsigned int m_program;
	int m_oldProgram;
	bool m_binding;

};


#endif // GL_SHADERPROGRAMOBJECT_H
