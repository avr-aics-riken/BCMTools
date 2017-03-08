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

#ifndef GL_TEXTUREOBJECT_H
#define GL_TEXTUREOBJECT_H

/*! @brief TextureObject�N���X
 *
 * @note OpenGL Texture. �e�N�X�`���ւ̃C���[�W�������݂��s����
 */
class TextureObject
{
public:
	/*! @brief �R���X�g���N�^
	 */
	TextureObject();

	/*! @brief �T�C�Y�w���t�R���X�g���N�^
	 *  @param[in] width ��
	 *  @param[in] height ����
	 *  @param[in] colorComponent �F�^�C�v
	 *  @param[in] colorBit �F�r�b�g��
	 *  @param[in] depthBit �f�v�X�r�b�g��
	 */
	TextureObject(int width, int height, int colorComponent, int colorBit = 32, int depthBit = -1);

	/*! @brief �f�X�g���N�^
	 */
	~TextureObject();

	/*! @brief �e�N�X�`���o�C���h
	 */
	void BindTexture();

	/*! @brief �e�N�X�`���A���o�C���h
	 */
	void UnbindTexture();

	/*! @brief �e�N�X�`��id�̎擾
	 *  @return �e�N�X�`��id
	 */
	unsigned int GetTexture() const;

	/*! @brief �e�N�X�`���̕��擾
	 *  @return ��
	 */
	int Width() const;

	/*! @brief �e�N�X�`���̍����擾
	 *  @return ����
	 */
	int Height() const;

	/*! @brief �e�N�X�`���̑傫�����ύX
	 *  @param[in] w ��
	 *  @param[in] h ����
	 *  @param[in] colorComponent �F�^�C�v
	 *  @param[in] colorBit �F�r�b�g��
	 *  @param[in] depthBit �f�v�X�r�b�g��
	 *  @return void
	 */
	void Resize(int w, int h, int colorComponent, int colorBit = 32, int depthBit = -1);

	/*! @brief �e�N�X�`���̑傫���݂̂��ύX
	 *  @param[in] w ��
	 *  @param[in] h ����
	 *  @return void
	 */
	void Resize(int w, int h);

	/*! @brief �e�N�X�`���ւ̃C���[�W��������
	 *  @param[in] pImage �C���[�W�̃|�C���^
	 *  @return ���������ꍇtrue, ���s�����ꍇfalse
	 */
	bool WriteImage(const float* pImage);

	/*! @brief �e�N�X�`���ւ̃C���[�W��������
	 *  @param[in] pImage �C���[�W�̃|�C���^
	 *  @return ���������ꍇtrue, ���s�����ꍇfalse
	 */
	bool WriteImage(const int* pImage);

	bool WriteImage(const int*pImage, const int wigth, const int height);

private:
	unsigned int m_texname;
	unsigned int m_savedTex;
	int m_width;
	int m_height;
	int m_color_component;
	int m_color_bit;
	int m_depth_bit;
	void release();
	void create();
};

#endif // GL_TEXTUREOBJECT_H
