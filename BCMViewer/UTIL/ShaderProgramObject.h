/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef GL_SHADERPROGRAMOBJECT_H
#define GL_SHADERPROGRAMOBJECT_H

#include <string>

/*! @brief ShaderObjectクラス
 *
 * @note GLSL用ShaderObjectクラス.VertexShader, FragmentShader共通
 * GLSLプログラムのロードとコンパイルを行う
 */
class ShaderObject
{
public:
	enum SHADERTYPE
	{
		FRAGMENT_SHADER = 0x8B30, // GL_FRAGMENT_SHADER
		VERTEX_SHADER   = 0x8B31  // GL_VERTEX_SHADER
	};

	/*! @brief コンストラクタ
	 */
	ShaderObject();

	/*! @brief デストラクタ
	 */
	~ShaderObject();

	/*! @brief ソースファイルロード処理
	 *  @param[in] filename ファイル名
	 *  @param[in] shaderType シェーダタイプ
	 *  @return 成功した場合true, 失敗した場合false
	 */
	bool LoadFromFile(const std::string& filename, SHADERTYPE shaderType);

	/*! @brief ソースファイルロード処理
	 *  @param[in] programSource ソースプログラム
	 *  @param[in] shaderType シェーダタイプ
	 *  @return 成功した場合true, 失敗した場合false
	 */
	bool LoadFromMemory(const std::string& programSource, SHADERTYPE shaderType);

	/*! @brief 解放処理
	 *  @return void
	 */
	void Release();

	/*! @brief シェーダid取得
	 *  @return シェーダid
	 */
	unsigned int GetShader() const
	{
		return m_shader;
	}

private:
	unsigned int m_shader;
	std::string	m_source;
};

/*! @brief ProgramObjectクラス
 *
 * @note GLSL用ProgramObjectクラス.
 * ShaderObjectクラスを利用しGLSLプログラムをリンクする
 */

class ProgramObject
{
public:

	/*! @brief コンストラクタ
	 */
	ProgramObject();

	/*! @brief コンストラクタ（リンク機能つき）
	 *  @param[in] vertexShader vertexプログラム
	 *  @param[in] fragmentShader fragmentプログラム
	 */
	ProgramObject(const ShaderObject& vertexShader, const ShaderObject& fragmentShader);

	/*! @brief ShaderObjectのリンク
	 *  @param[in] vertexShader vertexプログラム
	 *  @param[in] fragmentShader fragmentプログラム
	 *  @return 成功した場合true, 失敗した場合false
	 */
	bool Link(const ShaderObject& vertexShader, const ShaderObject& fragmentShader);

	/*! @brief ProgramObjectのバインド
	 *  @return void
	 */
	void Bind();

	/*! @brief ProgramObjectのアンバインド
	 *  @return void
	 */
	void Unbind();

	/*! @brief ProgramObjectの解放
	 *  @return void
	 */
	void Release();

	/*! @brief Uniform変数の設定
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

