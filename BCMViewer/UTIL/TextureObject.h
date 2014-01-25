/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef GL_TEXTUREOBJECT_H
#define GL_TEXTUREOBJECT_H

/*! @brief TextureObjectクラス
 *
 * @note OpenGL Texture. テクスチャへのイメージ書き込みが行える
 */
class TextureObject
{
public:
	/*! @brief コンストラクタ
	 */
	TextureObject();

	/*! @brief サイズ指定付コンストラクタ
	 *  @param[in] width 幅
	 *  @param[in] height 高さ
	 *  @param[in] colorComponent 色タイプ
	 *  @param[in] colorBit 色ビット数
	 *  @param[in] depthBit デプスビット数
	 */
	TextureObject(int width, int height, int colorComponent, int colorBit = 32, int depthBit = -1);

	/*! @brief デストラクタ
	 */
	~TextureObject();

	/*! @brief テクスチャバインド
	 */
	void BindTexture();

	/*! @brief テクスチャアンバインド
	 */
	void UnbindTexture();

	/*! @brief テクスチャidの取得
	 *  @return テクスチャid
	 */
	unsigned int GetTexture() const;

	/*! @brief テクスチャの幅取得
	 *  @return 幅
	 */
	int Width() const;

	/*! @brief テクスチャの高さ取得
	 *  @return 高さ
	 */
	int Height() const;

	/*! @brief テクスチャの大きさを変更
	 *  @param[in] w 幅
	 *  @param[in] h 高さ
	 *  @param[in] colorComponent 色タイプ
	 *  @param[in] colorBit 色ビット数
	 *  @param[in] depthBit デプスビット数
	 *  @return void
	 */
	void Resize(int w, int h, int colorComponent, int colorBit = 32, int depthBit = -1);

	/*! @brief テクスチャの大きさのみを変更
	 *  @param[in] w 幅
	 *  @param[in] h 高さ
	 *  @return void
	 */
	void Resize(int w, int h);

	/*! @brief テクスチャへのイメージ書き込み
	 *  @param[in] pImage イメージのポインタ
	 *  @return 成功した場合true, 失敗した場合false
	 */
	bool WriteImage(const float* pImage);

	/*! @brief テクスチャへのイメージ書き込み
	 *  @param[in] pImage イメージのポインタ
	 *  @return 成功した場合true, 失敗した場合false
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

