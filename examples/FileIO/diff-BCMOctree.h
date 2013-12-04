26,28c26,29
<     Z,        ///< Z(Morton)オーダリング
<     HILBERT,  ///< ヒルベルトオーダリング
<     RANDOM,   ///< ランダムシャッフル
---
>     Z,             ///< Z(Morton)オーダリング
>     HILBERT,       ///< ヒルベルトオーダリング
>     RANDOM,        ///< ランダムシャッフル
> 	PEDIGREELIST,  ///< Pedigreeリスト順(ファイルロード用)
57a59,68
> 
>   /// コンストラクタ.
>   ///
>   ///  @param[in] rootGrid ルートノード配置情報
>   ///  @param[in] pedigrees ペディグリリスト
>   ///
>   ///  @note rootGridは、デストラクタにより解放される．
>   ///
>   BCMOctree(RootGrid* rootGrid, const std::vector<Pedigree>& pedigrees);
> 
76a88,90
> 
>   /// ルートグリッドを取得
>   const RootGrid* getRootGrid() const { return rootGrid; }
