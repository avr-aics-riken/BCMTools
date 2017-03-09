# BCMTools


## REVISION HISTORY

---
- 2017-03-09  Version 1.0.3
  - build libBCMconfig.a

---
- 2017-03-09  Version 1.0.2
  - modify implementation of 'real_type'
  - build examples, add silo, cutlib, and hdf5
  - Change Vec3d >> Vec3r
    - include/BCMPolylib.h, BlockBase.h
    - src/BCMPolylib.cpp
    - examples/ParallelMeshGeneration/BlockBoundingBox.h, Config.h, main.cpp, SiloWriter.h
    - Utils/*.h, *.cpp except SiloWriter.h
  - double >> REAL_TYPE
    - examples/ParallelMeshGeneration/Config.h, main.cpp
    - Utils/*.h, *.cpp except SiloWriter.h
    -include/*.h
    - src/BlockManager.cpp (cout realted)

---
- 2017-03-08  Version 1.0.1
  - bug fix : src/CMakeLists.txt
  - add flag to switch float <=> double

---
- 2017-03-08  Version 1.0.0
  - cmake version
  - Tested

|Compiler|OpenMP|MPI |example|
|:--|:--:|:--:|:--:|
|Intel 17.0.1 |ok|ok|||
|GNU 6.2.0    |ok|ok|||
|fx10         |||||
|K            |||||
|fx100        |||||


---
- 2014-05-30  Version 0.9.6
  - update request pulls #6 - #12

---
- 2014-03-23  Version 0.9.5
  - update Vec3.h


---
- 2014-03-17  Version 0.9.4
  - common namespace of Vec3class


---
- 2014-01-25  Version 0.9.3
  - r10 move FileIO to HDMlib [ver 0.9.3]
  - Change year from 2013 to 2014 in copyright


---
- 2014-01-20
  - r9 add examples/ParallelMeshGeneration


---
- 2013-12-16  Version 0.9.2
  - r8 add Utils/include/*.h to install-dir/include
  - remove copyright from ConfigFile.h/cpp

---
- 2013-12-15  Version 0.9.1
  - r7 build only libBCM.a

---
- 2013-12-04
  - r6 BCMTools 20130709
    - insert copyright

  - r5 BCMTools 20130417

  - r4 BCMTools 20120916

  - r3 BCMTools 20120607

  - r2 add BCMviewer

  - r1 Initial commit


---
- 2013-07-09 (jonishi@iis.u-tokyo.ac.jp)
  - BlockManager.cpp におけるループのスレッド並列化
  - 各ブロックに割当てるメモリを連続にするデータクラスの生成・登録の実装
  - Vec3の曖昧さ回避
  - make_setting.{fx,gnu,ibm.intel} -> Makefile.inc への変更
  - SolverTest のアップデート
  - examples の各Makefileを修正
  - examples/FileIO -> FileIO への移動
  - examples/util -> Utils への移動
  - その他（バグなど）


---
- 2013-04-17 (jonishi@iis.u-tokyo.ac.jp)
  - FileIO の追加
  - SolverTest の追加
  - Makefile.inc -> make_setting.{fx,gnu,ibm,intel} への変更
  - DataClass.list の廃止
  - データ同期方法の追加
  - その他（バグなど）
