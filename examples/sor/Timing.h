#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h> // gettimeofday
#include <ctime> // clock
#include <cstdio>

#ifdef TIMING
#define TimingStart(item) Timing::start(item)
#define TimingStop(item) Timing::stop(item)
#else
#define TimingStart(item)
#define TimingStop(item)
#endif


/// 時間測定識別番号.
enum TimingItem {
    TOTAL,
    INIT,
    RUN,
    SOR,
    VCUPDATE,
    BC,
    FINAL,
    N_TIMING_ITEM,  ///< アイテム数
};


/// 計算時間測定クラス.
class Timing {

    double cStart;
    double wStart;
    double cTime;
    double wTime;

public:

    /// コンストラクタ.
    Timing() { cTime = wTime = 0.0; }

    /// 測定開始.
    void start() {
        cStart = getCTime();
        wStart = getWTime();
    }

    /// 測定停止.
    void stop() {
        cTime += getCTime() - cStart;
        wTime += getWTime() - wStart;
    }

    /// 結果出力(CPU時間のみ).
    void print(const char* label) {
        printf("%s: %f sec\n", label, wTime);
    }

private:

    /// CPU時間取得.
    double getCTime() {
        return (double)clock() / CLOCKS_PER_SEC;
    }

    /// 経過時間取得.
    double getWTime() {
        struct timeval tv;
        gettimeofday(&tv, 0);
        return (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
    }

// static
private:

    /// 計算時間測定クラス配列.
    static Timing timer[];

public:

    /// アイテムを指定して測定開始.
    static void start(TimingItem item) { timer[item].start(); }

    /// アイテムを指定して測定停止.
    static void stop(TimingItem item) { timer[item].stop(); }

    /// アイテムを指定して結果出力.
    static void print(TimingItem item, const char* label) { timer[item].print(label); }
};


#endif // TIMING_H
