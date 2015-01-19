/******************************************************************************/
/*! @file       MyDataFunction.h
    @brief      一次元波データに対する処理をまとめたクラス、定義群
*******************************************************************************
    ファイルの詳細な説明
*******************************************************************************
    @date       作成日(2014/10/24)
    @author		YumaNagai
    @par        History
    - 20XX/10/24 YumaNagai
      -# Initial Version
	  -# 作成
******************************************************************************/


#ifndef _MYDATAFUNCTION_H
#define _MYDATAFUNCTION_H

#define _USE_MATH_DEFINES	// for C++
#include <math.h>
#include <iostream>
#include <algorithm>		// swap
#include <complex>
#include <vector>

//typedef宣言
typedef std::complex<double> Complex;
typedef unsigned int uint;

#pragma region 列挙体定義

/*!-----------------------------
	@enum 窓関数列挙体
	@brief 窓関数の種類
---------------------------------*/
enum WaveBandType
{
	Delta = 0,		///< @brief δ波(0.1～3.0Hz)
	Theta,			///< @brief θ波(4.0～7.0Hz)
	Alpha,			///< @brief α波(8.0～12.0Hz)
	Beta,			///< @brief β波(13.0～30.0Hz)
	Gamma,			///< @brief γ波(30.0～100.0Hz)
	WaveBandMax,	///< @brief 番兵
};

/*!-----------------------------
	@enum 窓関数列挙体
	@brief 窓関数の種類
---------------------------------*/
enum WindowType
{
	Rectangle = 0,	///< @brief 矩形窓
	Hanning,		///< @brief ハニング窓
	Hamming,		///< @brief ハミング窓
	Blackman,		///< @brief ブラックマン窓
	WinTypeMax		///< @brief 番兵
};

#pragma endregion

#pragma region 高速フーリエ変換を行うクラス

/*! 
 @class FftTransration
 @brief 高速フーリエ変換を行うクラス
*/
class MyFourierTransform
{
public:


private:
	typedef unsigned int uint;

	//メンバー変数
	static Complex one;		///<@brief 要素(1.0,0.0)の複素数
	static Complex ione;		///<@brief 要素(0.0,1.0)の複素数

	Complex *fftArray;			///<@brief 複素数の配列
	size_t size;							///<@brief fftArrayのサイズ
	size_t use;								///<@brief コンストラクタで与えたサイズ
	bool verbose;							///<@brief 配列ダンプのモード指定に使う

	//メンバー関数
	
	/*! @brief ビット判定を行う関数  */
	void BitReverse();

	/*!--------------------------------------------- 
		@brief s以上の最小の2のべき乗を返す関数
		@param [in] s サンプル長
		@return s以上の最小の2のべき乗
	-------------------------------------------------*/
	size_t NextPow2(size_t s);
	
	
	/*!-------------------------------------------------
		@brief コンストラクタ
		@param s [in] 信号の長さ(サンプル長) 
	---------------------------------------------------*/
	MyFourierTransform(size_t s) : use(s), size(NextPow2(s)), verbose(false)
	{
		one = std::complex<double>(1.0,0.0);

		fftArray = new Complex[size];
		for(int ii=0;ii<size;ii++) fftArray[ii] = std::complex<double>(0.0,0.0);
	}

	/*!-------------------------------------------------
		@brief 波形データを配列にセットする関数
		@param [in] wavedata 波データの配列(窓関数をかけた波)
		@param [in] 窓関数の種類
	----------------------------------------------------*/
	void Setwavedata(double *wavedata,WindowType window);
	
public:
	
	/*!-----------------------------------------------------
		@brief 生成メソッド
		@param [in] sample 信号の長さ(サンプル長) 
		@param [in] wavedata 波データの配列(窓関数をかけた波)
		@param [in] window 窓関数の種類
		@return 生成初期化したFFTTransformクラスの実体
	---------------------------------------------------------*/
	static MyFourierTransform* Load(size_t sample,double *wavedata,WindowType window)
	{
		MyFourierTransform* ffttrans = new MyFourierTransform(sample);	//データ生成
		ffttrans->Setwavedata(wavedata,window);				//波データセット
		ffttrans->FFT();									//高速フーリエ変換を行う
		return ffttrans;									//高速フーリエ変換を行ったクラスデータを返す
	}

	///< @brief デストラクタ
	~MyFourierTransform()
	{
		delete[] fftArray;
	};

	///< @brief 逆変換を行うかどうかのフラグ
	///< @param[in] b true:逆フーリエ変換モード false:通常拘束フーリエ変換モード
	bool setVerbose(bool b){  return verbose = b; };
  
	void Dump();									///< @brief データダンプ(データ表示)
	void FFT(bool isReverse=false);					///< @brief 高速フーリエ変換
	void IFFT();									///< @brief 逆高速フーリエ変換
	size_t GetfftArraySize(){ return size;}			///< @brief FFTで使用する配列サイズの取得

	
	///< @brief 配列要素にアクセスするための演算子
	///< @param[in] index 配列添字番号
	///< @return 添字番号の配列要素
	const Complex& operator[](int index)const
	{
		if (index < 0) return fftArray[use - index];
		return fftArray[index];
	};

	///< @brief 配列のアドレスを返す関数
	///< @param index 配列添字番号
	///< @return 添字番号の配列アドレス
	const  Complex* ReturnPtr(int index)
	{
		if (index < 0) return (fftArray - (use - index));
		return (fftArray+index);
	}
};

#pragma endregion

#pragma region ファジー解析クラス

///< @struct MeberSipFunctionData
///< @brief  メンバーシップ関数の値を保存しておく構造体
struct MeberSipFunctionData
{
public:

	#pragma region メンバー変数

	double sstart;		//Smallスタート位置
	double send;		//Small終了位置
	double lstart;		//Largeスタート位置
	double lend;		//Largeスタート位置

	double sdegree;		///<@brief Smallの時の傾き
	double sintercept;	///<@brief Smallの時の切片
	double ldegree;		///<@brief Largeの時の傾き
	double lintercept;	///<@brief Largeの時の切片
	
	#pragma endregion

	#pragma region 関数群

	///<@brief コンストラクタ(初期化)
	MeberSipFunctionData()
	{
		sstart = send = lstart = lend = 0;
		sdegree = sintercept = ldegree = lintercept = 0.0;
	}

	///<@brief メンバーシップ関数作成メソッド
	///<@param[in] ls Smallスタート位置
	///<@param[in] le Small終了位置
	///<@param[in] hs Largeスタート位置
	///<@param[in] he Large終了位置
	///<@param[in] lh Large or Small
	void makeMenberSipFunction(double ss,double se,double ls,double le)
	{
		sstart = ss;
		send = se;
		lstart = ls;
		lend = le;

		//式 
		//(x1,y1) (x2,y2) から傾きを求める a=(y2-y1)/(x2-x1) 
		// b = y1 - a * x1

		//'Small'の時
		sdegree =  -1.0 / (send - sstart);
		sintercept = 1.0 - sdegree * sstart;		

		//'Large'の時
		ldegree = 1.0 / (lend - lstart);			
		lintercept = -ldegree * lstart;	
	}

	///<@brief 適合度を計算し返す関数
	///<@@aram[in] value 調べる値(包含率)
	///<@param[in] lh 'S'か'L'の文字列
	double fuzzyset(double value,char sl)
	{	
		double grade = 0.0;

		//'Small'の時
		if(sl=='S')
		{
			// sstart以下の時
			if(sstart > value) grade = 1.0;

			// sstart <= value <= send
			else if(sstart <= value && send >= value) grade = (sdegree * value) + sintercept;

			// send以上の時
			else if(send < value) grade = 0.0;
		}
	
		//Large
		else
		{
			// lstart以下の時
			if(lstart > value) grade = 0.0;

			// lstart <= value <= lend
			else if(lstart <= value && lend >= value) grade = (ldegree * value) + lintercept;
			
			// lend以上の時
			else if(lend < value) grade = 1.0;
		}
		return grade;
	}

	#pragma endregion

};

///< @struct FuzzyrlueTable
///< @brief  1組のファジールールを表す構造体
struct FuzzyRlue
{
public:
	char deltaState;
	char thetaState;
	char alphaState;
	char betaState;
	char gammaState;
	char output;
	double noWeight;
	double yesWeight;

	//コンストラクタ
	FuzzyRlue()
	{
		deltaState = ' ';
		thetaState = ' ';
		alphaState = ' ';
		betaState = ' ';
		gammaState = ' ';
		output = ' '; 
		noWeight = yesWeight = 0;
	}
};

#pragma endregion

#pragma region ヘルパーメソッド群

/*!-------------------------------------------------------------
	@brief  各波の包含率を取得するメソッド(0.1～100Hzまでを対象)
	@param [in] FFT解析クラスポインタ
	@return 指定した波周波数帯の包含率
---------------------------------------------------------------*/
double getWaveProbability(MyFourierTransform* fft,WaveBandType wbtype);

/*!------------------------------------------------------------------------
	@brief 文字列分割関数
	@param[in]  str		走査文字列
	@param[in]	delim	分割する文字
	@param[out] outlist	分割した文字列の位置を示したポインタの配列
	@param[out] count	1行で分割した文字列の個数
	@return		1行で分割した文字列の個数
----------------------------------------------------------------------------*/
int split(char *str,const char *delim,char *outlist[],int &count);

#pragma endregion

#endif