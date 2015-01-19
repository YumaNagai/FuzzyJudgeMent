#include "MyDataFunction.h"

#pragma region 高速フーリエ変換クラスMyFourierTransformに関する処理

//static要素変数の初期化
Complex MyFourierTransform::one = std::complex<double>(1.0,0.0);
Complex MyFourierTransform::ione = std::complex<double>(0.0,1.0);

#pragma region ヘルパーメソッド

/*!--------------------------------------------- 
	@brief s以上の最小の2のべき乗を返す関数
	@param [in] s サンプル長
	@return s以上の最小の2のべき乗
-------------------------------------------------*/
size_t MyFourierTransform::NextPow2(size_t s)
{
	size_t n = 1;
  
	//2の左1ビットシフトは 元の数 * 2
	while (n < s)	n <<= 1;
	return n;
}

//ビットを反転させる
void MyFourierTransform::BitReverse()
{
	uint k, b, a;
  
	for (uint ii = 0; ii < size; ii++)
	{
		k = 0;
		b = size >> 1;
		a = 1;
		
		while (b >= a)
		{
			if (b & ii) k |= a;
			if (a & ii) k |= b;
			b >>= 1;				// 1/2倍
			a <<= 1;				// 2倍
		}

		//i番目の要素とk番目の要素を入れ替える
		if (ii < k) swap(fftArray[ii], fftArray[k]);
	}
}

//ダンプ
void MyFourierTransform::Dump()
{
	//ダンプモードのフラグが立っている場合
	//配列データサイズを代入し、立っていない場合コンストラクタで与えたデータサイズを代入する
	uint end = verbose ? size : use;
  
	for (uint ii = 0; ii < end; ii++) std::cout << fftArray[ii] << " " <<std:: endl;
}


/*!-------------------------------------------------
	@brief 波形データを配列にセットする関数
	@param [in] wavedata 波データの配列(窓関数をかけた波)
	@param [in] window 窓関数の種類
----------------------------------------------------*/
void MyFourierTransform::Setwavedata(double *wavedata,WindowType window)
{
	for(int ii=0;ii<size;ii++)
	{
		Complex data(0.0,0.0);
	
		//波データの範囲内の場合
		if(ii < use)
		{
			switch(window)
			{
				//矩形窓
				case Rectangle:
				{
					data.real(wavedata[ii]);
					data.imag(0.0);
				}break;

				//ハニング窓
				case Hanning:
				{
					data.real(wavedata[ii] * (0.5 - 0.5 * cos(2.0 * M_PI * ii/(size-1))));
					data.imag(0.0);
				
				}break;
		
				//ハミング窓
				case Hamming:
				{
					data.real(wavedata[ii] * (0.54 - 0.46 * cos(2.0 * M_PI * ii/(size-1))));
					data.imag(0.0);
			
				}break;

				case WinTypeMax: break;
			}
		}

		fftArray[ii] = data;
	}
}

#pragma endregion

#pragma region 高速フーリエ変換

//高速フーリエ変換
void MyFourierTransform::FFT(bool isReverse)
{
	#pragma region 旧ソースコード

	/*Complex *outdata =  new Complex[size];
	Complex *table = new Complex[size/2];

	for(int ii=0;ii<size/2;ii++)
	{
		double arg = -2.0 * M_PI * (double)ii / (double)size;
		table[ii] = Complex(cos(arg),sin(arg));
	}

	int jj = 0;
	int n1=	0;
	int n2 = size/2;
	int k=0;
	
	outdata[0] = fftArray[0];
	outdata[size-1] = fftArray[size-1];

	for(int ii=1;ii<size-1;ii++)
	{
		n1 = n2;
		
		while( jj >= n1)
		{
			jj-=n1;
			n1/=2;
		}

		jj+=n1;
		
		if(ii < jj)
		{
			outdata[ii] = fftArray[jj];
			outdata[jj] = fftArray[ii];
		}

		else if( ii == jj)
		{
			outdata[ii] = fftArray[ii];
		}
	}

	n1 = 0;
	n2 = 1;

	int m = 2;
	int sign = (size >= 0 ? 1 : (-1));
	while( m <= size)
	{
		n1 = n2;
		n2 += n2;

		for(int jj=0;jj<n1;jj++)
		{
			for(int kk= jj; kk<size;kk+=n2)
			{
				int p = (jj * (size/2/n1));
				double t1 = table[p].real() * outdata[kk+n1].real() -  sign * table[p].imag() * outdata[kk+n1].imag();
				double t2 = sign * table[p].imag() * outdata[kk+n1].real() +  table[p].imag() * outdata[kk+n1].imag();
				Complex pp = Complex(t1,t2);	
				outdata[kk+n1] = outdata[k] - pp;
				outdata[kk] += pp;	
			}
		}
		m*=2;
	}

	for(int ii=0;ii<size;ii++)
	{
		fftArray[ii] = outdata[ii];
	}*/

	#pragma endregion

	BitReverse();
	size_t m = 2;
	Complex w, ww, t;

	//バタフライ演算
	while (m <= size)
	{
		double arg = -2.0 * M_PI / m;
		w = Complex(cos(arg), sin(arg));
		
		//-1乗 -(-2.0*PI/size) = 2.0*PI/size
		if (isReverse) w = one / w; 
    
		for(uint i=0;i<size; i+=m)
		{
			ww = 1.0;
			
			for(uint j=0;j<m/2;j++)
			{
				int a = i + j;
				int b = i + j + m/2;
	
				t = ww * fftArray[b];
	
				fftArray[b] = fftArray[a] - t;
				fftArray[a] = fftArray[a] + t;
	
				ww *= w;
			}
		}
		m *= 2;
	}
}

//逆高速フーリエ変換
void MyFourierTransform::IFFT()
{
	FFT(true);
	float s = (float)size;
	for (uint ii = 0; ii < size; ii++) fftArray[ii] /= s;
}

#pragma endregion

#pragma endregion

#pragma region ヘルパーメソッド群

/*!-------------------------------------------------------------
	@brief  各波の包含率を取得するメソッド(0.1～100Hzまでを対象)
	@param [in] FFT解析クラスポインタ
	@return 指定した波周波数帯の包含率
---------------------------------------------------------------*/
double getWaveProbability(MyFourierTransform* fft,WaveBandType wbtype)
{
	/*
	double sum = 0;
	double bandpower = 0;
	double power = 0;
	double bandpower2 = 0;
	int ii=0;*/
	
	const std::complex<double> *fftArray = fft->ReturnPtr(0);
	double spectrumwidth = 512.0 / fft->GetfftArraySize();
	double minmamband = 0.5;
    double bandpower = 0;
    double sum = 0;
    double power = 0;
    int ii = 1;                 //0が直流成分

    //分解能が足りない場合の補正
    if(spectrumwidth > minmamband)
    {
        ii = 1; //直流を含まない最初のインデクス
    }

    //分解能が足りている場合
    else
    {
        double dd = spectrumwidth;

        //分解能が足りる場合該当インデクスまで進める
        while(minmamband >= dd)
        {
            dd +=  spectrumwidth;
            ii++;
        }
    }

    for (; ii < (int)(100.0 / spectrumwidth); ii++)
    {
        power = (abs(fftArray[ii]) * abs(fftArray[ii]));

        switch (wbtype)
        {
            //δ波
			case Delta:
            {
                //bool bandmin = (ii >  (int)(minmamband / spectrumwidth));
                bool bandmax = (ii <= (int)(3.0 / spectrumwidth));

                if (bandmax == true) bandpower += power;

            } break;

            //θ波
			case Theta:
            {
                bool bandmin = (ii > (int)(4.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(7.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;
            } break;

            //α波
			case Alpha:
            {
                bool bandmin = (ii > (int)(8.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(12.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;

            //β波
            case Beta:
            {
                bool bandmin = (ii > (int)(13.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(30.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;

            //γ波
            case Gamma:
            {
                bool bandmin = (ii > (int)(30.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(100.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;
        }

        /*NeuroSky周波数定義では3.0～4.0となど周波数間が空く場合が
        あるためその分の周波数配列要素に関しては加算を除外する*/
        //δ-θ間
        bool aa = (ii > (int)(3.0 / spectrumwidth));
        bool bb = (ii <= (int)(4.0 / spectrumwidth));

        //θ-α間
        bool cc = (ii > (int)(7.0 / spectrumwidth));
        bool dd = (ii <= (int)(8.0 / spectrumwidth));

        //α-β間
        bool ee = (ii > (int)(12.0 / spectrumwidth));
        bool ff = (ii <= (int)(13.0 / spectrumwidth));

        //空文実行
        if ((aa && bb) || (cc && dd) || (ee && ff))
        {
            ;
        }

        //その他
        else
        {
            sum += power;
        }
    }

    return (bandpower/sum);
}

///<@brief 文字列分割関数
///<@param[in]  str		走査文字列
///<@param[in]	delim	分割する文字
///<@param[out] outlist	分割した文字列の位置を示したポインタの配列
///<@param[out] count	1行で分割した文字列の個数
///<@return		1行で分割した文字列の個数
int split(char *str,const char *delim,char *outlist[],int &count)
{
	char *tk;
	count = 0;
    tk = strtok( str, delim );
    
	while( tk != NULL) 
	{
        outlist[count++] = tk;
        tk = strtok( NULL, delim );
    }
    return count;
}


#pragma endregion

