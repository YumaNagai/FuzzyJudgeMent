#include "MyDataFunction.h"

#pragma region �����t�[���G�ϊ��N���XMyFourierTransform�Ɋւ��鏈��

//static�v�f�ϐ��̏�����
Complex MyFourierTransform::one = std::complex<double>(1.0,0.0);
Complex MyFourierTransform::ione = std::complex<double>(0.0,1.0);

#pragma region �w���p�[���\�b�h

/*!--------------------------------------------- 
	@brief s�ȏ�̍ŏ���2�ׂ̂����Ԃ��֐�
	@param [in] s �T���v����
	@return s�ȏ�̍ŏ���2�ׂ̂���
-------------------------------------------------*/
size_t MyFourierTransform::NextPow2(size_t s)
{
	size_t n = 1;
  
	//2�̍�1�r�b�g�V�t�g�� ���̐� * 2
	while (n < s)	n <<= 1;
	return n;
}

//�r�b�g�𔽓]������
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
			b >>= 1;				// 1/2�{
			a <<= 1;				// 2�{
		}

		//i�Ԗڂ̗v�f��k�Ԗڂ̗v�f�����ւ���
		if (ii < k) swap(fftArray[ii], fftArray[k]);
	}
}

//�_���v
void MyFourierTransform::Dump()
{
	//�_���v���[�h�̃t���O�������Ă���ꍇ
	//�z��f�[�^�T�C�Y�������A�����Ă��Ȃ��ꍇ�R���X�g���N�^�ŗ^�����f�[�^�T�C�Y��������
	uint end = verbose ? size : use;
  
	for (uint ii = 0; ii < end; ii++) std::cout << fftArray[ii] << " " <<std:: endl;
}


/*!-------------------------------------------------
	@brief �g�`�f�[�^��z��ɃZ�b�g����֐�
	@param [in] wavedata �g�f�[�^�̔z��(���֐����������g)
	@param [in] window ���֐��̎��
----------------------------------------------------*/
void MyFourierTransform::Setwavedata(double *wavedata,WindowType window)
{
	for(int ii=0;ii<size;ii++)
	{
		Complex data(0.0,0.0);
	
		//�g�f�[�^�͈͓̔��̏ꍇ
		if(ii < use)
		{
			switch(window)
			{
				//��`��
				case Rectangle:
				{
					data.real(wavedata[ii]);
					data.imag(0.0);
				}break;

				//�n�j���O��
				case Hanning:
				{
					data.real(wavedata[ii] * (0.5 - 0.5 * cos(2.0 * M_PI * ii/(size-1))));
					data.imag(0.0);
				
				}break;
		
				//�n�~���O��
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

#pragma region �����t�[���G�ϊ�

//�����t�[���G�ϊ�
void MyFourierTransform::FFT(bool isReverse)
{
	#pragma region ���\�[�X�R�[�h

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

	//�o�^�t���C���Z
	while (m <= size)
	{
		double arg = -2.0 * M_PI / m;
		w = Complex(cos(arg), sin(arg));
		
		//-1�� -(-2.0*PI/size) = 2.0*PI/size
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

//�t�����t�[���G�ϊ�
void MyFourierTransform::IFFT()
{
	FFT(true);
	float s = (float)size;
	for (uint ii = 0; ii < size; ii++) fftArray[ii] /= s;
}

#pragma endregion

#pragma endregion

#pragma region �w���p�[���\�b�h�Q

/*!-------------------------------------------------------------
	@brief  �e�g�̕�ܗ����擾���郁�\�b�h(0.1�`100Hz�܂ł�Ώ�)
	@param [in] FFT��̓N���X�|�C���^
	@return �w�肵���g���g���т̕�ܗ�
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
    int ii = 1;                 //0����������

    //����\������Ȃ��ꍇ�̕␳
    if(spectrumwidth > minmamband)
    {
        ii = 1; //�������܂܂Ȃ��ŏ��̃C���f�N�X
    }

    //����\������Ă���ꍇ
    else
    {
        double dd = spectrumwidth;

        //����\�������ꍇ�Y���C���f�N�X�܂Ői�߂�
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
            //�g
			case Delta:
            {
                //bool bandmin = (ii >  (int)(minmamband / spectrumwidth));
                bool bandmax = (ii <= (int)(3.0 / spectrumwidth));

                if (bandmax == true) bandpower += power;

            } break;

            //�Ɣg
			case Theta:
            {
                bool bandmin = (ii > (int)(4.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(7.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;
            } break;

            //���g
			case Alpha:
            {
                bool bandmin = (ii > (int)(8.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(12.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;

            //���g
            case Beta:
            {
                bool bandmin = (ii > (int)(13.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(30.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;

            //���g
            case Gamma:
            {
                bool bandmin = (ii > (int)(30.0 / spectrumwidth));
                bool bandmax = (ii <= (int)(100.0 / spectrumwidth));

                if (bandmin == true && bandmax == true) bandpower += power;

            } break;
        }

        /*NeuroSky���g����`�ł�3.0�`4.0�ƂȂǎ��g���Ԃ��󂭏ꍇ��
        ���邽�߂��̕��̎��g���z��v�f�Ɋւ��Ă͉��Z�����O����*/
        //��-�Ɗ�
        bool aa = (ii > (int)(3.0 / spectrumwidth));
        bool bb = (ii <= (int)(4.0 / spectrumwidth));

        //��-����
        bool cc = (ii > (int)(7.0 / spectrumwidth));
        bool dd = (ii <= (int)(8.0 / spectrumwidth));

        //��-����
        bool ee = (ii > (int)(12.0 / spectrumwidth));
        bool ff = (ii <= (int)(13.0 / spectrumwidth));

        //�󕶎��s
        if ((aa && bb) || (cc && dd) || (ee && ff))
        {
            ;
        }

        //���̑�
        else
        {
            sum += power;
        }
    }

    return (bandpower/sum);
}

///<@brief �����񕪊��֐�
///<@param[in]  str		����������
///<@param[in]	delim	�������镶��
///<@param[out] outlist	��������������̈ʒu���������|�C���^�̔z��
///<@param[out] count	1�s�ŕ�������������̌�
///<@return		1�s�ŕ�������������̌�
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

