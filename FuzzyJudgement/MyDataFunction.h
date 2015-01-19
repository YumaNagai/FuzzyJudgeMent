/******************************************************************************/
/*! @file       MyDataFunction.h
    @brief      �ꎟ���g�f�[�^�ɑ΂��鏈�����܂Ƃ߂��N���X�A��`�Q
*******************************************************************************
    �t�@�C���̏ڍׂȐ���
*******************************************************************************
    @date       �쐬��(2014/10/24)
    @author		YumaNagai
    @par        History
    - 20XX/10/24 YumaNagai
      -# Initial Version
	  -# �쐬
******************************************************************************/


#ifndef _MYDATAFUNCTION_H
#define _MYDATAFUNCTION_H

#define _USE_MATH_DEFINES	// for C++
#include <math.h>
#include <iostream>
#include <algorithm>		// swap
#include <complex>
#include <vector>

//typedef�錾
typedef std::complex<double> Complex;
typedef unsigned int uint;

#pragma region �񋓑̒�`

/*!-----------------------------
	@enum ���֐��񋓑�
	@brief ���֐��̎��
---------------------------------*/
enum WaveBandType
{
	Delta = 0,		///< @brief �g(0.1�`3.0Hz)
	Theta,			///< @brief �Ɣg(4.0�`7.0Hz)
	Alpha,			///< @brief ���g(8.0�`12.0Hz)
	Beta,			///< @brief ���g(13.0�`30.0Hz)
	Gamma,			///< @brief ���g(30.0�`100.0Hz)
	WaveBandMax,	///< @brief �ԕ�
};

/*!-----------------------------
	@enum ���֐��񋓑�
	@brief ���֐��̎��
---------------------------------*/
enum WindowType
{
	Rectangle = 0,	///< @brief ��`��
	Hanning,		///< @brief �n�j���O��
	Hamming,		///< @brief �n�~���O��
	Blackman,		///< @brief �u���b�N�}����
	WinTypeMax		///< @brief �ԕ�
};

#pragma endregion

#pragma region �����t�[���G�ϊ����s���N���X

/*! 
 @class FftTransration
 @brief �����t�[���G�ϊ����s���N���X
*/
class MyFourierTransform
{
public:


private:
	typedef unsigned int uint;

	//�����o�[�ϐ�
	static Complex one;		///<@brief �v�f(1.0,0.0)�̕��f��
	static Complex ione;		///<@brief �v�f(0.0,1.0)�̕��f��

	Complex *fftArray;			///<@brief ���f���̔z��
	size_t size;							///<@brief fftArray�̃T�C�Y
	size_t use;								///<@brief �R���X�g���N�^�ŗ^�����T�C�Y
	bool verbose;							///<@brief �z��_���v�̃��[�h�w��Ɏg��

	//�����o�[�֐�
	
	/*! @brief �r�b�g������s���֐�  */
	void BitReverse();

	/*!--------------------------------------------- 
		@brief s�ȏ�̍ŏ���2�ׂ̂����Ԃ��֐�
		@param [in] s �T���v����
		@return s�ȏ�̍ŏ���2�ׂ̂���
	-------------------------------------------------*/
	size_t NextPow2(size_t s);
	
	
	/*!-------------------------------------------------
		@brief �R���X�g���N�^
		@param s [in] �M���̒���(�T���v����) 
	---------------------------------------------------*/
	MyFourierTransform(size_t s) : use(s), size(NextPow2(s)), verbose(false)
	{
		one = std::complex<double>(1.0,0.0);

		fftArray = new Complex[size];
		for(int ii=0;ii<size;ii++) fftArray[ii] = std::complex<double>(0.0,0.0);
	}

	/*!-------------------------------------------------
		@brief �g�`�f�[�^��z��ɃZ�b�g����֐�
		@param [in] wavedata �g�f�[�^�̔z��(���֐����������g)
		@param [in] ���֐��̎��
	----------------------------------------------------*/
	void Setwavedata(double *wavedata,WindowType window);
	
public:
	
	/*!-----------------------------------------------------
		@brief �������\�b�h
		@param [in] sample �M���̒���(�T���v����) 
		@param [in] wavedata �g�f�[�^�̔z��(���֐����������g)
		@param [in] window ���֐��̎��
		@return ��������������FFTTransform�N���X�̎���
	---------------------------------------------------------*/
	static MyFourierTransform* Load(size_t sample,double *wavedata,WindowType window)
	{
		MyFourierTransform* ffttrans = new MyFourierTransform(sample);	//�f�[�^����
		ffttrans->Setwavedata(wavedata,window);				//�g�f�[�^�Z�b�g
		ffttrans->FFT();									//�����t�[���G�ϊ����s��
		return ffttrans;									//�����t�[���G�ϊ����s�����N���X�f�[�^��Ԃ�
	}

	///< @brief �f�X�g���N�^
	~MyFourierTransform()
	{
		delete[] fftArray;
	};

	///< @brief �t�ϊ����s�����ǂ����̃t���O
	///< @param[in] b true:�t�t�[���G�ϊ����[�h false:�ʏ�S���t�[���G�ϊ����[�h
	bool setVerbose(bool b){  return verbose = b; };
  
	void Dump();									///< @brief �f�[�^�_���v(�f�[�^�\��)
	void FFT(bool isReverse=false);					///< @brief �����t�[���G�ϊ�
	void IFFT();									///< @brief �t�����t�[���G�ϊ�
	size_t GetfftArraySize(){ return size;}			///< @brief FFT�Ŏg�p����z��T�C�Y�̎擾

	
	///< @brief �z��v�f�ɃA�N�Z�X���邽�߂̉��Z�q
	///< @param[in] index �z��Y���ԍ�
	///< @return �Y���ԍ��̔z��v�f
	const Complex& operator[](int index)const
	{
		if (index < 0) return fftArray[use - index];
		return fftArray[index];
	};

	///< @brief �z��̃A�h���X��Ԃ��֐�
	///< @param index �z��Y���ԍ�
	///< @return �Y���ԍ��̔z��A�h���X
	const  Complex* ReturnPtr(int index)
	{
		if (index < 0) return (fftArray - (use - index));
		return (fftArray+index);
	}
};

#pragma endregion

#pragma region �t�@�W�[��̓N���X

///< @struct MeberSipFunctionData
///< @brief  �����o�[�V�b�v�֐��̒l��ۑ����Ă����\����
struct MeberSipFunctionData
{
public:

	#pragma region �����o�[�ϐ�

	double sstart;		//Small�X�^�[�g�ʒu
	double send;		//Small�I���ʒu
	double lstart;		//Large�X�^�[�g�ʒu
	double lend;		//Large�X�^�[�g�ʒu

	double sdegree;		///<@brief Small�̎��̌X��
	double sintercept;	///<@brief Small�̎��̐ؕ�
	double ldegree;		///<@brief Large�̎��̌X��
	double lintercept;	///<@brief Large�̎��̐ؕ�
	
	#pragma endregion

	#pragma region �֐��Q

	///<@brief �R���X�g���N�^(������)
	MeberSipFunctionData()
	{
		sstart = send = lstart = lend = 0;
		sdegree = sintercept = ldegree = lintercept = 0.0;
	}

	///<@brief �����o�[�V�b�v�֐��쐬���\�b�h
	///<@param[in] ls Small�X�^�[�g�ʒu
	///<@param[in] le Small�I���ʒu
	///<@param[in] hs Large�X�^�[�g�ʒu
	///<@param[in] he Large�I���ʒu
	///<@param[in] lh Large or Small
	void makeMenberSipFunction(double ss,double se,double ls,double le)
	{
		sstart = ss;
		send = se;
		lstart = ls;
		lend = le;

		//�� 
		//(x1,y1) (x2,y2) ����X�������߂� a=(y2-y1)/(x2-x1) 
		// b = y1 - a * x1

		//'Small'�̎�
		sdegree =  -1.0 / (send - sstart);
		sintercept = 1.0 - sdegree * sstart;		

		//'Large'�̎�
		ldegree = 1.0 / (lend - lstart);			
		lintercept = -ldegree * lstart;	
	}

	///<@brief �K���x���v�Z���Ԃ��֐�
	///<@@aram[in] value ���ׂ�l(��ܗ�)
	///<@param[in] lh 'S'��'L'�̕�����
	double fuzzyset(double value,char sl)
	{	
		double grade = 0.0;

		//'Small'�̎�
		if(sl=='S')
		{
			// sstart�ȉ��̎�
			if(sstart > value) grade = 1.0;

			// sstart <= value <= send
			else if(sstart <= value && send >= value) grade = (sdegree * value) + sintercept;

			// send�ȏ�̎�
			else if(send < value) grade = 0.0;
		}
	
		//Large
		else
		{
			// lstart�ȉ��̎�
			if(lstart > value) grade = 0.0;

			// lstart <= value <= lend
			else if(lstart <= value && lend >= value) grade = (ldegree * value) + lintercept;
			
			// lend�ȏ�̎�
			else if(lend < value) grade = 1.0;
		}
		return grade;
	}

	#pragma endregion

};

///< @struct FuzzyrlueTable
///< @brief  1�g�̃t�@�W�[���[����\���\����
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

	//�R���X�g���N�^
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

#pragma region �w���p�[���\�b�h�Q

/*!-------------------------------------------------------------
	@brief  �e�g�̕�ܗ����擾���郁�\�b�h(0.1�`100Hz�܂ł�Ώ�)
	@param [in] FFT��̓N���X�|�C���^
	@return �w�肵���g���g���т̕�ܗ�
---------------------------------------------------------------*/
double getWaveProbability(MyFourierTransform* fft,WaveBandType wbtype);

/*!------------------------------------------------------------------------
	@brief �����񕪊��֐�
	@param[in]  str		����������
	@param[in]	delim	�������镶��
	@param[out] outlist	��������������̈ʒu���������|�C���^�̔z��
	@param[out] count	1�s�ŕ�������������̌�
	@return		1�s�ŕ�������������̌�
----------------------------------------------------------------------------*/
int split(char *str,const char *delim,char *outlist[],int &count);

#pragma endregion

#endif