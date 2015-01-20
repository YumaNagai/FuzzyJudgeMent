#include <vector>
#include "MyDataFunction.h"

using namespace std;

#pragma region (�P��)�ڕW�l-�K���x�Z�o�֐�

void calcGoodOfFit(vector<double>& gof,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable,double *bandp,
				   double prob,int start,int end,int width,WaveBandType type,bool diffFlag)
{
	double Nwgrade = 0;		//No�d�ݕt���K���x
	double Ngrade  = 0;		//No�K���x
	double Ywgrade = 0;		//Yes�d�ݕt���K���x
	double Ygrade  = 0;		//Yes�K���x
	int rulecount = ruletable.size();
	
	int index =0;

	switch(type)
	{
		case Delta: index = 0; break;
		case Theta: index = 1; break;
		case Alpha: index = 2; break;
		case Beta:  index = 3; break;
		case Gamma: index = 4; break;
		default: break;
	}

	for(int ii=0; ii<end; ii+=width)
	{
		for(int jj=start+ii; jj<=end ; jj+=width)
		{
			Ygrade = Ywgrade = Ngrade = Nwgrade = 0;
			data[index].makeMenberSipFunction(ii,jj,ii,jj);
		
			//���[���e�[�u���ɂ��K���x�v�Z
			for(int kk=0;kk<rulecount;kk++)
			{								
				double min = DBL_MAX;
					
				for(int aa=0;aa<5;aa++)
				{
					char state;
						
					switch(aa)
					{
						case 0: state = ruletable[kk].deltaState; break;
						case 1: state = ruletable[kk].thetaState; break;
						case 2: state = ruletable[kk].alphaState; break;
						case 3: state = ruletable[kk].betaState;  break;
						case 4: state = ruletable[kk].gammaState; break;
						default: break;
					}

					double ppp = data[aa].fuzzyset(bandp[aa],state);
					if(min >= ppp) min = ppp;
				}

				//Yes�o��
				Ygrade += min;
				Ywgrade += (min * ruletable[kk].yesWeight);
										
				//No�o��
				Ngrade += min;
				Nwgrade += (min * ruletable[kk].noWeight);
			}

			double Ytgrade,Ntgrade;

			Ytgrade = (Ygrade == 0) ? 0 : Ywgrade/Ygrade;
			Ntgrade = (Ngrade == 0) ? 0 : Nwgrade/Ngrade;

			//deltatmgof.push_back(probability[ii] - Ntgrade);
			if(diffFlag)	gof.push_back(prob - Ntgrade);
			else 			gof.push_back(Ntgrade);
		}
	}
}

#pragma endregion

#pragma region �P��K���x�Z�o�֐�(�������ϓ� ���ώZ�o)

void averageGoodOfFit(vector<vector<double>>& avrgof,vector<MyFourierTransform*>& fftArray,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable,
					  vector<double>& prob,int start,int end,int width,WaveBandType type,bool diffFlag)
{
	//double Nwgrade = 0;		//No�d�ݕt���K���x
	//double Ngrade  = 0;		//No�K���x
	//double Ywgrade = 0;		//Yes�d�ݕt���K���x
	//double Ygrade  = 0;		//Yes�K���x
	//int rulecount = ruletable.size();
	int index =0;

	switch(type)
	{
		case Delta: index = 0; break;
		case Theta: index = 1; break;
		case Alpha: index = 2; break;
		case Beta:  index = 3; break;
		case Gamma: index = 4; break;
		default: break;
	}

	int datacount = fftArray.size();

	for(int kk=0; kk<end; kk+=width)
	{
		vector<vector<double>> datavector;
		vector<double> averVector;
		
		for(int ll=start+kk; ll<=end ; ll+=width)
		{
			for(int dd=0;dd<5;dd++)
			{
				if(dd!=index) data[dd].makeMenberSipFunction(kk,ll,kk,ll);
			}

			for(int ii=0;ii<datacount;ii++)
			{
				double bandp[5]={0};
				vector<double> ff;

				bandp[0] = getWaveProbability(fftArray[ii],Delta) * 100 ;
				bandp[1] = getWaveProbability(fftArray[ii],Theta) * 100 ;
				bandp[2] = getWaveProbability(fftArray[ii],Alpha) * 100 ;
				bandp[3] = getWaveProbability(fftArray[ii],Beta)  * 100 ;
				bandp[4] = getWaveProbability(fftArray[ii],Gamma) * 100 ;

				calcGoodOfFit(ff,data,ruletable,bandp,prob[ii],start,end,width,type,diffFlag);
				
				datavector.push_back(ff);
			}

			for(int jj=0;jj<datavector[0].size();jj++)
			{
				double aa = 0;
				
				for(int ii=0;ii<datavector.size();ii++)
				{
					aa += datavector[ii][jj];
				}
				averVector.push_back((double)(aa/datacount));
			}
			avrgof.push_back(averVector);
		}
	}
}

#pragma endregion

#pragma region �K���x�Z�o+�t�@�C���o��

void writegof(FILE* fp,vector<MyFourierTransform*>& fftArray,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable)
{
	double Nwgrade = 0;		//No�d�ݕt���K���x
	double Ngrade  = 0;		//No�K���x
	double Ywgrade = 0;		//Yes�d�ݕt���K���x
	double Ygrade  = 0;		//Yes�K���x
	int rulecount = ruletable.size();
	int datacount = fftArray.size();

	for(int ii=0;datacount;ii++)
	{
		Ygrade = Ywgrade = Ngrade = Nwgrade = 0;
		double bandp[5]={0};
		double deltap = getWaveProbability(fftArray[ii],Delta) * 100 ;
		double thetap = getWaveProbability(fftArray[ii],Theta) * 100 ;
		double alphap = getWaveProbability(fftArray[ii],Alpha) * 100 ;
		double betap  = getWaveProbability(fftArray[ii],Beta)  * 100 ;
		double gammap = getWaveProbability(fftArray[ii],Gamma) * 100 ;

		//���[���e�[�u���ɂ��K���x�v�Z
		for(int kk=0;kk<rulecount;kk++)
		{								
			double min = DBL_MAX;
					
			for(int aa=0;aa<5;aa++)
			{
				char state;
						
				switch(aa)
				{
					case 0: state = ruletable[kk].deltaState; break;
					case 1: state = ruletable[kk].thetaState; break;
					case 2: state = ruletable[kk].alphaState; break;
					case 3: state = ruletable[kk].betaState;  break;
					case 4: state = ruletable[kk].gammaState; break;
					default: break;
				}

				double ppp = data[aa].fuzzyset(bandp[aa],state);
				if(min >= ppp) min = ppp;
			}

			//Yes�o��
			Ygrade += min;
			Ywgrade += (min * ruletable[kk].yesWeight);
										
			//No�o��
			Ngrade += min;
			Nwgrade += (min * ruletable[kk].noWeight);
		}

		double Ytgrade,Ntgrade;

		Ytgrade = (Ygrade == 0) ? 0 : Ywgrade/Ygrade;
		Ntgrade = (Ngrade == 0) ? 0 : Nwgrade/Ngrade;

		fprintf(fp,"%f,%f,%f,%f,%f,%f,%f,\n",bandp[0],bandp[1],bandp[2],bandp[3],bandp[4],Ytgrade,Ntgrade);
	}
}




#pragma endregion

#pragma region ���ʏ����o��

//����1
void result_write(FILE* fp,char *fname,vector<vector<double>>& data,int start,int width,int end)
{
	string buffer=",";
	int datacount = data.size();

	for(int ii=1;ii<=datacount;ii++)
	{
		char tmp[50];
		sprintf(tmp,"data%d,",ii);
		string t = tmp;
		buffer += t;
	}

	fp = fopen(fname,"w");
	int pp = 0;
	string disp; 
	fprintf(fp,"%s\n",buffer.c_str());

	for(int kk=0; kk<end; kk+=width)
	{
		for(int ll=start+kk; ll<=end ; ll+=width)
		{
			char dd[64];
			sprintf(dd,"%d=%d",kk,ll);
			disp = dd;
			string result;
			double ave = 0;
			for(int ii=0;ii<data.size();ii++)
			{
				ave += data[ii][pp];
				char tt[20];
				sprintf(tt,"%f,",data[ii][pp]);
				string kk = tt;
				result += kk;
			}
			fprintf(fp,"%s,%s%f\n",disp.c_str(),result.c_str(),double(ave/data.size()));
			pp++;
		}
	}
	fclose(fp);
}

//����1
void result_write_2(FILE* fp,char *fname,vector<vector<double>>& avergof,int start,int width,int end)
{
	fp = fopen(fname,"w");
	string buffer;
	
	for(int ii=0;ii<end;ii+=width)
	{
		for(int jj=ii+start;jj<=end;jj+=width)
		{
			char buffer2[32];
			sprintf(buffer2,"%d=%d,",ii,jj);
			string bb = buffer2;
			buffer += bb;
		}
	}
	
	fprintf(fp,",%s\n",buffer.c_str());
	int pp=0;

	for(int kk=0; kk<end; kk+=width)
	{
		for(int ll=start+kk; ll<=end ; ll+=width)
		{
			char dd[64];
			string disp;
			sprintf(dd,"%d=%d",kk,ll);
			disp = dd;
			string result;

			for(int ii=0;ii<avergof.size();ii++)
			{
				char tt[20];
				sprintf(tt,"%f,",avergof[ii][pp]);
				string kk = tt;
				result += kk;
			}
			fprintf(fp,"%s,%s\n",disp.c_str(),result.c_str());
			pp++;
		}
	}
	fclose(fp);
}


#pragma endregion

#pragma region ���C���֐�

//�G���g���|�C���g
int main(int argc,char *argv[])
{
	vector<vector<double>> datatmp;		//�ۑ��̂��߂̕ϐ�
	int datasize = 0;					//�T���v�������f�[�^�̐�
	int datacount = 0;							//�f�[�^���̂̐�

	vector<double*> data;			//�g�`�f�[�p�̂̃e�[�u��
	vector<double> probability;		//�m���ۑ��p�e�[�u��
	vector<FuzzyRlue> ruletable;	//���[���e�[�u��
	char buf[1024];	
	
	#pragma region �f�[�^���͕�

	//�f�[�^�̓ǂݍ���-------------------------------------------------
	FILE* fp = fopen("NormalData.csv","r");
	while(fgets(buf,sizeof(buf),fp)!=NULL)
	{
		char *outputlist[1024];		//������
		vector<double> aa;			//1�s�ǂݍ��ނ��߂̔z��
		double kk;					//�ꎞ�ϐ�

		//�J���}�ŕ���
		split(buf,",",outputlist,datacount);
		
		for(int ii=0;ii<datacount;ii++)
		{			
			sscanf(outputlist[ii],"%lf",&kk);
			aa.push_back(kk);
		}
		datatmp.push_back(aa);
	}
	datasize = datatmp.size();
	for(int ii=0;ii<datacount;ii++)
	{ 	
		double* realdata = new double[datasize];
		for(int jj=0;jj<datasize;jj++)
		{
			realdata[jj] = datatmp[jj][ii];
		}
		data.push_back(realdata);
	}
	fclose(fp);
	datatmp.clear();
	//-----------------------------------------------------------------

	//�m���e�[�u���ǂݍ���----------------------------------------------
	fp = fopen("NormalDataProbability.csv","r");
	while(fgets(buf,sizeof(buf),fp)!=NULL)
	{
		double dd;
		int cc= 0;
		char *outputlist[1024];

		split(buf,",",outputlist,cc);
		
		for(int ii=0;ii<cc;ii++)
		{
			sscanf(outputlist[ii],"%lf",&dd);
			probability.push_back(dd);
		}
	}
	fclose(fp);
	//------------------------------------------------------------------

	//�t�@�W�[���[���ǂݍ���--------------------------------------------
	fp = fopen("FuzzyRuletable.csv","r");
	while(fgets(buf,sizeof(buf),fp)!=NULL)
	{
		FuzzyRlue table;
		
		sscanf(buf,"%c,%c,%c,%c,%c,%c,%lf,%lf",
			&table.deltaState,&table.thetaState,&table.alphaState,&table.betaState,&table.gammaState,&table.output,
			&table.yesWeight,&table.noWeight);
		
		ruletable.push_back(table);
	}
	fclose(fp);
	//------------------------------------------------------------------

	#pragma endregion	

	#pragma region FFT�ϊ�

	//�f�[�^�ϊ���(FFT�{���֐�)
	vector<MyFourierTransform*> fftArray;
	for(int ii=0;ii<datacount;ii++)
	{
		MyFourierTransform* fft = MyFourierTransform::Load(datasize,data[ii],Hamming);
		fftArray.push_back(fft);
	}

	#pragma endregion

	//[0]:�� [1]:�� [2]:�� [3]:�� [4]:��
	MeberSipFunctionData bandfgydata[5];
	double Nwgrade = 0;		//No�d�ݕt���K���x
	double Ngrade  = 0;		//No�K���x
	double Ywgrade = 0;		//Yes�d�ݕt���K���x
	double Ygrade  = 0;		//Yes�K���x

	bool diffFlag = true; //true:�K���x - �ڕW�l�@or false:�K���x�̂� 

	//�ǉ�
	vector<vector<double>> deltagof;
	vector<vector<double>> thetagof;
	vector<vector<double>> alphagof;
	vector<vector<double>> betagof;
	vector<vector<double>> gammagof;

	int start = 10;		//�J�n�X�^�[�g�_
	int end = 100;		//�I��
	int width = 10;		//��

	int dst=40,ded=50;
	int tst=20,ted=30;
	int ast=20,aed=40;
	int bst=20,bed=30;
	int gst=20,ged=30;
		
	bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	bandfgydata[1].makeMenberSipFunction(tst,ted,tst,ted);
	bandfgydata[2].makeMenberSipFunction(ast,aed,ast,aed);			
	bandfgydata[3].makeMenberSipFunction(bst,bed,bst,bed);
	bandfgydata[4].makeMenberSipFunction(gst,ged,gst,ged);


	#pragma region �t�@�W�[�v�Z��

	////�f�[�^�����J��Ԃ�
	//for(int ii=0;ii<datacount;ii++)
	//{
	//	double bandp[5]={0};
	//	double deltap = getWaveProbability(fftArray[ii],Delta) * 100 ;
	//	double thetap = getWaveProbability(fftArray[ii],Theta) * 100 ;
	//	double alphap = getWaveProbability(fftArray[ii],Alpha) * 100 ;
	//	double betap  = getWaveProbability(fftArray[ii],Beta)  * 100 ;
	//	double gammap = getWaveProbability(fftArray[ii],Gamma) * 100 ;
	//	
	//	bandp[0] = deltap;
	//	bandp[1] = thetap;
	//	bandp[2] = alphap;
	//	bandp[3] = betap;
	//	bandp[4] = gammap;

	//	//fprintf(fp,"%f,%f,%f,%f,%f,\n",deltap,thetap,alphap,betap,gammap);

	//	//�t�@�C���ɋL�^
	//	char fname[1024];
	//	int rulecount = ruletable.size();
	//	
	//	#pragma region 1���̕ϓ�������

	//	
	//	int sk=30,ek=40;

	//	//�ɂ���--------------------------------------------------------------
	//	vector<double> deltatmgof;
	//	
	//	//�ȊO�Œ�
	//	bandfgydata[1].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[2].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);			
	//	calcGoodOfFit(deltatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Delta,diffFlag);
	//	deltagof.push_back(deltatmgof);

	//	//----------------------------------------------------------------------------------------
	//	
	//	//�Ƃɂ���--------------------------------------------------------------	
	//	//vector<double> thetatmgof;
	//	//
	//	////�ƈȊO�Œ�
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[2].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//calcGoodOfFit(thetatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Theta,diffFlag);
	//	//thetagof.push_back(thetatmgof);
	//	
	//	//---------------------------------------------------------------------------

	//	
	//	//���ɂ���--------------------------------------------------------------
	//	//vector<double> alphatmgof;
	//	//
	//	////���ȊO�Œ�
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,sts,ted);			
	//	//bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);
	//	//calcGoodOfFit(alphatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Alpha,diffFlag);
	//	//alphagof.push_back(alphatmgof);
	//	//-----------------------------------------------------------------------

	//	//���ɂ���--------------------------------------------------------------
	//	//vector<double> betatmgof;

	//	////���ȊO�Œ�
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,sts,ted);
	//	//bandfgydata[2].makeMenberSipFunction(ast,aed,ast,aed);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);
	//	//calcGoodOfFit(betatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Beta,diffFlag);
	//	//betagof.push_back(betatmgof);
	//	//-------------------------------------------------------------------------

	//	//���ɂ���--------------------------------------------------------------
	//	//vector<double> gammatmgof;
	//	
	//	////���ȊO�Œ�
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,tst,ted);
	//	//bandfgydata[2].makeMenberSipFunction(ast,aed,ast,aed);			
	//	//bandfgydata[3].makeMenberSipFunction(bst,bed,bst,bed);
	//	//calcGoodOfFit(gammatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Gamma,diffFlag);
	//	//gammagof.push_back(gammatmgof);
	//	

	//	#pragma endregion

	//	#pragma region (���\�[�X�R�[�h) �ϓ��{�K���x�Z�o

	//	////�ϓ�
	//	//for(int dp1=0; dp1<end; dp1+=width)
	//	//{
	//	//	for(int dp2=start+dp1; dp2<=end; dp2+=width)
	//	//	{
	//	//		#pragma region �ƕϓ���

	//	//		//�ƕϓ�
	//	//		for(int tp1=0; tp1<end; tp1+=width)
	//	//		{
	//	//			for(int tp2=start+tp1; tp2<=end; tp2+=width)
	//	//			{
	//	//				#pragma region ���ϓ���

	//	//				//���ϓ�
	//	//				for(int ap1=0; ap1<end; ap1+=width)
	//	//				{
	//	//					for(int ap2=start+ap1; ap2<=end; ap2+=width)
	//	//					{
	//	//						#pragma region ���ϓ���

	//	//						//���ϓ�
	//	//						for(int bp1=0; bp1<end; bp1+=width)
	//	//						{
	//	//							for(int bp2=start+bp1; bp2<=end; bp2+=width)
	//	//							{
	//	//								#pragma region ���ϓ���

	//	//								//���ϓ�
	//	//								for(int gp1=0; gp1<end; gp1+=width)
	//	//								{
	//	//									for(int gp2=start+gp1; gp2<=end; gp2+=width)
	//	//									{
	//	//										#pragma region ���[���e�[�u���ɂ��K���x�Z�o

	//	//										Ygrade = Ywgrade = Ngrade = Nwgrade = 0;

	//	//										bandfgydata[0].makeMenberSipFunction(dp1,dp2,dp1,dp2);
	//	//										bandfgydata[1].makeMenberSipFunction(tp1,tp2,tp1,tp2);
	//	//										bandfgydata[2].makeMenberSipFunction(ap1,ap2,ap1,ap2);
	//	//										bandfgydata[3].makeMenberSipFunction(bp1,bp2,bp1,bp2);
	//	//										bandfgydata[4].makeMenberSipFunction(gp1,gp2,gp1,gp2);

	//	//										//���[���e�[�u���ɂ��K���x�v�Z
	//	//										for(int jj=0;jj<rulecount;jj++)
	//	//										{
	//	//											double datamin[5]={0};
	//	//												
	//	//											datamin[0] = bandfgydata[0].fuzzyset(deltap,ruletable[jj].deltaState);
	//	//											datamin[1] = bandfgydata[1].fuzzyset(thetap,ruletable[jj].thetaState);
	//	//											datamin[2] = bandfgydata[2].fuzzyset(alphap,ruletable[jj].alphaState);
	//	//											datamin[3] = bandfgydata[3].fuzzyset(betap,ruletable[jj].betaState);
	//	//											datamin[4] = bandfgydata[4].fuzzyset(gammap,ruletable[jj].gammaState);
	//	//						
	//	//											//And�����̂��ߍŏ��l�����߂�
	//	//											double min = DBL_MAX;
	//	//								
	//	//											for(int kk = 0; kk<5; kk++)
	//	//											{ 
	//	//												if( min >= datamin[kk]) min = datamin[kk];
	//	//											}

	//	//						
	//	//											//Yes�o��
	//	//											Ygrade += min;
	//	//											Ywgrade += (min * ruletable[jj].yesWeight);
	//	//								
	//	//											//No�o��
	//	//											Ngrade += min;
	//	//											Nwgrade += (min * ruletable[jj].noWeight);

	//	//										}	

	//	//										#pragma endregion

	//	//										double Ytgrade,Ntgrade;

	//	//										Ytgrade = (Ygrade == 0) ? 0 : Ywgrade/Ygrade;
	//	//										Ntgrade = (Ngrade == 0) ? 0 : Nwgrade/Ngrade;

	//	//										fprintf(fp,"%d=%d,%d=%d,%d=%d,%d=%d,%d=%d,%f,%f,%f,\n",
	//	//											    dp1,dp2,tp1,tp2,ap1,ap2,bp1,bp2,gp1,gp2,Ytgrade,Ntgrade,probability[ii]-Ntgrade);
	//	//									
	//	//										
	//	//									}
	//	//								}	

	//	//								#pragma endregion
	//	//							}
	//	//						}	

	//	//						#pragma endregion
	//	//					}
	//	//				}

	//	//				#pragma endregion
	//	//			}
	//	//		}

	//	//		#pragma endregion
	//	//	}
	//	//}

	//	#pragma endregion
	//	
	//	//fprintf(fp,"deltaP:%f,thetaP:%f,alphaP:%f,betaP:%f,gammaP:%f,�ڕW�l%f\n",deltap,thetap,alphap,betap,gammap,probability[ii]);

	//	//fclose(fp);
	//}

	//fclose(fp);

	#pragma endregion

	#pragma region �t�@�C���o�͕�
	
	//result_write_2(fp,"DeltaAverage.csv",avergof,start,width,end);	

	////��------------------------------------------------------------
	//result_write(fp,"deltaT-GOF.csv",deltagof,start,width,end);
	//result_write(fp,"thetaT-GOF.csv",thetagof,start,width,end);
	//result_write(fp,"alphaT-GOF.csv",alphagof,start,width,end);
	//result_write(fp,"betaT-GOF.csv" ,betagof ,start,width,end);
	//result_write(fp,"gammaT-GOF.csv",gammagof,start,width,end);


	#pragma endregion
	
	#pragma region �f�[�^�J����

	//�f�[�^�J��
	for(int ii=0;ii<datacount;ii++)
	{
		if(fftArray[ii]!=NULL) delete fftArray[ii];
	}

	//�f�[�^�J��
	for(int ii=0;ii<datacount;ii++)
	{
		if(data[ii]!=NULL) delete[] data[ii];
	}

	fftArray.clear();
	data.clear();
	ruletable.clear();
	probability.clear();
	deltagof.clear();
	thetagof.clear();
	alphagof.clear();
	betagof.clear();
	gammagof.clear();

	#pragma endregion

	return 0;


}

#pragma endregion

