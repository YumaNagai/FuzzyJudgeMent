#include <vector>
#include "MyDataFunction.h"

using namespace std;

#pragma region (単一)目標値-適合度算出関数

void calcGoodOfFit(vector<double>& gof,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable,double *bandp,
				   double prob,int start,int end,int width,WaveBandType type,bool diffFlag)
{
	double Nwgrade = 0;		//No重み付き適合度
	double Ngrade  = 0;		//No適合度
	double Ywgrade = 0;		//Yes重み付き適合度
	double Ygrade  = 0;		//Yes適合度
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
		
			//ルールテーブルによる適合度計算
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

				//Yes出力
				Ygrade += min;
				Ywgrade += (min * ruletable[kk].yesWeight);
										
				//No出力
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

#pragma region 単一適合度算出関数(他成分変動 平均算出)

void averageGoodOfFit(vector<vector<double>>& avrgof,vector<MyFourierTransform*>& fftArray,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable,
					  vector<double>& prob,int start,int end,int width,WaveBandType type,bool diffFlag)
{
	//double Nwgrade = 0;		//No重み付き適合度
	//double Ngrade  = 0;		//No適合度
	//double Ywgrade = 0;		//Yes重み付き適合度
	//double Ygrade  = 0;		//Yes適合度
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

#pragma region 適合度算出+ファイル出力

void writegof(FILE* fp,vector<MyFourierTransform*>& fftArray,MeberSipFunctionData *data,vector<FuzzyRlue>& ruletable)
{
	double Nwgrade = 0;		//No重み付き適合度
	double Ngrade  = 0;		//No適合度
	double Ywgrade = 0;		//Yes重み付き適合度
	double Ygrade  = 0;		//Yes適合度
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

		//ルールテーブルによる適合度計算
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

			//Yes出力
			Ygrade += min;
			Ywgrade += (min * ruletable[kk].yesWeight);
										
			//No出力
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

#pragma region 結果書き出し

//その1
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

//その1
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

#pragma region メイン関数

//エントリポイント
int main(int argc,char *argv[])
{
	vector<vector<double>> datatmp;		//保存のための変数
	int datasize = 0;					//サンプルしたデータの数
	int datacount = 0;							//データ自体の数

	vector<double*> data;			//波形デー用ののテーブル
	vector<double> probability;		//確率保存用テーブル
	vector<FuzzyRlue> ruletable;	//ルールテーブル
	char buf[1024];	
	
	#pragma region データ入力部

	//データの読み込み-------------------------------------------------
	FILE* fp = fopen("NormalData.csv","r");
	while(fgets(buf,sizeof(buf),fp)!=NULL)
	{
		char *outputlist[1024];		//仮決め
		vector<double> aa;			//1行読み込むための配列
		double kk;					//一時変数

		//カンマで分割
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

	//確率テーブル読み込み----------------------------------------------
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

	//ファジールール読み込み--------------------------------------------
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

	#pragma region FFT変換

	//データ変換部(FFT＋窓関数)
	vector<MyFourierTransform*> fftArray;
	for(int ii=0;ii<datacount;ii++)
	{
		MyFourierTransform* fft = MyFourierTransform::Load(datasize,data[ii],Hamming);
		fftArray.push_back(fft);
	}

	#pragma endregion

	//[0]:δ [1]:θ [2]:α [3]:β [4]:γ
	MeberSipFunctionData bandfgydata[5];
	double Nwgrade = 0;		//No重み付き適合度
	double Ngrade  = 0;		//No適合度
	double Ywgrade = 0;		//Yes重み付き適合度
	double Ygrade  = 0;		//Yes適合度

	bool diffFlag = true; //true:適合度 - 目標値　or false:適合度のみ 

	//追加
	vector<vector<double>> deltagof;
	vector<vector<double>> thetagof;
	vector<vector<double>> alphagof;
	vector<vector<double>> betagof;
	vector<vector<double>> gammagof;

	int start = 10;		//開始スタート点
	int end = 100;		//終了
	int width = 10;		//幅

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


	#pragma region ファジー計算部

	////データ個数分繰り返し
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

	//	//ファイルに記録
	//	char fname[1024];
	//	int rulecount = ruletable.size();
	//	
	//	#pragma region 1つずつの変動を見る

	//	
	//	int sk=30,ek=40;

	//	//δについて--------------------------------------------------------------
	//	vector<double> deltatmgof;
	//	
	//	//δ以外固定
	//	bandfgydata[1].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[2].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);			
	//	calcGoodOfFit(deltatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Delta,diffFlag);
	//	deltagof.push_back(deltatmgof);

	//	//----------------------------------------------------------------------------------------
	//	
	//	//θについて--------------------------------------------------------------	
	//	//vector<double> thetatmgof;
	//	//
	//	////θ以外固定
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[2].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//calcGoodOfFit(thetatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Theta,diffFlag);
	//	//thetagof.push_back(thetatmgof);
	//	
	//	//---------------------------------------------------------------------------

	//	
	//	//αについて--------------------------------------------------------------
	//	//vector<double> alphatmgof;
	//	//
	//	////α以外固定
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,sts,ted);			
	//	//bandfgydata[3].makeMenberSipFunction(sk,ek,sk,ek);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);
	//	//calcGoodOfFit(alphatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Alpha,diffFlag);
	//	//alphagof.push_back(alphatmgof);
	//	//-----------------------------------------------------------------------

	//	//βについて--------------------------------------------------------------
	//	//vector<double> betatmgof;

	//	////β以外固定
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,sts,ted);
	//	//bandfgydata[2].makeMenberSipFunction(ast,aed,ast,aed);			
	//	//bandfgydata[4].makeMenberSipFunction(sk,ek,sk,ek);
	//	//calcGoodOfFit(betatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Beta,diffFlag);
	//	//betagof.push_back(betatmgof);
	//	//-------------------------------------------------------------------------

	//	//γについて--------------------------------------------------------------
	//	//vector<double> gammatmgof;
	//	
	//	////γ以外固定
	//	//bandfgydata[0].makeMenberSipFunction(dst,ded,dst,ded);			
	//	//bandfgydata[1].makeMenberSipFunction(tst,ted,tst,ted);
	//	//bandfgydata[2].makeMenberSipFunction(ast,aed,ast,aed);			
	//	//bandfgydata[3].makeMenberSipFunction(bst,bed,bst,bed);
	//	//calcGoodOfFit(gammatmgof,bandfgydata,ruletable,bandp,probability[ii],start,end,width,Gamma,diffFlag);
	//	//gammagof.push_back(gammatmgof);
	//	

	//	#pragma endregion

	//	#pragma region (旧ソースコード) 変動＋適合度算出

	//	////δ変動
	//	//for(int dp1=0; dp1<end; dp1+=width)
	//	//{
	//	//	for(int dp2=start+dp1; dp2<=end; dp2+=width)
	//	//	{
	//	//		#pragma region θ変動部

	//	//		//θ変動
	//	//		for(int tp1=0; tp1<end; tp1+=width)
	//	//		{
	//	//			for(int tp2=start+tp1; tp2<=end; tp2+=width)
	//	//			{
	//	//				#pragma region α変動部

	//	//				//α変動
	//	//				for(int ap1=0; ap1<end; ap1+=width)
	//	//				{
	//	//					for(int ap2=start+ap1; ap2<=end; ap2+=width)
	//	//					{
	//	//						#pragma region β変動分

	//	//						//β変動
	//	//						for(int bp1=0; bp1<end; bp1+=width)
	//	//						{
	//	//							for(int bp2=start+bp1; bp2<=end; bp2+=width)
	//	//							{
	//	//								#pragma region γ変動分

	//	//								//γ変動
	//	//								for(int gp1=0; gp1<end; gp1+=width)
	//	//								{
	//	//									for(int gp2=start+gp1; gp2<=end; gp2+=width)
	//	//									{
	//	//										#pragma region ルールテーブルによる適合度算出

	//	//										Ygrade = Ywgrade = Ngrade = Nwgrade = 0;

	//	//										bandfgydata[0].makeMenberSipFunction(dp1,dp2,dp1,dp2);
	//	//										bandfgydata[1].makeMenberSipFunction(tp1,tp2,tp1,tp2);
	//	//										bandfgydata[2].makeMenberSipFunction(ap1,ap2,ap1,ap2);
	//	//										bandfgydata[3].makeMenberSipFunction(bp1,bp2,bp1,bp2);
	//	//										bandfgydata[4].makeMenberSipFunction(gp1,gp2,gp1,gp2);

	//	//										//ルールテーブルによる適合度計算
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
	//	//											//And結合のため最小値を求める
	//	//											double min = DBL_MAX;
	//	//								
	//	//											for(int kk = 0; kk<5; kk++)
	//	//											{ 
	//	//												if( min >= datamin[kk]) min = datamin[kk];
	//	//											}

	//	//						
	//	//											//Yes出力
	//	//											Ygrade += min;
	//	//											Ywgrade += (min * ruletable[jj].yesWeight);
	//	//								
	//	//											//No出力
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
	//	//fprintf(fp,"deltaP:%f,thetaP:%f,alphaP:%f,betaP:%f,gammaP:%f,目標値%f\n",deltap,thetap,alphap,betap,gammap,probability[ii]);

	//	//fclose(fp);
	//}

	//fclose(fp);

	#pragma endregion

	#pragma region ファイル出力部
	
	//result_write_2(fp,"DeltaAverage.csv",avergof,start,width,end);	

	////δ------------------------------------------------------------
	//result_write(fp,"deltaT-GOF.csv",deltagof,start,width,end);
	//result_write(fp,"thetaT-GOF.csv",thetagof,start,width,end);
	//result_write(fp,"alphaT-GOF.csv",alphagof,start,width,end);
	//result_write(fp,"betaT-GOF.csv" ,betagof ,start,width,end);
	//result_write(fp,"gammaT-GOF.csv",gammagof,start,width,end);


	#pragma endregion
	
	#pragma region データ開放部

	//データ開放
	for(int ii=0;ii<datacount;ii++)
	{
		if(fftArray[ii]!=NULL) delete fftArray[ii];
	}

	//データ開放
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

