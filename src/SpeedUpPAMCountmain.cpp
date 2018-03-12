// SpeedUpPAM Count �㷨������

// �������ݸ�ʽҪ��
// �ڶ���Ϊ��������
// �������ݾ��Զ���Ϊ�ָ���

#include"SpeedUpPAMCount.h"

int main()
{
	cout<<"Speed Up PAM Count Running Start..."<<endl;
/*	cout<<"Data Filename:"<<endl;
	string filename;
	cin>>filename;
	cout<<"Data Size:"<<endl;
	int nrobject;
	cin>>nrobject;
	cout<<"Features Number:"<<endl;
	int nrfeatures;
	cin>>nrfeatures;  */

	//���ļ��ж�ȡ������Ϣ
	cout<<"Input Data Filename:"<<endl;
	string filename;
	cin>>filename;
	// Data Size
	int nrobject=0;
	// Features Number
	int nrfeatures=0;
	// Cluster Nunmber
	int N=0;
	// Max Iterations
	int MaxIter=0;

	//�������ݼ�����
	DATA DataSet(filename, nrobject, nrfeatures, N, MaxIter);

	cout<<endl;
	cout<<"Data Filename: "<<filename<<endl;
	cout<<"Data Size: "<<nrobject<<endl;
	cout<<"Feature Number: "<<nrfeatures<<endl;
	cout<<"Cluster Number: "<<N<<endl;
	cout<<"Max Iterations: "<<MaxIter<<endl;
	// ��ȡ�ļ����������ļ���չ��
	string outfile=SubString(filename)+"ClusterResult.csv";
	cout<<"Output Filename: "<<outfile<<endl;
	

/*	cout<<"Data Filename: bezdekIris.csv"<<endl;
	string filename="bezdekIris.csv";
	cout<<"Data Size: 150"<<endl;
	int nrobject=150;
	cout<<"Feature Number: 4"<<endl;
	int nrfeatures=4;
	cout<<"Cluster number: 3"<<endl;
	int N=3;
	cout<<"Max Iterations: 20"<<endl;
	int MaxIter=20;
	cout<<"Output Filename: bezdekIrisClusterResult.csv"<<endl;
	string outfile="bezdekIrisClusterResult.csv";  */

/*	cout<<"Data Filename: authentication.csv"<<endl;
	string filename="authentication.csv";
	cout<<"Data Size: 1372"<<endl;
	int nrobject=1372;
	cout<<"Feature Number: 5"<<endl;
	int nrfeatures=5;
	cout<<"Cluster Number: 5"<<endl;
	int N=5;
	cout<<"Max Iterations: 20"<<endl;
	int MaxIter=20;
	cout<<"Output Filename: authenticationClusterResult.csv"<<endl;
	string outfile="authenticationClusterResult.csv";    */
	
	//����ʱ��������
	clock_t start_time=clock();

	//��ӡ���ݼ�
	//DataSet.print();

	//���ĵ�����
	// ��������Ϊ�ر�ţ�����ֵΪ���ĵ������ݼ��еı��
	int *Medoids=new int[N];
	//��ʼ�����ĵ�
	DataSet.RandomMedoids(Medoids, N);

	//������ĵ�
/*	cout<<"��ʼ�����ĵ㣺"<<endl;
	for(int i=0; i<N; i++)
	{
		cout<<Medoids[i]+1<<"\t";
		DataSet.print(cout, Medoids[i]);
		cout<<endl;
	} 
	cout<<endl;  */


	// ��ʼ����
	DataSet.ClusterFirst(Medoids,N);

/*	// ���������
	cout<<"��������"<<endl;
	DataSet.print_label();
	cout<<endl;  */

	// �����ܴ���
	//double total_cost=DataSet.TotalCost();

	//cout<<"�����ܴ��ۣ�"<<endl;
	//cout<<total_cost<<endl;

	// ����Ѱ�����ž���
	int Iter=0;
	//cout<<"������ʼ��"<<endl;
	while( Iter< MaxIter )
	{
		// ѡ�������滻�����滻���ĵ�
		pair<int,int> Replace=DataSet.SelectO(Medoids,N);
		if( Replace.first > -1 )    // �и����滻��
		{
			Iter++;
			int index=Medoids[Replace.first];
			Medoids[Replace.first]=Replace.second;    // �������ĵ㣨��ţ�
			// ���·���
			//( ���º�����ĵ����飬���ĵ㱻�滻�������������滻���ĵ���ţ��ظ�����
			DataSet.Cluster(Medoids,Replace.first,index ,N);
			//cout<<"total cost : "<< DataSet.TotalCost()<<endl;
		}
		else
			break;
	}

	//������ĵ�
/*	cout<<"�������ĵ㣺"<<endl;
	for(int i=0; i<N; i++)
	{
		cout<<Medoids[i]<<"\t";
		DataSet.print(cout, Medoids[i]);
		cout<<endl;
	} 
	cout<<endl;

	// ���������
	cout<<"��������"<<endl;
	DataSet.print_label();
	cout<<endl;  */

	// ���վ����ܴ���
	double total_cost=DataSet.TotalCost();

	//cout<<"���վ����ܴ��ۣ�"<<endl;
	//cout<<total_cost<<endl;

	// ����ʱ��������
	clock_t end_time=clock();
	double dtime=running_time(start_time,end_time);

	//�������ļ������
	ofstream fout;
	open_outfile(fout,outfile);
	fout<<"Speed Up PAM Count Algorithm to "<<filename<<" Data File Classficate into "<<N<<" Cluster:"<<endl;
	fout<<"Speed Up PAM Count Running time: "<<","<<dtime<<endl;
	fout<<"Cluster Total Cost: "<<","<<total_cost<<endl;
	fout<<"Iterations: "<<","<<Iter<<endl;
	fout<<"Number of distance calculation: "<<","<<CountNum<<endl;
	fout<<endl;
	fout<<"Optimal Medoids: "<<endl;
	for(int i=0; i<N; i++)
	{
		fout<<Medoids[i]+1<<",";
		DataSet.print(fout,Medoids[i]);
		fout<<endl;
	}
	fout<<endl;
	fout<<"Cluster Results:"<<endl;
	DataSet.print_label(fout);


	cout<<"Speed Up PAM Count Running Over!"<<endl;
	

}
