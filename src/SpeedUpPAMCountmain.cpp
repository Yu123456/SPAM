// SpeedUpPAM Count 算法主程序

// 输入数据格式要求
// 第二行为属性名称
// 所有数据均以逗号为分隔符

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

	//从文件中读取下述信息
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

	//创建数据集对象
	DATA DataSet(filename, nrobject, nrfeatures, N, MaxIter);

	cout<<endl;
	cout<<"Data Filename: "<<filename<<endl;
	cout<<"Data Size: "<<nrobject<<endl;
	cout<<"Feature Number: "<<nrfeatures<<endl;
	cout<<"Cluster Number: "<<N<<endl;
	cout<<"Max Iterations: "<<MaxIter<<endl;
	// 提取文件名，不含文件扩展名
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
	
	//运行时间计算起点
	clock_t start_time=clock();

	//打印数据集
	//DataSet.print();

	//中心点数组
	// 数组索引为簇编号，索引值为中心点在数据集中的编号
	int *Medoids=new int[N];
	//初始化中心点
	DataSet.RandomMedoids(Medoids, N);

	//输出中心点
/*	cout<<"初始化中心点："<<endl;
	for(int i=0; i<N; i++)
	{
		cout<<Medoids[i]+1<<"\t";
		DataSet.print(cout, Medoids[i]);
		cout<<endl;
	} 
	cout<<endl;  */


	// 初始分派
	DataSet.ClusterFirst(Medoids,N);

/*	// 输出聚类结果
	cout<<"聚类结果："<<endl;
	DataSet.print_label();
	cout<<endl;  */

	// 聚类总代价
	//double total_cost=DataSet.TotalCost();

	//cout<<"聚类总代价："<<endl;
	//cout<<total_cost<<endl;

	// 迭代寻求最优聚类
	int Iter=0;
	//cout<<"迭代开始："<<endl;
	while( Iter< MaxIter )
	{
		// 选出最优替换及被替换中心点
		pair<int,int> Replace=DataSet.SelectO(Medoids,N);
		if( Replace.first > -1 )    // 有更优替换点
		{
			Iter++;
			int index=Medoids[Replace.first];
			Medoids[Replace.first]=Replace.second;    // 交换中心点（序号）
			// 重新分派
			//( 更新后的中心点数组，中心点被替换对象索引，被替换中心点序号，簇个数）
			DataSet.Cluster(Medoids,Replace.first,index ,N);
			//cout<<"total cost : "<< DataSet.TotalCost()<<endl;
		}
		else
			break;
	}

	//输出中心点
/*	cout<<"最终中心点："<<endl;
	for(int i=0; i<N; i++)
	{
		cout<<Medoids[i]<<"\t";
		DataSet.print(cout, Medoids[i]);
		cout<<endl;
	} 
	cout<<endl;

	// 输出聚类结果
	cout<<"聚类结果："<<endl;
	DataSet.print_label();
	cout<<endl;  */

	// 最终聚类总代价
	double total_cost=DataSet.TotalCost();

	//cout<<"最终聚类总代价："<<endl;
	//cout<<total_cost<<endl;

	// 运行时间计算结束
	clock_t end_time=clock();
	double dtime=running_time(start_time,end_time);

	//聚类结果文件流输出
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
