// Speed Up PAM 聚类算法
// 本程序编写于2015年11月29日
// 实现 k-medoids 算法中主流算法 PAM 
// 参考文献：
// 手稿，2015年11月26日 及 PAM 算法

#ifndef SPEEDUPPAMCOUNT_H
#define SPEEDUPPAMCOUNT_H

#include<iostream>
#include<fstream>
#include<sstream>   //字符串流
#include<ctime>     //运行时间
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<map>
#include<set>

using namespace std;

static unsigned long int CountNum=0;

//类的前向声明
class Point;      //样本点类
class DATA;       // 数据集类

//非类成员函数声明
// 运行时函数
double running_time(const clock_t &, const clock_t &);
//打开输入文件流函数
ifstream& open_file(ifstream &, const string &);
//打开输出文件流函数
ofstream& open_outfile(ofstream &, const string &);
// Point 对象关系操作符重载函数
//由于需要使用Point 私有成员，设置成 Point 类的友元函数
inline bool operator!=( const Point &, const Point &);
// Point 对象关系操作符重载函数
inline bool operator==( const Point &, const Point &);
//字符串转换函数
// ( 目标字符串，被替换字符或字符串，替换字符或字符串）
void string_replace(string &, const string &,const string &);
//提取从逗号开始至末尾的整型
//( 目标字符串）
int SubStringToInt(string&);
// 提取从开始到点号前的字符串
//( 目标字符串）
string SubString(const string&);



// 非类成员函数定义

// 运行时函数
double running_time(const clock_t &s, const clock_t &e)
{
	return static_cast<double>(e-s)/CLOCKS_PER_SEC;
}

//打开输入文件流函数
ifstream& open_file(ifstream &in, const string &file)
{
	in.close();
	in.clear();
	in.open(file.c_str());
	return in;
}

//打开输出文件流函数
ofstream& open_outfile(ofstream &out, const string &file)
{
	out.close();
	out.clear();
	out.open(file.c_str(), ofstream::trunc);   //清空原有文件
	return out;
}

//字符串转换函数
void string_replace(string &st, const string &pre_str,const string &post_str)
{
	string::size_type pos=0;
	string::size_type n=pre_str.size();
	while((pos=st.find_first_of(pre_str,pos))!=string::npos)
	{
		st.replace(pos,n,post_str);
		pos++;
	}
}

//提取从逗号开始至末尾的整型数
//( 目标字符串）
int SubStringToInt(string &str)
{
	string::size_type pos_start=str.find_first_of(",",0);
	++pos_start;
	string::size_type pos_end=str.find_first_of(",",pos_start);
	// 夹在两个逗号之间的数字提取出来
	string st=str.substr(pos_start,pos_end-pos_start);
	istringstream istr(st);
	int integer=0;
	istr>>integer;
	return integer;
}

// 提取从开始到点号前的字符串
//( 目标字符串）
string SubString(const string &str)
{
	return str.substr(0,str.size()-4);
}



// Point 类定义
class Point
{
	//友元函数
	friend bool operator!=( const Point &, const Point &);

private:
	int NrFeatures;     //特征个数/样本维数
	double *p;          // 样本点指针，指向动态数组
	void del_p()
	{
		if( p != NULL )
			delete [] p;
	}

public:
	// 默认构造函数
	Point():NrFeatures(0),p(NULL){}
	// 含参构造函数
	Point( istringstream&, int );   
	// 复制构造函数
	Point(const Point&);
	// 赋值操作符
	Point& operator=(const Point &);
	// 析构函数
	~Point()
	{
		del_p();
	}
	// 打印函数
	void print(ostream&);
	//打印函数
	void print(ofstream&);
	// 返回特征个数
	int get_NrF()
	{
		return NrFeatures;
	}
	// 距离计算函数
	double dist_square( Point &);
};

// 复制构造函数定义
// 定义成值型指针复制
Point::Point( const Point &point):p(NULL),NrFeatures(point.NrFeatures)
{
	//cout<<"Point 复制构造函数运行……"<<endl;
	if ( point.p != NULL)
	{
		p=new double[NrFeatures];
		for( int i=0; i<NrFeatures; i++)
			p[i]=point.p[i];
	}
}
// 赋值操作符
Point& Point::operator=(const Point &point)
{
	//cout<<"Point 赋值操作运行……"<<endl;
	//释放左边对象
	del_p();

	NrFeatures=point.NrFeatures;
	if( point.p != NULL )
	{
		p=new double[NrFeatures];
		for(int i=0; i<NrFeatures; i++)
			p[i]=point.p[i];
	}
	else 
	{
		p=NULL;
	}
	return *this;
}

// 含参构造函数定义
Point::Point( istringstream &istr, int nrf ):NrFeatures(nrf)
{
	p=new double[nrf];
	double dstr;
	for(int i=0; i<nrf; i++)
	{
		istr>>dstr;
		p[i]=dstr;
	}
}

//打印函数
void Point::print(ostream &os)
{
	for(int i=0; i<NrFeatures; i++)
		os<<p[i]<<"\t";
}
//打印函数
void Point::print(ofstream &fout)
{
	for(int i=0; i<NrFeatures; i++)
		fout<<p[i]<<",";
}

// 距离计算函数
double Point::dist_square( Point &rp)
{
	CountNum++;
	double sum=0.0;
	for(int i=0; i<NrFeatures; i++)
		sum+=(p[i]-rp.p[i])*(p[i]-rp.p[i]);
	return sum;
}



// 数据集类定义
class DATA
{
private:
	int NrObject;    // 数据集大小
	string *DimStr;   // 特征名称
	Point *data;     // 数据集
	int *label;      // 样本类别标记，数组索引值为簇编号
	double *DistSquare;    // 样本与簇中心点距离平方
	DATA(const DATA&);     //复制构造函数，只声明不定义用于禁止复制

public:
	DATA(const string&, int &, int &, int &, int &);   //构造函数
	~DATA()
	{
		if( DimStr != NULL )
			delete [] DimStr;
		if( data != NULL )
			delete [] data;
		if( label != NULL )
			delete [] label;
		if( DistSquare != NULL )
			delete [] DistSquare;
	}
	//打印函数
	void print();
	// 重载打印函数
	void print(ostream &, int);
	// 中心点文件流输出
	void print(ofstream &, int);
	// 打印类别编号
	void print_label();
	// 文件流输出函数
	void print_label(ofstream&);
	// 比较两个 Point 对象是否相同
	bool compare(const int, const int) const;
	// 比较两个 Point 对象是否相同
	bool compare(const int, const int);
	// 数据与中心点比较函数 (确保重复样本也可以执行)
	// 如果数据与中心点组中某一个相同，返回 true, 否则返回 false
	// ( 数据编号，中心点数组 , 簇个数)
	bool compare(int , int *, int );
	// 返回数据集大小
	int get_nro()
	{
		return NrObject;
	}
	// 计算两个 Point 对象的距离平方
	double dist_square(int num1, int num2)
	{
		return data[num1].dist_square(data[num2]);

	}
	// 设置簇标号
	void set_label(int index, int n_label)
	{
		label[index]=n_label;
	}
	// 设置距离平方
	void set_DistSquare(int index, double dist)
	{
		DistSquare[index]=dist;
	}
	// 获取簇标号
	int get_label(int index)
	{
		return label[index];
	}
	// 聚类完成后，总代价计算函数
	double TotalCost();
	//初始化中心点
	// ( 中心点数组，簇个数)
	void RandomMedoids(int *, int nc);
	// 初始分派函数
	// ( 中心点数组指针，簇个数）
	void ClusterFirst(int *, int);
	// 交换函数, 返回交换代价
	double SwapOH(int*,int,int,int);
	// 选出最优替换及被替换中心点，返回其序号
	pair<int, int> SelectO(int*, int);
	// 去除原簇中心点后，返回第二小距离值平方
	// ( 中心点组，计算点编号，所属簇索引,簇个数） 
	double SecondDist( int *, int , int , int );
	// 重新分派聚类函数
	// ( 更新后的中心点数组，中心点被替换对象索引，被替换中心点序号，簇个数）
	void Cluster(int *,int , int ,int );
	// 计算最小距离平方及对应簇标号
	// 返回值第一个为最小值簇标号，第二只最小值时距离平方
	//(中心点数组，待分派样本序号，簇个数）
	pair<int,double> MinDistLabel(int*,int,int);
};

// DATA 类成员函数定义
// 构造函数定义
// 形参表 ( 文件名，数据集大小，特征个数, 簇个数，最大迭代次数)
DATA::DATA(const string &file, int &nro, int &nrf, int &nc, int &maxiter):NrObject(nro),
	DimStr(NULL),data(NULL),label(NULL),DistSquare(NULL)
{
	// 文件流
	ifstream in;
	open_file(in,file);
	string in_str;
	getline(in,in_str);
	nro=SubStringToInt(in_str);
	getline(in,in_str);
	nrf=SubStringToInt(in_str);
	getline(in,in_str);
	nc=SubStringToInt(in_str);
	getline(in,in_str);
	maxiter=SubStringToInt(in_str);

	NrObject=nro;

	DimStr=new string[nrf];
	data=new Point[nro];
	label=new int[nro];
	DistSquare=new double[nro];
	
	// 读入特征名称
	getline(in,in_str);
	string_replace(in_str,",","\t");
	istringstream str(in_str);
	string out_str;
	for(int i=0; i<nrf; i++)
	{
		str>>out_str;
		DimStr[i]=out_str;
	}
	for(int j=0; j<nro; j++)
	{
		getline(in,in_str);
		string_replace(in_str,",","\t");
		istringstream str(in_str);
		data[j]=Point(str,nrf);
	}
	//关闭文件流
	in.close();
	in.clear();

	return;
}

//打印函数
void DATA::print()
{
	int nrf=data[0].get_NrF();  // 特征个数
	cout<<"序号"<<"\t";
	for(int i=0; i<nrf; i++)
		cout<<DimStr[i]<<"\t";
	cout<<endl;
	for(int i=0; i<NrObject; i++)
	{
		cout<<i+1<<"\t";
		data[i].print(cout);
		cout<<endl;
	}
}
// 重载打印函数
// ( 输出流，需要打印样本编号)
void DATA::print(ostream &os, int num)
{
	data[num].print(os);
}
// 中心点文件流输出
void DATA::print(ofstream &fout, int num)
{
	data[num].print(fout);
}
// 打印类别编号
void DATA::print_label()
{
	int nrf=data[0].get_NrF();   //特征个数
	cout<<"序号"<<"\t";
	for(int i=0; i<nrf; i++)
		cout<<DimStr[i]<<"\t";
	cout<<"类别编号"<<endl;
	for(int i=0; i<NrObject; i++)
	{
		cout<<i+1<<"\t";
		data[i].print(cout);
		cout<<label[i]+1<<endl;
	}
}

// 文件流输出函数
void DATA::print_label( ofstream &fout)
{
	int nrf=data[0].get_NrF();   //特征个数

	//数据按簇输出，相同簇在一起
	multimap<int,int> data_map;
	//簇信息
	set<int> set_label;
	// 按照簇类别保存 < 簇序号，样本序号>
	for(int i=0; i<NrObject; i++)
	{
		set_label.insert(label[i]);
		data_map.insert(make_pair(label[i],i));
	}

	// 重定义
	typedef multimap<int,int>::size_type sz_type;

	// 按照簇类别输出
	for( set<int>::iterator iter=set_label.begin(); iter != set_label.end(); iter++)
	{
		fout<<"Cluster "<<(*iter)+1<<endl;
		fout<<"Number"<<",";
		for( int i=0; i<nrf; i++)
			fout<<DimStr[i]<<",";
		fout<<endl;
		sz_type entries=data_map.count(*iter);
		multimap<int,int>::iterator miter=data_map.find(*iter);
		for(sz_type cnt=0; cnt != entries; ++cnt, ++miter)
		{
			int k=miter->second;
			fout<<k+1<<",";
			data[k].print(fout);
			fout<<endl;
		}
		fout<<endl;
	}
}

// 比较两个 Point 对象是否相同
// 相同返回 true, 否则返回 false
bool DATA::compare(const int lp, const int rp) const
{
	return !(data[lp] != data[rp] );
}

// 重载比较两个 Point 对象是否相同
// 相同返回 true, 否则返回 false
bool DATA::compare(const int lp, const int rp) 
{
	return !(data[lp] != data[rp] );
}

// 数据与中心点比较函数 (确保重复样本也可以执行)
// 如果数据与中心点组中某一个相同，返回 true, 否则返回 false
// ( 数据编号，中心点数组 , 簇个数)
bool DATA::compare(int nu, int *medoids, int nc)
{
	for( int i=0; i<nc; i++)
	{
		if( data[nu] == data[medoids[i]] )
			return true;
	}
	return false;
}

// 聚类完成后，总代价计算函数
double DATA::TotalCost()
{
	double sum=0.0;
	for( int i=0; i<NrObject; i++)
	{
		sum+=sqrt(DistSquare[i]);
	}
	return sum;
}


// 初始化中心点函数
//  ( 中心点指针，簇个数 )
void DATA::RandomMedoids( int *medoids, int nc)
{
	srand(unsigned(time(NULL)));
	medoids[0]=rand()%NrObject;
	int i=1,j;
	while( i<nc)
	{
		medoids[i]=rand()%NrObject;
		for(j=0; j<i; j++ )
		{
			// 比较是否重复
			if( compare(medoids[i],medoids[j]) )
			{
				break;
			}
		}
		// 如果不相同，i++
		if( j==i )
			i++;
	}
}

// 初始分派函数
// （中心点数组指针，簇个数）
void DATA::ClusterFirst( int *medoids, int nc)
{
	for(int i=0; i<NrObject; i++)
	{
		DistSquare[i]=0.0;
		for( int j=0; j<nc; j++)
		{
			if( j==0 )  // 初始化 label=0
			{
				DistSquare[i]=dist_square(i,medoids[j]);
				// 设置簇标号
				set_label(i,j);
			}
			else
			{
				double present_dist=dist_square(i,medoids[j]);
				if( present_dist < DistSquare[i] )
				{
					DistSquare[i]=present_dist;
					set_label(i,j);
				}
			}
		}
	}
}

// 中心点交换函数，返回交换代价
// ( 中心点数组，中心点被替换对象索引，替换点数据序号，簇个数）
double DATA::SwapOH(int *medoids, int oj, int oh, int nc)
{
	// oh 替换 oj
	double *DistSquareOhArray=new double[nc];
	double *DistSquareOiArray=new double[nc];
	bool *CosTheta=new bool[nc];
	double dist_jh=dist_square(medoids[oj],oh);  // d(O_j,O_j')^2
	for(int i=0; i<nc; i++)
	{
		double dist_ih=dist_square(medoids[i],oh);   // d(O_i, O_j')^2
		double dist_ij=dist_square(medoids[i],medoids[oj]);   // d(O_i,O_j)^2
		double dist_sqrt_ij=sqrt(dist_ij);     // d(O_i,O_j)
		double cos_d=(dist_ij+dist_jh-dist_ih)/(2.0*dist_sqrt_ij);  // cos * d(O_j,O_j')
		                                                     // 当 i=j 时，分母为0，出现非数
		                                                     // 但是，此种情形后面不会被调用
		DistSquareOhArray[i]=dist_ih/4.0;
		double dist=dist_sqrt_ij/2.0-cos_d;
		DistSquareOiArray[i]=dist*dist;
		CosTheta[i]=cos_d<0 ? true:false;
	}
	// 交换代价
	double TC=0.0;
	for( int i=0; i<NrObject; i++)
	{
		int class_label=get_label(i);   // i 所属簇的编号
		if( class_label != oj )  //  P belongs to the cluster not represented by O_j
		{
			if( CosTheta[class_label] )   // theta > pi/2
			{
				// 如果不满足条件5，需要计算交换代价，否则交换代价为0
				if( DistSquare[i] > DistSquareOiArray[class_label] )
				{
					double dist2=dist_square(i,oh);   // d(P, O_j')^2
					// d(P, O_i)^2 > d(P, O_j')^2  计算交换代价
					if( DistSquare[i] > dist2 )
						TC+=sqrt(dist2)-sqrt(DistSquare[i]);
				}
			}
			else     // theta <=pi/2
			{
				// 如果不满足条件4，需要计算交换代价，否则交换代价为0
				if( DistSquare[i] > DistSquareOhArray[class_label] )
				{
					double dist2=dist_square(i,oh);   // d(P,O_j')^2;
					// d(P,O_i)^2 > d(P, O_j')^2  计算交换代价
					if( DistSquare[i] > dist2 )
						TC+=sqrt(dist2)-sqrt(DistSquare[i]);
				}
			}
		}
		else    // P belongs to the cluster represented by O_j
		{
			// 此时无论属于 C_j 还是属于其他类别，交换代价都发生变化
			double dist2=dist_square(i,oh);   // d(P,O_j')^2
			// d(P,O_j')^2 > d(P,O_j)^2, 寻找第二小距离平方
			if( dist2 > DistSquare[i] )
			{
				// 第二小距离平方
				double djj2=SecondDist(medoids,i,class_label,nc);
				if( dist2 < djj2 )   // d(P,O_j') < d(P, O_l)
					TC+=sqrt(dist2)-sqrt(DistSquare[i]);
				else
					TC+=sqrt(djj2)-sqrt(DistSquare[i]);
			}
			else  // d(P,O_j') <= d(P,O_j)
			{
				TC+=sqrt(dist2)-sqrt(DistSquare[i]);
			}
		}
	}
	//释放动态内存
	delete [] DistSquareOhArray;
	delete [] DistSquareOiArray;
	delete [] CosTheta;

	return TC;
}

// 选出最优替换及被替换中心点，返回其序号
// 返回 pair<int, int> ,第一个为被替换中心点索引，第二个为替换中心点数据序号
// 当该值小于0时，说明没有更优替换点
// (中心点数组，簇个数）
pair<int, int> DATA::SelectO(int *medoids, int nc)
{
	int best_Oh=-1;        // 替换中心点数据编号
	int best_Oi=0;         // 被替换中心点索引
	double TC=0.0;         // 交换代价
	bool firsttime=true;       // 标记计算首次 TC
	for( int i=0; i<NrObject; i++ )
	{
		// 如果该样本就是中心点或者与中心点相同，则跳过
		if( !compare(i, medoids,nc) )
		{
			for(int j=0; j<nc; j++)
			{
				// 该迭代首次计算，需要先保存 TC 值，以备后续比较
				if(firsttime)
				{
					// 交换代价
					TC=SwapOH(medoids,j,i,nc);
					best_Oh=i;
					best_Oi=j;
					firsttime=false;
				}
				else
				{
					double newTC=SwapOH(medoids,j,i,nc);
					if( newTC < TC )
					{
						TC=newTC;
						best_Oh=i;
						best_Oi=j;
					}
				}
			}
		}
	}
	if( TC<0.0 )
		return make_pair(best_Oi,best_Oh);
	else
		return make_pair(-1,-1);
}

// 去除原簇中心点后，返回最小距离值
// (中心点组，计算点编号，所属簇索引,簇个数） 
double DATA::SecondDist( int *medoids, int oj, int oi, int nc)
{
	bool firsttime=true;
	double dji=0.0;
	for(int i=0; i<nc; i++)
	{
		if( i != oi )
		{
			if( firsttime )
			{
				dji=dist_square(medoids[i],oj);
				firsttime=false;
			}
			else
			{
				double newdji=dist_square(medoids[i],oj);
				if( newdji < dji )
					dji=newdji;
			}
		}
	}
	return dji;
}

// 重新分派聚类函数
// ( 更新后的中心点数组，中心点被替换对象索引，被替换中心点序号，簇个数）
void DATA::Cluster(int *medoids,int oi, int index,int nc)
{
	// medoids[oi] 替换 data[index]
	double *DistSquareOhArray=new double[nc];
	double *DistSquareOiArray=new double[nc];
	bool *CosTheta=new bool[nc];
	double dist_jh=dist_square(index,medoids[oi]);  // d(O_j,O_j')^2
	for(int i=0; i<nc; i++)
	{
		double dist_ih=dist_square(medoids[i],medoids[oi]);   // d(O_i, O_j')^2
		double dist_ij=dist_square(medoids[i],index);   // d(O_i,O_j)^2
		double dist_sqrt_ij=sqrt(dist_ij);     // d(O_i,O_j)
		double cos_d=(dist_ij+dist_jh-dist_ih)/(2.0*dist_sqrt_ij);  // cos * d(O_j,O_j')
		                                                     // 当 i=j 时，分母为0，出现非数
		                                                     // 但是，此种情形后面不会被调用
		DistSquareOhArray[i]=dist_ih/4.0;
		double dist=dist_sqrt_ij/2.0-cos_d;
		DistSquareOiArray[i]=dist*dist;
		CosTheta[i]=cos_d<-0.5 ? true:false;
	}
	// 分派簇类别/簇编号
	for( int i=0; i<NrObject; i++)
	{
		int class_label=get_label(i);   // i 所属簇的编号
		if( class_label != oi )  //  P belongs to the cluster not represented by O_j
		{
			if( CosTheta[class_label] )   // theta > pi/2
			{
				// 如果满足条件5，P类别编号无需改变，否则进一步判断
				if( DistSquare[i] > DistSquareOiArray[class_label] )
				{
					double dist2=dist_square(i,medoids[oi]);   // d(P, O_j')^2
					// d(P, O_i)^2 > d(P, O_j')^2  更新簇标号及对应距离平方
					if( DistSquare[i] > dist2 )
					{
						set_label(i,oi);
						set_DistSquare(i,dist2);
					}
				}
			}
			else     // theta <=pi/2
			{
				// 如果满足条件4，P类别编号无需改变，否则进一步判断
				if( DistSquare[i] > DistSquareOhArray[class_label] )
				{
					double dist2=dist_square(i,medoids[oi]);   // d(P,O_j')^2;
					// d(P,O_i)^2 > d(P, O_j')^2  更新簇标号及对应的距离平方
					if( DistSquare[i] > dist2 )
					{
						set_label(i,oi);
						set_DistSquare(i,dist2);
					}
				}
			}
		}
		else    // P belongs to the cluster represented by O_j
		{
			double dist2=dist_square(i,medoids[oi]);   // d(P,O_j')^2
			// d(P,O_j')^2 > d(P,O_j)^2, 寻找第二小距离平方及对应的簇标号
			if( dist2 > DistSquare[i] )
			{
				// 最小距离平方及对应簇标号
				pair<int,double> pairdistlabel=MinDistLabel(medoids,i,nc);
				set_DistSquare(i,pairdistlabel.second);
				set_label(i,pairdistlabel.first);
			}
			else  // d(P,O_j') <= d(P,O_j) 只需更新距离平方，簇标号无需更新
			{
				set_DistSquare(i,dist2);
			}
		}
	}
	//释放动态内存
	delete [] DistSquareOhArray;
	delete [] DistSquareOiArray;
	delete [] CosTheta;
}

// 计算最小距离平方及对应簇标号
// 返回值第一个为最小值簇标号，第二只最小值时距离平方
//(中心点数组，待分派样本序号，簇个数）
pair<int,double> DATA::MinDistLabel(int *medoids,int nu,int nc)
{
	double dist=0.0;
	int newlabel=0;
	bool firsttime=true;
	for( int i=0; i<nc; i++)
	{
		if( firsttime )
		{
			dist=dist_square(medoids[i],nu);
			newlabel=i;
			firsttime=false;
		}
		else
		{
			double newdist=dist_square(medoids[i],nu);
			if( newdist < dist )
			{
				dist=newdist;
				newlabel=i;
			}
		}
	}
	return make_pair(newlabel,dist);
}






// 类非成员函数定义

// Point 对象关系操作符重载函数
//由于需要使用Point 私有成员，设置成 Point 类的友元函数
bool operator!=( const Point &p1, const Point &p2)
{
	for(int i=0; i<p1.NrFeatures; i++)
	{
		if( p1.p[i] != p2.p[i] )
			return true;
	}
	return false;
}
// Point 对象关系操作符重载函数
//由于需要使用Point 私有成员，设置成 Point 类的友元函数
bool operator==( const Point &p1, const Point &p2)
{
	return !(p1 != p2);
}




#endif