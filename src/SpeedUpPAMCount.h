// Speed Up PAM �����㷨
// �������д��2015��11��29��
// ʵ�� k-medoids �㷨�������㷨 PAM 
// �ο����ף�
// �ָ壬2015��11��26�� �� PAM �㷨

#ifndef SPEEDUPPAMCOUNT_H
#define SPEEDUPPAMCOUNT_H

#include<iostream>
#include<fstream>
#include<sstream>   //�ַ�����
#include<ctime>     //����ʱ��
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<map>
#include<set>

using namespace std;

static unsigned long int CountNum=0;

//���ǰ������
class Point;      //��������
class DATA;       // ���ݼ���

//�����Ա��������
// ����ʱ����
double running_time(const clock_t &, const clock_t &);
//�������ļ�������
ifstream& open_file(ifstream &, const string &);
//������ļ�������
ofstream& open_outfile(ofstream &, const string &);
// Point �����ϵ���������غ���
//������Ҫʹ��Point ˽�г�Ա�����ó� Point �����Ԫ����
inline bool operator!=( const Point &, const Point &);
// Point �����ϵ���������غ���
inline bool operator==( const Point &, const Point &);
//�ַ���ת������
// ( Ŀ���ַ��������滻�ַ����ַ������滻�ַ����ַ�����
void string_replace(string &, const string &,const string &);
//��ȡ�Ӷ��ſ�ʼ��ĩβ������
//( Ŀ���ַ�����
int SubStringToInt(string&);
// ��ȡ�ӿ�ʼ�����ǰ���ַ���
//( Ŀ���ַ�����
string SubString(const string&);



// �����Ա��������

// ����ʱ����
double running_time(const clock_t &s, const clock_t &e)
{
	return static_cast<double>(e-s)/CLOCKS_PER_SEC;
}

//�������ļ�������
ifstream& open_file(ifstream &in, const string &file)
{
	in.close();
	in.clear();
	in.open(file.c_str());
	return in;
}

//������ļ�������
ofstream& open_outfile(ofstream &out, const string &file)
{
	out.close();
	out.clear();
	out.open(file.c_str(), ofstream::trunc);   //���ԭ���ļ�
	return out;
}

//�ַ���ת������
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

//��ȡ�Ӷ��ſ�ʼ��ĩβ��������
//( Ŀ���ַ�����
int SubStringToInt(string &str)
{
	string::size_type pos_start=str.find_first_of(",",0);
	++pos_start;
	string::size_type pos_end=str.find_first_of(",",pos_start);
	// ������������֮���������ȡ����
	string st=str.substr(pos_start,pos_end-pos_start);
	istringstream istr(st);
	int integer=0;
	istr>>integer;
	return integer;
}

// ��ȡ�ӿ�ʼ�����ǰ���ַ���
//( Ŀ���ַ�����
string SubString(const string &str)
{
	return str.substr(0,str.size()-4);
}



// Point �ඨ��
class Point
{
	//��Ԫ����
	friend bool operator!=( const Point &, const Point &);

private:
	int NrFeatures;     //��������/����ά��
	double *p;          // ������ָ�룬ָ��̬����
	void del_p()
	{
		if( p != NULL )
			delete [] p;
	}

public:
	// Ĭ�Ϲ��캯��
	Point():NrFeatures(0),p(NULL){}
	// ���ι��캯��
	Point( istringstream&, int );   
	// ���ƹ��캯��
	Point(const Point&);
	// ��ֵ������
	Point& operator=(const Point &);
	// ��������
	~Point()
	{
		del_p();
	}
	// ��ӡ����
	void print(ostream&);
	//��ӡ����
	void print(ofstream&);
	// ������������
	int get_NrF()
	{
		return NrFeatures;
	}
	// ������㺯��
	double dist_square( Point &);
};

// ���ƹ��캯������
// �����ֵ��ָ�븴��
Point::Point( const Point &point):p(NULL),NrFeatures(point.NrFeatures)
{
	//cout<<"Point ���ƹ��캯�����С���"<<endl;
	if ( point.p != NULL)
	{
		p=new double[NrFeatures];
		for( int i=0; i<NrFeatures; i++)
			p[i]=point.p[i];
	}
}
// ��ֵ������
Point& Point::operator=(const Point &point)
{
	//cout<<"Point ��ֵ�������С���"<<endl;
	//�ͷ���߶���
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

// ���ι��캯������
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

//��ӡ����
void Point::print(ostream &os)
{
	for(int i=0; i<NrFeatures; i++)
		os<<p[i]<<"\t";
}
//��ӡ����
void Point::print(ofstream &fout)
{
	for(int i=0; i<NrFeatures; i++)
		fout<<p[i]<<",";
}

// ������㺯��
double Point::dist_square( Point &rp)
{
	CountNum++;
	double sum=0.0;
	for(int i=0; i<NrFeatures; i++)
		sum+=(p[i]-rp.p[i])*(p[i]-rp.p[i]);
	return sum;
}



// ���ݼ��ඨ��
class DATA
{
private:
	int NrObject;    // ���ݼ���С
	string *DimStr;   // ��������
	Point *data;     // ���ݼ�
	int *label;      // ��������ǣ���������ֵΪ�ر��
	double *DistSquare;    // ����������ĵ����ƽ��
	DATA(const DATA&);     //���ƹ��캯����ֻ�������������ڽ�ֹ����

public:
	DATA(const string&, int &, int &, int &, int &);   //���캯��
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
	//��ӡ����
	void print();
	// ���ش�ӡ����
	void print(ostream &, int);
	// ���ĵ��ļ������
	void print(ofstream &, int);
	// ��ӡ�����
	void print_label();
	// �ļ����������
	void print_label(ofstream&);
	// �Ƚ����� Point �����Ƿ���ͬ
	bool compare(const int, const int) const;
	// �Ƚ����� Point �����Ƿ���ͬ
	bool compare(const int, const int);
	// ���������ĵ�ȽϺ��� (ȷ���ظ�����Ҳ����ִ��)
	// ������������ĵ�����ĳһ����ͬ������ true, ���򷵻� false
	// ( ���ݱ�ţ����ĵ����� , �ظ���)
	bool compare(int , int *, int );
	// �������ݼ���С
	int get_nro()
	{
		return NrObject;
	}
	// �������� Point ����ľ���ƽ��
	double dist_square(int num1, int num2)
	{
		return data[num1].dist_square(data[num2]);

	}
	// ���ôر��
	void set_label(int index, int n_label)
	{
		label[index]=n_label;
	}
	// ���þ���ƽ��
	void set_DistSquare(int index, double dist)
	{
		DistSquare[index]=dist;
	}
	// ��ȡ�ر��
	int get_label(int index)
	{
		return label[index];
	}
	// ������ɺ��ܴ��ۼ��㺯��
	double TotalCost();
	//��ʼ�����ĵ�
	// ( ���ĵ����飬�ظ���)
	void RandomMedoids(int *, int nc);
	// ��ʼ���ɺ���
	// ( ���ĵ�����ָ�룬�ظ�����
	void ClusterFirst(int *, int);
	// ��������, ���ؽ�������
	double SwapOH(int*,int,int,int);
	// ѡ�������滻�����滻���ĵ㣬���������
	pair<int, int> SelectO(int*, int);
	// ȥ��ԭ�����ĵ�󣬷��صڶ�С����ֵƽ��
	// ( ���ĵ��飬������ţ�����������,�ظ����� 
	double SecondDist( int *, int , int , int );
	// ���·��ɾ��ຯ��
	// ( ���º�����ĵ����飬���ĵ㱻�滻�������������滻���ĵ���ţ��ظ�����
	void Cluster(int *,int , int ,int );
	// ������С����ƽ������Ӧ�ر��
	// ����ֵ��һ��Ϊ��Сֵ�ر�ţ��ڶ�ֻ��Сֵʱ����ƽ��
	//(���ĵ����飬������������ţ��ظ�����
	pair<int,double> MinDistLabel(int*,int,int);
};

// DATA ���Ա��������
// ���캯������
// �βα� ( �ļ��������ݼ���С����������, �ظ���������������)
DATA::DATA(const string &file, int &nro, int &nrf, int &nc, int &maxiter):NrObject(nro),
	DimStr(NULL),data(NULL),label(NULL),DistSquare(NULL)
{
	// �ļ���
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
	
	// ������������
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
	//�ر��ļ���
	in.close();
	in.clear();

	return;
}

//��ӡ����
void DATA::print()
{
	int nrf=data[0].get_NrF();  // ��������
	cout<<"���"<<"\t";
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
// ���ش�ӡ����
// ( ���������Ҫ��ӡ�������)
void DATA::print(ostream &os, int num)
{
	data[num].print(os);
}
// ���ĵ��ļ������
void DATA::print(ofstream &fout, int num)
{
	data[num].print(fout);
}
// ��ӡ�����
void DATA::print_label()
{
	int nrf=data[0].get_NrF();   //��������
	cout<<"���"<<"\t";
	for(int i=0; i<nrf; i++)
		cout<<DimStr[i]<<"\t";
	cout<<"�����"<<endl;
	for(int i=0; i<NrObject; i++)
	{
		cout<<i+1<<"\t";
		data[i].print(cout);
		cout<<label[i]+1<<endl;
	}
}

// �ļ����������
void DATA::print_label( ofstream &fout)
{
	int nrf=data[0].get_NrF();   //��������

	//���ݰ����������ͬ����һ��
	multimap<int,int> data_map;
	//����Ϣ
	set<int> set_label;
	// ���մ���𱣴� < ����ţ��������>
	for(int i=0; i<NrObject; i++)
	{
		set_label.insert(label[i]);
		data_map.insert(make_pair(label[i],i));
	}

	// �ض���
	typedef multimap<int,int>::size_type sz_type;

	// ���մ�������
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

// �Ƚ����� Point �����Ƿ���ͬ
// ��ͬ���� true, ���򷵻� false
bool DATA::compare(const int lp, const int rp) const
{
	return !(data[lp] != data[rp] );
}

// ���رȽ����� Point �����Ƿ���ͬ
// ��ͬ���� true, ���򷵻� false
bool DATA::compare(const int lp, const int rp) 
{
	return !(data[lp] != data[rp] );
}

// ���������ĵ�ȽϺ��� (ȷ���ظ�����Ҳ����ִ��)
// ������������ĵ�����ĳһ����ͬ������ true, ���򷵻� false
// ( ���ݱ�ţ����ĵ����� , �ظ���)
bool DATA::compare(int nu, int *medoids, int nc)
{
	for( int i=0; i<nc; i++)
	{
		if( data[nu] == data[medoids[i]] )
			return true;
	}
	return false;
}

// ������ɺ��ܴ��ۼ��㺯��
double DATA::TotalCost()
{
	double sum=0.0;
	for( int i=0; i<NrObject; i++)
	{
		sum+=sqrt(DistSquare[i]);
	}
	return sum;
}


// ��ʼ�����ĵ㺯��
//  ( ���ĵ�ָ�룬�ظ��� )
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
			// �Ƚ��Ƿ��ظ�
			if( compare(medoids[i],medoids[j]) )
			{
				break;
			}
		}
		// �������ͬ��i++
		if( j==i )
			i++;
	}
}

// ��ʼ���ɺ���
// �����ĵ�����ָ�룬�ظ�����
void DATA::ClusterFirst( int *medoids, int nc)
{
	for(int i=0; i<NrObject; i++)
	{
		DistSquare[i]=0.0;
		for( int j=0; j<nc; j++)
		{
			if( j==0 )  // ��ʼ�� label=0
			{
				DistSquare[i]=dist_square(i,medoids[j]);
				// ���ôر��
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

// ���ĵ㽻�����������ؽ�������
// ( ���ĵ����飬���ĵ㱻�滻�����������滻��������ţ��ظ�����
double DATA::SwapOH(int *medoids, int oj, int oh, int nc)
{
	// oh �滻 oj
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
		                                                     // �� i=j ʱ����ĸΪ0�����ַ���
		                                                     // ���ǣ��������κ��治�ᱻ����
		DistSquareOhArray[i]=dist_ih/4.0;
		double dist=dist_sqrt_ij/2.0-cos_d;
		DistSquareOiArray[i]=dist*dist;
		CosTheta[i]=cos_d<0 ? true:false;
	}
	// ��������
	double TC=0.0;
	for( int i=0; i<NrObject; i++)
	{
		int class_label=get_label(i);   // i �����صı��
		if( class_label != oj )  //  P belongs to the cluster not represented by O_j
		{
			if( CosTheta[class_label] )   // theta > pi/2
			{
				// �������������5����Ҫ���㽻�����ۣ����򽻻�����Ϊ0
				if( DistSquare[i] > DistSquareOiArray[class_label] )
				{
					double dist2=dist_square(i,oh);   // d(P, O_j')^2
					// d(P, O_i)^2 > d(P, O_j')^2  ���㽻������
					if( DistSquare[i] > dist2 )
						TC+=sqrt(dist2)-sqrt(DistSquare[i]);
				}
			}
			else     // theta <=pi/2
			{
				// �������������4����Ҫ���㽻�����ۣ����򽻻�����Ϊ0
				if( DistSquare[i] > DistSquareOhArray[class_label] )
				{
					double dist2=dist_square(i,oh);   // d(P,O_j')^2;
					// d(P,O_i)^2 > d(P, O_j')^2  ���㽻������
					if( DistSquare[i] > dist2 )
						TC+=sqrt(dist2)-sqrt(DistSquare[i]);
				}
			}
		}
		else    // P belongs to the cluster represented by O_j
		{
			// ��ʱ�������� C_j ��������������𣬽������۶������仯
			double dist2=dist_square(i,oh);   // d(P,O_j')^2
			// d(P,O_j')^2 > d(P,O_j)^2, Ѱ�ҵڶ�С����ƽ��
			if( dist2 > DistSquare[i] )
			{
				// �ڶ�С����ƽ��
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
	//�ͷŶ�̬�ڴ�
	delete [] DistSquareOhArray;
	delete [] DistSquareOiArray;
	delete [] CosTheta;

	return TC;
}

// ѡ�������滻�����滻���ĵ㣬���������
// ���� pair<int, int> ,��һ��Ϊ���滻���ĵ��������ڶ���Ϊ�滻���ĵ��������
// ����ֵС��0ʱ��˵��û�и����滻��
// (���ĵ����飬�ظ�����
pair<int, int> DATA::SelectO(int *medoids, int nc)
{
	int best_Oh=-1;        // �滻���ĵ����ݱ��
	int best_Oi=0;         // ���滻���ĵ�����
	double TC=0.0;         // ��������
	bool firsttime=true;       // ��Ǽ����״� TC
	for( int i=0; i<NrObject; i++ )
	{
		// ����������������ĵ���������ĵ���ͬ��������
		if( !compare(i, medoids,nc) )
		{
			for(int j=0; j<nc; j++)
			{
				// �õ����״μ��㣬��Ҫ�ȱ��� TC ֵ���Ա������Ƚ�
				if(firsttime)
				{
					// ��������
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

// ȥ��ԭ�����ĵ�󣬷�����С����ֵ
// (���ĵ��飬������ţ�����������,�ظ����� 
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

// ���·��ɾ��ຯ��
// ( ���º�����ĵ����飬���ĵ㱻�滻�������������滻���ĵ���ţ��ظ�����
void DATA::Cluster(int *medoids,int oi, int index,int nc)
{
	// medoids[oi] �滻 data[index]
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
		                                                     // �� i=j ʱ����ĸΪ0�����ַ���
		                                                     // ���ǣ��������κ��治�ᱻ����
		DistSquareOhArray[i]=dist_ih/4.0;
		double dist=dist_sqrt_ij/2.0-cos_d;
		DistSquareOiArray[i]=dist*dist;
		CosTheta[i]=cos_d<-0.5 ? true:false;
	}
	// ���ɴ����/�ر��
	for( int i=0; i<NrObject; i++)
	{
		int class_label=get_label(i);   // i �����صı��
		if( class_label != oi )  //  P belongs to the cluster not represented by O_j
		{
			if( CosTheta[class_label] )   // theta > pi/2
			{
				// �����������5��P���������ı䣬�����һ���ж�
				if( DistSquare[i] > DistSquareOiArray[class_label] )
				{
					double dist2=dist_square(i,medoids[oi]);   // d(P, O_j')^2
					// d(P, O_i)^2 > d(P, O_j')^2  ���´ر�ż���Ӧ����ƽ��
					if( DistSquare[i] > dist2 )
					{
						set_label(i,oi);
						set_DistSquare(i,dist2);
					}
				}
			}
			else     // theta <=pi/2
			{
				// �����������4��P���������ı䣬�����һ���ж�
				if( DistSquare[i] > DistSquareOhArray[class_label] )
				{
					double dist2=dist_square(i,medoids[oi]);   // d(P,O_j')^2;
					// d(P,O_i)^2 > d(P, O_j')^2  ���´ر�ż���Ӧ�ľ���ƽ��
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
			// d(P,O_j')^2 > d(P,O_j)^2, Ѱ�ҵڶ�С����ƽ������Ӧ�Ĵر��
			if( dist2 > DistSquare[i] )
			{
				// ��С����ƽ������Ӧ�ر��
				pair<int,double> pairdistlabel=MinDistLabel(medoids,i,nc);
				set_DistSquare(i,pairdistlabel.second);
				set_label(i,pairdistlabel.first);
			}
			else  // d(P,O_j') <= d(P,O_j) ֻ����¾���ƽ�����ر���������
			{
				set_DistSquare(i,dist2);
			}
		}
	}
	//�ͷŶ�̬�ڴ�
	delete [] DistSquareOhArray;
	delete [] DistSquareOiArray;
	delete [] CosTheta;
}

// ������С����ƽ������Ӧ�ر��
// ����ֵ��һ��Ϊ��Сֵ�ر�ţ��ڶ�ֻ��Сֵʱ����ƽ��
//(���ĵ����飬������������ţ��ظ�����
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






// ��ǳ�Ա��������

// Point �����ϵ���������غ���
//������Ҫʹ��Point ˽�г�Ա�����ó� Point �����Ԫ����
bool operator!=( const Point &p1, const Point &p2)
{
	for(int i=0; i<p1.NrFeatures; i++)
	{
		if( p1.p[i] != p2.p[i] )
			return true;
	}
	return false;
}
// Point �����ϵ���������غ���
//������Ҫʹ��Point ˽�г�Ա�����ó� Point �����Ԫ����
bool operator==( const Point &p1, const Point &p2)
{
	return !(p1 != p2);
}




#endif