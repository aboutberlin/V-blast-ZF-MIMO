#include <cstdlib>
#include <limits>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

double uniform_rand(double lowerBndr, double upperBndr)	//用来产生标准正态分布（上网找的），原理不懂
{
	return lowerBndr + ((double) rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
}

double GGSN(double mean, double sigma)	//产生标准正态分布，mean平均值，sigma方差
{
	double x, y, r2;
	do 
	{
		x = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
		y = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
		r2 = x * x + y * y;
	} while (r2 > 1.0 || r2 == 0.0);
	return mean + sigma * y * sqrt(-2.0 * log(r2) / r2);
}


	/***复数类***/
class Complex
{
public:
	Complex(double r=0,double i=0):real(r),imag(i){}	//带默认参数的构造函数，默认为0+0i
	Complex daoshu(){return Complex(real/(real*real+imag*imag),-imag/(real*real+imag*imag));} //复数的倒数
	Complex operator+ (const Complex &c2){return Complex (real+c2.real,imag+c2.imag);}		  //复数的加法
	Complex operator- (const Complex &c2){return Complex (real-c2.real,imag-c2.imag);}		  //复数的减法
	Complex operator* (const Complex &c2){return Complex (real*c2.real-imag*c2.imag,imag*c2.real+real*c2.imag);} //复数的乘法
	Complex operator/ (const Complex &c2){return Complex ((real*c2.real+imag*c2.imag)/(c2.real*c2.real+c2.imag*c2.imag),
		(imag*c2.real-real*c2.imag)/(c2.real*c2.real+c2.imag*c2.imag));}					  //复数的除法
	Complex operator! (){return Complex(real,-imag);}					 //复数的共轭
	Complex operator* (const double f){return Complex (f*real,f*imag);}	 //复数乘常数
	double mochang(){return sqrt(real*real+imag*imag);}					 //复数的模长
	void change(double c,double i){real=c;imag=i;}						 //改变复数的实部和虚部
	double r_find(){return real;}										 //返回复数的实部
	double i_find(){return imag;}										 //返回复数的虚部

	bool operator != (const Complex&c)									 //判断两个复数是否相等
	{
		if (real==c.real && imag==c.imag)
			return false;
		return true;
	}

			/***定义复数的输出***/
	friend ostream& operator << (ostream& output,Complex c)	
	{
		output << "(" << setw(13) << c.real;  //先输出("实部"
		if (c.imag>=0)		 //判断虚部正负，输出"+"或"-"
			output << "+";
		else {output << "-"; }
		output << setw(13) << fabs(c.imag) << "i)" << "   ";  //最后输出"|虚部|i")
		return output;
	}

			/***定义复数的输入（没用到）***/
	friend istream& operator >> (istream& input,Complex& c)	
	{
		input >> c.real >> c.imag;
		return input;
	}

private:
	double real,imag; //私有数据：复数的实部和虚部
};

		/***矩阵类***/
class Matrix
{
private:
	int row,col;		//矩阵的行和列
	Complex p[20][20];  //矩阵的元素为复数
public:
	Matrix(int r=1,int c=1):row(r),col(c){}		//带默认参数的构造函数，默认为1*1矩阵

	Matrix operator* (const Matrix &);			//矩阵乘法
	Matrix operator- (const Matrix &);			//矩阵减法	
	Matrix operator+ (const Matrix &);			//矩阵加法
	Matrix modified(Matrix &);					//用改进格拉姆-施密特正交化对矩阵进行QR分解
	Complex qietiao(int ,int);					//QPSK解调
	Matrix outl(Matrix,int);					//消去矩阵的某一列
	Matrix outr(Matrix,int);					//消去矩阵的某一行
	Matrix ni_LM(Matrix);						//计算下三角矩阵的逆矩阵
	Matrix ni_UM(Matrix);						//计算上三角矩阵的逆矩阵
	double SQRD(Matrix,Matrix,Matrix);			//SQRD算法
	Matrix ZF(Matrix ,Matrix ,Matrix );		    //ZF算法
	Matrix MMSE(Matrix ,Matrix ,Matrix );		//MMSE算法
	Matrix MMSE_QR(Matrix ,Matrix ,Matrix );	//MMSE-QR算法

	Matrix operator* (const double f)			//矩阵乘常数
	{
		Matrix m(row,col);
		for (int i=0;i<row;i++)
			for (int j=0;j<col;j++)
				m.p[i][j]=p[i][j]*f;
		return m;
	}

	Matrix operator* (const Complex f)			//矩阵乘复数
	{
		Matrix m(row,col);
		for (int i=0;i<row;i++)
			for (int j=0;j<col;j++)
				m.p[i][j]=p[i][j]*f;
		return m;
	}

	void Transmit()								//产生发射信号c
	{
		for (int i=0;i<row;i++)
			for (int j=0;j<col;j++)
			{
				if (GGSN(0,1)>=0)
					p[i][j].change(pow(-1,int(rand())),0);
				else
					p[i][j].change(0,pow(-1,int(rand())));
			}
	}

	double Two_norm();								//计算向量的内积
	Matrix slice(int ,int ,int ,int );				//相当于python的切片
	void str_m(int ,int ,int ,int ,Matrix );		//将一个小矩阵赋值给一个大矩阵的某一个位置
	Matrix swap(Matrix ,int , int );				//交换矩阵某两列的向量
	void back(Matrix&,Matrix,Matrix );				//后向代入法（改）
	double QR(Matrix ,Matrix ,Matrix );				//QR算法
	Matrix MMSE_SQRD(Matrix H,Matrix x,Matrix c);	//MMSE-SQRD算法
	void GaoSi();									//产生信道矩阵
	double V_Blast(Matrix,Matrix,Matrix);			//V-Blast算法
	Matrix ML(Matrix ,Matrix , Matrix );


			/***定义矩阵的转置***/
	Matrix Matrix::operator !()
	{
		Matrix m(col,row);
		for (int i = 0; i < col; i++)//第二层循环：将所有行转置对应列
			for (int j = 0; j < row; j++)//第一层循环：将某一行转置成对应列
				m.p[i][j] = !p[j][i];
		return m;
	}

			/***定义矩阵输入***/
	friend istream& operator >> (istream& input,Matrix& m)  
	{
	int i=0,j=0;
	for (i=0;i<m.row;i++)
		for (j=0;j<m.col;j++)
		{
			cout << "Input number in the " << i+1 << " row and the " << j+1 << " column: ";
			input >> m.p[i][j];  //使用复数重载的提取运算符
		}
	cout << endl;
	return input;
	}

			/***定义矩阵输出***/
	friend ostream& operator << (ostream& output , Matrix m)
	{
		for (int i=0;i<m.row;i++)
		{
			for (int j=0;j<m.col;j++)
				output << m.p[i][j];
			output << endl;
		}
		output << endl;
		return output;
	}
};

			/***消去某一列***/
Matrix Matrix::outl(Matrix H,int c)		//形参列表：H为待消去矩阵，c表示消去第c列向量
										//返回：消去第c列后的矩阵
{
	Matrix H_(row,col-1);
	H_.str_m(0,row-1,0,c-1,H.slice(0,row-1,0,c-1));		   //将原矩阵第0列到第c-1列的矩阵利用str_m函数赋给返回矩阵
	H_.str_m(0,row-1,c,col-1,H.slice(0,row-1,c+1,col-1));  //将原矩阵第c+1列到最后一列的矩阵利用str_m函数赋给返回矩阵
	return H_;
}

			/***消去某一行***/
Matrix Matrix::outr(Matrix H,int r)		//形参列表：H为待消去矩阵，r表示消去第r行向量
										//返回：消去第r行后的矩阵
{
	H=!H;					//先转置
	H=H.outl(H,r);			//然后利用消去某一列的函数
	return (!H);			//再转置回来
}

			/***QPSK解调***/
Complex Matrix::qietiao(int r,int col)	//形参列表：r，col表示解调矩阵的第r行第col列的元素
										//返回：解调后的坐标
{
	Complex c[4];			//定义星座图上的四个点
	c[0].change(1,0);
	c[1].change(0,1);
	c[2].change(-1,0);
	c[3].change(0,-1);
			
		/***判断离哪一个点最近***/
	double temp=(p[r][col]-c[0]).mochang();int n=0;
	for (int i=1;i<4;i++)
		if (temp>(p[r][col]-c[i]).mochang())
		{
			temp=(p[r][col]-c[i]).mochang();
			n=i;
		}
	return c[n];		//返回距离最近的点的坐标（即解调后的坐标）
}

			/***定义矩阵乘法***/
Matrix Matrix::operator *(const Matrix &m2)
{
	int i,j;
	Matrix m(row,m2.col);
	if (col==m2.row)  //判断能否相乘
	{
		for (i = 0; i<row; i++)//第三层循环：计算相乘后矩阵所有的元素
			for (j = 0; j<m2.col;j++)//第二层循环：计算相乘后矩阵某一行的元素
				for (int k = 0; k<col;k++)//第一层循环：计算相乘后矩阵某一行某一列的元素
					m.p[i][j] = m.p[i][j] + p[i][k] * m2.p[k][j];
		return m;
	}
	cout << "Two matrices can not be multiplied!" << endl;  //不能相乘告知user
	return m;
}
		
		/***定义矩阵减法***/
Matrix Matrix::operator -(const Matrix &m2)  //与矩阵加法同理
{
	Matrix m(row,col);
	if (row==m2.row && col==m2.col)
	{
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				m.p[i][j]=p[i][j]-m2.p[i][j];
		return m;
	}
	cout << "Two matrices can not be subtracted!" << endl;
	return m;
}

		/***定义矩阵加法***/
Matrix Matrix::operator +(const Matrix &m2)
{
	Matrix m(row,col);
	if (row==m2.row && col==m2.col)  //判断两个矩阵能否进行加法运算
	{
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				m.p[i][j]=p[i][j]+m2.p[i][j];  //使用复数"+"的重载
		return m;
	}
	cout << "Two matrices can not be added!" << endl;  //不能相加跳出警告
	return m;
}

		/***产生符合标准正态分布的矩阵***/
void Matrix::GaoSi()
{
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			p[i][j].change(GGSN(0,1),GGSN(0,1));		//用复数的change函数和产生标准正态分布的GGSN函数
}

		/***后向代入法（改）***/
void Matrix::back(Matrix&x,Matrix A,Matrix b)	//x为接受矩阵，A为上三角矩阵，b为解向量
{
	for (int k=0;k<A.row;k++)					//A需要满足行数多于列数，否则运行失败（肯定满足的）
	{
		x.p[b.row-1-k][0]=b.p[b.row-1-k][0]/A.p[A.col-1-k][A.col-1-k];	//计算解的向量
		x.p[b.row-1-k][0]=x.qietiao(b.row-1-k,0);						//进行解调
		if (k!=A.col-1)
			for (int i=b.row-1-k;i<b.row;i++)		
				b.p[b.row-2-k][0]=b.p[b.row-2-k][0]-x.p[i][0]*A.p[b.row-2-k][i];	//修正右端向量
	}
}

		/***改进格拉姆-施密特正交化的QR分解***/
Matrix Matrix::modified(Matrix &CGS)	//CGS为需要进行QR分解的初始矩阵且为值传递
										//返回R
{
	Complex rjk;						//rjk表示QR分解中R矩阵第j行第k列的元素
	Matrix R(CGS.col,CGS.col);			//QR分解中R矩阵
	Matrix qk(CGS.row,1),b(CGS.row,1);	//qk表示需要CGS中需要正交的列向量，b表示qk之后的每一个列向量
	for (int i=0;i<CGS.col;i++)
	{
		qk=CGS.slice(0,CGS.row-1,i,i);		//用切片得到第i列未正交的列向量
		rjk.change(sqrt(qk.Two_norm()),0);	//计算模长
		qk=qk*(1/rjk.r_find());				//单位化
		R.p[i][i]=rjk;						
		CGS.str_m(0,CGS.row,i,i,qk);
			/***将后面所有列减去前一列的分量***/
		for (int j=i+1;j<CGS.col;j++)
		{
			b=CGS.slice(0,CGS.row-1,j,j);	//求解b和rjk
			rjk=(!qk*b).p[0][0];
			R.p[i][j]=rjk;
			b=b-qk*rjk;
			CGS.str_m(0,CGS.row,j,j,b);		//CGS最后变成Q
		}
	}
	return R;
}

			/***QR算法***/
double Matrix::QR(Matrix H,Matrix x,Matrix c)	//形参列表：H为信道矩阵，x为接受矩阵，c为发射矩阵
												//返回：误码率
{
	double n=double(2*c.row),w=0.0;				//n表示总共发射的比特数，w用来计数错误的比特数，两者相除为误码率
	Matrix R(H.col,H.col),y(H.col,1);			//R为QR分解后的R矩阵
	R=H.modified(H);							//使用modified函数后H变为Q
	y=(!H)*x;											
	x.back(x,R,y);								//根据QR算法得到x
		/***计算误码率***/
	for (int i=0;i<c.row;i++)
		if (c.p[i][0]!=x.p[i][0])
		{
			if (c.p[i][0].r_find()!=x.p[i][0].r_find())
				w=w+1.00;
			if (c.p[i][0].i_find()!=x.p[i][0].i_find())
				w=w+1.00;
		}
	return w/n;
}

		/***交换矩阵某两列***/
Matrix Matrix::swap(Matrix m,int c1, int c2)	//m表示待交换的矩阵，c1和c2表示待交换的列的位置
												//返回交换后的矩阵
{
	Matrix one=m.slice(0,m.row-1,c1,c1);
	m.str_m(0,m.row-1,c1,c1,m.slice(0,m.row-1,c2,c2));
	m.str_m(0,m.row-1,c2,c2,one);
	return m;
}

			/***SQRD算法***/
double Matrix::SQRD(Matrix CGS,Matrix x,Matrix c)	//形参列表：CGS为信道矩阵，x为接受矩阵，c为发射矩阵
													//返回：误码率
{
//	cout << CGS;
	Complex rjk;		
	Matrix R(CGS.col,CGS.col),S(CGS.col,CGS.col);	//S为单位矩阵
	Matrix qk(CGS.row,1),b(CGS.row,1);				//rjk，qk，b，R的含义和QR算法中的完全一样	

	for (int e=0;e<CGS.col;e++)
		S.p[e][e].change(1,0);						//将S变为单位矩阵

	for (int k=0;k<CGS.col;k++)						
	{	
			/***首先进行排序***/
		int f=k;
		double temp=CGS.slice(0,CGS.row-1,k,k).Two_norm();
		for (int v=k+1;v<CGS.col;v++)
		{
			if (temp>CGS.slice(0,CGS.row-1,v,v).Two_norm())
			{
				f=v;
				temp=CGS.slice(0,CGS.row-1,v,v).Two_norm();
			}
		}
		CGS=CGS.swap(CGS,k,f);
		R=R.swap(R,k,f);
		S=S.swap(S,k,f);

				/***剩下步骤和QR算法一样***/
		qk=CGS.slice(0,CGS.row-1,k,k);			//切片得到第i列未正交的列向量
		rjk.change(sqrt(qk.Two_norm()),0);	
		qk=qk*(1/rjk.r_find());
		R.p[k][k]=rjk;
		CGS.str_m(0,CGS.row,k,k,qk);
		for (int j=k+1;j<CGS.col;j++)
		{
			b=CGS.slice(0,CGS.row-1,j,j);		//求解b，rjk
			rjk=(!qk*b).p[0][0];
			R.p[k][j]=rjk;
			b=b-qk*rjk;
			CGS.str_m(0,CGS.row,j,j,b);
		}
	}

			/***计算误码率***/
	double n=double(2*c.row),w=0.0;
	Matrix y(CGS.col,1);
	y=(!CGS)*x;
	x.back(x,R,y);
	x.str_m(0,c.row-1,0,0,S*x.slice(0,c.row-1,0,0));
//	cout << c << x;
	for (int i=0;i<c.row;i++)
		if (c.p[i][0]!=x.p[i][0])
		{
			if (c.p[i][0].r_find()!=x.p[i][0].r_find())
				w=w+1.00;
			if (c.p[i][0].i_find()!=x.p[i][0].i_find())
				w=w+1.00;
		}
	return w/n;
}

		/***向量内积***/
double Matrix::Two_norm()
{
	Complex num(0,0);
	for (int i=0; i<row; i++)
		for (int j=0; j<col; j++)
			num=num+(!p[i][j])*p[i][j];
	return num.mochang();
}

		/***矩阵切片***/
Matrix Matrix::slice(int row1,int row2,int col1,int col2)	//需要切下矩阵的第row1行到第row2行，第col1列到第col2列
															//矩阵的行和列从0开始算
{	
	int u=0,v=0;
	Matrix a(row2-row1+1,col2-col1+1);
	for (int i=row1; i<row2+1; i++)
	{
		v=0;
		for (int j=col1; j<col2+1; j++)
		{
			a.p[u][v]=p[i][j];	//将对应元素赋给a
			v+=1;
		}
		u+=1;
	}
	return a;
}

		/***将小矩阵赋值给大矩阵的指定位置***/
void Matrix::str_m(int row1,int row2,int col1,int col2,Matrix x)//将矩阵x插入到大矩阵第row1到row2行，第col1到col2列
{
	int u=0,v=0;
	for (int i=row1; i<row2+1; i++)
	{
		v=0;
		for (int j=col1; j<col2+1; j++)
		{
			p[i][j]=x.p[u][v];	//将小矩阵对应元素赋给大矩阵
			v+=1;
		}
		u+=1;
	}
}

		/***下三角求逆矩阵***/
Matrix Matrix::ni_LM(Matrix L)		//L为需要求解的矩阵
{
	Matrix L_=L;
	for(int i=0;i<row;i++)
		for (int j=0;j<col;j++)
		{
			if (i==j)
				L.p[i][j]=L_.p[i][j].daoshu();	//对角元素为原来元素的倒数
			else
			{
				Complex sum;
				for (int k=0;k<i;k++)
					sum=sum-(L_.p[i][k]*L.p[k][j]);
				L.p[i][j]=sum/(L_.p[i][i]);	//下三角元素计算公式
			}
		}
	return L;
}

		/***上三角求逆矩阵***/
Matrix Matrix::ni_UM(Matrix L)	//利用下三角求逆矩阵解决上三角逆矩阵
{
	L=!L;
	L=L.ni_LM(L);
	return (!L);
}

		/***V_Blast算法***/
double Matrix::V_Blast(Matrix H,Matrix x,Matrix c)	//形参列表：CGS为信道矩阵，x为接受矩阵，c为发射矩阵
													//返回：误码率
{
	Matrix y(1,1), R(H.col,H.col);
	Complex c_;	double w=0.00, n=double(2*c.row),temp;
	int f=0;
	for (int i=0;i<H.col;)
	{
		f=0;
			/***求H的伪逆***/
		Matrix G(H.col,H.row);//H的伪逆
		G=(!H)*H;
		R=G.modified(G);
		R=G.ni_UM(R);
		G=(R*(!G))*(!H);

		/***找G中范数最小的行向量***/
		temp=G.slice(0,0,0,G.col-1).Two_norm();
		for (int j=1;j<G.row;j++)
			if (temp>G.slice(j,j,0,G.col-1).Two_norm())
			{
				f=j;
				temp=G.slice(j,j,0,G.col-1).Two_norm();
			}
		y=G.slice(f,f,0,G.col-1)*x;
		c_=y.qietiao(0,0);
		if (c.p[f][0]!=c_)
		{
			if (c.p[f][0].r_find()!=c_.r_find())
				w=w+1.00;
			if (c.p[f][0].i_find()!=c_.i_find())
				w=w+1.00;
		}
		x=x-H.slice(0,H.row-1,f,f)*c_;
		H=H.outl(H,f);
		c=c.outr(c,f);
	}
	return w/n;
}


Matrix Matrix::ZF(Matrix H,Matrix x,Matrix c)
{
	double w=0.00;
	Matrix G(H.col,H.row),Szf(H.row,1),R(H.col,H.col);
	
	/***求H的伪逆***/
	G=(!H)*H;
	R=G.modified(G);
	R=G.ni_UM(R);
	G=(R*(!G))*(!H);

	Szf=G*x;
	for (int i=0;i<H.row;i++)
	{
		Szf.p[i][0]=Szf.qietiao(i,0);
		if (c.p[i][0]!=Szf.p[i][0])
		{
			if (c.p[i][0].r_find()!=Szf.p[i][0].r_find())
				w=w+1.00;
			if (c.p[i][0].i_find()!=Szf.p[i][0].i_find())
				w=w+1.00;
		}
	}
	cout << w/(2*c.row) << "    ";
	return Szf;
}

Matrix Matrix::MMSE(Matrix H,Matrix x,Matrix c)
{
	double w=0.00;
	Matrix H_(H.row+H.col,H.col),x_(H.row+H.col,1),I(H.col,H.col),ZERO(H.col,1);
	for (int i=0;i<H.col;i++)
		I.p[i][i]=1;
	H_.str_m(0,H.row-1,0,H.col-1,H);
	H_.str_m(H.row,H_.row-1,0,H.col-1,I);
	x_.str_m(0,x.row-1,0,0,x);
	x_.str_m(x.row,x_.row-1,0,0,ZERO);
	Matrix G(H_.col,H_.row),Smmse(H.row,1),R(H_.col,H_.col);
	G=(!H_)*H_;
	R=G.modified(G);
	R=G.ni_UM(R);
	G=(R*(!G))*(!H_);
	Smmse=G*x_;
	for (int j=0;j<H.row;j++)
	{
		Smmse.p[j][0]=Smmse.qietiao(j,0);
		if (c.p[j][0]!=Smmse.p[j][0])
		{
			if (c.p[i][0].r_find()!=Smmse.p[j][0].r_find())
				w=w+1.00;
			if (c.p[i][0].i_find()!=Smmse.p[j][0].i_find())
				w=w+1.00;
		}
	}
	cout << w/(2*c.row) << endl;
	return Smmse;
}

Matrix Matrix::MMSE_QR(Matrix H,Matrix x,Matrix c)
{
	cout << H;
	Matrix R(H.col,H.col),Smmse_QR(H.row,1);
	R=H.modified(H);
	cout << H*R;
	H=!H;
	Smmse_QR=H*x;
	for (int j=0;j<H.row;j++)
		Smmse_QR.p[j][0]=Smmse_QR.qietiao(j,0);
	return Smmse_QR;
}

Matrix Matrix::MMSE_SQRD(Matrix H,Matrix x,Matrix c)
{
	Complex rjk;
	Matrix one(H.row+H.col,1),two(H.col,1),Smmse_SQRD(H.col,1);
	Matrix Q(H.row+H.col,H.col),I(H.col,H.col),R(H.col,H.col),norm(H.col,1),x_(H.row+H.col,1),ZERO(H.col,1);
	for (int i=0;i<H.col;i++)
		I.p[i][i]=1;
	Q.str_m(0,H.row-1,0,H.col-1,H);
	Q.str_m(H.row,Q.row-1,0,H.col-1,I);
	x_.str_m(0,x.row-1,0,0,x);
	x_.str_m(x.row,x_.row-1,0,0,ZERO);
	for (int j=0;j<H.col;j++)
		norm.p[j][0].change(Q.slice(0,Q.row-1,j,j).Two_norm(),0);
	for (int h=0;h<H.col;h++)
	{
		double temp=norm.p[h][0].r_find();
		int k=h;
		for(int f=h+1;f<H.col;f++)
			if (temp>norm.p[f][0].r_find())
			{
				temp=norm.p[f][0].r_find();
				k=f;
			}
		Matrix Q_=Q.slice(0,H.row+h,0,Q.col-1);
		Q_=Q_.swap(Q_,k,h);
		Q.str_m(0,H.row+h,0,Q.col-1,Q_);
		R=R.swap(R,k,h);
		I=I.swap(I,k,h);

		norm.p[k][0].change(norm.p[h][0].r_find(),0);
		norm.p[h][0].change(temp,0);
		Matrix b(Q.row,1),qk(Q.row,1);
		qk=Q.slice(0,Q.row-1,h,h);			//切片得到第h列未正交的列向量
		rjk.change(sqrt(norm.p[h][0].r_find()),0);	
		qk=qk*(1/rjk.r_find());
		R.p[h][h]=rjk;
		Q.str_m(0,Q.row-1,k,k,qk);
		for (int hh=h+1;hh<Q.col;hh++)
		{
			b=Q.slice(0,Q.row-1,hh,hh);		//求解b，rjk进行qk的正交化
			rjk=(!qk*b).p[0][0];
			R.p[h][hh]=rjk;
			b=b-qk*rjk;
			Q.str_m(0,Q.row-1,hh,hh,b);
			norm.p[hh][0]=norm.p[hh][0]-rjk*rjk;
		}
	}
//	cout << R;
//	cout << Q;
	Smmse_SQRD=I*(!Q*x_);
	return Smmse_SQRD;
}

Matrix Matrix::ML(Matrix H,Matrix x, Matrix c)
{
	Matrix c__(H.row,4);
	Matrix c_(H.col,4);
	Matrix c2(H.col,1);
	Matrix c3(H.col,1);
	Matrix c4(H.col,1);
	for (int i=0;i<H.col;i++)
		for (int k=0;k<4;k++)
		{
			if (k==0)
				c_.p[i][k].change(1,0);
			if (k==1)
				c_.p[i][k].change(0,1);
			if (k==2)
				c_.p[i][k].change(-1,0);
			if (k==3)
				c_.p[i][k].change(0,-1);
		}

	c__.str_m(0,H.row-1,0,0,H*c_.slice(0,H.col-1,0,0));
	c__.str_m(0,H.row-1,1,1,H*c_.slice(0,H.col-1,1,1));
	c__.str_m(0,H.row-1,2,2,H*c_.slice(0,H.col-1,2,2));
	c__.str_m(0,H.row-1,3,3,H*c_.slice(0,H.col-1,3,3));
	c__.str_m(0,H.row-1,0,0,x-c_.slice(0,H.col-1,0,0));
	c__.str_m(0,H.row-1,1,1,x-c_.slice(0,H.col-1,1,1));
	c__.str_m(0,H.row-1,2,2,x-c_.slice(0,H.col-1,2,2));
	c__.str_m(0,H.row-1,3,3,x-c_.slice(0,H.col-1,3,3));
	for (int j=0;j<H.row;j++)
	{
		double temp=c__.p[j][0].mochang();
		int f=0;
		for (int r=1;r<4;r++)
		{
			if (temp>c__.p[j][r].mochang())
			{
				temp=c__.p[j][r].mochang();
				f=r;
			}
		}
		x.p[j][0]=c_.p[0][r];
	}
	return x;
}
			
	/***主函数***/
int main()
{
	srand((int)time(0));
	int Hrow=5,Hcol=4;
	double ber1=0,ber2=0,ber3=0;
	for (double j=0;j<0.1;)
	{
		ber1=0; ber2=0;ber3=0;
	for (int i=0;i<1;i++)
	{
		Matrix H(Hrow,Hcol),v(Hrow,1),c(Hcol,1),x(Hrow,1);
		double k1=sqrt(pow(10,j)/(double(Hrow*Hcol)));
		H.GaoSi();
		v.GaoSi();
		c.Transmit();
		H=H*k1;
		v=v*(1/sqrt(2));
		x=H*c+v;
//		c.ZF(H,x,c);
//		c.MMSE(H,x,c);
//		c.MMSE_QR(H,x,c);
//		c.MMSE_SQRD(H,x,c);
//		ber1=ber1+x.V_Blast(H,x,c);
//		ber2=ber2+x.QR(H,x,c);
//		ber3=ber3+x.SQRD(H,x,c);
	}
//	cout << ber1/i << "   " << ber2/i << "  " << ber3/i << endl;
	j=j+0.1;
	}
	return 0;
}
