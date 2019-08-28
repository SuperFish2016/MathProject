

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
// 本程序实现matrix类
// 对matrix类的定义
#include "matrix.h"

matrix::matrix(buffer * b): // 缺省构造函数，产生0行0列空矩阵
	rownum(0),colnum(0),istrans(0),isneg(0),issym(0)
{
		if(b){ // 如用户给出b，则使用b作为数据缓存
			buf=b;
			buf->alloc(0);
			}
		else // 否则，产生一个新的缓存，类别是当前缺省的缓存类
			buf = getnewbuffer(0);
}

matrix::matrix(size_t n, buffer * b): // 产生n阶单位方阵
		rownum(n),colnum(1),istrans(0),isneg(0),issym(0)
{
		if(b){ // 如用户给出b，则使用b作为数据缓存
			buf=b;
			buf->alloc(n);
			}
		else // 否则，产生一个新的缓存，类别是当前缺省的缓存类
			buf = getnewbuffer(n);
}

matrix unit(size_t n)	// 产生n阶单位矩阵
{
	if(n==0) throw TMESSAGE("n must larger than 0\n");
	matrix m(n,n);
	for(size_t i=0; i<n; i++)
	for(size_t j=0; j<n; j++)
		if(i==j) m.set(i,j,1.0);
		else m.set(i,j,0.0);
	m.issym = 1;
	return m;
}


matrix::matrix(size_t r, size_t c, buffer * b):
		rownum(r),colnum(c),istrans(0),isneg(0),issym(0)
{
		if(b){ // 如用户给出b，则使用b作为数据缓存
			buf=b;
			buf->alloc(r*c);
			}
		else // 否则，产生一个新的缓存，类别是当前缺省的缓存类
			buf = getnewbuffer(r*c);
}

matrix::matrix(const matrix& m) // 拷贝构造函数
{
	rownum = m.rownum;	// 拷贝行数与列数
	colnum = m.colnum;
	istrans = m.istrans;	// 拷贝转置标记
	isneg = m.isneg;		// 拷贝负值标记
	issym = m.issym;		// 拷贝对称标记
	buf = m.buf;		// 设置指向相同缓存的指针
	buf->refnum++;	// 缓存的引用数增1
}

matrix::matrix(const char * filename, buffer * b): // 从数据文件构造矩阵
	istrans(0), isneg(0), issym(0)
{
	char label[10];
	ifstream in(filename); // 打开文件流
	in >> label; // 文件开始必须有matrix关键词
	if(strcmp(label, "matrix")!=0) throw TMESSAGE("format error!");
	in >> rownum >> colnum; // 读取行数和列数
	if(!in.good()) throw TMESSAGE("read file error!");
	// 申请或产生缓存
	if(b) { buf=b;
		buf->alloc(rownum*colnum);
	}
	else buf = getnewbuffer(rownum*colnum);
	size_t line;
	for(size_t i=0; i<rownum; i++) { // 依次按行读取
		in >> line; // 读行号
		if(line != i+1) throw TMESSAGE("format error!");
		in.width(sizeof(label));
		in >> label; // 行号后跟冒号
		if(label[0] != ':') throw TMESSAGE("format error!");
		DOUBLE a;
		for(size_t j=0; j<colnum; j++) { // 读取一行数据
			in >> a;
			set(i,j,a);
		}
		if(!in.good()) throw TMESSAGE("read file error!");
	}
	checksym(); // 检查是否为对称阵
}

matrix::matrix(void * data, size_t r, size_t c, buffer * b):
		rownum(r),colnum(c),istrans(0),isneg(0),issym(0) // 数据构造函数
{
	if(b){
		buf=b;
		buf->alloc(r*c);
		}
	else
		buf = getnewbuffer(r*c);
	DOUBLE * d = (DOUBLE *)data;
	for(size_t i=0; i<r*c; i++)  // 这里进行数据拷贝，因此原数据的内存是可以释放的
		buf->set(i,d[i]);
	checksym();	// 检查是否为对称阵
}

DOUBLE matrix::operator()(size_t r, size_t c)
{
	if(r>= rownum || c>= colnum)
		throw TMESSAGE("Out range!");
	return value(r,c);
}

matrix& matrix::operator=(const matrix& m)  // 赋值重载
{
	rownum = m.rownum; //  行数和列数的拷贝
	colnum = m.colnum;
	istrans = m.istrans; // 转置标志的拷贝
	isneg = m.isneg;	// 取负标志的拷贝
	issym = m.issym;	// 对称标志的拷贝
	if(buf == m.buf)  // 如果原缓存与m的缓存一样，则返回
		return (*this);
	buf->refnum--; // 原缓存不同，则原缓存的引用数减1
	if(!buf->refnum)delete buf; // 减1后的原缓存如果引用数为0，则删除原缓存
	buf = m.buf; // 将原缓存指针指向m的缓存
	buf->refnum++; // 新缓存的引用数增1
	checksym();	// 检查是否为对称阵
	return (*this); // 返回自己的引用
}

matrix& matrix::operator=(DOUBLE a) // 通过赋值运算符将矩阵所有元素设为a
{
	if(rownum == 0 || colnum == 0) return (*this);
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
		set(i,j,a);
	if(rownum == colnum)issym = 1;
	return (*this);
}

matrix matrix::operator-() // 矩阵求负，产生负矩阵
{
	matrix mm(*this);
	mm.neg();
	return mm;
}

ostream& operator<<(ostream& o, matrix& m) // 流输出运算符
{
	// 先输出关键字matrix,然后是行号，制表符，列号，换行
	o << "matrix " << m.rownum << '\t' << m.colnum << endl;
	for(size_t i=0; i<m.rownum; i++) { // 依次输出各行
		o<< (i+1) <<':'; // 行号后跟冒号。注意输出的行号是从1开始
				// 是内部编程的行号加1
		size_t k=8; // 后面输入一行的内容，每次一行数字超过八个时
			// 就换一行
		for(size_t j=0; j<m.colnum; j++) {
			o<<'\t'<<m.value(i,j);
			if(--k==0) {
				k=8;
				o<<endl;
			}
		}
		o<<endl; // 每行结束时也换行
	}
	return o;
}

istream& operator>>(istream& in, matrix& m) // 流输入运算符
{
	char label[10];
	in.width(sizeof(label));
	in >> label; // 输入matrix关键字
	if(strcmp(label, "matrix")!=0)
		throw TMESSAGE("format error!\n");
	in >> m.rownum >> m.colnum; // 输入行号和列号
	if(!in.good()) throw TMESSAGE("read file error!");
	m.buf->refnum--; // 原缓存引用数减1
	if(!m.buf->refnum) delete m.buf; // 如原缓存引用数为0，则删去原缓存
	m.isneg = m.istrans = 0; // 转置和取负标志清零
	m.buf = getnewbuffer(m.rownum*m.colnum); // 按缺省子类产生新缓存
	size_t line; // 按矩阵行输入
	for(size_t i=0; i<m.rownum; i++) {
		in >> line; // 先输入行号
		if(line != i+1) throw TMESSAGE("format error!\n");
		in.width(sizeof(label)); // 行号后应跟冒号
		in >> label;
		if(label[0] != ':') throw TMESSAGE("format error!\n");
		DOUBLE a; // 随后是本行的各个数值
		for(size_t j=0; j<m.colnum; j++) {
			in >> a;
			m.set(i,j,a);
		}
		if(!in.good()) throw TMESSAGE("read file error!");
	}
	m.checksym();	// 检查是否为对称阵
	return in;
}

matrix& matrix::operator*=(DOUBLE a) // 矩阵数乘常数a，结果放在原矩阵
{
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
		set(i,j,a*value(i,j));
	return (*this);
}

matrix matrix::operator*(DOUBLE a) // 矩阵数乘常数a，原矩阵内容不变，返回一新矩阵
{
	matrix m(rownum, colnum);
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
		m.set(i,j,a*value(i,j));
	return m;
}

matrix matrix::operator+(const matrix& m) // 矩阵相加，产生一新的矩阵并返回它
{
	if(rownum != m.rownum || colnum != m.colnum) // 对应行列必须相同
		throw TMESSAGE("can not do add of matrix\n");
	matrix mm(rownum, colnum); // 产生一同自己同形的矩阵
	DOUBLE a;
    matrix temp = m;
	for(size_t i=0; i<rownum; i++) // 求和
	for(size_t j=0; j<colnum; j++)
	{
        a = value(i,j)+temp.value(i,j);
		mm.set(i,j,a);
	}
	mm.checksym();	// 检查是否为对称阵
	return mm;
}

matrix& matrix::operator+=(const matrix &m) // 矩阵求和，自己内容改变为和
{
	DOUBLE a;
    matrix temp = m;
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
	{
        a = value(i,j)+temp.value(i,j);
		set(i,j,a);
	}
	checksym();	// 检查是否为对称阵
	return (*this);
}

matrix matrix::operator+(DOUBLE a)	// 矩阵加常数，指每一元素加一固定的常数，产生
											// 新矩阵，原矩阵不变
{
	matrix m(rownum, colnum);
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
		m.set(i,j,a+value(i,j));
	return m;
}

matrix& matrix::operator+=(DOUBLE a)	// 矩阵自身加常数，自身内容改变
{
	for(size_t i=0; i<rownum; i++)
	for(size_t j=0; j<colnum; j++)
		set(i,j,a+value(i,j));
	return (*this);
}

matrix operator-(DOUBLE a, matrix& m) { // 常数减矩阵，产生新的矩阵
	return (-m)+a;
};


matrix matrix::operator-(const matrix& m) // 矩阵相减，产生新的矩阵
{
	matrix mm(*this); // 产生一同自己同形的矩阵
    matrix t = m;
    mm += (-t);	// 加上相应的负矩阵
	return mm;
}

matrix& matrix::operator-=(const matrix& m) // 矩阵相减，结果修改原矩阵
{
    matrix mm = m;
    (*this) += (-mm);
	return (*this);
}

matrix matrix::operator*(const matrix& m) // 矩阵相乘，原矩阵内容不变，产生一新矩阵
{
	if(colnum != m.rownum) // 必须满足相乘条件
		throw TMESSAGE("can not multiply!");
	matrix mm(rownum,m.colnum); // 计算并产生一合要求的矩阵放乘积
	DOUBLE a;
    matrix tem = m;
	for(size_t i=0; i<rownum; i++) // 计算乘积
	for(size_t j=0; j<m.colnum; j++){
		a = 0.0;
		for(size_t k=0; k<colnum; k++)
            a += value(i,k)*tem.value(k,j);
		mm.set(i,j,a);
	}
	mm.checksym();	// 检查是否为对称阵
	return mm; // 返回乘积
}

matrix& matrix::operator*=(const matrix& m) // 矩阵相乘，自己修改成积矩阵
{
	(*this) = (*this)*m;
	return (*this);
}

matrix matrix::t()  // 矩阵转置，产生新的矩阵
{
	matrix mm(*this);
	mm.trans();
	return mm;
}

int matrix::isnear(const matrix& m, double e) // 检查二矩阵是否近似相等
{
	if(rownum != m.rownum || colnum != m.colnum) return 0;
    matrix temp = m;
	for(size_t i=0; i< rownum; i++)
	for(size_t j=0; j< colnum; j++)
        if(fabs(value(i,j)-temp.value(i,j)) > e) return 0;
	return 1;
}

int matrix::isnearunit(double e) // 检查矩阵是否近似为单位矩阵
{
	if(rownum != colnum) return 0;
	return isnear(unit(rownum), e);
}

matrix matrix::row(size_t r) // 提取第r行行向量
{
	matrix mm(1, colnum);
	for(size_t i=0; i< colnum; i++)
		mm.set(0, i, value(r,i));
	return mm;
}

matrix matrix::col(size_t c) // 提取第c列列向量
{
	matrix mm(rownum, 1);
	for(size_t i=0; i< rownum; i++)
		mm.set(i, value(i, c));
	return mm;
}

void matrix::swapr(size_t r1, size_t r2, size_t k) // 交换矩阵r1和r2两行
{
	DOUBLE a;
	for(size_t i=k; i<colnum; i++) {
		a = value(r1, i);
		set(r1, i, value(r2, i));
		set(r2, i, a);
	}
}

void matrix::swapc(size_t c1, size_t c2, size_t k) // 交换c1和c2两列
{
	DOUBLE a;
	for(size_t i=k; i<colnum; i++) {
		a = value(i, c1);
		set(i, c1, value(i, c2));
		set(i, c2, a);
	}
}

DOUBLE matrix::maxabs(size_t &r, size_t &c, size_t k) // 求第k行和第k列后的主元及位置
{
	DOUBLE a=0.0;
	for(size_t i=k;i<rownum;i++)
	for(size_t j=k;j<colnum;j++)
		if(a < fabs(value(i,j))) {
			r=i;c=j;a=fabs(value(i,j));
		}
	return a;
}

size_t matrix::zgsxy(matrix & m, int fn) // 进行主高斯消元运算，fn为参数，缺省为0
 /* 本矩阵其实是常数阵，而矩阵m必须是方阵
	运算过程其实是对本矩阵和m同时作行初等变换，
	运算结果m的对角线相乘将是行列式，而本矩阵则变换成
	自己的原矩阵被m的逆阵左乘，m的秩被返回，如果秩等于阶数
	则本矩阵中的内容已经是唯一解
 */
{
	if(rownum != m.rownum || m.rownum != m.colnum) // 本矩阵行数必须与m相等
						// 且m必须是方阵
		throw TMESSAGE("can not divide!");
	lbuffer * bb = getnewlbuffer(rownum); // 产生一维数为行数的长整数缓存区
	lbuffer & b = (*bb); // 用引用的办法使下面的程序容易懂
	size_t is;
	DOUBLE a;
    size_t i,j,rank,k=0;
    for(k=0; k<rownum; k++) { // 从第0行到第k行进行主高斯消元
		if(m.maxabs(is, i, k)==0) // 求m中第k级主元，主元所在的行，列存在is,i中
			break; // 如果主元为零，则m不可逆，运算失败
		rank = k+1; // rank存放当前的阶数
		b.retrieve(k) = i;  // 将长整数缓存区的第k个值设为i
		if(i != k)
			m.swapc(k, i); // 交换m中i,k两列
		if(is != k) {
			m.swapr(k, is, k); // 交换m中i,k两行,从k列以后交换
			swapr(k, is); // 交换本矩阵中i,k两行
		}
		a = m.value(k,k);  // 取出主元元素
		for (j=k+1;j<rownum;j++) // 本意是将m的第k行除以主元
				// 但只需把第k行的第k+1列以上除以主元即可
				// 这样还保留了主元作行列式运算用
			m.set(k,j,m.value(k,j)/a);
		for (j=0;j<colnum;j++) // 将本矩阵的第k行除以主元
			set(k,j,value(k,j)/a);
		// 上面两步相当于将m和本矩阵构成的增广矩阵第k行除以主元
		// 下面对增广矩阵作行基本初等变换使第k行的其余列均为零
		// 但0值无必要计算，因此从第k+1列开始计算
		for(j=k+1;j<rownum;j++) // j代表列，本矩阵的行数就是m的列数
		for(i=0;i<rownum;i++) //i代表行，依次对各行计算，k行除外
			if(i!=k)
				m.set(i,j,m.value(i,j)-m.value(i,k)*m.value(k,j));
		// 对本矩阵亦作同样的计算
		for(j=0;j<colnum;j++)
		for(i=0;i<rownum;i++)
			if(i!=k)
				set(i,j,value(i,j)-m.value(i,k)*value(k,j));
	} // 主高斯消元循环k结束
	if(fn == 1) {
		for(j=0; j<rank; j++)
		for(i=0; i<rownum; i++)
			if(i==j) m.set(i,i,1.0);
			else
				m.set(i,j,0.0);
		for(k = rank; k>0; k--)
			m.swapc(k-1,(size_t)b[k-1]);
	}
	for(k = rank; k>0; k--) // 将本矩阵中的各行按b中内容进行调节
		if((size_t)b[k-1] != k-1)
			swapr(k-1,(size_t)b[k-1]); // 行交换
	delete bb; // 释放长整数缓存
	return rank; // 返回mm的秩
}

matrix& matrix::operator/=(const matrix m) // 利用重载的除法符号/=来解方程
	// 本矩阵设为d,则方程为mx=d,考虑解写成x=d/m的形式，
	// 而方程的解也存放在d中，则实际编程时写d/=m
{
    matrix mm = m;
    if(zgsxy(mm)!=rownum) // 如秩不等于m的阶数，则方程无解
		throw TMESSAGE("can not divide!");
	return *this;
}

matrix matrix::operator/(matrix m) // 左乘m的逆矩阵产生新矩阵
{
	m.inv();	// m的逆矩阵
	return (*this)*m;
}

matrix& matrix::inv()		// 用全选主元高斯-约当法求逆矩阵
{
	if(rownum != colnum || rownum == 0)
		throw TMESSAGE("Can not calculate inverse");
	size_t i,j,k;
	 DOUBLE d,p;
	lbuffer * isp = getnewlbuffer(rownum); // 产生一维数为行数的长整数缓存区
	lbuffer * jsp = getnewlbuffer(rownum); // 产生一维数为行数的长整数缓存区
	lbuffer& is = *isp;	// 使用引用使程序看起来方便
	lbuffer& js = *jsp;
	for(k=0; k<rownum; k++)
		{
			d = maxabs(i, j, k); // 全主元的位置和值
			is[k] = i;
			js[k] = j;
			if(d==0.0) {
				delete isp;
				delete jsp;
				throw TMESSAGE("can not inverse");
			}
		  if ((size_t)is[k] != k) swapr(k,(size_t)is[k]);
		  if ((size_t)js[k] != k) swapc(k,(size_t)js[k]);
			p = 1.0/value(k,k);
			set(k,k,p);
		  for (j=0; j<rownum; j++)
			 if (j!=k) set(k,j,value(k,j)*p);
		  for (i=0; i<rownum; i++)
			 if (i!=k)
				for (j=0; j<rownum; j++)
				  if (j!=k) set(i,j,value(i,j)-value(i,k)*value(k,j));
		  for (i=0; i<rownum; i++)
			 if (i!=k) set(i,k,-value(i,k)*p);
		} // end for k
	 for (k=rownum; k>0; k--)
		{ if ((size_t)js[k-1]!=k-1) swapr((size_t)js[k-1], k-1);
		  if ((size_t)is[k-1]!=k-1) swapc((size_t)is[k-1], k-1);
		}
  delete isp;
  delete jsp;
	checksym();	// 检查是否为对称阵
  return (*this);
}

matrix matrix::operator~()	// 求逆矩阵，但产生新矩阵
{
	matrix m(*this);
	m.inv();
	return m;
}

matrix operator/(DOUBLE a, matrix& m) // 求逆矩阵再乘常数
{
	matrix mm(m);
	mm.inv();
	if(a != 1.0) mm*=a;
	return mm;
}

matrix& matrix::operator/=(DOUBLE a) // 所有元素乘a的倒数，自身改变
{
	return operator*=(1/a);
}

matrix matrix::operator/(DOUBLE a) // 所有元素乘a的倒数，产生新的矩阵
{
	matrix m(*this);
	m/=a;
	return m;
}

DOUBLE matrix::det(DOUBLE err)		// 求行列式的值
{
	if(rownum != colnum || rownum == 0)
		throw TMESSAGE("Can not calculate det");
	matrix m(*this);
	size_t rank;
	return m.detandrank(rank, err);
}

size_t matrix::rank(DOUBLE err)	// 求矩阵的秩
{
	if(rownum==0 || colnum==0)return 0;
	size_t k;
	k = rownum > colnum ? colnum : rownum;
	matrix m(k,k);		// 产生k阶方阵
	for(size_t i=0; i<k; i++)
	for(size_t j=0; j<k; j++)
		m.set(i,j,value(i,j));
	size_t r;
	m.detandrank(r, err);
	return r;
}

DOUBLE matrix::detandrank(size_t & rank, DOUBLE err)	// 求方阵的行列式和秩
{
	if(rownum != colnum || rownum == 0)
		throw TMESSAGE("calculate det and rank error!");
	size_t i,j,k,is,js;
	double f,detv,q,d;
	f=1.0; detv=1.0;
	rank = 0;
	 for (k=0; k<rownum-1; k++)
		{
			q = maxabs(is, js, k);
			if(q<err) return 0.0;	// 如主元太小，则行列式的值被认为是0
			rank++; // 秩增1
			if(is!=k) { f=-f; swapr(k,is,k); }
			if(js!=k) { f=-f; swapc(k,js,k); }
			q = value(k,k);
			detv *= q;
		  for (i=k+1; i<rownum; i++)
			 {
				d=value(i,k)/q;
				for (j=k+1; j<rownum; j++)
					set(i,j,value(i,j)-d*value(k,j));
			 }
		} // end loop k
	 q = value(rownum-1,rownum-1);
	 if(q != 0.0 ) rank++;
	return f*detv*q;
}

void matrix::checksym()	// 检查和尝试调整矩阵到对称矩阵
{
	issym = 0;	// 先假设矩阵非对称
	if(rownum != colnum) return;	// 行列不等当然不是对称矩阵
	DOUBLE a,b;
	for(size_t i=1; i<rownum; i++) // 从第二行开始检查
	for(size_t j=0; j<i; j++) // 从第一列到第i-1列
	{
		a = value(i,j);
		b = value(j,i);
		if( fabs(a-b) >= defaulterr ) return; // 发现不对称，返回
		if(a!=b)set(i,j,b); // 差别很小就进行微调
	}
	issym = 1;	// 符合对称阵标准
}

void matrix::house(buffer & b, buffer & c)// 用豪斯荷尔德变换将对称阵变为对称三对
											// 角阵，b返回主对角线元素，c返回次对角线元素
{
	if(!issym) throw TMESSAGE("not symatry");
	size_t i,j,k;
	DOUBLE h,f,g,h2,a;
	 for (i=rownum-1; i>0; i--)
		{ h=0.0;
		  if (i>1)
			 for (k=0; k<i; k++)
				{ a = value(i,k); h += a*a;}
		  if (h == 0.0)
			 { c[i] = 0.0;
				if (i==1) c[i] = value(i,i-1);
				b[i] = 0.0;
			 }
		  else
			 { c[i] = sqrt(h);
				a = value(i,i-1);
				if (a > 0.0) c[i] = -c[i];
				h -= a*c[i];
				set(i,i-1,a-c[i]);
				f=0.0;
				for (j=0; j<i; j++)
				  { set(j,i,value(i,j)/h);
					 g=0.0;
					 for (k=0; k<=j; k++)
						g += value(j,k)*value(i,k);
					 if(j<=i-2)
						for (k=j+1; k<i; k++)
							g += value(k,j)*value(i,k);
					 c[j] = g/h;
					 f += g*value(j,i);
				  }
				h2=f/(2*h);
				for (j=0; j<i; j++)
				  { f=value(i,j);
					 g=c[j] - h2*f;
					 c[j] = g;
					 for (k=0; k<=j; k++)
						set(j,k, value(j,k)-f*c[k]-g*value(i,k) );
				  }
				b[i] = h;
			 }
		}
	 for (i=0; i<=rownum-2; i++) c[i] = c[i+1];
	 c[rownum-1] = 0.0;
	 b[0] = 0.0;
	 for (i=0; i<rownum; i++)
		{ if ((b[i]!=0.0)&&(i>=1))
			 for (j=0; j<i; j++)
				{ g=0.0;
				  for (k=0; k<i; k++)
					 g=g+value(i,k)*value(k,j);
				  for (k=0; k<i; k++)
					 set(k,j,value(k,j)-g*value(k,i));
				}
		  b[i] = value(i,i);
		  set(i,i,1.0);
		  if (i>=1)
			 for (j=0; j<=i-1; j++)
				{ set(i,j,0.0);
				  set(j,i,0.0); }
		}
	 return;
}

void matrix::trieigen(buffer& b, buffer& c, size_t l, DOUBLE eps)
						// 计算三对角阵的全部特征值与特征向量
{	size_t i,j,k,m,it;
	double d,f,h,g,p,r,e,s;
	c[rownum-1]=0.0; d=0.0; f=0.0;
	for (j=0; j<rownum; j++)
		{ it=0;
		  h=eps*(fabs(b[j])+fabs(c[j]));
		  if (h>d) d=h;
		  m=j;
		  while ((m<rownum)&&(fabs(c[m])>d)) m+=1;
		  if (m!=j)
			 { do
				  { if (it==l) throw TMESSAGE("fial to calculate eigen value");
					 it += 1;
					 g=b[j];
					 p=(b[j+1]-g)/(2.0*c[j]);
					 r=sqrt(p*p+1.0);
					 if (p>=0.0) b[j]=c[j]/(p+r);
					 else b[j]=c[j]/(p-r);
					 h=g-b[j];
					 for (i=j+1; i<rownum; i++)
						b[i]-=h;
					 f=f+h; p=b[m]; e=1.0; s=0.0;
					 for (i=m-1; i>=j; i--)
						{ g=e*c[i]; h=e*p;
						  if (fabs(p)>=fabs(c[i]))
							 { e=c[i]/p; r=sqrt(e*e+1.0);
								c[i+1]=s*p*r; s=e/r; e=1.0/r;
							 }
						  else
				{ e=p/c[i]; r=sqrt(e*e+1.0);
								c[i+1]=s*c[i]*r;
								s=1.0/r; e=e/r;
							 }
						  p=e*b[i]-s*g;
						  b[i+1]=h+s*(e*g+s*b[i]);
						  for (k=0; k<rownum; k++)
							 {
								h=value(k,i+1);
								set(k,i+1, s*value(k,i)+e*h);;
								set(k,i,e*value(k,i)-s*h);
							 }
						  if(i==0) break;
						}
					 c[j]=s*p; b[j]=e*p;
				  }
				while (fabs(c[j])>d);
			 }
		  b[j]+=f;
		}
	 for (i=0; i<=rownum; i++)
		{ k=i; p=b[i];
		  if (i+1<rownum)
			 { j=i+1;
				while ((j<rownum)&&(b[j]<=p))
				  { k=j; p=b[j]; j++;}
			 }
		  if (k!=i)
			 { b[k]=b[i]; b[i]=p;
				for (j=0; j<rownum; j++)
				  { p=value(j,i);
					 set(j,i,value(j,k));
					 set(j,k,p);
				  }
			 }
		}
}

matrix matrix::eigen(matrix & eigenvalue, DOUBLE eps, size_t l)
	 // 计算对称阵的全部特征向量和特征值
		// 返回按列排放的特征向量，而eigenvalue将返回一维矩阵，为所有特征值
		// 组成的列向量
{
	if(!issym) throw TMESSAGE("it is not symetic matrix");
	eigenvalue = matrix(rownum); // 产生n行1列向量准备放置特征值
	matrix m(*this); // 复制自己产生一新矩阵
	if(rownum == 1) {	// 如果只有1X1矩阵，则特征向量为1，特征值为value(0,0)
		m.set(0,0,1.0);
		eigenvalue.set(0,value(0,0));
		return m;
	}
	buffer * b, *c;
	b = getnewbuffer(rownum);
	c = getnewbuffer(rownum);
	m.house(*b,*c);	// 转换成三对角阵
	m.trieigen(*b,*c,l,eps); // 算出特征向量和特征值
	for(size_t i=0; i<rownum; i++) // 复制b的内容到eigenvalue中
		eigenvalue.set(i,(*b)[i]);
	return m;
}

void matrix::hessenberg()	// 将一般实矩阵约化为赫申伯格矩阵
{
	 size_t i,j,k;
	 double d,t;
	 for (k=1; k<rownum-1; k++)
		{ d=0.0;
		  for (j=k; j<rownum; j++)
			 { t=value(j,k-1);
				if (fabs(t)>fabs(d))
				  { d=t; i=j;}
			 }
		  if (fabs(d)!=0.0)
			 { if (i!=k)
				  { for (j=k-1; j<rownum; j++)
						{
						  t = value(i,j);
						  set(i,j,value(k,j));
						  set(k,j,t);
						}
					 for (j=0; j<rownum; j++)
						{
						  t = value(j,i);
						  set(j,i,value(j,k));
						  set(j,k,t);
						}
				  }
				for (i=k+1; i<rownum; i++)
				  {
					 t = value(i,k-1)/d;
					 set(i,k-1,0.0);
					 for (j=k; j<rownum; j++)
						  set(i,j,value(i,j)-t*value(k,j));
					 for (j=0; j<rownum; j++)
						  set(j,k,value(j,k)+t*value(j,i));
				  }
			 }
		}
}

void matrix::qreigen(matrix & u1, matrix & v1, size_t jt, DOUBLE eps)
 // 求一般实矩阵的所有特征根
// a和b均返回rownum行一列的列向量矩阵，返回所有特征根的实部和虚部
{
	matrix a(*this);
	a.hessenberg();	// 先算出赫申伯格矩阵
	u1 = matrix(rownum);
	v1 = matrix(rownum);
	buffer * uu = getnewbuffer(rownum);
	buffer * vv = getnewbuffer(rownum);
	buffer &u = *uu;
	buffer &v = *vv;
	 size_t m,it,i,j,k,l;
	 size_t iir,iic,jjr,jjc,kkr,kkc,llr,llc;
	 DOUBLE b,c,w,g,xy,p,q,r,x,s,e,f,z,y;
	 it=0; m=rownum;
	 while (m!=0)
		{ l=m-1;
		  while ( l>0 && (fabs(a.value(l,l-1))>eps*
			(fabs(a.value(l-1,l-1))+fabs(a.value(l,l))))) l--;
		  iir = m-1; iic = m-1;
		  jjr = m-1; jjc = m-2;
		  kkr = m-2; kkc = m-1;
		  llr = m-2; llc = m-2;
		  if (l==m-1)
			 { u[m-1]=a.value(m-1,m-1); v[m-1]=0.0;
				m--; it=0;
			 }
		  else if (l==m-2)
			 { b=-(a.value(iir,iic)+a.value(llr,llc));
				c=a.value(iir,iic)*a.value(llr,llc)-
					a.value(jjr,jjc)*a.value(kkr,kkc);
				w=b*b-4.0*c;
				y=sqrt(fabs(w));
				if (w>0.0)
				  { xy=1.0;
					 if (b<0.0) xy=-1.0;
					 u[m-1]=(-b-xy*y)/2.0;
					 u[m-2]=c/u[m-1];
					 v[m-1]=0.0; v[m-2]=0.0;
				  }
				else
				  { u[m-1]=-b/2.0; u[m-2]=u[m-1];
					 v[m-1]=y/2.0; v[m-2]=-v[m-1];
				  }
				m=m-2; it=0;
			 }
		  else
			 {
			 if (it>=jt) {
				delete uu;
				delete vv;
				throw TMESSAGE("fail to calculate eigenvalue");
			 }
				it++;
				for (j=l+2; j<m; j++)
					a.set(j,j-2,0.0);
				for (j=l+3; j<m; j++)
					a.set(j,j-3,0.0);
				for (k=l; k+1<m; k++)
				  { if (k!=l)
						{ p=a.value(k,k-1); q=a.value(k+1,k-1);
						  r=0.0;
						  if (k!=m-2) r=a.value(k+2,k-1);
						}
					 else
						{
						  x=a.value(iir,iic)+a.value(llr,llc);
						  y=a.value(llr,llc)*a.value(iir,iic)-
								a.value(kkr,kkc)*a.value(jjr,jjc);
						  iir = l; iic = l;
						  jjr = l; jjc = l+1;
						  kkr = l+1; kkc = l;
						  llr = l+1; llc = l+1;
						  p=a.value(iir,iic)*(a.value(iir,iic)-x)
								+a.value(jjr,jjc)*a.value(kkr,kkc)+y;
						  q=a.value(kkr,kkc)*(a.value(iir,iic)+a.value(llr,llc)-x);
						  r=a.value(kkr,kkc)*a.value(l+2,l+1);
						}
					 if ((fabs(p)+fabs(q)+fabs(r))!=0.0)
						{ xy=1.0;
						  if (p<0.0) xy=-1.0;
						  s=xy*sqrt(p*p+q*q+r*r);
						  if (k!=l) a.set(k,k-1,-s);
						  e=-q/s; f=-r/s; x=-p/s;
						  y=-x-f*r/(p+s);
						  g=e*r/(p+s);
						  z=-x-e*q/(p+s);
						  for (j=k; j<m; j++)
							 {
								iir = k; iic = j;
								jjr = k+1; jjc = j;
								p=x*a.value(iir,iic)+e*a.value(jjr,jjc);
								q=e*a.value(iir,iic)+y*a.value(jjr,jjc);
								r=f*a.value(iir,iic)+g*a.value(jjr,jjc);
								if (k!=m-2)
								  { kkr = k+2; kkc = j;
									 p=p+f*a.value(kkr,kkc);
									 q=q+g*a.value(kkr,kkc);
									 r=r+z*a.value(kkr,kkc);
									 a.set(kkr,kkc,r);
								  }
								a.set(jjr,jjc,q);
								a.set(iir,iic,p);
							 }
						  j=k+3;
						  if (j>=m-1) j=m-1;
						  for (i=l; i<=j; i++)
							 { iir = i; iic = k;
								jjr = i; jjc = k+1;
								p=x*a.value(iir,iic)+e*a.value(jjr,jjc);
								q=e*a.value(iir,iic)+y*a.value(jjr,jjc);
								r=f*a.value(iir,iic)+g*a.value(jjr,jjc);
								if (k!=m-2)
								  { kkr = i; kkc = k+2;
									 p=p+f*a.value(kkr,kkc);
									 q=q+g*a.value(kkr,kkc);
									 r=r+z*a.value(kkr,kkc);
									 a.set(kkr,kkc,r);
								  }
								a.set(jjr,jjc,q);
								a.set(iir,iic,p);
							 }
						}
				  }
			 }
		}
	for(i=0;i<rownum;i++) {
		u1.set(i,u[i]);
		v1.set(i,v[i]);
	}
	delete uu;
	delete vv;
}

DOUBLE gassrand(int rr)	// 返回一零均值单位方差的正态分布随机数
{
	static DOUBLE r=3.0;	// 种子
	if(rr) r = rr;
	int i,m;
	DOUBLE s,w,v,t;
	s=65536.0; w=2053.0; v=13849.0;
	t=0.0;
	for (i=1; i<=12; i++)
		{ r=r*w+v; m=(int)(r/s);
		  r-=m*s; t+=r/s;
		}
	t-=6.0;
	return(t);
}

gassvector::gassvector(matrix & r):	//r必须是正定对称阵，为正态随机向量的协方差
	a(r)
{
	if(r.rownum != r.colnum) throw TMESSAGE("must be a sqare matrix");
	if(!r.issym) throw TMESSAGE("must be a symetic matrix");
	matrix evalue;
	a = a.eigen(evalue);
	for(size_t i=0; i<a.colnum; i++) {
		DOUBLE e = sqrt(evalue(i));
		for(size_t j=0; j<a.rownum; j++)
			a.set(j,i,a.value(j,i)*e);
	}
}

matrix gassvector::operator()(matrix & r) // 返回给定协方差矩阵的正态随机向量
{
	a = r;
	matrix evalue;
	a = a.eigen(evalue);
	size_t i;
	for(i=0; i<a.colnum; i++) {
		DOUBLE e = sqrt(evalue(i));
		for(size_t j=0; j<a.rownum; j++)
			a.set(j,i,a.value(j,i)*e);
	}
	matrix rr(a.rownum);	// 产生列向量
	for(i=0; i<a.rownum; i++)
		rr.set(i,gassrand());
	return a*rr;
}

matrix gassvector::operator()()	// 返回已设定的协方差矩阵的正态随机向量
{
	matrix rr(a.rownum);
	for(size_t i=0; i<a.rownum; i++)
		rr.set(i,gassrand());
	return a*rr;
}

void gassvector::setdata(matrix & r) // 根据协方差矩阵设置增益矩阵
{
	if(!r.issym) throw TMESSAGE("r must be symetic!");
	a = r;
	matrix evalue;
	a = a.eigen(evalue);
	for(size_t i=0; i<a.colnum; i++) {
   	if(evalue(i)<0.0) throw TMESSAGE("var matrix not positive!");
		DOUBLE e = sqrt(evalue(i));
		for(size_t j=0; j<a.rownum; j++)
			a.set(j,i,a.value(j,i)*e);
	}
}

int matrix::ispositive()		// 判定矩阵是否对称非负定阵，如是返回1，否则返回0
{
	if(!issym) return 0;
	matrix ev;
	eigen(ev);
	for(size_t i=0; i<rownum; i++)
		if(ev(i)<0) return 0;
	return 1;
}

matrix matrix::cholesky(matrix& dd)	// 用乔里斯基分解法求对称正定阵的线性
		// 方程组ax=d返回方程组的解，本身为a，改变为分解阵u,d不变
{
	if(!issym) throw TMESSAGE("not symetic!");
	if(dd.rownum != colnum) throw TMESSAGE("dd's rownum not right!");
	matrix md(dd);
	size_t i,j,k,u,v;
	if(value(0,0)<=0.0) throw TMESSAGE("not positive");
	set(0,0,sqrt(value(0,0))); //	 a[0]=sqrt(a[0]);
	buffer& a = (*buf);
	buffer& d = (*(md.buf));
	size_t n = rownum;
	size_t m = dd.colnum;
	for (j=1; j<n; j++) a[j]=a[j]/a[0];
	for (i=1; i<n; i++)
		{ u=i*n+i;
		  for (j=1; j<=i; j++)
			 { v=(j-1)*n+i;
				a[u]=a[u]-a[v]*a[v];
			 }
		  if (a[u]<=0.0) throw TMESSAGE("not positive");
		  a[u]=sqrt(a[u]);
		  if (i!=(n-1))
			 { for (j=i+1; j<n; j++)
				  { v=i*n+j;
					 for (k=1; k<=i; k++)
						a[v]=a[v]-a[(k-1)*n+i]*a[(k-1)*n+j];
					 a[v]=a[v]/a[u];
				  }
			 }
		}
	for (j=0; j<m; j++)
		{ d[j]=d[j]/a[0];
		  for (i=1; i<=n-1; i++)
			 { u=i*n+i; v=i*m+j;
				for (k=1; k<=i; k++)
				  d[v]=d[v]-a[(k-1)*n+i]*d[(k-1)*m+j];
				d[v]=d[v]/a[u];
			 }
		}
	for (j=0; j<=m-1; j++)
		{ u=(n-1)*m+j;
		  d[u]=d[u]/a[n*n-1];
		  for (k=n-1; k>=1; k--)
			 { u=(k-1)*m+j;
				for (i=k; i<=n-1; i++)
				  { v=(k-1)*n+i;
					 d[u]=d[u]-a[v]*d[i*m+j];
				  }
				v=(k-1)*n+k-1;
				d[u]=d[u]/a[v];
			 }
		}
	if(n>1)
	for(j=1; j<n; j++)
	for(i=0; i<j; i++)
		a[i+j*n]=0.0;
	return md;
}

DOUBLE lineopt(matrix& aa,matrix& bb, matrix& cc, matrix & xx)
 // 线性规划最优点寻找程序，aa为mXn不等式约束条件左端系数矩阵，bb为不等式约束
 // 条件的右端项，为m维向量，cc为目标函数系数，n维向量，xx返回极小点，为n维向量
{
	if(aa.rownum != bb.rownum || aa.colnum != cc.rownum ||
		aa.colnum != xx.rownum) throw TMESSAGE("dimenstion not right!");
	size_t n=aa.colnum, m=aa.rownum;
	size_t i,mn,k,j;
	matrix a(m,n+m);
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++)
			a.set(i,j,aa(i,j));
		for(j=n;j<n+m;j++)
			if(j-n == i) a.set(i,j,1.0);
			else a.set(i,j,0.0);
	}
	matrix c(m+n);
	c = 0.0;
	for(i=0;i<m;i++)
		c.set(i,cc(i));
	lbuffer* jjs = getnewlbuffer(m);
	lbuffer& js = (*jjs);
	DOUBLE s,z,dd,y; //,*p,*d;

	for (i=0; i<m; i++) js[i]=n+i;
	matrix p(m,m);
	matrix d;
	mn=m+n; s=0.0;
	matrix x(mn);
	while (1)
		{ for (i=0; i<m; i++)
			 for (j=0; j<m; j++)
				p.set(i,j,a(i,(size_t)js[j]));
		  p.inv();
			d = p*a;
			x = 0.0;
		  for (i=0; i<m; i++)
			 { s=0.0;
				for (j=0; j<=m-1; j++)
					s+=p(i,j)*bb(j);
				x.set((size_t)js[i],s);
			 }
		  k=mn; dd=1.0e-35;
		  for (j=0; j<mn; j++)
			 { z=0.0;
				for (i=0; i<m; i++)
					z+=c((size_t)js[i])*d(i,j);
				z-=c(j);
				if (z>dd) { dd=z; k=j;}
			 }
		  if (k==mn)
			 { s=0.0;
				for (j=0; j<n; j++) {
					xx.set(j,x(j));
					s+=c(j)*x(j);
				}
				delete jjs;
				return s;
			 }
		  j=m;
		  dd=1.0e+20;
		  for (i=0; i<=m-1; i++)
			 if (d(i,k)>=1.0e-20)   // [i*mn+k]>=1.0e-20)
				{ y=x(size_t(js[i]))/d(i,k);
				  if (y<dd) { dd=y; j=i;}
				}
		  if (j==m) { delete jjs;
							throw TMESSAGE("lineopt failed!");
						}
		  js[j]=k;
		}
}
