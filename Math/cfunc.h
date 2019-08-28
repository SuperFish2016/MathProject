#ifndef CFUNC_H
#define CFUNC_H

#include <math.h>

#include "matrix.h"
#include "cmatrix.h"
#include "func.h"

class calgo	// 复数算法类,专门生成自变量为实数，运算结果为复数的算法
{
 private:
	COMPLEX yfactor;		// 乘因子，初始化为1
	DOUBLE xfactor;		// x轴放大因子，初始化为1
	COMPLEX addconst;	// 加和，初始化为0
	DOUBLE xshift;		// x平移量，初始化为0
 public:
	unsigned refnum;		// 引用数,初始化为1
    calgo():refnum(1),yfactor(1.0),xfactor(1.0),addconst(0.0),xshift(0.0){}
	 // 构造函数,产生y=x线性函数
	calgo(DOUBLE xs, DOUBLE xf,COMPLEX adc=0.0, COMPLEX yf=1.0):refnum(1),
        yfactor(yf),addconst(adc),xshift(xs),xfactor(xf){}
	calgo(COMPLEX a):refnum(1),yfactor(0.0),xfactor(1.0),addconst(a),xshift(0.0)
        {}	 // 常函数的构造
    calgo(const calgo & alg):yfactor(alg.yfactor),xfactor(alg.xfactor),
        addconst(alg.addconst),xshift(alg.xshift),refnum(1){} // 拷贝构造函数
    virtual ~calgo(){} // 虚析构函数
	COMPLEX cal(DOUBLE x); // 计算算法值
	virtual COMPLEX calculate(DOUBLE x)
        {return x;} // 本身算法,将被继承子类改写
	virtual calgo * clone(); // 克隆自己，必须被继承子类改写
	calgo * mul(COMPLEX a);	// 乘a
	calgo * add(COMPLEX a);	// 加a
	calgo * neg(); // 取负
	calgo * setxfactor(DOUBLE x);		// 设置x轴因子
	calgo * xroom(DOUBLE x);	// 将xfactor扩大x倍
	calgo * setxshift(DOUBLE x); // 设置xshift的值
	calgo * xshiftas(DOUBLE x);	// 从当前开始平移x
};

class calgojoin : public calgo // 结合算法
{
 public:
	calgo * leftalgo;	// 左算法，初始化为0
	calgo * rightalgo; // 右算法，初始化为0
	method met;	// 指明算法
	calgojoin(calgo * l, calgo * r, method m):leftalgo(l),
		rightalgo(r), met(m)
		{ if(leftalgo)
				leftalgo->refnum++;
		  if(rightalgo)
				rightalgo->refnum++;
		};
	calgojoin(calgojoin& alg):calgo(alg),
		leftalgo(alg.leftalgo),rightalgo(alg.rightalgo),met(alg.met){
			if(leftalgo)
				leftalgo->refnum++;
			if(rightalgo)
				rightalgo->refnum++;};
			// 拷贝构造函数
	virtual ~calgojoin() {
		if(leftalgo) {	// 如左或者右算法已经没有被引用，则删除
			leftalgo->refnum--;
			if(!leftalgo->refnum) delete leftalgo;
		}
		if(rightalgo) {
			rightalgo->refnum--;
			if(!rightalgo->refnum) delete rightalgo;
		}
	};
	virtual calgo * clone(); // 克隆自己
	virtual COMPLEX calculate(DOUBLE x);	// 实施结合算法
};

class calgofun : public calgo	// 函数算法
{
 public:
	COMPLEX (*f)(DOUBLE);	// 函数指针
	calgofun(COMPLEX (*fun)(DOUBLE)):f(fun){};	// 用函数指针进行初始化
	calgofun(calgofun& alg):calgo(alg),f(alg.f){}; // 拷贝构造函数
	virtual COMPLEX calculate(DOUBLE x);	// 实施函数算法
	virtual calgo * clone(); // 克隆自己
};

class calgopoly : public calgo	// 复多项式
{
 public:
	cmatrix data;	//	n乘1矩阵，存放n-1次多项式的系数a(0)到a(n-1)
	calgopoly(cmatrix& d):data(d){};	// 用矩阵构造多项式
	calgopoly(calgopoly& alg):calgo(alg),data(alg.data){}; // 拷贝构造函数
	virtual COMPLEX calculate(DOUBLE x);	// 实施函数算法
	virtual calgo * clone();	// 克隆自己
};

class calgoenter : public calgo	// 等间隔插值算法
{
 public:
	algoenter * er;
	algoenter * ei;	// 两个插值算法构成结果的实部与虚部
    calgoenter(const matrix& ver, const matrix& vei, DOUBLE x0, DOUBLE h):
        er(new algoenter(ver,x0,h)), ei(new algoenter(vei,x0,h)) {}
    calgoenter(const calgoenter& alg):calgo(alg),er(alg.er),ei(alg.ei) {
		er->refnum++; ei->refnum++;} // 拷贝构造函数
	virtual ~calgoenter(){
		er->refnum--;
		if(!er->refnum) delete er;
		ei->refnum--;
        if(!ei->refnum) delete ei;}
	virtual COMPLEX calculate(DOUBLE x);	// 实施函数算法
	virtual calgo * clone();	// 克隆自己
};

class cfunc {	// 复函数类,计算自变量为实数而结果是复数的类
 public:
	calgo * alg;	// 决定函数的算法

	cfunc();	// 缺省构造函数
	cfunc(COMPLEX a);	// 常函数的构造函数
	cfunc(COMPLEX (*fun)(DOUBLE));	// 函数指针的构造函数
    cfunc(const cfunc & fn);	// 拷贝构造函数
	cfunc(calgo * a):alg(a){} // 算法构造函数，使用要小心，不能将一个算法产生
									// 两个函数，除非自己控制引用数的增加
	cfunc(cmatrix& m);	// 构造多项式，m为nX1矩阵，是n-1阶多项式系数，
				// 其中m(0,0)为常数项，m(n-1,0)为n-1次项。

	virtual ~cfunc() {		// 析构函数
		if(alg) {
			alg->refnum--;	// 引用数减一，如再无其它引用，则删除算法
			if(!alg->refnum)
				delete alg;
		}
	};

	COMPLEX operator()(DOUBLE x){return alg->cal(x);}; // 计算x的函数值
	cfunc& operator=(cfunc& fn);	// 赋值运算符
	cfunc& operator=(COMPLEX (*fn)(DOUBLE)); // 用函数指针的赋值运算符
	cfunc& operator=(COMPLEX a); // 常函数的赋值运算符

	cfunc& operator+=(cfunc& fn);	// 自身加一个函数
	cfunc& operator+=(COMPLEX a){alg=alg->add(a);return (*this);};
			//自身加一个常数
	cfunc& operator+=(COMPLEX (*f)(DOUBLE)); // 自身加一个函数指针
	cfunc operator+(cfunc& fn);	// 相加产生新函数
	cfunc operator+(COMPLEX a);	// 与常数相加产生新函数
	friend cfunc operator+(COMPLEX a, cfunc& f); // 同上但常数在前
	cfunc operator+(COMPLEX (*f)(DOUBLE)); // 加一个函数指针产生新函数
	friend cfunc operator+(COMPLEX (*f)(DOUBLE),cfunc& fn);
		// 同上但函数指针在前

	cfunc& neg(); // 自身取负
	cfunc operator-(); // 产生负函数

	cfunc& operator-=(cfunc& fn); // 自身减一个函数
	cfunc& operator-=(COMPLEX a){alg=alg->add(-a);return (*this);};
			//自身减一个常数
	cfunc& operator-=(COMPLEX (*f)(DOUBLE)); // 自身减一个函数指针
	cfunc operator-(cfunc& fn);	// 相减产生新函数
	cfunc operator-(COMPLEX a);	// 与常数相减产生新函数
	friend cfunc operator-(COMPLEX a, cfunc& f); // 同上但常数在前
	cfunc operator-(COMPLEX (*f)(DOUBLE)); // 减一个函数指针产生新函数
	friend cfunc operator-(COMPLEX (*f)(DOUBLE),cfunc& fn); // 函数指针减函数

	cfunc& operator*=(cfunc& fn);	// 自身乘一个函数
	cfunc& operator*=(COMPLEX a){alg=alg->mul(a);return (*this);};
		//自身乘一个常数
	cfunc& operator*=(COMPLEX (*f)(DOUBLE)); // 自身乘一个函数指针
	cfunc operator*(cfunc& fn);	// 相乘产生新函数
	cfunc operator*(COMPLEX a);	// 与常数相乘产生新函数
	friend cfunc operator*(COMPLEX a, cfunc& f); // 同上但常数在前
	cfunc operator*(COMPLEX (*f)(DOUBLE)); // 乘一个函数指针产生新函数
	friend cfunc operator*(COMPLEX (*f)(DOUBLE),cfunc& fn); // 函数指针乘函数

	cfunc& operator/=(cfunc& fn);	// 自身除以一个函数
	cfunc& operator/=(COMPLEX a){alg=alg->mul(1.0/a);return (*this);
			};//自身除以常数
	cfunc& operator/=(COMPLEX (*f)(DOUBLE)); // 自身除以一个函数指针
	cfunc operator/(cfunc& fn);	// 相除产生新函数
	cfunc operator/(COMPLEX a);	// 与常数相除产生新函数
	friend cfunc operator/(COMPLEX a, cfunc& f); // 常数除以函数
	cfunc operator/(COMPLEX (*f)(DOUBLE)); // 除以一个函数指针产生新函数
	friend cfunc operator/(COMPLEX (*f)(DOUBLE),cfunc& fn); // 函数指针除以函数

	void setxfactor(DOUBLE a);	// 设置x因子为a
	void xroom(DOUBLE a);	  // x方向扩大a倍
	void setxshift(DOUBLE a);	// 设置函数沿x轴平移a
	void shiftxas(DOUBLE a); // 函数沿x轴右移a

	cfunc& power(cfunc& f);	// 函数的f次乘幂，函数自身改变
	cfunc& power(COMPLEX a);	// 函数的a次幂，函数自身改变
	cfunc operator^(cfunc & fn);	// 函数的fn次乘幂，产生新函数，原函数不变
	cfunc operator^(COMPLEX a);  // 函数的a次幂，产生新函数，原函数不变
};

inline cfunc operator+(COMPLEX a, cfunc& f) // 常数加函数
{	return f+a; }

inline cfunc operator+(COMPLEX (*f)(DOUBLE),cfunc& fn) // 函数指针加函数
{	return fn+f;}

cfunc operator-(COMPLEX a, cfunc& f); // 常数减函数
cfunc operator-(COMPLEX (*f)(DOUBLE),cfunc& fn); // 函数指针减函数

inline cfunc operator*(COMPLEX a, cfunc& f) // 常数乘函数
{	return f*a; }
inline cfunc operator*(COMPLEX (*f)(DOUBLE),cfunc& fn) // 函数指针乘函数
{	return fn*f;}

cfunc operator/(COMPLEX a, cfunc& f); // 常数除以函数
cfunc operator/(COMPLEX (*f)(DOUBLE),cfunc& fn); // 函数指针除以函数

class cfuncenter: public cfunc {	// 插值复函数
 public:
	cfuncenter(cmatrix& s,DOUBLE t0, DOUBLE dt);
};

cfunc fourier(func& f,DOUBLE tb, DOUBLE te, DOUBLE dt, DOUBLE df);
	// 利用fft技术对函数f作傅里叶变换，其中tb为采样窗口的起始点，te为结束点，
	// 必须te>tb，dt为采样间隔，df为频率采样间隔，但返回的cfunc是作了插值的
	// 插值函数

#endif // CFUNC_H
