#ifndef BUFFER_H
#define BUFFER_H

// 功能：定义数据缓冲区各抽象类
#include <ostream>
using namespace std;

#ifndef DOUBLE

#define DOUBLE double
// #define DOUBLE long double  // 如double精度不够，则可将此行去注释，而将上行注释
				// 即可得到更精确的结果
#endif // DOUBLE

#include <stdio.h>

extern DOUBLE defaulterr;	// 缺省的误差标准，为十的负八次方
extern int doadjust;		// 逻辑变量，为1则作调整，否则不作，缺省为1
DOUBLE adjust(DOUBLE & a);  // 将实数调整为最靠近的小数点后二位数

class buffer {  // 定义抽象的双精度实数缓冲区类，此类以一维数组的形式专门存储大量的
		// 双精度实数，而存储方式则由子类实现，即可以用内存，也可以用磁盘
		// 空间，甚至服务器上的存储空间
 public:
	size_t refnum;	// 引用数，一个缓冲区可能有多个矩阵同时引用它，
			// 主要是在矩阵变量赋值导致多个矩阵具有同样的内容时
			// 防止数据的大量拷贝，而引用数就说明现在有多少个
			// 矩阵正在引用这个缓冲区
	buffer():refnum(1){};		// 构造函数，预设引用数为1
	virtual void alloc(size_t num)=0;  // 申请容量为num的存储空间，或将缓冲区容量
					//改为num，纯虚函数，子类必须实现
	virtual void release()=0;	// 释放占用的存储空间
	virtual DOUBLE& retrieve(size_t n)=0;  // 取出第n个双精度实型数据，n从0开始计算
	virtual buffer* clone()=0; // 克隆出一个和自己类别一样，内容也一样的新的缓冲区，
				// 为纯虚函数，子类必须实现
	virtual buffer* newbuf(size_t n=0)=0; // 产生一个新的和自己类别一样的具有指定
					// 尺寸的缓冲区，
	DOUBLE & operator[](size_t n){	// 重载数组操作符(), 检索第n个字节的内容
		return retrieve(n);};
	virtual size_t len()=0;		// 返回数组的长度，纯虚函数，必须由子类
					// 实现
	buffer* set(size_t n, DOUBLE d);
	// 将实数缓冲区的第n个实数的值设为d，并返回一个实数缓冲区的
	// 指针，此指针一般就是指向自己，而在一个缓冲区的引用数多于一个时
	// 此函数将先克隆出一个新的缓冲区，此新的缓冲区的引用数为1，并对此
	// 新的缓冲区的第n个实数的值设为d，然后返加新缓冲区的指针
	friend ostream& operator<<(ostream& o, buffer& b);
};

ostream& operator<<(ostream& o, buffer& b);
buffer * getnewbuffer(size_t n=0);  // 生成一个新的尺寸为n个实数的缓冲区，返回新产生的
				// 缓冲区的指针，新产生的缓冲区的类别与下面的
				// defbuffer指针指向的缓冲区的类别一致
extern buffer *defbuffer;	// 指向buffer子类的变量，作为产生新的缓冲区的缺省类别
void settomemory(); // 将defbuffer指向内存类，这以后再调用getnewbuffer将产生内存类的
		// 实数数组，即数组全部放在内存中
void settodisk(); // 将defbuffer指向磁盘类，这以后再调用getnewbuffer将产生数据存放在文件
		// 中的实数数组，即数组全部放在数据文件中

class membuffer : public buffer { // 一维的存放在内存中的双精度实数数组，请注意在DOS
			// 操作系统下将受到64k内存段容量限制，而如果在win95以上
			// 操作系统中用32位编译则没有这个限制。如果要在DOS下使用
			// 扩展内存或扩充内存，可以另外编写buffer的子类
 public:
	DOUBLE * buf;	// 指向存放双精度实数数据的内存缓冲区的指针
	size_t length;	// 数组的以实数个数计的长度
	membuffer():buffer(),buf(0),length(0){}; // 构造函数，产生一个容量为0的实数缓冲区
	membuffer(size_t lth):buffer(),length(lth){ // 构造函数，产生长度为lth的实数缓冲区
		if(lth > 0)
			buf = new DOUBLE[lth];
		else buf = 0; };
	~membuffer(){delete buf;}  // 析构函数，释放缓冲区
	void alloc(size_t num){  // 将缓冲区的容量改为num个实数，原来的内容全部删去
		delete []buf;
		length = num;
		buf = new DOUBLE[num]; };
	DOUBLE& retrieve(size_t n){ // 返回第n个字节，n的范围从0到length-1
		return buf[n]; };
	void release() { delete []buf; length = 0; };  // 释放内存，数组的长度变为0
	buffer* clone(); // 克隆自己
	buffer* newbuf(size_t n=0){ // 产生一个长度为n的与自己同类别的新缓冲区并返回
				// 相应的指针
		return new membuffer(n);};
	size_t len(){return length;}; // 返回数组的长度
};

class diskbuffer : public buffer { // 磁盘实数缓冲区类，将实数数组存放在临时文件中
 public:
	DOUBLE buf;	// 临时文件取出来的实数放在buf里
	size_t length;	// 缓冲区的长度，或者实数的个数
	size_t n;	// 当前的buf中存放的是第n个实数，也以此知道临时文件的当前指针位置
	FILE * tempfp;    // 指向临时文件的指针
	diskbuffer():buffer(),buf(0.0),length(0),tempfp(0),n(0){}; // 缺省构造函数，产生长度为0的
				// 磁盘实数缓冲区
	diskbuffer(size_t lth);	// 构造函数，产生长度为lth的磁盘实数缓冲区
	~diskbuffer();		// 析构函数，关闭并删除临时文件
	void alloc(size_t num);	// 将缓冲区的尺寸重定为num
	DOUBLE& retrieve(size_t n);	// 检索出文件的第n个实数放在buf中并返回buf的引用
				// 在检索之前先把当前的buf中的内容放回到文件中
	void release();		// 释放空间
	buffer* clone(); // 克隆自己
	buffer* newbuf(size_t n=0){return new diskbuffer(n);}; // 产生一个新的磁盘实数缓冲区
	size_t len(){return length;}; // 获得数组的长度
};

class lbuffer {	// 长整数缓冲区的抽象基类，用来存放许多个长整数
		// 继承的子类可以使用各种数据载体，如内存或者文件来存放长整数
 public:
	size_t refnum;	// 缓冲区引用数，原理与buffer的相同
	lbuffer():refnum(1){};	// 缺省构造函数
	virtual void alloc(size_t num)=0; // 将缓冲区尺寸修改为num,原内容全部删去
	virtual void release()=0; // 释放缓冲区
	virtual long& retrieve(size_t n)=0; // 检索第n个元素，并返回引用
	long & operator[](size_t n){ // 重载[]操作符，得到第n个元素，并返回引用
										// 但这样调用有危险，要注意引用数为1时才可当左值
		return retrieve(n);};
	virtual lbuffer* clone()=0; // 克隆一个类别和内容与自己完全一样的长整数缓冲区
	virtual lbuffer* newbuf(size_t n=0)=0; // 产生一个类别和自己一样的长整数缓冲区
	virtual size_t len()=0; // 返回长度
	lbuffer* set(size_t n, long v); // 设置第n个元素的值为v，n的范围是0到refnum-1
				// 如引用数大于1，则会产生新的缓冲区进行操作
				// 并返回新缓冲区指针
	friend ostream& operator<<(ostream& o, lbuffer& b);
};

ostream& operator<<(ostream& o, lbuffer& b);

class lmembuffer : public lbuffer { // 内存长整数缓冲区，DOS系统将受到64K段限制
 public:
	long * buf;	// 指向内存缓冲区的指针
	size_t length;	// 缓冲区长度
	lmembuffer():lbuffer(),buf(0),length(0){}; // 缺省构造函数，长度为0
	lmembuffer(size_t lth):lbuffer(),length(lth){ // 构造函数，lth为缓冲区长整数个数
		buf = new long[lth]; };
	~lmembuffer(){delete []buf;}	// 析构函数，释放buf指向的内存
	void alloc(size_t num){ // 将缓冲区长度修改为num，原内容消失
		delete []buf;
		length = num;
		buf = new long[num]; };
	long& retrieve(size_t n){ // 返回第n个长整数，n取值从0到length-1
		return buf[n]; };
	void release() { delete []buf; length = 0; }; // 释放内存，长度变为0
	lbuffer* clone(); // 克隆一个内容完全一样的长整数缓冲区
	lbuffer* newbuf(size_t n=0){return new lmembuffer(n);}; // 产生一个新的尺寸为n的缓冲区
	size_t len(){return length;}; // 返回缓冲区的长度
};

class ldiskbuffer : public lbuffer { // 长整数磁盘缓冲区子类，将多个长整数存放在临时文件内
 public:
	long buf;	// 当前的长整数
	size_t length; // 缓冲区长度，以长整数个数为单位
	size_t n;	// buf中的内容是第n个元素的内容，n也说明了临时文件当前的指向
	FILE * tempfp; // 临时文件指针
	ldiskbuffer():lbuffer(),buf(0),length(0),tempfp(0),n(0){}; // 缺省构造函数，长度为0
	ldiskbuffer(size_t lth); // 构造函数，lth为缓冲区长度
	~ldiskbuffer();	// 析构函数，将关闭和删除临时文件
	void alloc(size_t num); // 将缓冲区长度重新定为num
	long& retrieve(size_t nn); // 检索第nn个元素放到buf中，在此之前先存放buf的内容到
			// n的位置，再将文件的指针移到nn处取元素
	void release();	// 释放存储空间
	lbuffer* clone();	// 克隆一个内容完全一样的长整数磁盘缓冲区，但临时文件不同
	lbuffer* newbuf(size_t nn=0){return new ldiskbuffer(nn);}; // 产生一个长度为nn的
			// 长整数磁盘缓冲区子类变量，返回它的指针
	size_t len(){return length;}; // 返回缓冲区长度
};

lbuffer * getnewlbuffer(size_t n=0); // 产生一个新的长整数缓冲区，子类的类别与下面的
		// deflbuffer指向的缓冲区的类别相同
extern lbuffer *deflbuffer; // 缺省的长整数缓冲区指针

char * throwmessage(int l, char * f, char * c); // 产生错误信息字符串

#ifndef __LINE__
#define __LINE__ 0
#define __FILE__ "n"
#endif
	//  l为行号,f为文件名，c为出错信息，请不要直接调用，而是调用宏TMESSAGE
#define TMESSAGE(c) throwmessage(__LINE__,__FILE__,c)

#endif // BUFFER_H
