#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>

typedef int mode;
#define all 0
#define left 1
#define right 2
typedef int space;
#define space_col 0
#define space_row 1
#define space_nul 2
#define space_rank 3

typedef struct __fraction{
	int numerator;//分子
    int denominator;//分母
}fraction,*Fraction;

bool fraction_ifexist(Fraction a){
	/*作用：判断分数是否有效 
      若分母为0则返回0否则返回1*/
	return (a->denominator == 0) ? false : true;
};

int fraction_ifproper(Fraction a){
	/*作用：判断分数是否为真分数 
	  若分数是真分数则返回1否则返回0*/ 
	if(!fraction_ifexist(a)){
		printf("错误：分母为0");
		exit(-1);
	}
	return (a->denominator > a->numerator) ? true : false;
};

int gcd(int a,int b){
	//作用：求a与b最大公约数 
	b = abs(b);
	a = abs(a);
	if(b > a){
		int t = a;
		a = b;
		b = t;
	}
	return (b == 0) ? a : gcd(b,a%b);
}

void fraction_simplify(Fraction a){
	//作用：对a约分
	if(!fraction_ifexist(a)){
		printf("错误：分母为0");
		exit(-1);
	}
	int GCD = gcd(a->numerator,a->denominator);
	a->numerator /= GCD;
	a->denominator /= GCD;
	if(a->denominator < 0){
		a->numerator *= -1;
		a->denominator *= -1;
	}
	return;
}

int lcm(int a,int b){
	//作用：求a与b的最小公倍数 
	return abs(a)*abs(b)/gcd(a,b);
}

void fraction_complicated(Fraction a,Fraction b){
	//作用：a与b通分
	if(!fraction_ifexist(a) || !fraction_ifexist(b)){
		printf("错误：分母为0");
		exit(-1);
	}
	int LCM = lcm(a->denominator,b->denominator);
	a->numerator *= LCM/a->denominator;
	b->numerator *= LCM/b->denominator;
	a->denominator = b->denominator = LCM;
	return;
} 

Fraction fraction_add(Fraction a,Fraction b){
	//作用：加法a+b
	if(!fraction_ifexist(a) || !fraction_ifexist(b)){
		printf("错误：分母为0");
		exit(-1);
	}
	Fraction r = (Fraction)malloc(sizeof(fraction));
	fraction_complicated(a,b);
	r->numerator = a->numerator + b->numerator;
	r->denominator = a->denominator;
	fraction_simplify(r);
	return r;
}

Fraction fraction_sub(Fraction a,Fraction b){
	//作用：减法a-b 
	if(!fraction_ifexist(a) || !fraction_ifexist(b)){
		printf("错误：分母为0");
		exit(-1);
	}
	Fraction r = (Fraction)malloc(sizeof(fraction));
	fraction_complicated(a,b);
	r->numerator = a->numerator - b->numerator;
	r->denominator = a->denominator;
	fraction_simplify(r);
	return r;
}

Fraction fraction_mul(Fraction a,Fraction b){
	//作用：乘法a*b
	if(!fraction_ifexist(a) || !fraction_ifexist(b)){
		printf("错误：分母为0");
		exit(-1);
	}
	Fraction r = (Fraction)malloc(sizeof(fraction));
	r->numerator = a->numerator * b->numerator;
	r->denominator = a->denominator * b->denominator;
	fraction_simplify(r);
	return r;
}

Fraction fraction_div(Fraction a,Fraction b){
	//作用：除法a/b
	 if(!fraction_ifexist(a) || !fraction_ifexist(b)){
		printf("错误：分母为0");
		exit(-1);
	}
	Fraction r = (Fraction)malloc(sizeof(fraction));
	r->numerator = a->numerator * b->denominator;
	r->denominator = a->denominator * b->numerator;
	fraction_simplify(r);
	return r;
}

Fraction fraction_opposite(Fraction a){
	//作用：求分数a的相反数 
	fraction b;
	b.numerator = -1;
	b.denominator = 1;
	return fraction_mul(a,&b);
}

int fraction_abscmp(Fraction a,Fraction b){
	/*作用：比较分数a与b绝对值的大小
      若a < b返回-1
	  若a = b返回0
	  若a > b返回1*/
    fraction_complicated(a,b);
    if(fabs(a->numerator) < fabs(b->numerator))return -1;
    else if(fabs(a->numerator) == fabs(b->numerator))return 0;
    else return 1;
}

void fraction_print(Fraction a){
	//作用：对a约分并输出 
	 if(!fraction_ifexist(a)){
		printf("错误：分母为0");
		exit(-1);
	}
	fraction_simplify(a);
	if(a->denominator == 1 || a->numerator == 0){
		a->denominator = 1;
		printf("%+-8d",a->numerator);
	}else printf("%+d/%-8d",a->numerator,a->denominator); 
	return;
} 

Fraction fraction_create(int a,int b){
	//作用：返回一个以a为分子b为分母的分数
    Fraction r = malloc(sizeof(fraction));
    r->numerator = a;
    r->denominator = b;
    return r;
}

Fraction int_create(int a){
	//作用：返回一个a的整数(分母为1)
	Fraction r = malloc(sizeof(fraction));
    r->numerator = a;
    r->denominator = 1;
    return r;
}
///////////////////////////////////////////////////////
//////////////////////////////////////////////////////
/////////////////////////////////////////////////////
typedef struct __matrix{
	Fraction **base;//矩阵首元素
    int row;//矩阵行数
    int col;//矩阵列数
}matrix,*Matrix;

void matrix_null(Matrix a,int m,int n){
	//作用：创建m×n的零矩阵
	if(m < 0 || n < 0){
		printf("错误：负数");
		exit(-1);
	}
	a->base = (Fraction**)malloc(sizeof(Fraction*)*m);
	for(int i = 0;i < m;i++)
		a->base[i] = (Fraction*)malloc(sizeof(Fraction)*n);
	a->row = m;
	a->col = n;
	for(int i=0;i < a->row;i++)
		for(int j=0;j < a->col;j++)
			a->base[i][j] = fraction_create(0,1);
    return;
}

void matrix_identity(Matrix a,int n){
    //作用：创建n阶的单位矩阵
	if(n < 0){
		printf("错误：负数");
		exit(-1);
	}
    a->base = (Fraction**)malloc(sizeof(Fraction*)*n);
	for(int i = 0;i < n;i++)
		a->base[i] = (Fraction*)malloc(sizeof(Fraction)*n);
	a->row = n;
	a->col = n;
	for(int i=0;i < a->row;i++)
		for(int j=0;j < a->col;j++)
			if(i != j) a->base[i][j] = fraction_create(0,1);
			else a->base[i][j] = fraction_create(1,1);
	return;
}

void matrix_square(Matrix a,int n){
    //作用：创建n阶的零矩阵
	if(n < 0){
		printf("错误：负数");
		exit(-1);
	}
    a->base = (Fraction**)malloc(sizeof(Fraction*)*n);
	for(int i = 0;i < n;i++)
		a->base[i] = (Fraction*)malloc(sizeof(Fraction)*n);
	a->row = n;
	a->col = n;
	for(int i=0;i < a->row;i++)
		for(int j=0;j < a->col;j++)
			a->base[i][j] = fraction_create(0,1);
	return;
}

void matrix_edelete(matrix a,int i,int j){
    //作用：删除元素a[i][j]
    free(a.base[i][j]);
    return;
}

void matrix_set(matrix a,int i,int j,Fraction value){
	//作用：将a[i][j]修改为value
	if(!(((i >= 0) && (i < a.row)) && ((j >= 0) && (j < a.col)))){
		printf("错误：位置错误");
		exit(-1);
	}
	matrix_edelete(a,i,j);
	a.base[i][j] = value;
	return;	
}

void matrix_print(matrix a){
    //作用：输出矩阵
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			fraction_print(a.base[i][j]);
		}
		printf("\n");
	}
}

void matrix_delete(matrix a){
	//作用：删除矩阵 
	for(int i = 0;i < a.row;i++){
		free(a.base[i]);
	}
	free(a.base);
	a.base = NULL;
	a.col = a.row = 0;
	return;
}

bool matrix_ifrow(matrix a,int i){
	//作用：判断a矩阵是否有第i行 
	return (i >=0 && i < a.row) ? true : false;
}

bool matrix_ifcol(matrix a,int j){
	//作用：判断a矩阵是否有第j列 
	return (j >= 0 && j < a.col) ? true : false;
}

void matrix_radd(matrix a,int i1,int i2){
	//作用：i1行加i2行
	if(!matrix_ifrow(a,i1) || !matrix_ifrow(a,i2)){
		printf("错误：超出行数");
		exit(-1);
	} 
	for(int j = 0;j < a.col;j++){
		a.base[i1][j] = fraction_add(a.base[i1][j],a.base[i2][j]);
	}
	return;
}

void matrix_rreplace(matrix a,int i1,int i2){
	//作用：交换i1行、i2行
	if(!matrix_ifrow(a,i1) || !matrix_ifrow(a,i2)){
		printf("错误：超出行数");
		exit(-1);
	} 
	for(int j = 0;j < a.col;j++){
		Fraction t = a.base[i1][j];
		a.base[i1][j] = a.base[i2][j];
		a.base[i2][j] = t;
	}
	return;
}

void matrix_rmul(matrix a,int i,Fraction C){
	//作用：i行乘上分数C
	if(!matrix_ifrow(a,i)){
		printf("错误：超出行数");
		exit(-1);
	} 
	for(int j = 0;j < a.col;j++){
		a.base[i][j] = fraction_mul(a.base[i][j],C);
	}
	return;
}

void matrix_rcadd(matrix a,int i1,Fraction C,int i2){
	//作用：i1行加C倍的i2行
	if(!matrix_ifrow(a,i1) || !matrix_ifrow(a,i2)){
		printf("错误：超出行数");
		exit(-1);
	} 
	for(int j = 0;j < a.col;j++){
		a.base[i1][j] = fraction_add(a.base[i1][j],fraction_mul(a.base[i2][j],C));
	}
	return;
}

void matrix_cadd(matrix a,int j1,int j2){
	//作用：j1列加j2列
	if(!matrix_ifcol(a,j1) || !matrix_ifrow(a,j2)){
		printf("错误：超出列数");
		exit(-1);
	} 
	for(int i = 0;i < a.row;i++){
		a.base[i][j1] = fraction_add(a.base[i][j1],a.base[i][j2]);
	}
	return;
}

void matrix_creplace(matrix a,int j1,int j2){
	//作用：交换j1列、j2列
	if(!matrix_ifcol(a,j1) || !matrix_ifcol(a,j2)){
		printf("错误：超出列数");
		exit(-1);
	} 
	for(int i = 0;i < a.row;i++){
		Fraction t = a.base[i][j1];
		a.base[i][j1] = a.base[i][j2];
		a.base[i][j2] = t;
	}
	return;
}

void matrix_ccmul(matrix a,int j,Fraction C){
	//作用：j列乘上分数C
	if(!matrix_ifcol(a,j)){
		printf("错误：超出列数");
		exit(-1);
	} 
	for(int i = 0;i < a.row;i++){
		a.base[i][j] = fraction_mul(a.base[i][j],C);
	}
	return;
}

void matrix_ccadd(matrix a,int j1,Fraction C,int j2){
	//作用：j1列加C倍的j2列
	if(!matrix_ifcol(a,j1) || !matrix_ifcol(a,j2)){
		printf("错误：超出列数");
		exit(-1);
	} 
	for(int i = 0;i < a.row;i++){
		a.base[i][j1] = fraction_add(a.base[i][j1],fraction_mul(a.base[i][j2],C));
	}
	return;
}

void matrix_input(matrix a,int* n,int* d){
	/*作用：按行输入分数（以n为一组分子，d为一组分母）至a
	  注，这个函数并不稳定*/
    int count = 0;
    for(int i = 0;i < a.row;i++){
    	for(int j = 0;j < a.col;j++){
    		matrix_set(a,i,j,fraction_create(n[count],d[count]));
    		count++;
		}
	}
	return;
}

void matrix_inputint(matrix a,int* n){
	/*作用：按行输入整数（以n为一组整数，分母为1）至a
	  注，这个函数并不稳定*/
    int count = 0;
    for(int i = 0;i < a.row;i++){
    	for(int j = 0;j < a.col;j++){
    		matrix_set(a,i,j,fraction_create(n[count],1));
    		count++;
		}
	}
	return;
}

matrix matrix_rowcombine(matrix a,matrix b){
	//作用：按行合并矩阵a,b
	if(a.row != b.row){
		printf("错误：行数不一致");
		exit(-1);
	}
	matrix new;
	matrix_null(&new,a.row,a.col+b.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			matrix_set(new,i,j,a.base[i][j]);
		}
	}
	for(int i = 0;i < a.row;i++){
		for(int j = a.col;j < a.col+b.col;j++){
			matrix_set(new,i,j,b.base[i][j - a.col]);
		}
	}
	return new;
}

matrix matrix_colcombine(matrix a,matrix b){
	//作用：按列合并矩阵a,b
	if(a.col != b.col){
		printf("错误：列数不一致");
		exit(-1);
	}
	matrix new;
	matrix_null(&new,a.row+b.row,a.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			matrix_set(new,i,j,a.base[i][j]);
		}
	}
	for(int i = a.row;i < a.row+b.row;i++){
		for(int j = 0;j < a.col;j++){
			matrix_set(new,i,j,b.base[i - a.row][j]);
		}
	}
	return new;
}

matrix matrix_add(matrix a,matrix b){
	//作用：矩阵加法a+b
	if(a.col != b.col || a.row != b.row){
		printf("错误：维数不匹配");
		exit(-1);
	}
	matrix r;
	matrix_null(&r,a.row,a.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			r.base[i][j] = fraction_add(a.base[i][j],b.base[i][j]);
		}
	}
	return r;
}

matrix matrix_sub(matrix a,matrix b){
	//作用：矩阵减法a-b
	if(a.col != b.col || a.row != b.row){
		printf("错误：维数不匹配");
		exit(-1);
	}
	matrix r;
	matrix_null(&r,a.row,a.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			r.base[i][j] = fraction_sub(a.base[i][j],b.base[i][j]);
		}
	}
	return r;
}

matrix matrix_cmul(matrix a,Fraction C){
	//作用：标量乘法C*a
	matrix r;
	matrix_null(&r,a.row,a.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			r.base[i][j] = fraction_mul(a.base[i][j],C);
		}
	}
	return r;
}

matrix matrix_linear(Fraction C1,matrix a,Fraction C2,matrix b){
	//作用：矩阵线性组合C1*a+C2*b
	if(a.col != b.col || a.row != b.row){
		printf("错误：维数不匹配");
		exit(-1);
	}
	matrix r;
	matrix_null(&r,a.row,a.col);
	for(int i = 0;i < a.row;i++){
		for(int j = 0;j < a.col;j++){
			r.base[i][j] = fraction_add(fraction_mul(a.base[i][j],C1),fraction_mul(b.base[i][j],C2));
		}
	}
	return r;
}

matrix matrix_mul(matrix a,matrix b){
	//作用：矩阵乘法A*B
	if(a.col != b.row){
		printf("错误：维数不匹配");
		exit(-1);
	}
	matrix r;
	matrix_null(&r,a.row,b.col);
	for(int i = 0;i < r.row;i++){
		for(int j = 0;j < r.col;j++){
			for(int k = 0;k < a.col;k++)
			r.base[i][j] = fraction_add(r.base[i][j],fraction_mul(a.base[i][k],b.base[k][j]));
		}
	}
	return r;
}

matrix matrix_trans(matrix a){
	//作用：矩阵转置A^T
	matrix r;
	matrix_null(&r,a.col,a.row);
	for(int i = 0;i < r.row;i++){
		for(int j = 0;j < r.col;j++){
			r.base[i][j] = a.base[j][i];
		}
	}
	return r;
}

void matrix_mode(matrix* r,int col,mode code){
	if(code == all){
		return;
	}else if(code == left){
		r->col -= col;
		return;
	}else if(code == right){
		matrix new;
		matrix_null(&new,r->row,r->col-col);
		for(int i = 0;i < new.row;i++){
			for(int j = 0;j < new.col;j++){
				new.base[i][j] = r->base[i][j+col];
			}
		}
		*r = new;
		return;
	}else{
		printf("错误：无效模式码");
		exit(-1);
	}
}

bool matrix_REF_select(matrix a,int* posx,int* posy){
    //作用：REF算法的第一步——选择主元
	int maxx = *posx;
	for(int i = *posx;i < a.row;i++){
		if(fraction_abscmp(a.base[maxx][*posy],a.base[i][*posy]) < 0)maxx = i;
	}
    if(a.base[maxx][*posy]->numerator == 0){
        return false;
    }
	if(maxx != *posx)matrix_rreplace(a,*posx,maxx);
    return true;
}

void matrix_REF_eliminate(matrix a,int* posx,int* posy){
	//作用：REF算法中的第二步——主元列除主元以外其他变为0
	for(int i = *posx+1;i < a.row;i++){
		matrix_rcadd(a,i,fraction_div(a.base[i][*posy],fraction_opposite(a.base[*posx][*posy])),*posx);
    }
    return;
}

bool matrix_REF_poschange(matrix a,bool ifsuccess,int* posx,int* posy){
    //作用：REF算法的最后一步——移动位置
	(ifsuccess == true) ? (*posx)++,(*posy)++ : (*posy)++ ;
    if(*posx >= a.row || *posy > a.col-1)return false;
	else return true;
}

matrix matrix_REF(matrix a,matrix b,mode code){
    //作用：求增广矩阵的阶梯形矩阵
	matrix r = matrix_rowcombine(a,b);
	int posx = 0,posy = 0;
    bool flag = true;
	while(flag){
        bool issuccess = matrix_REF_select(r,&posx,&posy);
        if(issuccess)matrix_REF_eliminate(r,&posx,&posy);
        flag =  matrix_REF_poschange(r,issuccess,&posx,&posy);
    }
	matrix_mode(&r,a.col,code);
    return r;
}

void matrix_REFF_extra(matrix a,int* posx,int* posy){
	//REF算法额外的一步——当前主元变为1主元列其余元素化为0
    matrix_rmul(a,*posx,fraction_create(a.base[*posx][*posy]->denominator,a.base[*posx][*posy]->numerator));
    for(int i = (*posx)-1;i >= 0;i--){
    	matrix_rcadd(a,i,fraction_opposite(a.base[i][*posy]),*posx);
	}
}

matrix matrix_REFF(matrix a,matrix b,mode code){
    //作用：求增广矩阵的简化阶梯形
	matrix r = matrix_rowcombine(a,b);
	int posx = 0,posy = 0;
    bool flag = true;
	while(flag){
        bool issuccess = matrix_REF_select(r,&posx,&posy);
        if(issuccess){
			matrix_REF_eliminate(r,&posx,&posy);
			matrix_REFF_extra(r,&posx,&posy);
		}
        flag =  matrix_REF_poschange(r,issuccess,&posx,&posy);
    }
	matrix_mode(&r,a.col,code);
    return r;
}

Fraction matrix_det(matrix a){
	//作用：求矩阵a的行列式
	if(a.col != a.row){
		printf("错误：行列不匹配");
		exit(-1);
	}
	int count = 0;
	int posx = 0,posy = 0;
	bool flag = true;
	while(flag){
		bool issuccess;
		{
			int maxx = posx;
			for(int i = posx;i < a.row;i++){
				if(fraction_abscmp(a.base[maxx][posy],a.base[i][posy]) < 0)maxx = i;
			}
			if(a.base[maxx][posy]->numerator == 0){
				issuccess = false;
			}else{
				issuccess = true;
			}
			if(maxx != posx){
				matrix_rreplace(a,posx,maxx);
				count++;
			}
		}
		if(issuccess){
			matrix_REF_eliminate(a,&posx,&posy);
		}
		flag =  matrix_REF_poschange(a,issuccess,&posx,&posy);
	}
	Fraction r;
	r = (count%2) ? fraction_create(-1,1) : fraction_create(1,1);
	for(int i = 0;i < a.row;i++){
		r  = fraction_mul(r,a.base[i][i]);
	}
	return r;
}

matrix matrix_inverse(matrix a){
	if(a.col != a.row){
		printf("错误：行列不匹配");
		exit(-1);
	}
	if(matrix_det(a) == 0){
		printf("错误：奇异矩阵");
		exit(-1);
	}
	matrix i;
	matrix_identity(&i,a.row);
	matrix r = matrix_REFF(a,i,right);
	matrix_delete(i);
	return r;
}

matrix matrix_pow(matrix a,int k){
	//作用：矩阵乘幂A^k
	if(a.row != a.col){
		printf("错误：维数无法匹配");
		exit(-1);
	}
	matrix r;
	if(k < 0){
		r = matrix_inverse(a);
		k *= -1;
	}else{
		matrix_identity(&r,a.row);
	}
	while(k--){
		r = matrix_mul(r,a);
	}
	return r;
}

int matrix_dim(matrix a,space code){
	/*code若为col，计算矩阵a列空间维数
	  code若为row，计算矩阵a行空间维数
	  code若为nul，计算矩阵a零空间维数
	  code若为rank，计算矩阵a秩数*/
	//bug:使用dim函数后a会变成阶梯形
	int bind = 0;
	matrix r = a;

	int posx = 0,posy = 0;
    bool flag = true;
	while(flag){
        bool issuccess = matrix_REF_select(r,&posx,&posy);
        if(issuccess)matrix_REF_eliminate(r,&posx,&posy);
        (issuccess) ? posx++,posy++,bind++ : posy++;
    	if(posx >= r.row || posy > r.col-1)flag = false;
		else flag = true;
    }
	return (code == space_nul) ? r.col-bind : bind;
}
