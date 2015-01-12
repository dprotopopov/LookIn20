// Семинар 7. Решение линейных пространственно одномерных краевых задач.
//
//  (k(x)u')' - q(x) u = - f(x), xa < x < xb
//
//  u(xa) = ua, u(xb) = ub
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myprog.h"

#define IDX(i,j,n)					((n)*(j)+(i))
#define LIDX(i1,i2,nc1,nc2)			(IDX((i1),(i2),(nc1)))
#define GIDX(i1,i2,i11,i21,n1,n2)	(LIDX((i11+i1),(i21+i2),(n1+1),(n2+1)))

//static int debug=0x08+0x04+0x02+0x01;
static int debug=0;
static int np, mp, nl, ier, lp;
static int np1, np2, mp1, mp2;
static int mp_l, mp_r, mp_b, mp_t;
static char pname[MPI_MAX_PROCESSOR_NAME];
static char vname[48] = "ProNina2";
static char sname[48];
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3;

static FILE *Fi = NULL;
static FILE *Fo = NULL;

static int nx, n1, n2, n;
static double xa, xb, ua, ub, x0, a, b;

double k(double x);
double k(double x) {
	double s1 = (x-x0);
	double s2 = s1*s1;
	return a*(1.0 + s2);
}

// Первая производная
double k1(double x);
double k1(double x) {
	double s1 = (x-x0);
	return a*2.0*s1;
}

double q(double x);
double q(double x) {
	double c=pi*b;
	double s1 = (x-x0);
	double s2 = s1*s1;
	double s4 = s2*s2;
	return c*c*(1.0 + s4);
}

double u(double x);
double u(double x) {
	double c=pi*b;
	return ua*cos(c*x) + ub*sin(c*x);
}

// Первая производная
double u1(double x);
double u1(double x) {
	double c=pi*b;
	return c*(-ua*sin(c*x) + ub*cos(c*x));
}

// Вторая производная
double u2(double x);
double u2(double x) {
	double c=pi*b;
	double c2=c*c;
	return -c2*(ua*cos(c*x) + ub*sin(c*x));
}

double f(double x);
double f(double x) {
	return - k1(x)*u1(x) - k(x)*u2(x) + q(x)*u(x);
}

size_t round8bytes(size_t size);
size_t round8bytes(size_t size) {
	int r=size&0x7;
	int q=size&~0x7;
	if(r) return q+8;
	return q;
}

int main(int argc, char *argv[])
{
	int i, j, i1, i2, j1, j2;
	int i11, i12, i21, i22;
	int nc1, nc2, nc1m, nc2m, nc12, nc12m;
	int nccp, nccp1, nccp2, ncp, ncp1, ncp2, ncx, ncx1, ncx2, ncpx;
	double hx, hx2, s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
	double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al;
	double *yy1, *yy2, *yy3, *yy4;
	double *y0, *y1, *y2, *y3, *y4;
	double *aa0, *bb0, *cc0, *ff0;
	double *aa1, *bb1, *cc1, *ff1;
	int id1, id2, it;
	int ii11, ii12, ii21, ii22, nnc1, nnc2; // предыдущее вычисления
	double **yyy1; // предыдущее вычисления

	int ranks[128];
	MPI_Group gr0, gr1, gr2;
	MPI_Comm  cm0, cm1, cm2;

	MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

	fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	//sleep(1);


	sprintf(sname,"%s.p%02d",vname,mp);
	ier = fopen_m(&Fo,sname,"wt");
	if (ier!=0) mpierr("Protocol file not opened",1);

	if (mp==0) {
		sprintf(sname,"%s.d",vname);
		ier = fopen_m(&Fi,sname,"rt");
		if (ier!=0) mpierr("Data file not opened",2);
		i = fscanf(Fi,"xa=%le\n",&xa);
		i = fscanf(Fi,"xb=%le\n",&xb);
		i = fscanf(Fi,"x0=%le\n",&x0);
		i = fscanf(Fi,"a=%le\n",&a);
		i = fscanf(Fi,"b=%le\n",&b);
		i = fscanf(Fi,"ua=%le\n",&ua);
		i = fscanf(Fi,"ub=%le\n",&ub);
		i = fscanf(Fi,"n1=%d\n",&n1);
		i = fscanf(Fi,"n2=%d\n",&n2);
		i = fscanf(Fi,"n=%d\n",&n);
		i = fscanf(Fi,"lp=%d\n",&lp);
		fclose_m(&Fi);
		if (argc>1) sscanf(argv[1],"%d",&n1);
		if (argc>2) sscanf(argv[2],"%d",&n2);
		if (argc>3) sscanf(argv[3],"%d",&n);
	}

	if (np>1) {
		if (mp==0) {
			buf.ddata[0] = xa; 
			buf.ddata[1] = xb;
			buf.ddata[2] = x0;
			buf.ddata[3] = a; 
			buf.ddata[4] = b;
			buf.ddata[5] = ua; 
			buf.ddata[6] = ub;
			buf.idata[100] = n1; 
			buf.idata[101] = n2; 
			buf.idata[102] = n;
			buf.idata[103] = lp;
		}
		MPI_Bcast(buf.ddata,200,MPI_DOUBLE,0,MPI_COMM_WORLD);
		if (mp>0) {
			xa = buf.ddata[0]; 
			xb = buf.ddata[1];
			x0 = buf.ddata[2];
			a  = buf.ddata[3]; 
			b  = buf.ddata[4];
			ua = buf.ddata[5]; 
			ub = buf.ddata[6];
			n1 = buf.idata[100]; 
			n2 = buf.idata[101]; 
			n  = buf.idata[102]; // число итераций с увеличением числа ячеек
			lp = buf.idata[103];
		}
	}

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le nx=%d lp=%d\n",xa,xb,ua,ub,nx,lp);

	t1 = MPI_Wtime();

	My2DGrid(np,mp,n1,n2,&np1,&np2,&mp1,&mp2); // Распределение сетки между процессами
	//
	// mp = np1 * mp2 + mp1
	//
	if (mp1 ==     0) mp_l = -1; else mp_l = mp - 1;
	if (mp1 == np1-1) mp_r = -1; else mp_r = mp + 1;
	if (mp2 ==     0) mp_b = -1; else mp_b = mp - np1;
	if (mp2 == np2-1) mp_t = -1; else mp_t = mp + np1;

	// Создаём горизонтальтальные и вертикальные группы
	//
	// Base group:
	//
	MPI_Comm_group(MPI_COMM_WORLD,&gr0);
	cm0 = MPI_COMM_WORLD;
	//
	// Horizontal group:
	//
	for (m=0; m<np1; m++) ranks[m] = np1 * mp2 + m;
	MPI_Group_incl(gr0,np1,ranks,&gr1);
	MPI_Comm_create(MPI_COMM_WORLD,gr1,&cm1);
	//
	// Vertical group:
	//
	for (m=0; m<np2; m++) ranks[m] = np1 * m + mp1;
	MPI_Group_incl(gr0,np2,ranks,&gr2);
	MPI_Comm_create(MPI_COMM_WORLD,gr2,&cm2); 

	// Расчёты для максимального массива
	// Вычисляем максимальные размеры массивов и аллокируем их
	MyRange(np1,mp1,0,n1,&i11,&i12,&nc1); 
	MyRange(np2,mp2,0,n2,&i21,&i22,&nc2); 
	nc1m = (nc1-1)<<n; // Старший индекс в локальном массиве
	nc2m = (nc2-1)<<n; // Старший индекс в локальном массиве
	i11<<=n; // Младший индекс в глобальном массиве
	i12<<=n; // Старший индекс в глобальном массиве
	i21<<=n; // Младший индекс в глобальном массиве
	i22<<=n; // Старший индекс в глобальном массиве
	nc1=nc1m+1; // Размер локального массива
	nc2=nc2m+1; // Размер локального массива
	nc12 = nc1 * nc2; // Размер локального массива
	nc12m = nc12-1; // Старший индекс в локальном массиве

	ncp1 = 2*(np1-1); ncx1 = imax(nc1,ncp1);
	ncp2 = 2*(np2-1); ncx2 = imax(nc2,ncp2);
	ncp = imax(ncp1,ncp2); 
	ncx = imax(ncx1,ncx2);
	ncpx = imax(ncp,ncx);
	nccp1 = 2*(nc1-1); 
	nccp2 = 2*(nc2-1); 
	nccp = imax(nccp1,nccp2);

	xx = (double*)(malloc(round8bytes(sizeof(double)*nc1)));

	yyy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yyy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc12))); // предыдущее вычисления

	aa = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nccp2))));
	bb = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nccp2))));
	cc = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nccp2))));
	ff = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nccp2))));

	al = (double*)(malloc(round8bytes(sizeof(double)*imax(ncpx,nccp))));

	dd = (double*)(malloc(round8bytes(sizeof(double)*4*nccp)));
	ee = (double*)(malloc(round8bytes(sizeof(double)*4*nccp)));

	yy1 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	yy2 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	yy3 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	yy4 = (double*)(malloc(round8bytes(sizeof(double)*4*nccp)));

	al = (double*)(malloc(round8bytes(sizeof(double)*(ncpx,nccp))));
	y0 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	y2 = (double*)(malloc(round8bytes(sizeof(double)*imax(ncx,nccp))));
	y3 = (double*)(malloc(round8bytes(sizeof(double)*imax(ncx,nccp))));
	y4 = (double*)(malloc(round8bytes(sizeof(double)*9*ncp)));

	aa0 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	bb0 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	cc0 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	ff0 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));

	aa1 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	bb1 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	cc1 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));
	ff1 = (double*)(malloc(round8bytes(sizeof(double)*imax(nc1,nc2))));

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {

		MyRange(np1,mp1,0,n1,&i11,&i12,&nc1);
		MyRange(np2,mp2,0,n2,&i21,&i22,&nc2);

		nc1m = (nc1-1)<<it; // Старший индекс в локальном массиве
		nc2m = (nc2-1)<<it; // Старший индекс в локальном массиве
		i11<<=it; // Младший индекс в глобальном массиве
		i12<<=it; // Старший индекс в глобальном массиве
		i21<<=it; // Младший индекс в глобальном массиве
		i22<<=it; // Старший индекс в глобальном массиве
		nc1=nc1m+1; // Размер локального массива
		nc2=nc2m+1; // Размер локального массива
		nc12 = nc1 * nc2; // Размер локального массива
		nc12m = nc12-1; // Старший индекс в локальном массиве

		ncp1 = 2*(np1-1); ncx1 = imax(nc1,ncp1);
		ncp2 = 2*(np2-1); ncx2 = imax(nc2,ncp2);
		ncp = imax(ncp1,ncp2); ncx = imax(ncx1,ncx2);
		ncpx = imax(ncp,ncx);

		hx = (xb-xa)/(nx<<it); hx2 = hx * hx;

		// Глобальная матрица систоит из решётки локальных матриц
		// Глобальная строка образована последовательностью строк глобальной матрицы

		/////////////////////////////////////////////////////////////////////////////
		// Блок параллельных прогонок по строкам глобальной матрицы
		for(i2=0;i2<nc2;i2++){

			////////////////////////////////////////////////////////////////////////
			// Инициализация массива значений аргумента
			// сквозной нумерацией по строкам глобальной матрицы
			for (i1=0; i1<nc1; i1++) {
				xx[i1] = xa + hx * GIDX(i1,i2,i11,i21,n1,n2); // GIDX - индекс в глобальном массиве
			}
			////////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////////////////
			// Инициализация массива коэффициентов
			// см. Семинар 7. Решение линейных пространственно одномерных краевых задач.
			for (i1=0; i1<nc1; i1++) {
				s0 = k(xx[i1]); 
				s1 = k(xx[i1]-hx); 
				s2 = k(xx[i1]+hx);
				aa[i1] = 0.5 * (s0 + s1);
				bb[i1] = 0.5 * (s0 + s2);
				cc[i1] = hx2 * q(xx[i1]) + aa[i1] + bb[i1];
				ff[i1] = hx2 * f(xx[i1]);
			}

			if (GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,0,0,0,n1,n2)) { 
				// Первый элемент глобальной строки
				aa[0] = 0.0; 
				bb[0] = 0.0; 
				cc[0] = 1.0; 
				ff[0] = ua;
			}

			if (GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,n2,0,0,n1,n2)) { 
				// Последний элемент глобальной строки
				aa[nc12m] = 0.0; 
				bb[nc12m] = 0.0; 
				cc[nc12m] = 1.0; 
				ff[nc12m] = ub;
			}
			////////////////////////////////////////////////////////////////////////////


			if (lp>0) for (i1=0; i1<nc1; i1++) {
				fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
					GIDX(i1,i2,i11,i21,n1,n2),
					aa[i1],	bb[i1],	cc[i1], ff[i1]);
			}

			aa0[i2] = 0;
			bb0[i2] = 0;
			cc0[i2] = 0;
			ff0[i2] = 0;
			aa1[i2] = 0;
			bb1[i2] = 0;
			cc1[i2] = 0;
			ff1[i2] = 0;

			if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,i2,0,0,n1,n2)){ 
				// Сохраняем левый столбец глобальной матрицы
				// и присваиваем единичный коэффициент
				aa0[i2] = aa[0];
				bb0[i2] = bb[0];
				cc0[i2] = bb[0];
				ff0[i2] = ff[0];
				aa[0] = 0.0; 
				bb[0] = 0.0; 
				cc[0] = 1.0; 
				ff[0] = 0.0;
			}

			if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,i2,0,0,n1,n2)){ 
				// Сохраняем правый столбец глобальной матрицы
				// и присваиваем единичный коэффициент
				aa1[i2] = aa[nc1m];
				bb1[i2] = bb[nc1m];
				cc1[i2] = bb[nc1m];
				ff1[i2] = ff[nc1m];
				aa[nc1m] = 0.0; 
				bb[nc1m] = 0.0; 
				cc[nc1m] = 1.0; 
				ff[nc1m] = 0.0;
			}

			// Выполняем паралельную прогонку по строке глобальной матрицы
			// прогонки выполняются параллельно в строках решётки процессов
			ier = prog_rightpn(np,mp,cm1,nc1,0,
				aa, bb, cc, ff, al,
				&yy1[LIDX(0,i2,nc1,nc2)], // решение алгоритма прогонки
				y2, y3,	y4);

			if (ier!=0) mpierr("Bad solution 1",1);

			for (i1=0; i1<nc1; i1++) ff[LIDX(i1,i2,nc1,nc2)] = 0.0; 

			if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,i2,0,0,n1,n2)) ff[0] = 0.0;
			if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,i2,0,0,n1,n2)) ff[nc1m] = 1.0;

			if(GIDX(0,i2,i11,i21,n1,n2)>GIDX(n1,0,0,0,n1,n2)){ // Если не первая строка
				// Выполняем паралельную прогонку по строке глобальной матрицы
				// прогонки выполняются параллельно в строках решётки процессов
				ier = prog_rightpn(np,mp,cm1,nc1,0,
					aa, bb, cc, ff, al,
					&yy2[LIDX(0,i2,nc1,nc2)], // решение алгоритма прогонки
					y2, y3,	y4);
				if (ier!=0) mpierr("Bad solution 2",2);
			}

			if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,i2,0,0,n1,n2)) ff[0] = 1.0;
			if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,i2,0,0,n1,n2)) ff[nc1m] = 0.0;

			if(GIDX(nc1m,i2,i11,i21,n1,n2)<GIDX(0,n2,0,0,n1,n2)){ // Если не последняя строка
				// Выполняем паралельную прогонку по строке глобальной матрицы
				// прогонки выполняются параллельно в строках решётки процессов
				ier = prog_rightpn(np,mp,cm1,nc1,0,
					aa, bb, cc, ff, al,
					&yy3[LIDX(0,i2,nc1,nc2)], // решение алгоритма прогонки
					y2, y3,	y4);
				if (ier!=0) mpierr("Bad solution 3",3);
			}			
			ff[LIDX(0,i2,nc1,nc2)] = 0.0;
			ff[LIDX(nc1m,i2,nc1,nc2)] = 0.0;
		}
		// Конец блока параллельных прогонок по строкам глобальной матрицы
		/////////////////////////////////////////////////////////////////////////////

		for (i=0; i<8*nc2; i++) dd[i] = 0;
		for (i=0; i<8*nc2; i++) ee[i] = 0;

		////////////////////////////////////////////////////////////////////////
		// Подсчёт суммы сохранённых коэффициентов
		// и выполнение пересчёта в глобальной матрице
		j=0;
		for(i2=0;i2<nc2;i2++) {
			if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,0,0,0,n1,n2)) j+=4;
			else {
				cc0[i2] = cc0[i2] - bb0[i2] * yy3[LIDX(1,i2,nc1,nc2)];
				ff0[i2] = ff0[i2] + bb0[i2] * yy1[LIDX(1,i2,nc1,nc2)];
				bb0[i2] = 0.0;
				dd[j++] = aa0[i2];
				dd[j++] = bb0[i2];
				dd[j++] = cc0[i2];
				dd[j++] = ff0[i2];
			}
			if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,i2,0,0,n1,n2)) j+=4;
			else {
				cc1[i2] = cc1[i2] - aa1[i2] * yy2[LIDX(nc1m-1,i2,nc1,nc2)];
				ff1[i2] = ff1[i2] + aa1[i2] * yy1[LIDX(nc1m-1,i2,nc1,nc2)];
				aa1[i2] = 0.0;
				dd[j++] = aa1[i2];
				dd[j++] = bb1[i2];
				dd[j++] = cc1[i2];
				dd[j++] = ff1[i2];
			}
		}

		t2 = MPI_Wtime();
		MPI_Allreduce(dd,ee,8*nc2,MPI_DOUBLE,MPI_SUM,cm1);
		t2 = MPI_Wtime() - t2;

		// ee содержат одинаковые данные в процессах по горизонтали
		////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////
		// параллельная прогонка по столбцам решётки процессов

		j=0;
		for(i2=0;i2<2*nc2;i2++) {
			aa[i2] = ee[j++];   
			bb[i2] = ee[j++];
			cc[i2] = ee[j++]; 
			ff[i2] = ee[j++];
		}

		// агоритмы прогонок параллельно обрабатывают одинаковые данные

		j=0;
		i=8*nc2;
		if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,0,0,0,n1,n2)) j+=4;
		if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,n2,0,0,n1,n2)) i-=4;

		ier = prog_rightpn(np,mp,cm2,i-j,0,
			&aa[j], &bb[j], &cc[j], &ff[j], &al[j],
			&yy4[j], // решение алгоритма прогонки
			y2, y3,	y4);
		if (ier!=0) mpierr("Bad solution 4",4);

		if(GIDX(0,i2,i11,i21,n1,n2)==GIDX(0,0,0,0,n1,n2)) yy4[0]=0;
		if(GIDX(nc1m,i2,i11,i21,n1,n2)==GIDX(n1,n2,0,0,n1,n2)) yy4[2*nc2-1]=0;
		// пересчёт результатов

		for(i2=0;i2<nc2;i2++) {
			a1 = yy4[2*i2]; 
			b1 = yy4[2*i2+1];
			for (i1=0; i1<nc1; i1++)
				yy1[LIDX(i1,i2,nc1,nc2)] += a1 * yy3[LIDX(i1,i2,nc1,nc2)] + b1 * yy2[LIDX(i1,i2,nc1,nc2)];
		}
		////////////////////////////////////////////////////////////////////////

		t1 = MPI_Wtime() - t1;

		s0 = 0.0;
		for (i=0; i<nc; i++) {
			s1 = u(xx[i]); 
			s2 = dabs(s1-y1[i]); 
			s0 = dmax(s0,s2);
			if (lp>0)
				fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le d=%12le\n",
				i,xx[i],y1[i],s1,s2);
		}

		if (np>1) {
			s1 = s0; MPI_Allreduce(&s1,&s0,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		}

		if (mp==0) fprintf(stderr,"nx=%d t1=%le t2=%le dmax=%le\n",nx,t1,t2,s0);
		fprintf(Fo,"t1=%le t2=%le dmax=%le\n",t1,t2,s0);

		for(j=1;j<=it;j++){
			s0=0.0;
			for(i=0;i<=nc;i+=1<<j) 
				s0=dmax(s0,dabs(y1[i]-yyy1[it-j][i>>j]));

			MPI_Allreduce(&s0,&s1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			if (mp == 0) {
				fprintf(stderr,"Grid=%d nx=%d : nx=%d tv=%le s1=%le\n",
					np,nx<<it,(nx<<it)>>j,tv,s1);
				fflush(stderr);
			}
		}

		// Сохраняем в предыдущее значение
		for (i=0; i<nc; i++) yyy1[it][i] = y1[i];

		ii1=i1;
		ii2=i2;
		nnc=nc;
	}
	ier = fclose_m(&Fo);

	MPI_Finalize();
	return 0;
}
