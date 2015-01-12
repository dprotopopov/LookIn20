#define _CRT_SECURE_NO_WARNINGS

//
//  du/dt = (d/dx1)(k(x1,x2)du/dx1) + (d/dx2)(k(x1,x2)du/dx2) - u(x1,x2,t)
//
//  a1 < x1 < b1, a2 < x2 < b2, t>0
//
//  u(x1,x2,0) = g0(x1,x2)
//
//  u(x1,x2,0) = g11(x1,x2) = u0, 
//  u(b1,x2,t) = g12(t) = u1 - u10*exp(-omg1*t)
//  d/dx2 u(x1,a2,t) = g21(x1,t) = u0,
//  d/dx2 u(x1,b2,t) = g22(x1,t) = u0,
//
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <unistd.h>
#include <math.h>
#include "mycom.h"
#include "mynet.h"
#include "myio.h"
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
static char vname[48] = "LookIn20";
static char sname[48];
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3;

static FILE *Fi = NULL;
static FILE *Fo = NULL;

static int n1, n2, dn1, dn2, n, ntp, ntm, ntv;
static double a1, b1, a2, b2;
static double x10, x20, r0, q0, tau0;
static double x11, x12, x21, x22, x31, x32, k1, k2;
static double u0, u1, tau1, tmax, tq, epst;
static double tv, u10, omg0, omg1, gt;

double k(double x1, double x2);
double k(double x1, double x2) {
	double abcab[] = {x11, x12, x21, x22, x31, x32, x11, x12, x21, x22};
	double v0[] = {abcab[IDX(0,0,2)]-x1,abcab[IDX(1,0,2)]-x2};
	double v1[] = {abcab[IDX(0,1,2)]-x1,abcab[IDX(1,1,2)]-x2};
	double v2[] = {abcab[IDX(0,2,2)]-x1,abcab[IDX(1,2,2)]-x2};
	double w1[] = {abcab[IDX(0,1,2)]-abcab[IDX(0,0,2)],abcab[IDX(1,1,2)]-abcab[IDX(1,0,2)]};
	double w2[] = {abcab[IDX(0,2,2)]-abcab[IDX(0,0,2)],abcab[IDX(1,2,2)]-abcab[IDX(1,0,2)]};
	double sq1 = dabs(v0[0]*v1[1]-v0[1]*v1[0]);
	double sq2 = dabs(v1[0]*v2[1]-v1[1]*v2[0]);
	double sq3 = dabs(v0[0]*v2[1]-v0[1]*v2[0]);
	double sq = dabs(w1[0]*w2[1]-w1[1]*w2[0]);
	if(dabs(sq-sq1-sq2-sq3)<0.000001) return k2;
	return k1;
}

double f(double x1, double x2, double t);
double f(double x1, double x2, double t) {
	return u0;
}

double g0(double x1, double x2);
double g0(double x1, double x2) {
	return u0;
}

double g11(double t);
double g11(double t) {
	double s1 = omg1 * t;
	return u1 - u10 * dexp(-s1);
}

double g12(double t);
double g12(double t) {
	return u0;
}

double g21(double t);
double g21(double t) {
	return u0;
}

double g22(double t);
double g22(double t) {
	return u0;
}

size_t round8bytes(size_t size);
size_t round8bytes(size_t size){
	int r=size&0x7;
	int q=size&~0x7;
	if(r) return q+8;
	return q;
}

// наибольший общий делитель
int gcd(int n1, int n2);
int gcd(int n1, int n2){
	int x=imax(n1,n2);
	int y=imin(n1,n2);
	while(y!=0){
		int r=x%y;
		x=y;
		y=r;
	}
	return x;
}

int main(int argc, char *argv[])
{
	int m, i1, i2, j1, j2;
	int i11, i12, i21, i22;
	int nc1, nc2, nc1m, nc2m, nc12;
	int ncp, ncp1, ncp2, ncx, ncx1, ncx2, ncpx;
	double h1, h12, h2, h22, tau, gam1, gam2, s0, s1, s2, s3, s4, s5;
	double *xx1, *xx2, *aa1, *bb1, *aa2, *bb2, *yy0, *yy1, *yy2;
	double *aa, *bb, *cc, *ff, *al, *y1, *y2, *y3, *y4;
	double *ss_l, *rr_l, *ss_r, *rr_r, *ss_b, *rr_b, *ss_t, *rr_t; // Буферы для обмена с соседями
	int j, id1, id2, ncc1, ncc2, it, gcd1, gcd2;
	int ii11, ii12, ii21, ii22, nnc1, nnc2; // предыдущее вычисления
	double **yyy1; // предыдущее вычисления

	int ranks[128];
	MPI_Group gr0, gr1, gr2;
	MPI_Comm  cm0, cm1, cm2;

	MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

	//fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);fflush(stderr);
	//sleep(1);

	sprintf(sname,"%s.p%02d",vname,mp);
	ier = fopen_m(&Fo,sname,"wt");
	if (ier!=0) mpierr("Protocol file not opened",1);

	if (mp==0) {
		sprintf(sname,"%s.d",vname);
		ier = fopen_m(&Fi,sname,"rt");
		if (ier!=0) mpierr("Data file not opened",2);
		fscanf(Fi,"a1=%le\n",&a1);
		fscanf(Fi,"b1=%le\n",&b1);
		fscanf(Fi,"a2=%le\n",&a2);
		fscanf(Fi,"b2=%le\n",&b2);
		fscanf(Fi,"x10=%le\n",&x10);
		fscanf(Fi,"x20=%le\n",&x20);
		fscanf(Fi,"x11=%le\n",&x11); // координата треугольника
		fscanf(Fi,"x12=%le\n",&x12); // координата треугольника
		fscanf(Fi,"x21=%le\n",&x21); // координата треугольника
		fscanf(Fi,"x22=%le\n",&x22); // координата треугольника
		fscanf(Fi,"x31=%le\n",&x31); // координата треугольника
		fscanf(Fi,"x32=%le\n",&x32); // координата треугольника
		fscanf(Fi,"r0=%le\n",&r0);
		fscanf(Fi,"q0=%le\n",&q0);
		fscanf(Fi,"u0=%le\n",&u0); // Значение нулевого уровня
		fscanf(Fi,"u1=%le\n",&u1);
		fscanf(Fi,"k1=%le\n",&k1);
		fscanf(Fi,"k2=%le\n",&k2);
		fscanf(Fi,"omg0=%le\n",&omg0);
		fscanf(Fi,"omg1=%le\n",&omg1);
		fscanf(Fi,"tmax=%le\n",&tmax);
		fscanf(Fi,"tq=%le\n",&tq);
		fscanf(Fi,"epst=%le\n",&epst);
		fscanf(Fi,"n1=%d\n",&n1); // начальное число ячеек сетки
		fscanf(Fi,"n2=%d\n",&n2); // начальное число ячеек сетки
		fscanf(Fi,"dn1=%d\n",&dn1); // прирост числа ячеек сетки
		fscanf(Fi,"dn2=%d\n",&dn2); // прирост числа ячеек сетки
		fscanf(Fi,"n=%d\n",&n); // число итераций с увеличением числа ячеек
		fscanf(Fi,"ntp=%d\n",&ntp); // через сколько итераций выводить информацию
		fscanf(Fi,"ntm=%d\n",&ntm); // максимальное количество итераций
		fscanf(Fi,"lp=%d\n",&lp);
		fclose_m(&Fi);
		if (argc>1) sscanf(argv[1],"%d",&n1); // начальное число ячеек сетки
		if (argc>2) sscanf(argv[2],"%d",&n2); // начальное число ячеек сетки
		if (argc>3) sscanf(argv[3],"%d",&dn2); // прирост числа ячеек сетки
		if (argc>4) sscanf(argv[4],"%d",&dn2); // прирост числа ячеек сетки
		if (argc>5) sscanf(argv[5],"%d",&n); // число итераций с увеличением числа ячеек
		if (argc>6) sscanf(argv[6],"%d",&ntp); // через сколько итераций выводить информацию
		if (argc>7) sscanf(argv[7],"%d",&ntm); // максимальное количество итераций
	}

	if (np>1) { // Рассылка исходных данных по процессам из ведушего процесса
		if (mp==0) {
			buf.ddata[0]  = a1; // размеры сетки в условных единицах
			buf.ddata[1]  = b1; // размеры сетки в условных единицах
			buf.ddata[2]  = a2; // размеры сетки в условных единицах
			buf.ddata[3]  = b2; // размеры сетки в условных единицах
			buf.ddata[4]  = x10;
			buf.ddata[5]  = x20;
			buf.ddata[6]  = x11; // координата треугольника
			buf.ddata[7]  = x12; // координата треугольника
			buf.ddata[8]  = x21; // координата треугольника
			buf.ddata[9]  = x22; // координата треугольника
			buf.ddata[10]  = x31; // координата треугольника
			buf.ddata[11]  = x32; // координата треугольника
			buf.ddata[12] = r0;
			buf.ddata[13] = q0;
			buf.ddata[14] = u0; // Значение нулевого уровня
			buf.ddata[15] = u1;
			buf.ddata[16] = k1;
			buf.ddata[17] = k2;
			buf.ddata[18] = omg0;
			buf.ddata[19] = omg1;
			buf.ddata[20] = tmax;
			buf.ddata[21] = tq;
			buf.ddata[22] = epst;
			buf.idata[100] = n1; // количество ячеек сетки
			buf.idata[101] = n2; // количество ячеек сетки
			buf.idata[102] = dn1; // прирост числа ячеек сетки
			buf.idata[103] = dn2; // прирост числа ячеек сетки
			buf.idata[104] = n; // число итераций с увеличением числа ячеек
			buf.idata[105] = ntp; // через сколько итераций выводить информацию
			buf.idata[106] = ntm; // максимальное количество итераций
			buf.idata[107] = lp;
		}
		MPI_Bcast(buf.ddata,200,MPI_DOUBLE,0,MPI_COMM_WORLD); // Брэдкаст начальных параметров
		if (mp>0) {
			a1   = buf.ddata[0]; // размеры сетки в условных единицах
			b1   = buf.ddata[1]; // размеры сетки в условных единицах
			a2   = buf.ddata[2]; // размеры сетки в условных единицах
			b2   = buf.ddata[3]; // размеры сетки в условных единицах
			x10  = buf.ddata[4];
			x20  = buf.ddata[5];
			x11  = buf.ddata[6]; // координата треугольника
			x12  = buf.ddata[7]; // координата треугольника
			x21  = buf.ddata[8]; // координата треугольника
			x22  = buf.ddata[9]; // координата треугольника
			x31  = buf.ddata[10]; // координата треугольника
			x32  = buf.ddata[11]; // координата треугольника
			r0   = buf.ddata[12];
			q0   = buf.ddata[13];
			u0   = buf.ddata[14]; // Значение нулевого уровня
			u1   = buf.ddata[15];
			k1   = buf.ddata[16];
			k2   = buf.ddata[17];
			omg0 = buf.ddata[18];
			omg1 = buf.ddata[19];
			tmax = buf.ddata[20];
			tq   = buf.ddata[21];
			epst = buf.ddata[22];
			n1   = buf.idata[100]; // количество ячеек сетки
			n2   = buf.idata[101]; // количество ячеек сетки
			dn1  = buf.idata[102]; // прирост числа ячеек сетки
			dn2  = buf.idata[103]; // прирост числа ячеек сетки
			n    = buf.idata[104]; // число итераций с увеличением числа ячеек
			ntp  = buf.idata[105]; // через сколько итераций выводить информацию
			ntm  = buf.idata[106]; // максимальное количество итераций
			lp   = buf.idata[107];
		}
	}

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);fflush(Fo);

	fprintf(Fo,"a1=%le b1=%le a2=%le b2=%le\n",a1,b1,a2,b2);fflush(Fo);
	fprintf(Fo,"x10=%le x20=%le r0=%le q0=%le tau0=%le\n",x10,x20,r0,q0,tau0);fflush(Fo);
	fprintf(Fo,"x11=%le x12=%le x21=%le x22=%le x31=%le x32=%le k1=%le k2=%le\n",x11,x12,x21,x22,x31,x32,k1,k2);fflush(Fo);
	fprintf(Fo,"u0=%le u1=%le omg1=%le tmax=%le tq=%le epst=%le\n",u0,u1,omg1,tmax,tq,epst);fflush(Fo);
	fprintf(Fo,"n1=%d n2=%d ntp=%d ntm=%d lp=%d\n",n1,n2,ntp,ntm,lp);fflush(Fo);

	t1 = MPI_Wtime(); // Замер времени начала работы программы

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

	ncp1 = 2*(np1-1); ncx1 = imax(nc1,ncp1);
	ncp2 = 2*(np2-1); ncx2 = imax(nc2,ncp2);
	ncp = imax(ncp1,ncp2); ncx = imax(ncx1,ncx2);
	ncpx = imax(ncp,ncx);

	xx1 = (double*)(malloc(round8bytes(sizeof(double)*nc1)));
	xx2 = (double*)(malloc(round8bytes(sizeof(double)*nc2)));

	yy0 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	yy1 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	yy2 = (double*)(malloc(round8bytes(sizeof(double)*nc12))); // промежуточные вычисления

	yyy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yyy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc12))); // предыдущее вычисления

	aa1 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	bb1 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));

	aa2 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));
	bb2 = (double*)(malloc(round8bytes(sizeof(double)*nc12)));

	aa = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
	bb = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
	cc = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
	ff = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
	al = (double*)(malloc(round8bytes(sizeof(double)*ncpx)));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*ncx)));

	if (np>1) {
		y2 = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
		y3 = (double*)(malloc(round8bytes(sizeof(double)*ncx)));
		y4 = (double*)(malloc(round8bytes(sizeof(double)*9*ncp)));
	}

	// Буфферы для обмена с соседями
	rr_l = (double*)(malloc(round8bytes(sizeof(double)*nc2)));
	ss_l = (double*)(malloc(round8bytes(sizeof(double)*nc2)));
	rr_r = (double*)(malloc(round8bytes(sizeof(double)*nc2)));
	ss_r = (double*)(malloc(round8bytes(sizeof(double)*nc2)));
	rr_b = (double*)(malloc(round8bytes(sizeof(double)*nc1)));
	ss_b = (double*)(malloc(round8bytes(sizeof(double)*nc1)));
	rr_t = (double*)(malloc(round8bytes(sizeof(double)*nc1)));
	ss_t = (double*)(malloc(round8bytes(sizeof(double)*nc1)));

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {
		// k1 k3 q0 взяты из конспекта семинара. 
		// здесь они используются для рассчёта сходимости
		u10 = u1 - u0; /*omg0 = 1.0 / tau0; omg1 = 1.0 / tau1;*/
		h1 = (b1-a1)/(n1<<it); h12 = h1 * h1;
		h2 = (b2-a2)/(n2<<it); h22 = h2 * h2;
		tau = tmax/ntm;
		tau = dmin(tau,0.5 * dmin(h1,h2) / dmax(k1,k2)); 
		s0 = dmin(tmax/tau,1000000000.0); ntm = imax(ntm,(int)s0);

		fprintf(Fo,"u10=%le omg0=%le omg1=%le\n",u10,omg0,omg1);
		fprintf(Fo,"h1=%le h2=%le tau=%le ntm=%d\n",h1,h2,tau,ntm);

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

		ncp1 = 2*(np1-1); ncx1 = imax(nc1,ncp1);
		ncp2 = 2*(np2-1); ncx2 = imax(nc2,ncp2);
		ncp = imax(ncp1,ncp2); ncx = imax(ncx1,ncx2);
		ncpx = imax(ncp,ncx);

		fprintf(Fo,"Grid=%dx%d coord=(%d,%d)\n",np1,np2,mp1,mp2);fflush(Fo);
		fprintf(Fo,"i11=%d i12=%d nc1=%d\n",i11,i12,nc1);fflush(Fo);
		fprintf(Fo,"i21=%d i22=%d nc2=%d\n",i21,i22,nc2);fflush(Fo);
		fprintf(Fo,"ncp1=%d ncp2=%d ncp=%d\n",ncp1,ncp2,ncp);fflush(Fo);
		fprintf(Fo,"ncx1=%d ncx2=%d ncx=%d\n",ncx1,ncx2,ncx);fflush(Fo);

		fprintf(Fo,"n1=%d n2=%d h1=%le h2=%le tau=%le ntm=%d\n",
			n1<<it,n2<<it,h1,h2,tau,ntm);fflush(Fo);
		fprintf(Fo,"Grid=%dx%d\n",np1,np2);fflush(Fo);

		for (i1=0; i1<nc1; i1++) xx1[i1] = a1 + h1 * (i11 + i1); // grid for x1
		for (i2=0; i2<nc2; i2++) xx2[i2] = a2 + h2 * (i21 + i2); // grid for x2

		for (i2=0; i2<nc2; i2++){
			for (i1=0; i1<nc1; i1++) {
				m = LIDX(i1,i2,nc1,nc2);
				yy1[m] = k(xx1[i1],xx2[i2]);
			}
		}

		sprintf(sname,"%s_%02d_k.dat",vname,np);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,yy1);

		ntv = 0; tv = 0.0; gt = 1.0;

		for (i2=0; i2<nc2; i2++){
			for (i1=0; i1<nc1; i1++) {
				m = LIDX(i1,i2,nc1,nc2);
				yy1[m] = g0(xx1[i1],xx2[i2]);
			}
		}

		gam1 = tau / h12; 
		gam2 = tau / h22;

		for (i2=0; i2<nc2; i2++) {
			j2 = i21 + i2;
			for (i1=0; i1<nc1; i1++) {
				j1 = i11 + i1;
				m = LIDX(i1,i2,nc1,nc2);

				if ((j1==0) || (j1==(n1<<it))) {
					aa1[m] = 0.0; bb1[m] = 0.0;
				}
				else {
					s0 = k(xx1[i1],xx2[i2]);
					s1 = k(xx1[i1]-h1,xx2[i2]);
					s2 = k(xx1[i1]+h1,xx2[i2]);
					aa1[m] = gam1 * 2.0 * s0 * s1 / (s0 + s1);
					bb1[m] = gam1 * 2.0 * s0 * s2 / (s0 + s2);
				}

				if ((j2==0) || (j2==(n2<<it))) {
					aa2[m] = 0.0; bb2[m] = 0.0;
				}
				else {
					s0 = k(xx1[i1],xx2[i2]);
					s1 = k(xx1[i1],xx2[i2]-h2);
					s2 = k(xx1[i1],xx2[i2]+h2);
					aa2[m] = gam2 * 2.0 * s0 * s1 / (s0 + s1);
					bb2[m] = gam2 * 2.0 * s0 * s2 / (s0 + s2);
				}
			}
		}

		sprintf(sname,"%s_%02d_aa1.dat",vname,np);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,aa1);

		sprintf(sname,"%s_%02d_bb1.dat",vname,np);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,bb1);

		sprintf(sname,"%s_%02d_aa2.dat",vname,np);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,aa2);

		sprintf(sname,"%s_%02d_bb2.dat",vname,np);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,bb2);

		// Time loop:

		do {
			ntv++;

			// step 1:
			tv += tau; // Увеличение временного параметра (третье измерение)
			if (debug&0x01) { fprintf(Fo,"tv=%le\n",tv);fflush(Fo); }

			// НАЧАЛО АЛГОРИТМА ШАГА
			if (debug&0x08) { fprintf(Fo,"Begin algo\n");fflush(Fo); }

			// Задаём граничные условия по x1
			if (debug&0x04) { fprintf(Fo,"Begin let by x1\n");fflush(Fo); }

			if (mp_l<0) { // Если нет колонок левее
				for (i2=0; i2<nc2; i2++) { m = LIDX(0,i2,nc1,nc2); yy1[m]=g11(tv); }
			}
			if (mp_r<0) { // Если нет колонок правее
				for (i2=0; i2<nc2; i2++) { m = LIDX(nc1m,i2,nc1,nc2); yy1[m]=g12(tv); }
			}
			if (debug&0x04) { fprintf(Fo,"\t\tEnd let by x1\n");fflush(Fo); }

			// Берём дифференциал по оси x2
			// чтобы добавить крайние условия для производной по x2
			if (debug&0x04) { fprintf(Fo,"Begin diff by x2\n");fflush(Fo); }

			//// Вычисляем дифференциал по оси x2 следующее минус текущее
			//// Задаём граничные условия для производной по x2
			//// то есть присваиваем ноль дифференциалу по верхней границе
			//// Интегрируем обратно по оси i2 прогонкой по оси x2
			//// Это равносильно присвоению элементу верхней грацицы значения
			//// предыдущего элемента

			//// Вычисляем дифференциал по оси x2
			//// Следующее минус текущее

			//for (i1=0; i1<nc1; i1++) { m = LIDX(i1,0,nc1,nc2); ss_b[i1] = yy1[m]; }
			//for (i1=0; i1<nc1; i1++) { m = LIDX(i1,0,nc1,nc2); rr_b[i1] = yy1[m]; }
			//for (i1=0; i1<nc1; i1++) { m = LIDX(i1,nc2m,nc1,nc2); ss_t[i1] = yy1[m]; }
			//for (i1=0; i1<nc1; i1++) { m = LIDX(i1,nc2m,nc1,nc2); rr_t[i1] = yy1[m]; }

			//if(np2>1){
			//	// Обмениваемся копиями крайних строк фрагмента матрицы с соседями

			//	if (debug&0x02) { fprintf(Fo,"Begin BndAExch1D\n");fflush(Fo); }
			//	BndAExch1D(mp_b,nc1,ss_b,rr_b,
			//		mp_t,nc1,ss_t,rr_t);
			//	if (debug&0x02) { fprintf(Fo,"\t\tEnd BndAExch1D\n");fflush(Fo); }
			//}

			//for (m=0; m<nc12; m++) {
			//	i1=m%nc1;
			//	i2=m/nc1;
			//	j1 = i11 + i1;
			//	j2 = i21 + i2;
			//	id1 = nc1 * i2 + i1 - nc1; // Предыдущий элемент по столбцу
			//	id2 = nc1 * i2 + i1 + nc1; // Следующий элемент по столбцу
			//	s0=yy1[m];
			//	s1=(i2==0)?rr_b[i1]:yy1[id1];
			//	s2=(i2==nc2m)?rr_t[i1]:yy1[id2];
			//	yy2[m]=1.0*(s2-s0)+0.0*(s1-s0);
			//}

			//if (debug&0x04) { fprintf(Fo,"\t\tEnd diff by x2\n");fflush(Fo); }

			//// Задаём граничные условия для производной по x2
			//// то есть присваиваем ноль дифференциалу
			//if (debug&0x04) { fprintf(Fo,"Begin diff by x2 let zero\n");fflush(Fo); }

			//if (mp_b<0) { // Если нет строк ниже
			//	for (i1=0; i1<nc1; i1++) { m = LIDX(i1,0,nc1,nc2); yy2[m]=g21(tv)*h2; }
			//}
			//if (mp_t<0) { // Если нет строк ниже
			//	for (i1=0; i1<nc1; i1++) { m = LIDX(i1,nc2m,nc1,nc2); yy2[m]=g22(tv)*h2; }
			//}

			//if (debug&0x04) { fprintf(Fo,"\t\tEnd diff by x2 let zero\n");fflush(Fo); }

			//// Интегрируем обратно по оси i2 прогонкой по оси x2
			//// чтобы восстановить исходную функцию с которой сейчас работаем
			//if (debug&0x04) { fprintf(Fo,"Begin restore diff by x2\n");fflush(Fo); }

			//for (i1=0; i1<nc1; i1++) {
			//	j1 = i11 + i1;

			//	for (i2=0; i2<nc2; i2++) {
			//		j2 = i21 + i2;
			//		m = LIDX(i1,i2,nc1,nc2);
			//		//               c[i]*y[i]-b[i]*y[i+1]=f[i], i=0
			//		//  -a[i]*y[i-1]+c[i]*y[i]-b[i]*y[i+1]=f[i], 0<i<n-1
			//		//  -a[i]*y[i-1]+c[i]*y[i]            =f[i], i=n-1
			//		aa[i2] = 0.0; bb[i2] = 1.0; cc[i2] = 1.0; ff[i2] = yy2[m];
			//	}

			//	// Задание ведущего элемента
			//	if (debug&0x04) { fprintf(Fo,"Begin let main\n");fflush(Fo); }
			//	if(mp_b<0) { aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = yy1[0];}
			//	if(mp_t<0) { aa[nc2m] = 0.0; bb[nc2m] = 0.0; cc[nc2m] = 1.0; ff[nc2m] = yy1[nc2m];}
			//	if (debug&0x04) { fprintf(Fo,"\t\tEnd let main\n");fflush(Fo); }

			//	if (debug&0x02) { fprintf(Fo,"Begin prog_rightpn\n");fflush(Fo); }
			//	ier = prog_rightpn(np2,mp2,cm2,nc2,0,aa,bb,cc,ff,al,y1,y2,y3,y4);
			//	if (debug&0x01) { fprintf(Fo,"\tprog_rightpn returns %d\n",ier);fflush(Fo); }
			//	if (debug&0x02) { fprintf(Fo,"\t\tEnd prog_rightpn\n");fflush(Fo); }

			//	for (i2=0; i2<nc2; i2++) { m = LIDX(i1,i2,nc1,nc2); yy1[m] = y1[i2]; }
			//}

			if(ntv==1000){
				sprintf(sname,"%s_%02d_ntv.dat",vname,np);
				OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,yy1);
			}

			// Сохраняем предыдущее значение
			for (m=0; m<nc12; m++) yy0[m] = yy1[m];

			// Модифицированная формула из семинара 10

			//  (1+tau/2)*y[j+0](i1,i2) + B2(i1,i2)*(y[j+0](i1  ,i2+1)-y[j+0](i1  ,i2  )) - A2(i1,i2)*(y[j+0](i1  ,i2  )-y[j+0](i1  ,i2-1)) =
			//  = y[j+1/2](i1,i2) - B1(i1,i2)*(y[j+1/2](i1+1,i2  )-y[j+1/2](i1  ,i2  )) + A1(i1,i2)*(y[j+1/2](i1  ,i2  )-y[j+1/2](i1-1,i2  ))
			//
			//  (1+tau/2)*y[j+1/2](i1,i2) + B1(i1,i2)*(y[j+1/2](i1+1,i2  )-y[j+1/2](i1  ,i2  )) - A1(i1,i2)*(y[j+1/2](i1  ,i2  )-y[j+1/2](i1-1,i2  )) =
			//  = y[j+1](i1,i2) - B2(i1,i2)*(y[j+1](i1  ,i2+1)-y[j+1](i1  ,i2  )) + A2(i1,i2)*(y[j+1](i1  ,i2  )-y[j+1](i1  ,i2-1))

			// Начинаем вычислять вышеприведённую формулу

			// Делаем полшага по времени
			if (debug&0x04) { fprintf(Fo,"Begin first half by x2\n");fflush(Fo); }

			for (i1=0; i1<nc1; i1++) { m = LIDX(i1,0,nc1,nc2); ss_b[i1] = yy1[m]; }
			for (i1=0; i1<nc1; i1++) { m = LIDX(i1,0,nc1,nc2); rr_b[i1] = yy1[m]; }
			for (i1=0; i1<nc1; i1++) { m = LIDX(i1,nc2m,nc1,nc2); ss_t[i1] = yy1[m]; }
			for (i1=0; i1<nc1; i1++) { m = LIDX(i1,nc2m,nc1,nc2); rr_t[i1] = yy1[m]; }

			if(np2>1){
				// Обмениваемся копиями крайних строк фрагмента матрицы с соседями

				if (debug&0x02) { fprintf(Fo,"Begin BndAExch1D\n");fflush(Fo); }
				BndAExch1D(mp_b,nc1,ss_b,rr_b,
					mp_t,nc1,ss_t,rr_t);
				if (debug&0x02) { fprintf(Fo,"\t\tEnd BndAExch1D\n");fflush(Fo); }
			}

			// Вычисляем полшага
			// по столбцам x2
			//  (1+tau/2)*y[j+0](i1,i2) + B2(i1,i2)*(y[j+0](i1  ,i2+1)-y[j+0](i1  ,i2  )) - A2(i1,i2)*(y[j+0](i1  ,i2  )-y[j+0](i1  ,i2-1)) =

			for (m=0; m<nc12; m++) {
				i1=m%nc1;
				i2=m/nc1;
				j1 = i11 + i1;
				j2 = i21 + i2;
				id1=nc1*i2+i1-nc1; // Предыдущий элемент по столбцу
				id2=nc1*i2+i1+nc1; // Следующий элемент по столбцу
				s0=yy1[m];
				s1=(i2==0)?rr_b[i1]:yy1[id1];
				s2=(i2==nc2m)?rr_t[i1]:yy1[id2]; 

				yy2[m] = (1.0+tau/2)*s0+bb2[m]*(s2-s0)+aa2[m]*(s1-s0);
			}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd first half by x2\n");fflush(Fo); }

			// Восстанавливаем полшага
			// применяем алгоритм прогона
			// по строкам x1
			//  = y[j+1/2](i1,i2) - B1(i1,i2)*(y[j+1/2](i1+1,i2  )-y[j+1/2](i1  ,i2  )) + A1(i1,i2)*(y[j+1/2](i1  ,i2  )-y[j+1/2](i1-1,i2  ))

			if (debug&0x04) { fprintf(Fo,"Begin restore first half by x1\n");fflush(Fo); }

			for (i2=0; i2<nc2; i2++) {
				j2 = i21 + i2;
				for (i1=0; i1<nc1; i1++) {
					j1 = i11 + i1;

					m = LIDX(i1,i2,nc1,nc2);
					//               c[i]*y[i]-b[i]*y[i+1]=f[i], i=0
					//  -a[i]*y[i-1]+c[i]*y[i]-b[i]*y[i+1]=f[i], 0<i<n-1
					//  -a[i]*y[i-1]+c[i]*y[i]            =f[i], i=n-1
					aa[i1] = aa1[m]; bb[i1] = bb1[m]; cc[i1] = 1.0+aa1[m]+bb1[m]; ff[i1] = yy2[m];
				}

				if (debug&0x02) { fprintf(Fo,"Begin prog_rightpn\n");fflush(Fo); }
				ier = prog_rightpn(np1,mp1,cm1,nc1,0,aa,bb,cc,ff,al,y1,y2,y3,y4);
				if (debug&0x01) { fprintf(Fo,"\tprog_rightpn returns %d\n",ier);fflush(Fo); }
				if (debug&0x02) { fprintf(Fo,"\t\tEnd prog_rightpn\n");fflush(Fo); }

				for (i1=0; i1<nc1; i1++) { m = LIDX(i1,i2,nc1,nc2); yy1[m] = y1[i2]; }
			}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd restore first half by x1\n");fflush(Fo); }

			// Делаем вторые полшага по времени

			if (debug&0x04) { fprintf(Fo,"Begin second half by x1\n");fflush(Fo); }

			for (i2=0; i2<nc2; i2++) { m = LIDX(0,i2,nc1,nc2); ss_l[i2] = yy1[m]; }
			for (i2=0; i2<nc2; i2++) { m = LIDX(0,i2,nc1,nc2); rr_l[i2] = yy1[m]; }
			for (i2=0; i2<nc2; i2++) { m = LIDX(nc1m,i2,nc1,nc2); ss_r[i2] = yy1[m]; }
			for (i2=0; i2<nc2; i2++) { m = LIDX(nc1m,i2,nc1,nc2); rr_r[i2] = yy1[m]; }

			if(np1>1){
				// Обмениваемся копиями крайних столбцов фрагмента матрицы с соседями

				if (debug&0x02) { fprintf(Fo,"Begin BndAExch1D\n");fflush(Fo); }
				BndAExch1D(mp_l,nc2,ss_l,rr_l,
					mp_r,nc2,ss_r,rr_r);
				if (debug&0x02) { fprintf(Fo,"\t\tEnd BndAExch1D\n");fflush(Fo); }
			}

			// Вычисляем полшага
			// по строкам x1
			//  (1+tau/2)*y[j+1/2](i1,i2) + B1(i1,i2)*(y[j+1/2](i1+1,i2  )-y[j+1/2](i1  ,i2  )) - A1(i1,i2)*(y[j+1/2](i1  ,i2  )-y[j+1/2](i1-1,i2  )) =

			for (m=0; m<nc12; m++) {
				i1=m%nc1;
				i2=m/nc1;
				j1 = i11 + i1;
				j2 = i21 + i2;
				id1=nc1*i2+i1-1; // Предыдущий элемент по строке
				id2=nc1*i2+i1+1; // Следующий элемент по строке
				s0=yy1[m];
				s1=((i1==0)?rr_l[i2]:yy1[id1]);
				s2=((i1==nc1m)?rr_r[i2]:yy1[id2]); 

				yy2[m] = (1.0+tau/2)*s0+bb1[m]*(s2-s0)+aa1[m]*(s1-s0);
			}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd second half by x1\n");fflush(Fo); }

			// Восстанавливаем вторые полшага
			// применяем алгоритм прогона
			// по столбцам x2
			//  = y[j+1](i1,i2) - B2(i1,i2)*(y[j+1](i1  ,i2+1)-y[j+1](i1  ,i2  )) + A2(i1,i2)*(y[j+1](i1  ,i2  )-y[j+1](i1  ,i2-1))

			if (debug&0x04) { fprintf(Fo,"Begin restore second half by x2\n");fflush(Fo); }

			for (i1=0; i1<nc1; i1++) {
				j1 = i11 + i1;
				for (i2=0; i2<nc2; i2++) {
					j2 = i21 + i2;

					m = LIDX(i1,i2,nc1,nc2);
					//               c[i]*y[i]-b[i]*y[i+1]=f[i], i=0
					//  -a[i]*y[i-1]+c[i]*y[i]-b[i]*y[i+1]=f[i], 0<i<n-1
					//  -a[i]*y[i-1]+c[i]*y[i]            =f[i], i=n-1
					aa[i2] = aa2[m]; bb[i2] = bb2[m]; cc[i2] = 1.0+aa2[m]+bb2[m]; ff[i2] = yy2[m];
				}

				if (debug&0x02) { fprintf(Fo,"Begin prog_rightpn\n");fflush(Fo); }
				ier = prog_rightpn(np2,mp2,cm2,nc2,0,aa,bb,cc,ff,al,y1,y2,y3,y4);
				if (debug&0x01) { fprintf(Fo,"\tprog_rightpn returns %d\n",ier);fflush(Fo); }
				if (debug&0x02) { fprintf(Fo,"\t\tEnd prog_rightpn\n");fflush(Fo); }

				for (i2=0; i2<nc2; i2++) { m = LIDX(i1,i2,nc1,nc2); yy1[m] = y1[i2]; }
			}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd restore second half by x2\n");fflush(Fo); }
			if (debug&0x08) { fprintf(Fo,"\t\tEnd algo\n");fflush(Fo); }
			// КОНЕЦ АЛГОРИТМА ШАГА

			// step 2:
			// вывод различной информации

			if (ntv % ntp == 0) {
				gt = 0.0;
				for (m=0; m<nc12; m++) {
					s0 = yy1[m]-yy0[m]; 
					gt = dmax(gt,dabs(s0));
				}
				gt = gt / tau; 

				if (np>1) {
					s0 = gt; 
					MPI_Allreduce(&s0,&gt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
				}

				if (mp == 0) {
					t2 = MPI_Wtime() - t1;
					fprintf(stderr,"ntv=%d tv=%le tau=%le gt=%le tcpu=%le\n",ntv,tv,tau,gt,t2);
				}
			}

			if (lp>0) {
				fprintf(Fo,"ntv=%d tv=%le gt=%le\n",ntv,tv,gt);fflush(Fo);

				for (i2=0; i2<nc2; i2++) {
					j2 = i21 + i2;
					for (i1=0; i1<nc1; i1++) {
						j1 = i11 + i1;
						m = LIDX(i1,i2,nc1,nc2);
						fprintf(Fo,"i1=%8d i2=%8d x1=%12le x2=%12le y1=%12le\n",
							j1,j2,xx1[i1],xx2[i2],yy1[m]);fflush(Fo);
					}
				}
			}
		} while ((ntv<ntm) && (gt>epst));

		t1 = MPI_Wtime() - t1; // Замер времени от начала работы программы

		sprintf(sname,"%s_%02d_%02d.dat",vname,np,it);
		OutFun2DP(sname,np,mp,nc1,nc2,xx1,xx2,yy1);

		fprintf(Fo,"ntv=%d tv=%le gt=%le time=%le\n",ntv,tv,gt,t1);fflush(Fo);

		if (mp == 0) {
			fprintf(stderr,"Grid=%dx%d n1=%d n2=%d ntv=%d tv=%le gt=%le tcpu=%le\n",
				np1,np2,n1<<it,n2+it*dn2,ntv,tv,gt,t1);fflush(stderr);
		}

		for(j=1;j<=it;j++){
			s0=0.0;
			for(i1=0;i1<=nc1;i1+=1<<j) 
				for(i2=0;i2<=nc2;i2+=1<<j) 
					s0=dmax(s0,dabs(yy1[nc1*i2+i1]-yyy1[it-j][(nnc1*i2+i1)>>j]));

			MPI_Allreduce(&s0,&s1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			if (mp == 0) {
				fprintf(stderr,"Grid=%dx%d n1=%d n2=%d : n1=%d n2=%d tv=%le s1=%le\n",
					np1,np2,n1<<it,n2<<it,(n1<<it)>>j,(n2<<it)>>j,tv,s1);
				fflush(stderr);
			}
		}
		// Сохраняем в предыдущее значение
		for (m=0; m<nc12; m++) yyy1[it][m] = yy1[m];

		ii11=i11;
		ii12=i12;
		ii21=i21; 
		ii22=i22;
		nnc1=nc1;
		nnc2=nc2;
	}

	ier = fclose_m(&Fo);

	MPI_Finalize();
	return 0;
}
