//  Boundary problem:
//
//  (k(x)u')' - q(x) u = - f(x), xa < x < xb
//
//  u'(xa) = ua, u'(xb) = ub
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

static int np, mp, nl, ier, lp;
static char pname[MPI_MAX_PROCESSOR_NAME];
static char vname[48] = "ProNina2";
static char sname[48];
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3, gt, epst;
static double m11, m12, m21, m22;
static int n, ntp, ntm, ntv;

static FILE *Fi = NULL;
static FILE *Fo = NULL;

static int nx, n;
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

// Асинхронный циклический обмен данными между процессами
void CircleABndExch1D(int np, int mp, 
					  int nsl, int nrl, int nsr, int nrr,
					  double *bsl, double *brl, double *bsr, double *brr);
void CircleABndExch1D(int np, int mp, 
					  int nsl, int nrl, int nsr, int nrr,
					  double *bsl, double *brl, double *bsr, double *brr)
{
	int i, n, m;
	static MPI_Status Sta[4];
	static MPI_Request Req[4];
	MPI_Request *R;

	R = Req; m = 0;

	MPI_Irecv(brl,nrl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++;
	MPI_Irecv(brr,nrr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++;
	MPI_Isend(bsl,nsl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++;
	MPI_Isend(bsr,nsr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++;

	if (m>0) {
		MPI_Waitall(m,Req,Sta);

		n = 0;

		for (i=0; i<m; i++)
			if (Sta[i].MPI_ERROR != 0) n++;

		if (n>0) mpierr("Bad asynchronous exchange",-101);
	}
}

int main(int argc, char *argv[])
{
	int m, i, j, i1, i2, nc, ncm, ncp, ncx;
	double hx, hx2, s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
	double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al, *y0, *y1, *y2, *y3, *y4;
	int id1, id2, it;
	double *bsl, *brl, *bsr, *brr; // буферы для обмена данными между процессами
	int ii1, ii2, nnx; // предыдущее вычисления
	double **yy1; // предыдущее вычисления

	MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

	fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	//sleep(1);

	MyRange(np,mp,0,nx,&i1,&i2,&nc);
	ncm = (nc-1)<<n; 
	nc = ncm+1;
	ncp = 2*(np-1); 
	ncx = imax(nc,ncp);

	xx = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	aa = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	bb = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	cc = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	ff = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	al = (double*)(malloc(round8bytes(sizeof(double)*ncx)));

	y0 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y2 = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	if (np>1) {
		y3 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
		y4 = (double*)(malloc(round8bytes(sizeof(double)*9*ncp)));
	}

	// Буферы для обмена с соседями
	brl = (double*)(malloc(round8bytes(sizeof(double)*2)));
	bsl = (double*)(malloc(round8bytes(sizeof(double)*2)));
	brr = (double*)(malloc(round8bytes(sizeof(double)*2)));
	bsr = (double*)(malloc(round8bytes(sizeof(double)*2)));


	yy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc))); // предыдущее вычисления

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
		//i = fscanf(Fi,"ua=%le\n",&ua);
		//i = fscanf(Fi,"ub=%le\n",&ub);
		i = fscanf(Fi,"m11=%le\n",&m11);
		i = fscanf(Fi,"m12=%le\n",&m12);
		i = fscanf(Fi,"m21=%le\n",&m21);
		i = fscanf(Fi,"m22=%le\n",&m22);
		i = fscanf(Fi,"nx=%d\n",&nx);
		i = fscanf(Fi,"epst=%le\n",&epst);
		fscanf(Fi,"ntm=%d\n",&ntm); // максимальное количество итераций
		i = fscanf(Fi,"n=%d\n",&n);
		i = fscanf(Fi,"lp=%d\n",&lp);
		fclose_m(&Fi);
		if (argc>1) sscanf(argv[1],"%d",&nx);
		if (argc>2) sscanf(argv[2],"%d",&n);
	}

	if (np>1) {
		if (mp==0) {
			buf.ddata[0] = xa; 
			buf.ddata[1] = xb;
			buf.ddata[2] = x0;
			buf.ddata[3] = a; 
			buf.ddata[4] = b;
			//buf.ddata[5] = ua; 
			//buf.ddata[6] = ub;
			buf.ddata[7] = m11;
			buf.ddata[8] = m12;
			buf.ddata[9] = m21;
			buf.ddata[10] = m22;
			buf.ddata[14] = epst;
			buf.idata[101] = nx; 
			buf.idata[102] = n;
			buf.idata[103] = ntm;
			buf.idata[104] = lp;
		}
		MPI_Bcast(buf.ddata,200,MPI_DOUBLE,0,MPI_COMM_WORLD);
		if (mp>0) {
			xa = buf.ddata[0]; 
			xb = buf.ddata[1];
			x0 = buf.ddata[2];
			a  = buf.ddata[3]; 
			b  = buf.ddata[4];
			//ua = buf.ddata[5]; 
			//ub = buf.ddata[6];
			m11 = buf.ddata[7];
			m12 = buf.ddata[8];
			m21 = buf.ddata[9];
			m22 = buf.ddata[10];
			epst = buf.ddata[14];
			nx = buf.idata[101]; 
			n  = buf.idata[102]; // число итераций с увеличением числа ячеек
			ntm  = buf.idata[103];
			lp = buf.idata[104];
		}
	}

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le nx=%d lp=%d\n",xa,xb,ua,ub,nx,lp);

	t1 = MPI_Wtime();

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {
		// Канонический вид задачи
		// aa[i]y[i-1]-cc[i]y[i]+bb[i]y[i=1] = -ff[i]


		// Дополнительные условия

		// режим 1
		// u'(xa) = pi*b*u(xa)
		// u'(xb) = -pi*b*u(xb)

		// режим 2
		// u(xa) = u(xb)
		// u'(xa) = u'(xb)

		// (u,u')(xa) LR (u,u')(xb)

		// необходтмое условие сходимости
		// |LR|<=1

		MyRange(np,mp,0,nx,&i1,&i2,&nc);
		ncm = (nc-1)<<it; // Старший индекс в локальном массиве
		i1<<=it; // Младший индекс в глобальном массиве
		i2<<=it; // Старший индекс в глобальном массиве
		nc = ncm+1; // Размер локального массива
		ncp = 2*(np-1); 
		ncx = imax(nc,ncp);

		hx = (xb-xa)/(nx<<it); 
		hx2 = hx * hx;
		ua = u(xa); 
		ub = u(xb);

		fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

		for (i=0; i<nc; i++)
			xx[i] = xa + hx * (i1 + i);

		if (mp==0) {
			s0 = k(xx[0]); s2 = k(xx[1]);
			aa[0] = 0.0;
			bb[0] = 0.5 * (s0 + s2);
			cc[0] = 0.5 * hx2 * q(xx[0]) + bb[0];
		}
		else {
			s0 = k(xx[0]); s1 = k(xx[0]-hx); s2 = k(xx[0]+hx);
			aa[0] = 0.5 * (s0 + s1);
			bb[0] = 0.5 * (s0 + s2);
			cc[0] = hx2 * q(xx[0]) + aa[0] + bb[0];
		}

		for (i=1; i<ncm; i++) {
			s0 = k(xx[i]); s1 = k(xx[i-1]); s2 = k(xx[i+1]);
			aa[i] = 0.5 * (s0 + s1);
			bb[i] = 0.5 * (s0 + s2);
			cc[i] = hx2 * q(xx[i]) + aa[i] + bb[i];
		}

		if (mp==np-1) {
			s0 = k(xx[ncm]); s1 = k(xx[ncm-1]);
			aa[ncm] = 0.5 * (s0 + s1);
			bb[ncm] = 0.0;
			cc[ncm] = 0.5 * hx2 * q(xx[ncm]) + aa[ncm];
		}
		else {
			s0 = k(xx[ncm]); s1 = k(xx[ncm]-hx); s2 = k(xx[ncm]+hx);
			aa[ncm] = 0.5 * (s0 + s1);
			bb[ncm] = 0.5 * (s0 + s2);
			cc[ncm] = hx2 * q(xx[ncm]) + aa[ncm] + bb[ncm];
		}

		// (u[i],u'[i+1/2],u[i+1],u'[i+3/2],..
		// h*u'[i+1/2] = -u[i] + u[i+1]
		// u[i+1] = u[i]+u'[i+1/2]*h 

		//	u[i-1]	hu'[i-1/2]	u[i]	hu'[i+1/2]	u[i+1]
		//	-1		1			1
		//	                    1		-1			-1

		ntv = 0; gt = 1.0;

		// Задаём начальное значение функции
		for(i=0;i<ncm;i++) y1[i] = 0; 

		do{
			ntv++; 

			if(mp==0) y0[0] = ua;
			if(mp==np-1) y0[ncm] = ub;

			for(m=0;m<2*np;m++){
				// Повторяем цикл несколько раз, чтобы изменения дошли от левого края до правого и наоборот

				for(i=0;i<ncm;i++){
					y1[i  ] = y0[i+1] - y0[i  ]; // Прямые циклические вычисления
				}

				bsl[0] = y0[0  ];	bsl[1] = y1[0  ];
				bsr[0] = y0[ncm];	bsr[1] = y1[ncm];

				if(np>1) {
					// Циклический обмен данными между процессами
					CircleABndExch1D(np, mp, 2, 2, 2, 2, bsl, brl, bsr, brr);
				}
				else {
					// Если процесс один, то просто копируем содержимое буферов
					brr[0] = bsl[0];	brr[1] = bsl[1];
					brl[0] = bsr[0];	brl[1] = bsr[1];
				}

				if(mp==0) {
					// Для первого процесса накладываем дополнительные условия 
					//на первый отсчёт
					s1 = brl[0];
					s2 = brl[1];
					y0[0] = m11*s1 + pi*m12*s2*hx;
				}
				if(mp==np-1) {
					// Для последнего процесса накладываем дополнительные условия 
					// на последний отсчёт
					s1 = brr[0];
					s2 = brr[1];
					y1[ncm] = m21*s1/pi/hx + m22*s2;
				}

				if(0<mp && mp<np-1) {
					// Дообрабатываем то что не могли обработать из-за отсутствия данных
					// и насинаем встречные вычмсления
					y1[ncm] = brr[0] - y0[ncm];
					y0[ncm] = brr[0] - y1[ncm];
				}

				for(j=ncm;j>0;i++,j--){
					y0[j-1] = y0[j  ] - y1[j-1];	 // Встречные циклические вычисления		
				}
			}

			for(i=0;i<ncx;i++){
				s0 = y0[i];
				s1 = (i==0)?brl[0]:y0[i-1];
				s2 = (i==ncm)?brr[0]:y0[i+1];
				y2[i]=(s2-s0)+(s1-s0);
			}

			for(i=0;i<nc;i++){
				// Вычисляем значении для классической сжемы решения
				ff[i]=q(xx[i])*y0[i]*hx2-k1(xx[i])*y1[i]*hx-k(xx[i])*y2[i];			
			}

			// Накладываем начальные ограничения
			if(mp==0) ff[0] = ua;
			if(mp==np-1) ff[ncm] = ub;

			// Находим искомую функцию методом прогонки
			if (np<2) ier = prog_right(nc,aa,bb,cc,ff,al,y1);
			else      ier = prog_rightpm(np,mp,nc,0,aa,bb,cc,ff,al,y1,y2,y3,y4);

			t1 = MPI_Wtime() - t1;

			// Вычисляем отличие от реальной функции
			gt = 0.0;
			for (i=0; i<nc; i++) {
				s1 = u(xx[i]); s2 = dabs(s1-y1[i]); gt = dmax(gt,s2);
				if (lp>0)
					fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le gt=%12le\n",
					i,xx[i],y1[i],s1,gt);
			}

			if (np>1) {
				s1 = gt; 
				MPI_Allreduce(&s1,&gt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			}

			// повторяем если большое расхождение с текущим значением в качестве начального
		} while ((ntv<ntm) && (gt>epst));

		sprintf(sname,"%s_%02d_%02d.dat",vname,np,it);
		OutFun1DP(sname,np,mp,nc,xx,y1);

		if (mp==0) fprintf(stderr,"nx=%d t1=%le t2=%le dmax=%le\n",nx,t1,t2,s0);
		fprintf(Fo,"t1=%le t2=%le dmax=%le\n",t1,t2,s0);

		// Сохраняем в предыдущее значение
		for (i=0; i<nc; i++) yy1[it][i] = y1[i];

		ii1=i1;
		ii2=i2;
		nnx=nx;
	}
	ier = fclose_m(&Fo);

	MPI_Finalize();
	return 0;
}
