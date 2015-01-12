#define _CRT_SECURE_NO_WARNINGS

//
//  du/dt = (d/dx)(k(u)du/dx - r(u)u) - u, xa < x < xb, t>0
//  du/dt = (d/dx)(k(u)du/dx) - (d/dx)(r(u)u) - u, xa < x < xb, t>0
//
//  u(x,0) = g0(x), u(xa,t) = g1(t), u(xb,t) = g2(t) 
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
static int mp_l, mp_r;
static char pname[MPI_MAX_PROCESSOR_NAME];
static char vname[48] = "Dubowik12";
static char sname[48];
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3;

static FILE *Fi = NULL;
static FILE *Fo = NULL;

static int nx, dnx, n, ntp, ntm, ntv;
static double xa, xb, r0, q0, u0, u1, a, b;
static double k1, k2, tau0, tau1, tmax, tq, epst;
static double tv, u10, omg0, omg1, gt;

double k(double u);
double k(double u) {
	return u1+a*u*u;
}

double r(double u);
double r(double u) {
	double s=dsin(u);
	return u1+b*s*s;
}

double g0(double x);
double g0(double x) {
	return u0;
}

double g1(double t);
double g1(double t) {
	double s1 = omg1 * t;
	return u1 - u10 * exp(-s1);
}

double g2(double t);
double g2(double t) {
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
	int i, j, ii, i1, i2, nc, ncm, ncp, ncx;
	double hx, hx2, tau, gam, s0, s1, s2, s3;
	double *xx, *aa, *bb, *cc, *ff, *y0, *y1, *y2, *y3, *y4, *al;
	int id1, id2, it, gcdx;
	double *bsl, *brl, *bsr, *brr; // Буферы для обмена с соседями
	double **yy1; // предыдущее вычисления
	int ii1, ii2, nnc; // предыдущее вычисления

	MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);
	if (mp ==    0) mp_l = -1; else mp_l = mp - 1;
	if (mp == np-1) mp_r = -1; else mp_r = mp + 1;

	fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	//sleep(1);

	sprintf(sname,"%s.p%02d",vname,mp);
	ier = fopen_m(&Fo,sname,"wt");
	if (ier!=0) mpierr("Protocol file not opened",1);

	if (mp==0) {
		sprintf(sname,"%s.d",vname);
		ier = fopen_m(&Fi,sname,"rt");
		if (ier!=0) mpierr("Data file not opened",2);
		fscanf(Fi,"xa=%le\n",&xa);
		fscanf(Fi,"xb=%le\n",&xb);
		fscanf(Fi,"a=%le\n",&a);
		fscanf(Fi,"b=%le\n",&b);
		fscanf(Fi,"omg0=%le\n",&omg0);
		fscanf(Fi,"omg1=%le\n",&omg1);
		fscanf(Fi,"u0=%le\n",&u0);
		fscanf(Fi,"u1=%le\n",&u1);
		fscanf(Fi,"k1=%le\n",&k1);
		fscanf(Fi,"k2=%le\n",&k2);
		fscanf(Fi,"tmax=%le\n",&tmax);
		fscanf(Fi,"tq=%le\n",&tq);
		fscanf(Fi,"epst=%le\n",&epst);
		fscanf(Fi,"nx=%d\n",&nx); //  начальное число ячеек сетки
		//fscanf(Fi,"dnx=%d\n",&dnx); // прирост числа ячеек сетки
		fscanf(Fi,"n=%d\n",&n); // число итераций с увеличением числа ячеек
		fscanf(Fi,"ntp=%d\n",&ntp); // через сколько итераций выводить информацию
		fscanf(Fi,"ntm=%d\n",&ntm); // максимальное количество итераций
		fscanf(Fi,"lp=%d\n",&lp);
		fclose_m(&Fi);
		if (argc>1) sscanf(argv[1],"%d",&nx); // число ячеек сетки
		if (argc>2) sscanf(argv[2],"%d",&dnx); // прирост числа ячеек сетки
		if (argc>3) sscanf(argv[3],"%d",&n); // число итераций с увеличением числа ячеек
		if (argc>4) sscanf(argv[4],"%d",&ntp); // через сколько итераций выводить информацию
		if (argc>5) sscanf(argv[5],"%d",&ntm); // максимальное количество итераций
	}

	if (np>1) {
		if (mp==0) {
			buf.ddata[0]  = xa;
			buf.ddata[1]  = xb;
			buf.ddata[2]  = a;
			buf.ddata[3]  = b;
			buf.ddata[4]  = omg0;
			buf.ddata[5]  = omg1;
			buf.ddata[6]  = u0;
			buf.ddata[7]  = u1;
			buf.ddata[8]  = k1;
			buf.ddata[9]  = k2;
			buf.ddata[12] = tmax;
			buf.ddata[13] = tq;
			buf.ddata[14] = epst;
			buf.idata[100] = nx; // число ячеек сетки
			//buf.idata[101] = dnx; // прирост числа ячеек сетки
			buf.idata[102] = n;
			buf.idata[103] = ntp;
			buf.idata[104] = ntm;
			buf.idata[105] = lp;
		}
		MPI_Bcast(buf.ddata,200,MPI_DOUBLE,0,MPI_COMM_WORLD);
		if (mp>0) {
			xa   = buf.ddata[0];
			xb   = buf.ddata[1];
			a   = buf.ddata[2];
			b   = buf.ddata[3];
			omg0   = buf.ddata[4];
			omg1   = buf.ddata[5];
			u0   = buf.ddata[6];
			u1   = buf.ddata[7];
			k1   = buf.ddata[8];
			k2   = buf.ddata[9];
			tmax = buf.ddata[12];
			tq   = buf.ddata[13];
			epst = buf.ddata[14];
			nx   = buf.idata[100]; // число ячеек сетки
			//dnx  = buf.idata[101]; // прирост числа ячеек сетки
			n    = buf.idata[102]; // число итераций с увеличением числа ячеек
			ntp  = buf.idata[103];
			ntm  = buf.idata[104];
			lp   = buf.idata[105];
		}
	}

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
		np,mp,pname,tick);

	fprintf(Fo,"xa=%le xb=%le a=%le b=%le omg0=%le omg1=%le\n",xa,xb,a,b,omg0,omg1);
	fprintf(Fo,"u0=%le u1=%le k1=%le k2=%le\n",u0,u1,k1,k2);
	fprintf(Fo,"tmax=%le tq=%le epst=%le\n",tmax,tq,epst);
	fprintf(Fo,"nx=%d ntp=%d ntm=%d lp=%d\n",nx,ntp,ntm,lp);

	t1 = MPI_Wtime();

	// Рассчёты для максимального массива
	MyRange(np,mp,0,nx,&i1,&i2,&nc);
	ncm = (nc-1)<<n; // Старший индекс в локальном массиве
	i1<<=n; // Младший индекс в глобальном массиве
	i2<<=n; // Старший индекс в глобальном массиве
	nc = ncm+1; // Размер локального массива
	ncp = 2*(np-1); ncx = imax(nc,ncp);

	xx = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y0 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	aa = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	bb = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	cc = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	ff = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	al = (double*)(malloc(round8bytes(sizeof(double)*ncx)));

	if (np>1) {
		y2 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
		y3 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
		y4 = (double*)(malloc(round8bytes(sizeof(double)*9*ncp)));
	}

	// Буферы для обмена с соседями
	brl = (double*)(malloc(round8bytes(sizeof(double))));
	bsl = (double*)(malloc(round8bytes(sizeof(double))));
	brr = (double*)(malloc(round8bytes(sizeof(double))));
	bsr = (double*)(malloc(round8bytes(sizeof(double))));

	yy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc))); // предыдущее вычисления

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {
		// k1 k3 q0 взяты из конспекта семинара. 
		// здесь они используются для рассчёта сходимости
		u10 = u1 - u0; /*omg0 = 1.0 / tau0; omg1 = 1.0 / tau1;*/
		hx = (xb-xa)/(nx<<it); hx2 = hx * hx;
		tau = tmax/ntm;
		tau = dmin(tau,0.5 * hx2 / dmax(k1,k2));
		tau = dmin(tau,0.5 * hx / sqrt(dmax(k1,k2)));
		s0 = dmin(tmax/tau,1000000000.0); ntm = imax(ntm,(int)s0);

		MyRange(np,mp,0,nx,&i1,&i2,&nc);
		ncm = (nc-1)<<it; // Старший индекс в локальном массиве
		i1<<=it; // Младший индекс в глобальном массиве
		i2<<=it; // Старший индекс в глобальном массиве
		nc = ncm+1; // Размер локального массива
		ncp = 2*(np-1); 
		ncx = imax(nc,ncp);

		fprintf(Fo,"u10=%le omg0=%le omg1=%le\n",u10,omg0,omg1);
		fprintf(Fo,"hx=%le tau=%le ntm=%d\n",hx,tau,ntm);
		fprintf(Fo,"nx=%d hx=%le tau=%le ntm=%d\n",nx+it*dnx,hx,tau,ntm);

		fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

		for (i=0; i<nc; i++) xx[i] = xa + hx * (i1 + i); // grid

		ntv = 0; tv = 0.0; gt = 1.0;

		for (i=0; i<nc; i++) y1[i] = g0(xx[i]); // initial profile

		gam = tau / hx2;

		// Time loop:

		do {
			ntv++; 
			tv += tau;
			if (debug&0x01) { fprintf(Fo,"tv=%le\n",tv);fflush(Fo); }

			// НАЧАЛО АЛГОРИТМА ШАГА
			if (debug&0x08) { fprintf(Fo,"Begin algo\n");fflush(Fo); }

			// Задаём начальные условия и сохраняем значение искомой функции
			for (i=0; i<nc; i++) {
				ii = i1 + i;

				if (ii==0)       y0[i] = g1(tv);
				else if (ii==nx) y0[i] = g2(tv);
				else             y0[i] = y1[i];
			}

			// Модифицированная формула из семинара 9
			// (1+tau/2)*y[j+1](i) 
			//       - tau/hh*{k(y[j](i))k(y[j](i+1))/(k(y[j](i))+k(y[j](i+1)))*(y[j+1](i+1)-y[j+1](i))
			//              -k(y[j](i))k(y[j](i-1))/(k(y[j](i))+k(y[j](i-1)))*(y[j+1](i)-y[j+1](i-1))} 
			// = (1-tau/2)*y[j](i) 
			//       + tau/hh*{k(y[j](i))k(y[j](i+1))/(k(y[j](i))+k(y[j](i+1)))*(y[j](i+1)-y[j](i))
			//              -k(y[j](i))k(y[j](i-1))/(k(y[j](i))+k(y[j](i-1)))*(y[j](i)-y[j](i-1))} 
			//       + tau*(d/dx)(r(u)u)

			// Рассчитываем коэффициенты левой/правой части
			// на основе текущих значений u
			if (debug&0x04) { fprintf(Fo,"Begin calc aa bb cc\n");fflush(Fo); }

			bsl[0] = y0[0];
			bsr[0] = y0[ncm];
			brl[0] = y0[0];
			brr[0] = y0[ncm];

			if(np>1){
				// Обмениваемся копиями крайних ячеек фрагмента матрицы с соседями
				if (debug&0x02) { fprintf(Fo,"Begin BndExch1D\n");fflush(Fo); }
				BndExch1D(np, mp, 1, 1, 1, 1,
					bsl, brl, bsr, brr);
				if (debug&0x02) { fprintf(Fo,"\t\tEnd BndExch1D\n");fflush(Fo); }
			}

			for (i=0; i<nc; i++) {
				ii = i1 + i;
				id1 = i-1;
				id2 = i+1;
				//               c[i]*y[i]-b[i]*y[i+1]=f[i], i=0
				//  -a[i]*y[i-1]+c[i]*y[i]-b[i]*y[i+1]=f[i], 0<i<n-1
				//  -a[i]*y[i-1]+c[i]*y[i]            =f[i], i=n-1
				s0 = y0[i]; 
				s1 = (id1>=0)?y0[id1]:brl[0]; 
				s2 = (id1<=ncm)?y0[id2]:brr[0];
				s0 = k(s0);
				s1 = k(s1);
				s2 = k(s2);
				aa[i] = gam * s0 * s1 / (s0 + s1);
				bb[i] = gam * s0 * s2 / (s0 + s2);
				cc[i] = 1.0 + tau/2 + aa[i] + bb[i];
			}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd calc aa bb cc\n");fflush(Fo); }

			// Вычисляем значения правой части уравнения неявной схемы
			if (debug&0x04) { fprintf(Fo,"Begin calc ff\n");fflush(Fo); }

			for (i=0; i<nc; i++) {
				ii = i1 + i;
				id1 = i-1;
				id2 = i+1;
				s0 = y0[i]; 
				s1 = (id1>=0)?y0[id1]:brl[0]; 
				s2 = (id2<nc)?y0[id2]:brr[0];
				ff[i] = (1.0 - tau/2 - aa[i] - bb[i])*s0 - aa[i]*s1 - bb[i]*s2;
			}

			// Вычисляем tau*(d/dx)(r(u)u) и прибавляем к правой части уравнения неявной схемы
			if (debug&0x04) { fprintf(Fo,"Begin add tau*(d/dx)(r(u)u)\n");fflush(Fo); }

			for (i=0; i<nc; i++) y1[i]=r(y0[i])*y0[i];

			bsl[0] = y1[0];
			bsr[0] = y1[ncm];
			brl[0] = y1[0];
			brr[0] = y1[ncm];

			if(np>1){
				// Обмениваемся копиями крайних ячеек фрагмента матрицы с соседями
				if (debug&0x02) { fprintf(Fo,"Begin BndExch1D\n");fflush(Fo); }
				BndExch1D(np, mp, 1, 1, 1, 1,
					bsl, brl, bsr, brr);
				if (debug&0x02) { fprintf(Fo,"\t\tEnd BndExch1D\n");fflush(Fo); }
			}

			for (i=0; i<nc; i++) {
				ii = i1 + i;
				id1 = i-1;
				id2 = i+1;
				s0 = y1[i]; 
				s1 = (id1>=0)?y1[id1]:brl[0]; 
				s2 = (id2<nc)?y1[id2]:brr[0];
				ff[i] += tau*((s2-s0)-(s1-s0))/2/hx;
			}

			//if(ntv==100000) {
			//	sprintf(sname,"%s_%02d_ntv.dat",vname,np);
			//	OutFun1DP(sname,np,mp,nc,xx,ff);
			//}

			if (debug&0x04) { fprintf(Fo,"\t\tEnd add tau*(d/dx)(r(u)u)\n");fflush(Fo); }

			// Задание ведущего элемента
			// Поскольку задали два ведущих элемента, то система может быть несовместна
			// но поскольку прогонка справа то всё нормально?
			if (debug&0x04) { fprintf(Fo,"Begin let main\n");fflush(Fo); }

			if(mp_l<0) { aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = g1(tv+tau); }
			if(mp_r<0) { aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = g2(tv+tau); }

			if (debug&0x04) { fprintf(Fo,"\t\tEnd let main\n");fflush(Fo); }
			if (debug&0x04) { fprintf(Fo,"\t\tEnd calc ff\n");fflush(Fo); }

			// Находим y[j+1] алгоритмом прогона
			if (debug&0x04) { fprintf(Fo,"Begin restore\n");fflush(Fo); }

			if (np<2) ier = prog_right(nc,aa,bb,cc,ff,al,y1);
			else      ier = prog_rightpm(np,mp,nc,0,aa,bb,cc,ff,al,y1,y2,y3,y4);
			if (debug&0x01) { fprintf(Fo,"\tprog_rightpm returns %d\n",ier);fflush(Fo); }

			if (debug&0x04) { fprintf(Fo,"\t\tEnd restore\n");fflush(Fo); }
			if (debug&0x08) { fprintf(Fo,"\t\tEnd algo\n");fflush(Fo); }
			// КОНЕЦ АЛГОРИТМА ШАГА

			// if (ier!=0) mpierr("Bad solution",1);

			if (ntv % ntp == 0) {
				gt = 0.0;
				for (i=0; i<nc; i++) {
					s0 = y1[i]-y0[i]; 
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
				fprintf(Fo,"ntv=%d tv=%le gt=%le\n",ntv,tv,gt);
				for (i=0; i<nc; i++)
					fprintf(Fo,"i=%8d x=%12le y1=%12le\n",(i1+i),xx[i],y1[i]);
				fflush(Fo); 
			}

		} while ((ntv<ntm) && (gt>epst));

		t1 = MPI_Wtime() - t1;

		sprintf(sname,"%s_%02d_%02d.dat",vname,np,it);
		OutFun1DP(sname,np,mp,nc,xx,y1);

		fprintf(Fo,"ntv=%d tv=%le gt=%le time=%le\n",ntv,tv,gt,t1);
		if (mp == 0) fprintf(stderr,"ntv=%d tv=%le gt=%le tcpu=%le\n",ntv,tv,gt,t1);

		for(j=1;j<=it;j++){
			s0=0.0;
			for(i=0;i<=nc;i+=1<<j) 
				s0=dmax(s0,dabs(y1[i]-yy1[it-j][i>>j]));

			MPI_Allreduce(&s0,&s1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			if (mp == 0) {
				fprintf(stderr,"Grid=%d nx=%d : nx=%d tv=%le s1=%le\n",
					np,nx<<it,(nx<<it)>>j,tv,s1);
				fflush(stderr);
			}
		}
		// Сохраняем в предыдущее значение
		for (i=0; i<nc; i++) yy1[it][i] = y1[i];

		ii1=i1;
		ii2=i2;
		nnc=nc;
	}
	ier = fclose_m(&Fo);

	MPI_Finalize();
	return 0;
}
