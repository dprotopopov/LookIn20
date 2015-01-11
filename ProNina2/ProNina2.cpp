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

static int np, mp, nl, ier, lp;
static char pname[MPI_MAX_PROCESSOR_NAME];
static char sname[48] = "ProNina2.p00";
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3;

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

int main(int argc, char *argv[])
{
	int i, j, i1, i2, nc, ncm, ncp, ncx;
	double hx, hx2, s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
	double *xx, *aa, *bb, *cc, *dd, *ee, *ff, *al, *y0, *y1, *y2, *y3, *y4;
	int id1, id2, it;
	double *bsl, *brl, *bsr, *brr; // Буферы для обмена с соседями
	double **yy1; // предыдущее вычисления
	int ii1, ii2, nnc; // предыдущее вычисления

	MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

	fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	//sleep(1);

	sprintf(sname+7,"%02d",mp);
	ier = fopen_m(&Fo,sname,"wt");
	if (ier!=0) mpierr("Protocol file not opened",1);

	if (mp==0) {
		ier = fopen_m(&Fi,"ex11a.d","rt");
		if (ier!=0) mpierr("Data file not opened",2);
		i = fscanf(Fi,"xa=%le\n",&xa);
		i = fscanf(Fi,"xb=%le\n",&xb);
		i = fscanf(Fi,"x0=%le\n",&x0);
		i = fscanf(Fi,"a=%le\n",&a);
		i = fscanf(Fi,"b=%le\n",&b);
		i = fscanf(Fi,"ua=%le\n",&ua);
		i = fscanf(Fi,"ub=%le\n",&ub);
		i = fscanf(Fi,"nx=%d\n",&nx);
		i = fscanf(Fi,"n=%d\n",&n);
		i = fscanf(Fi,"lp=%d\n",&lp);
		fclose_m(&Fi);
		if (argc>1) sscanf(argv[1],"%d",&nx);
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
			buf.idata[100] = nx; 
			buf.idata[101] = n;
			buf.idata[102] = lp;
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
			nx = buf.idata[100]; 
			n    = buf.idata[101]; // число итераций с увеличением числа ячеек
			lp = buf.idata[102];
		}
	}

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	fprintf(Fo,"xa=%le xb=%le ua=%le ub=%le nx=%d lp=%d\n",xa,xb,ua,ub,nx,lp);

	t1 = MPI_Wtime();

	// Рассчёты для максимального массива
	MyRange(np,mp,0,nx,&i1,&i2,&nc);
	ncm = (nc-1)<<n; // Старший индекс в локальном массиве
	i1<<=n; // Младший индекс в глобальном массиве
	i2<<=n; // Старший индекс в глобальном массиве
	nc = ncm+1; // Размер локального массива
	ncp = 2*(np-1); ncx = imax(nc,ncp);

	fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

	xx = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y0 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	yy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc))); // предыдущее вычисления

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

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {
		MyRange(np,mp,0,nx,&i1,&i2,&nc);
		ncm = (nc-1)<<it; // Старший индекс в локальном массиве
		i1<<=it; // Младший индекс в глобальном массиве
		i2<<=it; // Старший индекс в глобальном массиве
		nc = ncm+1; // Размер локального массива
		ncp = 2*(np-1); 
		ncx = imax(nc,ncp);

		hx = (xb-xa)/(nx<<it); hx2 = hx * hx;

		for (i=0; i<nc; i++)
			xx[i] = xa + hx * (i1 + i);

		if (mp==0) {
			aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = ua;
		}
		else {
			s0 = k(xx[0]); s1 = k(xx[0]-hx); s2 = k(xx[0]+hx);
			aa[0] = 0.5 * (s0 + s1);
			bb[0] = 0.5 * (s0 + s2);
			cc[0] = hx2 * q(xx[0]) + aa[0] + bb[0];
			ff[0] = hx2 * f(xx[0]);
		}

		for (i=1; i<ncm; i++) {
			s0 = k(xx[i]); s1 = k(xx[i-1]); s2 = k(xx[i+1]);
			aa[i] = 0.5 * (s0 + s1);
			bb[i] = 0.5 * (s0 + s2);
			cc[i] = hx2 * q(xx[i]) + aa[i] + bb[i];
			ff[i] = hx2 * f(xx[i]);
		}

		if (mp==np-1) {
			aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = ub;
		}
		else {
			s0 = k(xx[ncm]); s1 = k(xx[ncm]-hx); s2 = k(xx[ncm]+hx);
			aa[ncm] = 0.5 * (s0 + s1);
			bb[ncm] = 0.5 * (s0 + s2);
			cc[ncm] = hx2 * q(xx[ncm]) + aa[ncm] + bb[ncm];
			ff[ncm] = hx2 * f(xx[ncm]);
		}

		if (lp>0)
			for (i=0; i<nc; i++)
				fprintf(Fo,"i=%8d a=%12le b=%12le c=%12le f=%12le\n",
				i,aa[i],bb[i],cc[i],ff[i]);

		if (np<2) {
			ier = prog_right(nc,aa,bb,cc,ff,al,y1);
			if (ier!=0) mpierr("Bad solution 1",1);
			t2 = 0.0;
		}
		else {
			y2 = (double*)(malloc(sizeof(double)*nc));
			y3 = (double*)(malloc(sizeof(double)*nc));
			y4 = (double*)(malloc(sizeof(double)*ncp));
			dd = (double*)(malloc(sizeof(double)*4*ncp));
			ee = (double*)(malloc(sizeof(double)*4*ncp));

			a0 = aa[0];   b0 = bb[0];   c0 = cc[0];   f0 = ff[0];
			a1 = aa[ncm]; b1 = bb[ncm]; c1 = cc[ncm]; f1 = ff[ncm];

			if (mp==0) {
				aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y1);
				if (ier!=0) mpierr("Bad solution 1",1);

				for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y2);
				if (ier!=0) mpierr("Bad solution 2",2);

				if (lp>0)
					for (i=0; i<nc; i++)
						fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le\n",
						i,xx[i],y1[i],y2[i]);
			}
			else if (mp<np-1) {
				aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
				aa[ncm] = 0.0; bb[ncm] = 0.0; cc[ncm] = 1.0; ff[ncm] = 0.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y1);
				if (ier!=0) mpierr("Bad solution 1",1);

				for (i=0; i<ncm; i++) ff[i] = 0.0; ff[ncm] = 1.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y2);
				if (ier!=0) mpierr("Bad solution 2",2);

				ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y3);
				if (ier!=0) mpierr("Bad solution 3",3);

				if (lp>0)
					for (i=0; i<nc; i++)
						fprintf(Fo,"i=%8d x=%12le y1=%12le y2=%12le y3=%12le\n",
						i,xx[i],y1[i],y2[i],y3[i]);
			}
			else {
				aa[0] = 0.0; bb[0] = 0.0; cc[0] = 1.0; ff[0] = 0.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y1);
				if (ier!=0) mpierr("Bad solution 1",1);

				ff[0] = 1.0; for (i=1; i<=ncm; i++) ff[i] = 0.0;
				ier = prog_right(nc,aa,bb,cc,ff,al,y3);
				if (ier!=0) mpierr("Bad solution 3",3);

				if (lp>0)
					for (i=0; i<nc; i++)
						fprintf(Fo,"i=%8d x=%12le y1=%12le y3=%12le\n",
						i,xx[i],y1[i],y3[i]);
			}

			for (i=0; i<4*ncp; i++) dd[i] = 0;
			for (i=0; i<4*ncp; i++) ee[i] = 0;

			if (mp==0) {
				c1 = c1 - a1 * y2[ncm-1];
				f1 = f1 + a1 * y1[ncm-1];
				a1 = 0.0;
				dd[0] = a1;
				dd[1] = b1;
				dd[2] = c1;
				dd[3] = f1;
			}
			else if (mp<np-1) {
				c0 = c0 - b0 * y3[1];
				f0 = f0 + b0 * y1[1];
				b0 = b0 * y2[1];
				c1 = c1 - a1 * y2[ncm-1];
				f1 = f1 + a1 * y1[ncm-1];
				a1 = a1 * y3[ncm-1];
				i = mp * 8 - 4;
				dd[i]   = a0;
				dd[i+1] = b0;
				dd[i+2] = c0;
				dd[i+3] = f0;
				dd[i+4] = a1;
				dd[i+5] = b1;
				dd[i+6] = c1;
				dd[i+7] = f1;
			}
			else {
				c0 = c0 - b0 * y3[1];
				f0 = f0 + b0 * y1[1];
				b0 = 0.0;
				i = mp * 8 - 4;
				dd[i]   = a0;
				dd[i+1] = b0;
				dd[i+2] = c0;
				dd[i+3] = f0;
			}

			t2 = MPI_Wtime();

			MPI_Allreduce(dd,ee,4*ncp,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

			t2 = MPI_Wtime() - t2;

			if (lp>0)
				for (i=0; i<4*ncp; i++)
					fprintf(Fo,"i=%8d d=%12le e=%12le\n",i,dd[i],ee[i]);

			for (i=0; i<ncp; i++) {
				j = 4*i;
				aa[i] = ee[j];   bb[i] = ee[j+1];
				cc[i] = ee[j+2]; ff[i] = ee[j+3];
			}

			ier = prog_right(ncp,aa,bb,cc,ff,al,y4);
			if (ier!=0) mpierr("Bad solution 4",4);

			if (mp==0){
				b1 = y4[0];
				for (i=0; i<nc; i++)
					y1[i] = y1[i] + b1 * y2[i];
			}
			else if (mp<np-1) {
				a1 = y4[2*mp-1]; b1 = y4[2*mp];
				for (i=0; i<nc; i++)
					y1[i] = y1[i] + a1 * y3[i] + b1 * y2[i];
			}
			else {
				a1 = y4[2*mp-1];
				for (i=0; i<nc; i++)
					y1[i] = y1[i] + a1 * y3[i];
			}
		}

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
