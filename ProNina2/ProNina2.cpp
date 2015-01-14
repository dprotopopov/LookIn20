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

#define LIDX(i,nc)		(i)
#define GIDX(i,i1,nx)	(LIDX((i1+i),(n1+1)))

static int np, mp, nl, ier, lp;
static char pname[MPI_MAX_PROCESSOR_NAME];
static char vname[48] = "ProNina2";
static char sname[48];
static MPI_Status status;
static union_t buf;
static double tick, t1, t2, t3, gt, epst;
//static double m11, m12, m21, m22;
static int n, ntp, ntm, ntv, mode;

static FILE *Fi = NULL;
static FILE *Fo = NULL;

static int nx, n;
static double xa, xb, ua, ub, x0, a, b, r;

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
	double c=pi*b/2;
	double s1 = (x-x0);
	double s2 = s1*s1;
	double s4 = s2*s2;
	return c*c*(1.0 + s4);
}

double u(double x);
double u(double x) {
	double c=pi*b/2;
	return dcos(c*x) + dsin(c*x);
}

// Первая производная
double u1(double x);
double u1(double x) {
	double c=pi*b/2;
	return c*(-dsin(c*x) + dcos(c*x));
}

// Вторая производная
double u2(double x);
double u2(double x) {
	double c=pi*b/2;
	double c2=c*c;
	return -c2*(dcos(c*x) + dsin(c*x));
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
	if(np<2) {
		// Если процесс один, то просто копируем содержимое буферов
		for(i=0;i<imin(nsl,nrr);i++) brr[i] = bsl[i];
		for(i=0;i<imin(nsr,nrl);i++) brl[i] = bsr[i];
	}
	else if(mp==0){
		if(nsr>0) { MPI_Isend(bsr,nsr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nrr>0) { MPI_Irecv(brr,nrr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nsl>0) { MPI_Isend(bsl,nsl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nrl>0) { MPI_Irecv(brl,nrl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
	} else {
		if(nrl>0) { MPI_Irecv(brl,nrl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nsl>0) { MPI_Isend(bsl,nsl,MPI_DOUBLE,(mp+np-1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nrr>0) { MPI_Irecv(brr,nrr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
		if(nsr>0) { MPI_Isend(bsr,nsr,MPI_DOUBLE,(mp+np+1)%np,MY_TAG,MPI_COMM_WORLD,R); R++; m++; }
	}

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
	int m, p, i, j, i1, i2, nc, ncm, ncp, ncx;
	double hx, hx2, s0, s1, s2, a0, b0, c0, f0, a1, b1, c1, f1;
	double *xx, *aa, *bb, *cc, *kk, *kk1, *qq, *ff, *uu, *y0, *y1, *y2, *_y2, *_y3, *_y4, *al;
	double *bsl, *brl, *bsr, *brr;
	int id1, id2, it;
	int ii1, ii2, nnx; // предыдущее вычисления
	double **yy1; // предыдущее вычисления
	double tau, gam;

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
		//i = fscanf(Fi,"ua=%le\n",&ua);
		//i = fscanf(Fi,"ub=%le\n",&ub);
		//i = fscanf(Fi,"m11=%le\n",&m11);
		//i = fscanf(Fi,"m12=%le\n",&m12);
		//i = fscanf(Fi,"m21=%le\n",&m21);
		//i = fscanf(Fi,"m22=%le\n",&m22);
		i = fscanf(Fi,"mode=%d\n",&mode);
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
			//buf.ddata[7] = m11;
			//buf.ddata[8] = m12;
			//buf.ddata[9] = m21;
			//buf.ddata[10] = m22;
			buf.ddata[14] = epst;
			buf.idata[100] = mode; 
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
			//m11 = buf.ddata[7];
			//m12 = buf.ddata[8];
			//m21 = buf.ddata[9];
			//m22 = buf.ddata[10];
			epst = buf.ddata[14];
			mode = buf.idata[100]; 
			nx = buf.idata[101]; 
			n  = buf.idata[102]; // число итераций с увеличением числа ячеек
			ntm  = buf.idata[103];
			lp = buf.idata[104];
		}
	}

	ua = u(xa); // Значения в крайних точках вычисляем чтобы не передавать в параметрах
	ub = u(xb); // Значения в крайних точках вычисляем чтобы не передавать в параметрах

	fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
	fprintf(Fo,"xa=%le xb=%le x0=%le ua=%le ub=%le a=%le b=%le nx=%d lp=%d mode=%d\n",
		xa,xb,x0,ua,ub,a,b,nx,lp,mode);
	//fprintf(Fo,"m11=%le m12=%le m21=%le m22=%le\n",	m11,m12,m21,m22);

	t1 = MPI_Wtime();

	MyRange(np,mp,0,nx,&i1,&i2,&nc);
	ncm = (nc-1)<<n; 
	nc = ncm+1;
	ncp = 2*(np-1); 
	ncx = imax(nc,ncp);

	bsl = (double*)(malloc(round8bytes(sizeof(double)*2)));
	brl = (double*)(malloc(round8bytes(sizeof(double)*2)));
	bsr = (double*)(malloc(round8bytes(sizeof(double)*2)));
	brr = (double*)(malloc(round8bytes(sizeof(double)*2)));

	xx = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	aa = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	bb = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	cc = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	ff = (double*)(malloc(round8bytes(sizeof(double)*(nc+1))));
	kk = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	kk1= (double*)(malloc(round8bytes(sizeof(double)*nc)));
	qq = (double*)(malloc(round8bytes(sizeof(double)*nc)));

	uu = (double*)(malloc(round8bytes(sizeof(double)*nc)));
	y0 = (double*)(malloc(round8bytes(sizeof(double)*(nc+1))));
	y1 = (double*)(malloc(round8bytes(sizeof(double)*(nc+1))));
	y2 = (double*)(malloc(round8bytes(sizeof(double)*(nc+1))));

	al = (double*)(malloc(round8bytes(sizeof(double)*ncx)));

	if (np>1) {
		_y2 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
		_y3 = (double*)(malloc(round8bytes(sizeof(double)*nc)));
		_y4 = (double*)(malloc(round8bytes(sizeof(double)*9*ncp)));
	}

	yy1 = (double**)(malloc(round8bytes(sizeof(double*)*(n+1)))); // предыдущее вычисления
	for(j=0;j<=n;j++) yy1[j] = (double*)(malloc(round8bytes(sizeof(double)*nc))); // предыдущее вычисления

	// Цикл с разными шагами сетки
	for(it=0;it<=n;it++) {
		// Канонический вид задачи
		// -aa[i]y[i-1]+cc[i]y[i]-bb[i]y[i+1] = ff[i]
		// aa[i]y[i-1] = ff[i]-cc[i]y[i]+bb[i]y[i+1]
		// bb[i]y[i+1] = ff[i]-cc[i]y[i]-aa[i]y[i-1]

		// Разностная схема
		// 1/2h*(k[i+1/2]*(y[i+1]-y[i])-k[i-1/2]*(y[i]-y[i-1])) - q[i]*y[i] = -f[i]


		// Цикл
		// y0[i] = uu[i]
		// y1[i] = uu[i+1/2]-uu[i-1/2]
		// y2[i] = uu[i+1]-uu[i+1]-2*uu[i]
		// ff[i] =-k1(xx[i])*y1[i]*hx-k(xx[i])*y2[i]+q(xx[i])*y0[i]*hx2

		// Дополнительные условия

		// режим 1
		// u'(xa) = pi*b*u(xa)
		// u'(xb) = -pi*b*u(xb)

		// режим 2
		// u(xa) = u(xb)
		// u'(xa) = u'(xb)

		// Оба режима сводятся к виду
		// (u,u')(xa) LR (u,u')(xb)


		MyRange(np,mp,0,nx,&i1,&i2,&nc);
		ncm = (nc-1)<<it; // Старший индекс в локальном массиве
		i1<<=it; // Младший индекс в глобальном массиве
		i2<<=it; // Старший индекс в глобальном массиве
		nc = ncm+1; // Размер локального массива
		ncp = 2*(np-1); 
		ncx = imax(nc,ncp);

		hx = (xb-xa)/(nx<<it); 
		hx2 = hx * hx;

		fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

		for (i=0; i<nc; i++) xx[i] = xa + hx * GIDX(i,i1,nx);
		for (i=0; i<nc; i++) kk[i] = k(xx[i]);
		for (i=0; i<nc; i++) kk1[i] = k1(xx[i]);
		for (i=0; i<nc; i++) qq[i] = q(xx[i]);

		tau = 1.0;

		s0 = dabs(kk[0]);for (i=0; i<nc; i++) s0 = dmax(s0,dabs(kk[i]));
		s1 = dabs(kk1[0]);for (i=0; i<nc; i++) s1 = dmax(s1,dabs(kk1[i]));
		s2 = dabs(qq[0]);for (i=0; i<nc; i++) s2 = dmax(s2,dabs(qq[i]));

		s2 = dmax(s0,s2);
		s2 = dmax(s1,s2);

		gt = s0; MPI_Allreduce(&gt,&s0,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		gt = s1; MPI_Allreduce(&gt,&s1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
		gt = s2; MPI_Allreduce(&gt,&s2,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

		s0 = dmax(s0,sqrt(s0));
		s1 = dmax(s1,sqrt(s1));
		s2 = dmax(s2,sqrt(s2));

		tau = 0.25*hx/s2;
		gam = tau/hx2;

		if (mp==0) {
			s0 = kk[0]; 
			s2 = kk[1];
			aa[0] = 0.0;
			bb[0] = gam * 2.0 * s0 * s2 / (s0 + s2);
			cc[0] = qq[0] + bb[0];
		}
		else {
			s0 = kk[0]; 
			s1 = k(xx[0]-hx); 
			s2 = k(xx[0]+hx);
			aa[0] = gam * 2.0 * s0 * s1 / (s0 + s1);
			bb[0] = gam * 2.0 * s0 * s2 / (s0 + s2);
			cc[0] = qq[0] + aa[0] + bb[0];
		}

		for (i=1; i<ncm; i++) {
			s0 = kk[i]; 
			s1 = kk[i-1]; 
			s2 = kk[i+1];
			aa[i] = gam * 2.0 * s0 * s1 / (s0 + s1);
			bb[i] = gam * 2.0 * s0 * s2 / (s0 + s2);
			cc[i] = qq[i] + aa[i] + bb[i];
		}

		if (mp==np-1) {
			s0 = kk[ncm]; 
			s1 = kk[ncm-1];
			aa[ncm] = gam * 2.0 * s0 * s1 / (s0 + s1);
			bb[ncm] = 0.0;
			cc[ncm] = qq[ncm] + aa[ncm];
		}
		else {
			s0 = kk[ncm]; 
			s1 = k(xx[ncm]-hx); 
			s2 = k(xx[ncm]+hx);
			aa[ncm] = gam * 2.0 * s0 * s1 / (s0 + s1);
			bb[ncm] = gam * 2.0 * s0 * s2 / (s0 + s2);
			cc[ncm] = qq[ncm] + aa[ncm] + bb[ncm];
		}


		sprintf(sname,"%s_%02d_kk.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,kk);
		sprintf(sname,"%s_%02d_qq.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,qq);
		sprintf(sname,"%s_%02d_aa.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,aa);
		sprintf(sname,"%s_%02d_bb.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,bb);
		sprintf(sname,"%s_%02d_cc.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,cc);
		// Рассчитываем функцию для самопроверки
		for(i=0;i<ncm;i++) uu[i] = u(xx[i]); 
		sprintf(sname,"%s_%02d_uu.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,uu);

		ntv = 0; gt = 1.0;

		// Задаём начальное значение функции
		// В качестве начального выбираем ненулевой вектор
		// поскольку ноль переходит в ноль
		for(i=0;i<nc;i++) uu[i] = 1; 
		for(i=0;i<nc+1;i++) y0[i] = 0; 
		for(i=0;i<nc+1;i++) y1[i] = 0; 
		for(i=0;i<nc+1;i++) y2[i] = 0; 

		do{
			ntv++; 

			if(mp==0) uu[0]=1;
			for(i=0;i<nc;i++) y0[i] = uu[i];
			y0[nc] = y0[nc-1]+y1[nc-1]*hx;

			if(mode==1&&mp==np-1){
				// u'(xa) = pi*b*u(xa)
				// u'(xb) = -pi*b*u(xb)
				s0=(y0[ncm]+y1[ncm]/pi/b)/2;
				s1=(pi*b*y0[ncm]+y1[ncm])/2;
				y0[ncm]=s0;
				y1[ncm]=s1;
			}

			if(mode==1&&mp==0){
				// u'(xa) = pi*b*u(xa)
				// u'(xb) = -pi*b*u(xb)
				s0=(y0[0]+y1[0]/pi/b)/2;
				s1=(pi*b*y0[0]+y1[0])/2;
				y0[0]=s0;
				y1[0]=s1;
			}

			bsl[0]=y0[0];
			bsl[1]=y1[0];
			bsr[0]=y0[ncm]; 
			bsr[1]=y1[ncm];
			CircleABndExch1D(np, mp, 2, 2, 2, 2, bsl, brl, bsr, brr);
			if(mode==2&&mp==0) y0[0]=(y0[0]+brl[0]); // u(xa)==u(xb)
			if(mode==2&&mp==0) y1[0]=(y1[0]+brl[1]); // u'(xa)==u'(xb)
			if(mode==2&&mp==np-1) y0[ncm]=(y0[ncm]+brr[0]); // u(xa)==u(xb)
			if(mode==2&&mp==np-1) y1[ncm]=(y1[ncm]+brr[1]); // u'(xa)==u'(xb)
			y0[nc]=brr[0];
			y1[nc]=brr[1];

			if(mode==1)for(i=nc;i-->0;) y1[i] = (y0[i+1]-y0[i])/hx;
			if(mode==2)for(i=0;i<nc+1;i++) y1[i] = (y0[i+1]-y0[i])/hx;

			for(i=0;i<nc;i++) y2[i] = (y1[i+1]-y1[i])/hx;
			for(i=0;i<nc;i++) ff[i] = qq[i]*y0[i] - tau*(kk1[i]*y1[i] + kk[i]*y2[i]);

			if (np<2) ier = prog_right(nc,aa,bb,cc,ff,al,uu);
			else      ier = prog_rightpm(np,mp,nc,0,aa,bb,cc,ff,al,uu,_y2,_y3,_y4);

			t2 = MPI_Wtime() - t1;

			// Вычисляем отличие от реальной функции
			gt = 0.0;
			for (i=0; i<nc; i++) {
				s1 = u(xx[i]); s2 = dabs(s1-uu[i]); gt = dmax(gt,s2);
				if (lp>0)
					fprintf(Fo,"i=%8d x=%12le y=%12le u=%12le gt=%12le\n",
					i,xx[i],uu[i],s1,gt);
			}

			if (np>1) {
				s1 = gt; 
				MPI_Allreduce(&s1,&gt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			}
			if (mp==0) fprintf(stderr,"tau=%le ntv=%d nx=%d t1=%le dmax=%le\n",tau,ntv,nx<<it,t2,gt);
			fprintf(Fo,"t1=%le dmax=%le\n",t1,gt);

			// повторяем если большое расхождение с текущим значением в качестве начального
		} while ((ntv<ntm) && (gt>epst));

		sprintf(sname,"%s_%02d_y0.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,y0);
		sprintf(sname,"%s_%02d_y1.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,y1);
		sprintf(sname,"%s_%02d_y2.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,y2);
		sprintf(sname,"%s_%02d_ff.dat",vname,np);
		OutFun1DP(sname,np,mp,nc,xx,ff);

		sprintf(sname,"%s_%02d_%02d.dat",vname,np,it);
		OutFun1DP(sname,np,mp,nc,xx,uu);

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
