#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;


//var 3 -(1/(x-3))u''+(1+x/2)u'+e^(x/2)u=2-x, u(-1)=u(1)=0 -> a=-1, b=1, alpha=beta=0
//p(x), q(x), r(x), f(x) -- соответсвующие функции u'', u', u и правой части соответственно
double p(double x) {
	return -(1.0/(x-3));
}

double q(double x) {
	return 1+x/2;
}

double r(double x) {
	return exp(x/2);
}

double f(double x) {
	return 2-x;
}

vector<vd> y_ex {	//Точное решение
	{-1,   -0.9, 	-0.8,   -0.7,   -0.6,   -0.5, 
		   -0.4, 	-0.3,   -0.2,   -0.1,      0, 
		  	0.1,  	 0.2,    0.3,    0.4,    0.5, 
		  	0.6,  	 0.7,    0.8,    0.9,     1},
	{0,  0.2916,  0.5356, 0.7361, 0.9016, 1.0331, 
		 1.1380,  1.2156, 1.2734, 1.3259, 1.3273,
		 1.3258,  1.3103, 1.2735, 1.2204, 1.1315,
		1.03330,  0.8847, 0.6800, 0.3930,      0}
};

vd tridiag_algo(double a, double b, int n) { //Основная функция решения
	double h=(b-a)/n;
	vector<double> A(n+1), B(n+1), C(n+1), G(n+1), s(n+1), t(n+1), y(n+1), x(n+1);
	y[0]=0, s[0]=0, t[0]=0, A[0]=0, C[0]=0, G[0]=0, x[0]=a; 
	y[n]=0, s[n]=0, t[n]=0, A[n]=0, C[n]=0, G[n]=0, x[n]=b;
	B[0]=-1, B[n]=-1;
	for (int i=1; i<n; i++) {  //Прогонка "вниз"
		x[i]=a+i*h;
		A[i]=-p(x[i])/(h*h)-q(x[i])/(2*h);
		B[i]=-2*p(x[i])/(h*h)-r(x[i]);
		C[i]=-p(x[i])/(h*h)+q(x[i])/(2*h);
		G[i]=f(x[i]);
		s[i]=C[i]/(B[i]-A[i]*s[i-1]);
		t[i]=(A[i]*t[i-1]-G[i])/(B[i]-A[i]*s[i-1]);
	}
	for (int i=n-1; i>0; i--) {	  //Прогонка "вверх" для нахождения искомых "y"
		y[i]=s[i]*y[i+1]+t[i];
	}
	cout << "\n-------------\nCase: n = " << n << ":\n"; //Вывод всех коэффициентов и игреков
	printf("   X          A          B          C        G        s        t         y\n");
	for (int i=0; i<n+1; i++) {
		printf("%.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n", x[i], A[i], B[i], C[i], G[i], s[i], t[i], y[i]);
	}
	return y;
}

vd Richard_extr (vd h_half_table, vd h_table, int s) {   //Экстраполяция по ричардсону, где half_table - игреки для разбиения 2*n
	vd res(h_table.size(), 0);
	for (int i=0; i<h_table.size(); i++) {
		res[i]=h_half_table[2*i]+(h_half_table[2*i]-h_table[i])/(pow(2, s)-1);
	}
	return res;
}

int main()
{
	double a=-1.0, b=1.0;
	//Случай n=20//
	vd y20 = tridiag_algo(a, b, 20);
	//------------------//
	//Случай n=10//
	vd y10 = tridiag_algo(a, b, 10);
	vd newy=Richard_extr(y20, y10, 2);	//Заметим, что мы берем s=2, так как у нас алгоритм O(h^2)
	printf("      X          Y_ex           Y_ut       Difference (modulo)\n");
	for (int i=0; i<11; i++) {
		printf("%.8f   %.8f   %.8f   %.8f\n", y_ex[0][2*i], y_ex[1][2*i], newy[i], abs(y_ex[1][2*i]-newy[i]));
	}
  	return 0;
}
