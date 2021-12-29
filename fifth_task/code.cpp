#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;


//y'(x)=1-sin(1.25*x+y)-(0.1*y)/(2+x), y(0)=0, h=0.05, N=10
//for 
int MAX = 1000;

double f(double x, double y) {
	return 1-sin(1.25*x+y)-(0.1*y)/(2+x);
}

vector<vd> y_table {
	{0, 	 0.05, 		0.10,	  0.15, 	0.20, 	  0.25, 
			 0.30, 		0.35, 	  0.40, 	0.45, 	  0.50,
			 0.55, 		0.60, 	  0.65, 	0.70, 	  0.75,
			 0.80, 		0.85, 	  0.90, 	0.95, 	  1.0},
	{0, 0.0472025, 0.0889605, 0.125664,  0.15769, 0.185419, 
		 0.209224,  0.229468, 0.246495, 0.260636,   0.2722,
		  0.28148,  0.288752, 0.294274, 0.298293,  0.30104,
		 0.302738,  0.303599, 0.303827, 0.303622, 0.303179}
};

vector<vd> Euler1(double x0, double h, int N) {
	vector<vd> table(2, vd(N+1, 0));
	int i0 = int(x0/h);
	table[0][0] = x0;
	table[1][0] = y_table[1][i0];
	double x_m, y_m;
	for (int i=1; i<N+1; i++) {
		x_m = table[0][i-1];
		y_m = table[1][i-1];
		table[0][i] = x0+i*h;
		table[1][i] = y_m + h*f(x_m, y_m);
	}
	return table; 
}

/*vector<vd> Euler2(double x0, double h, int N) {
	vector<vd> table(2, vd(N+1, 0));
	int i0 = int(x0/h);
	table[0][0] = x0;
	table[1][0] = y_table[1][i0];
	double x_m, y_m;
	for (int i=1; i<N+1; i++) {
		x_m = table[0][i-1];
		y_m = table[1][i-1];
		table[0][i] = x0+i*h;
		table[1][i] = y_m + h*f(x_m+h/2, y_m+h/2*f(x_m, y_m));
	}
	return table; 
}*/

vector<vd> Euler3(double x0, double h, int N) {
	vector<vd> table(2, vd(N+1, 0));
	int i0 = int(x0/h);
	table[0][0] = x0;
	table[1][0] = y_table[1][i0];
	double x_m, y_m, ynext;
	for (int i=1; i<N+1; i++) {
		x_m = table[0][i-1];
		y_m = table[1][i-1];
		table[0][i] = x0+i*h;
		ynext = y_m + h*f(x_m, y_m);
		table[1][i] = y_m + h/2*(f(x_m, y_m)+f(table[0][i], ynext));
	}
	return table; 
}

vector<vd> Richard_extr (vector<vd> h_half_table, vector<vd> h_table, int s) {
	vector<vd> res(2, vd(h_table[0].size(), 0));
	for (int i=0; i<h_table[0].size(); i++) {
		res[0][i]=h_table[0][i];
		res[1][i]=h_half_table[1][i]+(h_half_table[1][i]-h_table[1][i])/(pow(2, s)-1);
	}
	return res;
}

vector<vd> RungeKutt(double x0, double h, int N) {
	vector<vd> table(2, vd(N+1, 0));
	table[0][0] = 0;
	table[1][0] = y_table[1][x0];
	double x_m, y_m, k1, k2, k3, k4;
	for (int i=1; i<N+1; i++) {
		x_m = table[0][i-1];
		y_m = table[1][i-1];
		k1 = h*f(x_m, y_m);
		k2 = h*f(x_m+h/2, y_m+k1/2);
		k3 = h*f(x_m+h/2, y_m+k2/2);
		k4 = h*f(x_m+h, y_m+k3);
		table[0][i] = i*h;
		table[1][i] = y_m + 1.0/6*(k1+2*k2+2*k3+k4);
	}
	return table;
}

vector<vd> Adams(double x0, double xn, double h, vector<vd> runge) {
	int i0 = int(x0/h);
	int N=int((xn-x0)/h)+1;
	vector<vd> table(2, vd(N+1, 0));
	table[0][0] = x0;
	table[1][0] = runge[1][i0];
	//double dq[5][MAX];
	vector<vd> dq(5, vd(MAX, 0));
	//Записываем q0
	dq[0][0] = h*f(x0, runge[1][i0]);
	//Первые пять значений - точные
	//Также записываем q1, q2, q3, q4 в dq
	for (int i=1; i<5; i++) {
		table[0][i] = x0+i*h;
		table[1][i] = runge[1][i0+i];
		dq[0][i] = h*f(table[0][i], table[1][i]);
	}
	//По имеющимся q получаем (d)q (d^2)q (d^3)q (d^4)q
	for (int i=1; i<5; i++) {
		for (int j=0; j<5-i; j++) {
			dq[i][j] = dq[i-1][j+1]-dq[i-1][j];
		}
	}
	//Дальше вычисляем y_(m+1) по формуле и дополняем dq
	for (int i=5; i<N+1; i++) {
		table[0][i] = x0+i*h;
		table[1][i] = table[1][i-1] + dq[0][i-1] + 1.0/2*dq[1][i-2] 
				      + 5.0/12*dq[2][i-3] + 3.0/8*dq[3][i-4]
				      + 251.0/720*dq[4][i-5];
		dq[0][i] = h*f(table[0][i], table[1][i]);
		for (int j=1; j<5; j++) {
			dq[j][i-j] = dq[j-1][i-j+1]-dq[j-1][i-j]; 
		}
	}
	return table;
}

vector<vd> Adams_in (double x0, double h, int N, double eps) {
	vector<vd> table(2, vd(N+1, 0));
	int i0 = int(x0/h);
	double prev=0;
	table[0][0] = x0;
	table[1][0] = y_table[1][i0];
	//double dq[5][MAX];
	vector<vd> dq(5, vd(MAX, 0));
	//Записываем q0
	dq[0][0] = h*f(x0, y_table[1][i0]);
	//Первые пять значений - точные
	//Также записываем q1, q2, q3, q4 в dq
	for (int i=1; i<5; i++) {
		table[0][i] = x0+i*h;
		table[1][i] = y_table[1][i0+i];
		dq[0][i] = h*f(table[0][i], table[1][i]);
	}
	//По имеющимся q получаем (d)q (d^2)q (d^3)q (d^4)q
	for (int i=1; i<5; i++) {
		for (int j=0; j<5-i; j++) {
			dq[i][j] = dq[i-1][j+1]-dq[i-1][j];
		}
	}
	//Дальше вычисляем y_(m+1) по формуле и дополняем dq
	for (int i=5; i<N+1; i++) {
		table[0][i] = x0+i*h;
		table[1][i] = table[1][i-1] + dq[0][i-1] + 1.0/2*dq[1][i-2] 
				      + 5.0/12*dq[2][i-3] + 3.0/8*dq[3][i-4]
				      + 251.0/720*dq[4][i-5]; //y_m из метода интерполяции
		do {
			prev = table[1][i];
			dq[0][i] = h*f(table[0][i], table[1][i]);
			for (int j=1; j<5; j++) {
				dq[j][i-j] = dq[j-1][i-j+1]-dq[j-1][i-j]; 
			}
			table[1][i] = table[1][i-1] + dq[0][i] - 1.0/2*dq[1][i-1]
						  - 1.0/12*dq[2][i-2] - 1.0/24*dq[3][i-3]
						  - 19.0/720*dq[4][i-4]; //пересчитываем для интерполяции
		} while (abs(table[1][i]-prev)>=eps);
	}
	return table;
}

void printtable (vector<vd> table, vector<vd> y, double x0, double h, int N) {
	int i0 = int(x0/h);
	cout << "   x_i                y_i             Abs. diff.\n";
	for (int i=0; i<N+1; i++) {
		printf("%f   ->   %.12f     %.12f\n", table[0][i], table[1][i], abs(table[1][i]-y[1][i0+i]));
	} 
	cout << "--------------\n";
}

int main()
{
  //ios::sync_with_stdio(0);
  //cin.tie(0);
  double x = 0.0;
  double x0=0, h=0.05;
  int N = 10;
  string s="y";
  //-----------------------------------
  cout << "Test parameters:\n" << "h = " << h << ", N = " << N << "\n";
  vector<vd> Eul1 = Euler1(x0, h, N);
  //vector<vd> Eul2 = Euler2(x0, h, N);
  vector<vd> Eul3 = Euler3(x0, h, N);
  vector<vd> Runge = RungeKutt(x0, h, 2*N);
  vector<vd> Adam = Adams(4*h, 1, h, Runge);
  vector<vd> Adam_int = Adams_in(4*h, h, int((1-4*h)/h)+1, 0.00001);
  cout << "Euler:\n";
  printtable(Eul1, y_table, x0, h, N);
  /*cout << "Euler improved #1:\n";
  printtable(Eul2, y_table, x0, h, N);*/
  cout << "Euler improved:\n";
  printtable(Eul3, y_table, x0, h, N);
  cout << "Runge-Kutt:\n";
  printtable(Runge, y_table, x0, h, 2*N);
  cout << "Adams:\n";
  printtable(Adam, y_table, 4*h, h, int((1-4*h)/h));
  cout << "Adams interp.:\n";
  printtable(Adam_int, y_table, 4*h, h, int((1-4*h)/h));
  cout << "Summury: abs. diff. for y_N:\n";
  printf("Euler:             %.12f\n", abs(Eul1[1][N]-y_table[1][N]));
  //printf("Euler improved #1: %.12f\n", abs(Eul2[1][N]-y_table[1][N]));
  printf("Euler improved: %.12f\n", abs(Eul3[1][N]-y_table[1][N]));
  printf("Runge-Kutt:        %.12f\n", abs(Runge[1][2*N]-y_table[1][2*N]));
  printf("Adams:             %.12f\n", abs(Adam[1][int((1-4*h)/h)]-y_table[1][2*N]));
  printf("Adams interp:      %.12f\n", abs(Adam_int[1][int((1-4*h)/h)]-y_table[1][2*N]));
  return 0;
}
