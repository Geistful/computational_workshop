#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;

vector<vd> E_matr(int deg) { //Единичная матрица
	vector<vd> e(deg, vd (deg, 0));
	for (int i=0; i<deg; i++) {
		e[i][i]=1;
	}
	return e;
}

int sign (double temp) {	//Функция "знак" - нужна для подсчета константы "s"
	if (temp>0) return 1;
	if (temp<0) return -1;
	return 0;
}

double set_p_q_max(vector<vd> m, int &p, int &q) { //Позиция наибольшего по модулю числа
	double maxvalue=0;
	for (int i=0; i<m.size(); i++) {
		for (int j=i+1; j<m.size(); j++) {
			if (abs(m[i][j])>maxvalue) {
				p=i, q=j;
				maxvalue = abs(m[i][j]);
			}
		}
	}
	return maxvalue;
}

vector<vd> matr_x_matr (vector<vd> m1, vector<vd> m2) { //Произведение матриц
  vector<vd> mult(m1.size(), vd (m2[0].size(), 0));
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m2[0].size(); i++) {
		for (int j = 0; j < m1[0].size(); j++) {
			mult[k][i] += m1[k][j]*m2[j][i];
		}
	}
  }
  return mult;
} 

vd matr_x_vect (vector<vd> m1, vd m2) { //Произведение матрицы на вектор
  vd mult(m1.size(), 0);
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1[0].size(); i++) {
		mult[k] += m1[k][i]*m2[i];
	}
  }
  return mult;
}

vd dig_x_vect (double lambda, vd xi) {	//Число (дабл) умножается на вектор
	for (int i=0; i<xi.size(); i++) {
		xi[i]*=lambda;
	}
	return xi;
} 

//Подсчет констант d, с и s
double get_d (vector<vd> m, int p, int q) {
	return sqrt((m[p][p]-m[q][q])*(m[p][p]-m[q][q])+4*m[p][q]*m[p][q]);
}

double get_c (vector<vd> m, int p, int q) {
	return sqrt(0.5*(1+(abs(m[p][p]-m[q][q])/get_d(m, p, q))));
}

double get_s (vector<vd> m, int p, int q) {
	return sign(m[p][q]*(m[p][p]-m[q][q]))*
		   sqrt(0.5*(1-(abs(m[p][p]-m[q][q])/get_d(m, p, q))));
}

//Получение следующей итерации изначальной матрицы
vector<vd> a_get_next (vector<vd> a, int p, int q, double c, double s) {
	vector<vd> newa(a.size(), vd (a[0].size(), 0));
  	copy(begin(a), end(a), begin(newa));
	for (int i=0; i<a.size(); i++) {
		for (int j=0; j<a.size(); j++) {
			if (i!=p && i!=q && j!=p && j!=q)
				newa[i][j]=a[i][j];
			else if (i!=p && i!=q && j==p) {
				newa[i][j]=c*a[i][p]+s*a[i][q];   //change sign before s
				newa[j][i]=newa[i][j];
			}
			else if (i!=p && i!=q && j==q) {
				newa[i][j]=-s*a[i][p]+c*a[i][q];  //change sign before s
				newa[j][i]=newa[i][j];
			}
			else if (i==p && j==p)
				newa[i][j]=c*c*a[p][p]+2*c*s*a[p][q]+s*s*a[q][q];
			else if (i==q && j==q)
				newa[i][j]=s*s*a[p][p]-2*c*s*a[p][q]+c*c*a[q][q];
			else if (i==p && j==q) {
				newa[i][j]=0/*(c*c-s*s)*a[p][q]+c*s*(a[q][q]-a[p][p])*/;
				newa[j][i]=newa[i][j];
			}
		}
	}
	return newa;
}

void print_matr(vector<vd> m) { //Печать матрицы (обычная)
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m.size(); j++) {
			if (j==m.size()-1)
				cout << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

//Получаем следующую матрицу, в которой должны получиться собственные значения
vector<vd> x_get_next (vector<vd> x, int p, int q, double c, double s) {
	vector<vd> newx(x.size(), vd (x[0].size(), 0));
  	copy(begin(x), end(x), begin(newx));
	for (int i=0; i<newx.size(); i++) { //col
		for (int j=0; j<newx.size(); j++) { //row
			if (i!=p && i!=q)
				newx[j][i]=x[j][i];
			else if (i==p)
				newx[j][i]=c*x[j][p]+s*x[j][q];
			else if (i==q)
				newx[j][i]=-s*x[j][p]+c*x[j][q];
		}
	}
	return newx;
}

vector<vd> transpose_matr (vector<vd> m) { //Транспонирование матрицы
	for (int i = 0; i < m.size(); ++i) {
		for (int j=i+1; j<m.size(); j++) {
			double temp = m[i][j];
			m[i][j] = m[j][i];
			m[j][i] = temp;
		}	
	}
	return m;
}

void print_vect(vd v) { 	//Печать вектора
	for (int i=0; i<v.size(); i++) {
		if (i==0) 
			cout << "(" << v[i] << ", ";
		else if (i==v.size()-1)
			cout << v[i] << ")\n";
		else
			cout << v[i] << ", ";
	}
}

int main()
{
	vector<vd> matr {
  	{-1.48213, -0.03916, 1.08254},
	{-0.03916, 1.13958, 0.01617},
	{1.08254, 0.01617, -1.48271}
  	};
  	double eps = 0.00000000001, c, s;
  	int p, q;
  	vector<vd> x = E_matr(matr.size());
  	vector<vd> newmatr(matr.size(), vd (matr[0].size(), 0));
  	copy(begin(matr), end(matr), begin(newmatr));
  	cout << "Initial A:\n";
  	print_matr(matr);
  	cout << "Initial X:\n";
  	print_matr(x); 
  	while (set_p_q_max(newmatr, p, q)>=eps) {
  		c=get_c(newmatr, p, q);
		s=get_s(newmatr, p, q);
  		newmatr = a_get_next(newmatr, p, q, c, s);
  		x = x_get_next(x, p, q, c, s);
  	}
  	cout << "Final A matrix:\n";
  	print_matr(newmatr);
  	cout << "Eigenvalues:\n";
  	for (int i=0; i<newmatr.size(); i++) {
  		cout << "lambda " << i+1 << " = " << newmatr[i][i] << "\n";
  	}
  	cout << "Final X matrix:\n";
  	print_matr(x);
  	cout << "Eigenvectors:\n";
  	vector<vd> row_x = transpose_matr(x);
  	for (int i=0; i<x.size(); i++) {
  		cout << "E.v. " << i+1 << ": "; 
  		print_vect(row_x[i]);
  		double sum=0;
  		for (int j=0; j<row_x[i].size(); j++) {
  			sum+=(row_x[i][j])*(row_x[i][j]);
  		}
  		sum = sqrt(sum);
  		cout << "sum = " << sum <<"\n";
  		cout << " --> ";
  		print_vect(dig_x_vect(1.0/sum, row_x[i]));
  		cout << " (normalized)\n";
  	}
  	for (int i=0; i<matr.size(); i++) {
  		cout << "A["<< i+1 <<"]x:\n";
  		print_vect(matr_x_vect(matr, row_x[i]));
  		cout << "lambda["<< i+1 <<"]*x:\n";
  		print_vect(dig_x_vect(newmatr[i][i], row_x[i]));
  	}
	return 0;
}
