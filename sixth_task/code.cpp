#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;

//y1'=-125y1+123.45y2,
//y2'=123.45y1-123y2, y1(0)=1, y2(0)=1.
// а) Явный м. Э; б) Неявн. м. Э.; в) Интер. м. А. 3 пор.
// 1) [0, 0.5] точное реш. в точках t_i=ih, i=1..5 при h=0.1
// 2) [0, 0.5] прибл. реш --||-- при h=0.05  для а) б) в)
// 3)Исследовать на устойчивость
// 4)повторить п2 для h=0.001
double a=0, b=0.5;

vector<vd> y_ex {	//Точное решение [0, 0.5] h=0.1
	{0.008100, 0.003819, 0.003616, 0.003424, 0.003242, 0.003070},
	{0.008100, 0.950688, 0.900177, 0.852349, 0.807062, 0.764182}
};
//-------------From task #4----------//
vector<vd> E(int deg) {
	vector<vd> e(deg, vd (deg, 0));
	for (int i=0; i<deg; i++) {
		e[i][i]=1;
	}
	return e;
}

int sign (double temp) {
	if (temp>0) return 1;
	if (temp<0) return -1;
	return 0;
}

double set_p_q_max(vector<vd> m, int &p, int &q) {
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

vector<vd> matr_plus_matr (vector<vd> m1, vector<vd> m2) { //Поэлементная сумма матриц
  vector<vd> mult(m1.size(), vd (m1.size(), 0));
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1.size(); i++) {
		mult[k][i] = m1[k][i]+m2[k][i];
	}
  }
  return mult;
}

vector<vd> matr_x_matr (vector<vd> m1, vector<vd> m2) {
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

vd matr_x_vect (vector<vd> m1, vd m2) {
  vd mult(m1.size(), 0);
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1[0].size(); i++) {
		mult[k] += m1[k][i]*m2[i];
	}
  }
  return mult;
}

vd vect_sum (vd a, vd b) {		//Сумма векторов
	for (int i=0; i<a.size(); i++) {
		a[i]+=b[i];
	}
	return a;
}

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
//----------------------//

vector<vd> swap_row(vector<vd> &m, int p, int k) { //Меняем строки
	if (p==k) return m;
	double tmp;
	for (int i=0; i<m[0].size(); i++) {
		tmp = m[p][i];
		m[p][i] = m[k][i];
		m[k][i] = tmp;
	}
	return m;
}

vector<vd> swap_col(vector<vd> &m, int q, int k, vi &xorder) { //Меняем столбцы
	if (q==k) return m;
	double tmp;
	for (int i=0; i<m.size(); i++) {
		tmp = m[i][q];
		m[i][q] = m[i][k];
		m[i][k] = tmp;
	}
	tmp = xorder[q];
	xorder[q]=xorder[k];
	xorder[k]=tmp;
	return m;
}

vi gauss_forw_elim (vector<vd> &newmatr) {	//"Прямой ход" метода Гаусса
	vi xorder;
 	for (int i=0; i<newmatr.size(); i++) {
  		xorder.push_back(i);
  	}
  	int p, q;
  	set_p_q_max(newmatr, p, q);
  	double temp, maxval;
  	for (int k=0; k<newmatr.size(); k++) {
	  	maxval=0;
	  	swap_col(newmatr, q, k, xorder);
	  	swap_row(newmatr, p, k);
	  	temp=newmatr[k][k];
	  	for (int j=k; j<newmatr[0].size(); j++) {
	  		newmatr[k][j]/=temp;
	  	}
	  	for (int i=k+1; i<newmatr.size(); i++) {
	  		temp = newmatr[i][k];
	  		for (int j=k; j<newmatr[0].size(); j++) {
	  			newmatr[i][j]-=newmatr[k][j]*temp;
	  			if (abs(newmatr[i][j])>abs(maxval) && j!=newmatr[0].size()-1) {
	  				p=i, q=j;
	  				maxval = newmatr[i][j];
	  			}
	  		}
	  	}
	}
	return xorder;
}

vd jordan_method (vector<vd> &newmatr) {	//Метод Жордана - обратный ход оставляет только диагональ
	vi xorder = gauss_forw_elim(newmatr);
	vd coef;
	for (int k=1; k<newmatr.size(); k++) {
		for (int i=k; i<newmatr[0].size(); i++) {
			for (int j=0; j<k; j++) {
				if (i==k) {
					if (j==k-1) coef.push_back(newmatr[j][i]);
					else coef[j]=newmatr[j][i];
				}
				newmatr[j][i]-=newmatr[k][i]*coef[j];
			}
		}
	}
	vd x(newmatr.size(), 0);		//Здесь тоже сортируем корни
  	for (int i=0; i<newmatr.size(); i++) {
  		x[i]=newmatr[xorder[i]][newmatr[0].size()-1];
  	}
	return x;
}

vector<vd> unite_a_b (vector<vd> a, vd b) {		//Соединяем матрицы (A+B=A|B)
	vector<vd> newa(a.size(), vd(a.size()+1, 0));
	for (int i=0; i<newa.size(); i++) {
		for (int j=0; j<newa[0].size(); j++) {
			if (j==newa[0].size()-1)
				newa[i][j]=b[i];
			else newa[i][j]=a[i][j];
		}
	}
	return newa;
}

vector<vd> rev_matrix (vector<vd> matr) {		//Получение обратной матрицы
	vd b (matr.size(), 0);
	vector<vd> bigmatr = unite_a_b(matr, b);
	vector<vd> rev (bigmatr.size(), vd (bigmatr.size(), 0));
	for (int i=0; i<bigmatr.size(); i++) {
		vector<vd> temp = bigmatr;
		for (int j=0; j<temp.size(); j++) {
			if (j==i) temp[j][temp[0].size()-1]=1;
			else temp[j][temp[0].size()-1]=0;
		}
		vd x = jordan_method(temp);		//Используем метод Жодана для i=0..n
		for (int k=0; k<bigmatr.size(); k++) {
			rev[k][i] = x[k];
		}
	}
	return rev;
}

vector<vd> num_x_matr (double l, vector<vd> m) {
	for (int i=0; i<m.size(); i++) {
		for (int j=0; j<m[0].size(); j++) {
			m[i][j]*=l;
		}
	}
	return m;
}

vd double_x_vect (double num, vd v) {	//Умножение скаляра на массив
	vd res(v.size(), 0);
	for (int i=0; i<v.size(); i++) {
		res[i]=num*v[i];
	}
	return res;
}

vd get_eivals(vector<vd> matr, double eps) {
	vd eivals (matr.size(), 0);
	int p, q;
	double c, s;
  	vector<vd> x = E(matr.size());
  	vector<vd> newmatr(matr.size(), vd (matr[0].size(), 0));
  	copy(begin(matr), end(matr), begin(newmatr));
  	while (set_p_q_max(newmatr, p, q)>=eps) {
  		c=get_c(newmatr, p, q);
		s=get_s(newmatr, p, q);
  		newmatr = a_get_next(newmatr, p, q, c, s);
  		x = x_get_next(x, p, q, c, s);
  	}
  	for (int i=0; i<newmatr.size(); i++) {
  		eivals[i]=newmatr[i][i];
  		//cout << "lambda " << i+1 << " = " << newmatr[i][i] << "\n";
  	}
  	return eivals;
}

double max_elem_vect (vd vect) {
	double maxval=0;
	for (int i=0; i<vect.size(); i++) {
		maxval = max(maxval, abs(vect[i]));
	}
	return maxval;
}

vector<vd> euler_meth (vector<vd> m, double h) {
	vd y_cur={1, 1};
	int n=int((b-a)/h+1);
	vector<vd> eul(2, vd (n, 1));
	vector<vd> W = matr_plus_matr(E(m.size()), num_x_matr(h, m));
	for (int i=1; i<n; i++) {
		y_cur = matr_x_vect(W, y_cur);
		eul[0][i]=y_cur[0];
		eul[1][i]=y_cur[1];
	}
	return eul;
}

vector<vd> rev_euler_meth (vector<vd> m, double h) {
	vd y_cur={1, 1};
	int n=int((b-a)/h+1);
	vector<vd> rev_eul(2, vd (n, 1));
	vector<vd> W = rev_matrix(matr_plus_matr(E(m.size()), num_x_matr(-h, m)));
	for (int i=1; i<n; i++) {
		y_cur = matr_x_vect(W, y_cur);
		rev_eul[0][i]=y_cur[0];
		rev_eul[1][i]=y_cur[1];
	}
	return rev_eul;
}

vector<vd> adams_method (vector<vd> m, double h) {
	vd y_cur(2), y_prev(2), newy;
	int n=int((b-a)/h+1);
	vector<vd> adams(2, vd (n, 1));
	vector<vd> rev_sol = rev_euler_meth(m, h);
	adams[0][1]=rev_sol[0][1];
	adams[1][1]=rev_sol[1][1];
	vector<vd> w_first = rev_matrix(matr_plus_matr(E(m.size()), num_x_matr(-(5*h)/12, m)));
	vector<vd> W1=matr_x_matr(w_first, matr_plus_matr(E(m.size()), num_x_matr(2*h/3, m)));
	vector<vd> W2=matr_x_matr(w_first, num_x_matr(h/12, m));
 	for (int i=2; i<n; i++) {
 		y_prev[0] = adams[0][i-2], y_prev[1]=adams[1][i-2];
 		y_cur[0]=adams[0][i-1], y_cur[1]=adams[1][i-1];
		newy = vect_sum(matr_x_vect(W1, y_cur), double_x_vect(-1.0, matr_x_vect(W2, y_prev)));
		adams[0][i]=newy[0];
		adams[1][i]=newy[1];
	}
	return adams;
}

void print_for_h (vector<vd> m, double h) {
	int k = int(0.1/h);
	vector<vd> eul = euler_meth(m, h);
	vector<vd> reveul = rev_euler_meth(m, h);
	vector<vd> adams = adams_method(m, h);
	cout << "\n\n-----------------------------\nFOR h = " << h << "\n";
	cout << "Method  |   t=0.0                  t=0.1                t=0.2                 t=0.3              t=0.4                  t=0.5\n";
	cout << "Exact:   ";
	for (int i=0; i<6; i++) {
		printf("(%.5f,%.5f)   ", y_ex[0][i], y_ex[1][i]);
	}
	cout << "\nEuler:   ";
	for (int i=0; i<6; i++) {
		printf("(%.5f,%.5f)   ", eul[0][i*k], eul[1][i*k]);
	}
	cout << "\nRev. Euler: ";
	for (int i=0; i<6; i++) {
		printf("(%.5f,%.5f)   ", reveul[0][i*k], reveul[1][i*k]);
	}
	cout << "\nAdams: ";
	for (int i=0; i<6; i++) {
		printf("(%.5f,%.5f)   ", adams[0][i*k], adams[1][i*k]);
	}
	return;
}

int main()
{
	vector<vd> matr {
  	{-125, 123.45},
	{123.45, -123}
  	};
  	double eps = 0.001, c, s, maxeival=0;
  	vd eigenvals = get_eivals(matr, eps);
  	cout << "Eigenval #1 = " << eigenvals[0] << "\n";
  	cout << "Eigencal #2 = " << eigenvals[1] << "\n";
  	maxeival = max_elem_vect(eigenvals);
  	cout << "|Max| = " << maxeival << "\n";
  	cout << "So Euler method is stable for h < " << 2.0/maxeival <<"\n";
  	cout << "Rev. Euler is stable for every h\n";
  	cout << "Adams method is stable for h < " << 6.0/maxeival <<"\n";
  	print_for_h(matr, 0.05);
  	cout <<"\n------------\n";
  	print_for_h(matr, 0.001);
  	return 0;
}
