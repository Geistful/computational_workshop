#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;

vector<vd> E_matr(int deg) {	//Единичная матрица
	vector<vd> e(deg, vd (deg, 0));
	for (int i=0; i<deg; i++) {
		e[i][i]=1;
	}
	return e;
}

vector<vd> D_matr (vector<vd> m) {	//Диагональная матрица
	vector<vd> d(m.size(), vd(m.size(), 0));
	for (int i=0; i<m.size(); i++) {
		d[i][i] = m[i][i];
	}
	return d;
} 

vector<vd> L_matr (vector<vd> m) {	//Левая(нижняя) матрица
	vector<vd> d(m.size(), vd(m.size(), 0));
	for (int i=0; i<m.size(); i++) {
		for (int j=0; j<i; j++) {
			d[i][j] = m[i][j];
		}
	}
	return d;
} 

vector<vd> R_matr (vector<vd> m) {	//Правая(верхняя) матрица
	vector<vd> d(m.size(), vd(m.size(), 0));
	for (int i=0; i<m.size(); i++) {
		for (int j=i+1; j<m.size(); j++) {
			d[i][j] = m[i][j];
		}
	}
	return d;
} 

double inf_norm (vector<vd> m) {	//Норма ||.||_inf для матриц
	double maxvalue=0;
	for (int i=0; i<m.size(); i++) {
		double colsum=0;
		for (int j=0; j<m[0].size(); j++) {
			colsum+=abs(m[i][j]);
		}
		maxvalue=max(maxvalue, colsum);
	}
	return maxvalue;
}
//--------------Для поиска собственных чисел--------------// 
//------используются следующие функции из задания №3------//
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

int sign (double temp) {
	if (temp>0) return 1;
	if (temp<0) return -1;
	return 0;
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
//-----------------------------------------//

vd get_eivals (vector<vd> m, double eps) {	//Получаем массив собственных чисел
	double c, s;
  	int p, q;
  	vector<vd> x = E_matr(m.size());
  	while (set_p_q_max(m, p, q)>=eps) {
  		c=get_c(m, p, q);
		s=get_s(m, p, q);
  		m = a_get_next(m, p, q, c, s);
  		x = x_get_next(x, p, q, c, s);
  	}
  	vd eigenval (m.size(), 0);
  	for (int i=0; i<m.size(); i++) {
  		eigenval[i]=m[i][i];
  	}
  	return eigenval;
}

double vect_inf_norm (vd v) {	//Соответствующая ||.||_inf для массива
	double maxv=0;
	for (int j=0; j<v.size(); j++) {
  		maxv=max(maxv, abs(v[j]));
  	}
  	return maxv;
}

//--------Некоторые функции-операции для матриц(векторов)
vector<vd> matr_plus_matr (vector<vd> m1, vector<vd> m2) { //Поэлементная сумма матриц
  vector<vd> mult(m1.size(), vd (m1.size(), 0));
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1.size(); i++) {
		mult[k][i] = m1[k][i]+m2[k][i];
	}
  }
  return mult;
}

vector<vd> matr_minus_matr (vector<vd> m1, vector<vd> m2) { //Соответственно разность
  vector<vd> mult(m1.size(), vd (m1.size(), 0));
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1.size(); i++) {
		mult[k][i] = m1[k][i]-m2[k][i];
	}
  }
  return mult;
}

vd vect_plus_vect (vd v1, vd v2) {	//Сумма массивов
	vd v3(v1.size(), 0);
	for (int i=0; i<v1.size(); i++) {
		v3[i]=v1[i]+v2[i];
	}
	return v3;
}

vd double_x_vect (double num, vd v) {	//Умножение скаляра на массив
	vd res(v.size(), 0);
	for (int i=0; i<v.size(); i++) {
		res[i]=num*v[i];
	}
	return res;
}

vector<vd> matr_x_matr (vector<vd> m1, vector<vd> m2) {	//Произведение матриц
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

vd matr_x_vect (vector<vd> m1, vd m2) { 	//Произведение матрицы и массива
  vd mult(m1.size(), 0);
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1[0].size(); i++) {
		mult[k] += m1[k][i]*m2[i];
	}
  }
  return mult;
}
//-------------------------//

//----Функции из первой задачи для решения методом Гаусса--------//
vector<vd> swap_row(vector<vd> &m, int p, int k) {	
	if (p==k) return m;
	double tmp;
	for (int i=0; i<m[0].size(); i++) {
		tmp = m[p][i];
		m[p][i] = m[k][i];
		m[k][i] = tmp;
	}
	return m;
}

vector<vd> swap_col(vector<vd> &m, int q, int k, vi &xorder) {
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

vi gauss_forw_elim (vector<vd> &newmatr) {
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

vd gauss_back_subs (vector<vd> newmatr, vi xorder) {
	vd x(newmatr.size(), 0);
  	x[newmatr.size()-1]=newmatr[newmatr.size()-1][newmatr[0].size()-1];
  	for (int i=newmatr.size()-2; i>=0; i--) {
  		x[i]=newmatr[i][newmatr[0].size()-1];
  		for (int j=i+1; j<newmatr[0].size()-1; j++) {
  			x[i]-=newmatr[i][j]*x[j];
  		}	
  	}
  	vd x_ordered(newmatr.size(), 0);
  	for (int i=0; i<newmatr.size(); i++) {
  		x_ordered[i]=x[xorder[i]];
  	}
  	return x_ordered;
}

vd gauss_method (vector<vd> &newmatr) {
	vi xorder = gauss_forw_elim(newmatr);
	vd x = gauss_back_subs(newmatr, xorder);
	return x;
}

vector<vd> unite_a_b (vector<vd> a, vd b) {
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

vector<vd> rev_diag (vector<vd> matr) { //Матрица, обратная диагональной (для матрицы D)
	vector<vd> rev (matr.size(), vd (matr.size(), 0));
	for (int i=0; i<matr.size(); i++) {
		for (int j=0; j<matr.size(); j++) {
			if (i==j) rev[j][i] = 1.0/matr[i][j];
			else rev[j][i]=0; 
		}
	}
	return rev;
}

void print_long_matr(vector<vd> m) {
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m[0].size(); j++) {
			if (j==m[0].size()-1)
				cout << "| " << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

void print_matr(vector<vd> m) {
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m.size(); j++) {
			if (j==m.size()-1)
				cout << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

void print_vect(vd v) {
	for (int i=0; i<v.size(); i++) {
		if (i==0) 
			printf("(%.16f, ", v[i]);
		else if (i==v.size()-1)
			printf("%.16f)\n", v[i]);
		else
			printf("%.16f, ", v[i]);
	}
}
//----------------------------//

vd zeid_method (vector<vd> h, vd g, vd x_k) {	//Метод Зейдля
	vd next(x_k.size(), 0);
	for (int i=0; i<h.size(); i++) {
		for (int j=0; j<h.size(); j++) {
			if (j<i)
				next[i] += h[i][j]*next[j];
			else
				next[i] += h[i][j]*x_k[j];
		}
		next[i]+=g[i];
	}
	return next;
}

vector<vd> num_x_matr (double l, vector<vd> m) {	//Поэлементное умножение числа на матрицу
	for (int i=0; i<m.size(); i++) {
		for (int j=0; j<m[0].size(); j++) {
			m[i][j]*=l;
		}
	}
	return m;
} 

double Aposteric_val (double Hd_inf, vd last_diff_vect) {	//Оценка апостериорная
	return (Hd_inf/(1-Hd_inf))*vect_inf_norm(last_diff_vect);
}

double Aprioral_val (double Hd_inf, int cnt, double x_inf, vd gd) {	//Оценка априорная
	return pow(Hd_inf, cnt)*x_inf+(pow(Hd_inf, cnt)/
  		  (1-Hd_inf))*vect_inf_norm(gd);
}

vd relax_method (vector<vd> h, vd g, double q, vd x_k) {	//Метод релаксации
	vd next(x_k.size(), 0);
	for (int i=0; i<h.size(); i++) {
		double sum_i = 0;
		for (int j=0; j<h.size(); j++) {
			if (j<i)
				sum_i += h[i][j]*next[j];
			else if (j==i)
				sum_i -= x_k[i];
			else
				sum_i += h[i][j]*x_k[j];
		}
		sum_i+=g[i];
		next[i]=x_k[i]+q*sum_i;
	}
	return next;
}

double vect_max_num_modulo (vd values) {	//Получение максимального по модулю числа из массива
	double maxval = 0;
	for (int i=0; i<values.size(); i++) {
		maxval=max(maxval, abs(values[i]));
	}
	return maxval;
}

double vect_metric (vd v1, vd v2) {		//Аналог "расстояния": ||v1-v2||_inf
	return vect_inf_norm(vect_plus_vect(v1, double_x_vect(-1.0, v2)));
}

int main()
{
	vector<vd> a {
  	{8.29381, 0.995516, -0.560617},
	{0.995516, 6.298198, 0.595772 },
	{-0.560617, 0.595772, 4.997407}
    };
    vd b {0.766522, 3.844422, 5.239231};
    double eps = 0.00000000001;
    vector<vd> A = unite_a_b(a, b);
    cout << "A|B:\n";
    print_long_matr(A);
  	vd x = gauss_method(A);
  	cout << "Solution:\n";
  	print_vect(x);
  	cout << "Matrix H_{D}:\n";
  	vector<vd> Ea=E_matr(a.size()), Da=D_matr(a),
  			   Drev=rev_diag(Da);
  	vector<vd> Hd = matr_minus_matr(Ea, matr_x_matr(Drev, a));
  	print_matr(Hd);
  	vd gd = matr_x_vect(Drev, b);
  	cout << "Vector g_{D}:\n";
  	print_vect(gd);
  	double Hd_inf = inf_norm(Hd);
  	cout << "Inf. norm of H_{D} = " << Hd_inf << "\n";
  	int k=0;
  	double x_inf = vect_inf_norm(x);
  	while (Aprioral_val(Hd_inf, k, x_inf, gd)>=eps) {
  		k++;
  	}
  	cout << "Aprioral num. of iterations:\n";
  	cout << k << "\n";
  	vd v1(a.size(), 0);
  	vd temp_prev(a.size(), 0);
  	int cnt = 0;
  	while (vect_metric(x, v1)>=eps) {
  		temp_prev = v1;
  		v1 = vect_plus_vect(matr_x_vect(Hd, v1), gd);
  		cnt++;
  	}
  	cout << "Actual num. of iterations:\n";
  	cout << cnt << "\n";
  	print_vect(v1);
	cout << "Factual error:\n";
	double fact = vect_metric(x, v1);
	cout << fact << "\n";
	cout << "For k = " << cnt << "\n";
	double apr = Aprioral_val(Hd_inf, cnt, x_inf, gd);
	vd last_diff_vect = vect_plus_vect(v1, double_x_vect(-1.0, temp_prev));
	double apo = Aposteric_val(Hd_inf, last_diff_vect);
	cout << "Aprioric error is " << apr << "\n";
	cout << "Diff (Aprioric-factual): " << apr-fact << "\n";
	cout << "Aposteric error is " << apo << "\n";
	cout << "Diff (aposteric-factual): " << apo-fact << "\n";
	//--------------Part from task#3 for eigenvalues-------//
	vd eigen = get_eivals(a, 0.001);
	print_vect(eigen);
	double max_eival=vect_max_num_modulo(eigen);
	cout << "Lusternick approximation = ";
	vd luster_vect = vect_plus_vect(temp_prev, double_x_vect(1.0/(1-max_eival), last_diff_vect));
	print_vect(luster_vect);
	cout << "Factual diff. for approximation = " << vect_metric(x, luster_vect) << "\n";
	cout << "Zeidel method:\n";
	vd x_first_zeidel(a.size(), 0);
	cnt=0;
	while (vect_metric(x, x_first_zeidel)>=eps) {
		temp_prev = x_first_zeidel;
  		x_first_zeidel = zeid_method(Hd, gd, x_first_zeidel);
  		cnt++;
  	}
  	cout << "Zeidel sol:\n";
  	print_vect(x_first_zeidel);
  	cout << "Factual diff. for approximation = " << vect_metric(x, x_first_zeidel) << "\n";
  	cout << "Number of interations for Zeidel:\n";
  	cout << cnt << "\n";
  	vd x_first_relax(a.size(), 0);
  	vd h_eigen = get_eivals(Hd, 0.001);
  	double h_radius = vect_max_num_modulo(h_eigen);
  	double q = 2.0/(1+sqrt(1-h_radius*h_radius));
	cnt=0;
  	while (vect_metric(x, x_first_relax)>=eps) {
  		x_first_relax = relax_method(Hd, gd, q, x_first_relax);
  		cnt++;
  	}
  	cout << "Relax sol:\n";
  	print_vect(x_first_relax);
  	cout << "Factual diff. for approximation = " << vect_metric(x, x_first_relax) << "\n";
  	cout << "Number of interations for Relax:\n";
  	cout << cnt << "\n";
	return 0;
}
