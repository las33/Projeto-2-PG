#include <iostream>
using std::cerr;
using std::endl;
#include <stdlib.h>
#include <GL/glut.h> 
#include <windows.h>
#include <cmath>
#include <gl/Gl.h>
#include <gl/Glu.h>
#include <vector>

//Leitura de arquivo
#include <fstream>
#include <string>

using namespace std;

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

class Point {
public:
	float x, y, z;
	void setxyz(float x2, float y2, float z2) { x = x2; y = y2; z = z2; }

};

class Point_2D {
public:
	float x, y;
	void setxyz(int x2, int y2) { x = x2; y = y2; }

};

class Triangulo {
public:
	int  v1, v2, v3;
	void setvs(int a, int b, int c) { v1 = a; v2 = b; v3 = c; }
};

bool eq(float a, float b) {
	if (abs(a - b)<(10e-9))
		return true;
	return false;
}

class Linha {
public:

	float a, b, c, d;
	Linha() :a(0), b(0), c(0), d(0) {}
	Linha(float a, float b, float c, float d) :a(a), b(b), c(c), d(d) {}
	Linha(const Linha &l) :a(l.a), b(l.b), c(l.c), d(l.d) {}
	Linha operator /(float x) const {
		if (eq(x, 0))
			x = 1.0f;
		return Linha(a / x, b / x, c / x, d / x);
	}
	Linha operator %(const Linha &x) const {
		return Linha(a - x.a*a, b - x.b*a, c - x.c*a, d - x.d*a);
	}
	Linha operator ^(const Linha &x) const {
		return Linha(a, b - x.b*b, c - x.c*b, d - x.d*b);
	}
	Linha operator +(const Linha &x) const {
		return Linha(a, b, c - x.c*c, d - x.d*c);
	}
	string to_string() {
		char res[150];
		sprintf(res, "(%.03f, %.03f, %.03f, %.03f)", a, b, c, d);
		return res;
	}
};

class Escalona {
public:

	Linha l1, l2, l3;
	Escalona() :l1(0, 0, 0, 0), l2(0, 0, 0, 0), l3(0, 0, 0, 0) {}
	Escalona(const Linha &l1, const Linha &l2, const Linha &l3) :l1(l1), l2(l2), l3(l3) {}
	Escalona(const Escalona &l) :l1(l.l1), l2(l.l2), l3(l.l3) {}

	Linha esc() {
		l1 = l1 / l1.a;
		l2 = l2%l1;
		l3 = l3%l1;
		l2 = l2 / l2.b;
		l1 = l1^l2;
		l3 = l3^l2;
		l3 = l3 / l3.c;
		l1 = l1 + l3;
		l2 = l2 + l3;
		return Linha(l1.d, l2.d, l3.d, 0);
	}
};

Point U, N, C, V, Pl, Ia, Il, Od;
vector<Point> points;
vector<Triangulo> triangulos;
vector<Point> normaisTriangulos;
vector<Point> normaisPoints;
vector<Point_2D> points_tela;

float z_buffer[SCREEN_WIDTH + 1][SCREEN_HEIGHT + 1];
float d, hx, hy, ka, kd, n, ks;
bool gr = false;

vector<string> split(string str) {
	vector<string> ret;
	while (!str.empty())
	{
		ret.push_back(str.substr(0, str.find(' ')));
		str.erase(0, str.find(' '));
		while (!(str.empty()) && (str[0] == ' ')) { str.erase(0, 1); }
	}

	return ret;
}

void calcular_normal_triangulos();
void coordenadasToVista();
void FindU();
void drawTriangle(Triangulo t, Point_2D p1, Point_2D p2, Point_2D p3);
void inicializa_z_buffer();
void drawLine(Triangulo t, Point_2D p1, Point_2D p2, Point_2D p3);
void calculoGeral(int xscan, int scanlineY, Triangulo t, bool line);

void myInit() {
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glOrtho(0.0, float(SCREEN_WIDTH), float(SCREEN_HEIGHT), 0.0, -5.0, 5.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//camera
	string vetor, numero;
	ifstream camera("calice2.cfg");//arquivo deve estar em \bin do projeto :)
	if (camera.is_open())
	{
		getline(camera, vetor);
		C.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(camera, vetor);
		N.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(camera, vetor);
		V.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(camera, numero);
		d = atof(split(numero)[0].c_str());
		hx = atof(split(numero)[1].c_str());
		hy = atof(split(numero)[2].c_str());

		camera.close();
	}

	//iluminacao
	ifstream iluminacao("iluminacao.txt");//arquivo deve estar em \bin do projeto :)
	if (iluminacao.is_open())
	{
		getline(iluminacao, vetor);
		Pl.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(iluminacao, numero);
		ka = atof(numero.c_str());

		getline(iluminacao, vetor);
		Ia.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(iluminacao, numero);
		kd = atof(split(numero)[0].c_str());

		getline(iluminacao, vetor);
		Od.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(iluminacao, numero);
		ks = atof(split(numero)[0].c_str());

		getline(iluminacao, vetor);
		Il.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));

		getline(iluminacao, numero);
		n = atof(split(numero)[0].c_str());

		iluminacao.close();
	}

	//objetos
	ifstream objetos("calice2.byu");//arquivo deve estar em \bin do projeto :)
	int pts, t;
	if (objetos.is_open())
	{
		getline(objetos, vetor);
		pts = atoi(split(vetor)[0].c_str());
		t = atoi(split(vetor)[1].c_str());

		Point ponto;
		ponto.x = 0.0;
		ponto.y = 0.0;
		ponto.z = 0.0;
		points.push_back(ponto);
		for (int i = 1; i <= pts; i++) {
			getline(objetos, vetor);
			ponto.setxyz(atof(split(vetor)[0].c_str()), atof(split(vetor)[1].c_str()), atof(split(vetor)[2].c_str()));
			points.push_back(ponto);
		}

		Triangulo triangulo;
		triangulo.v1 = 0;
		triangulo.v2 = 0;
		triangulo.v3 = 0;
		triangulos.push_back(triangulo);
		for (int i = 1; i <= t; i++) {
			getline(objetos, vetor);
			triangulo.setvs(atoi(split(vetor)[0].c_str()), atoi(split(vetor)[1].c_str()), atoi(split(vetor)[2].c_str()));
			triangulos.push_back(triangulo);
		}
		objetos.close();
	}
}

void start() {
	FindU();
	coordenadasToVista();
	calcular_normal_triangulos();
	inicializa_z_buffer();
	int r = 0, t = 0;
	for (int i = 1; i < triangulos.size(); i++) {
		drawTriangle(triangulos[i], points_tela[triangulos[i].v1], points_tela[triangulos[i].v2], points_tela[triangulos[i].v3]);
		t++;
	}
}

void drawDot(int x, int y) {
	glVertex2f(x, y);
}

Point subtracaoVetor(Point p1, Point p2) {
	p1.setxyz(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
	return p1;
}

Point produtoVetorEscalar(float k, Point p1) {
	p1.setxyz(p1.x * k, p1.y * k, p1.z * k);
	return p1;
}

float produtoInterno(Point p1, Point p2) {
	float valor;
	valor = p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
	return valor;
}

Point projecao(Point v, Point n) {
	float produtoInternoVN = produtoInterno(v, n);
	float produtoInternoN = produtoInterno(n, n);

	float valor = produtoInternoVN / produtoInternoN;

	n.x *= valor;
	n.y *= valor;
	n.z *= valor;

	return n;
}

Point normaliza(Point p) {
	float pInterno = produtoInterno(p, p);
	float norma = sqrt(pInterno);
	p.x /= norma;
	p.y /= norma;
	p.z /= norma;

	return p;
}

Point produtoVetorial(Point p1, Point p2) {
	float i = p1.y*p2.z - p1.z*p2.y;
	float j = p1.z*p2.x - p1.x*p2.z;
	float k = p1.x*p2.y - p1.y*p2.x;
	Point p;
	p.setxyz(i, j, k);

	return p;
}

void FindU() {
	N = normaliza(N);
	Point projec = projecao(V, N);

	V.x = V.x - projec.x;
	V.y = V.y - projec.y;
	V.z = V.z - projec.z;

	V = normaliza(V);
	U = produtoVetorial(N, V);
}

Point convertToVista(Point P_mundo) {
	P_mundo.x -= C.x;
	P_mundo.y -= C.y;
	P_mundo.z -= C.z;

	Point p_vista;
	p_vista.x = P_mundo.x *U.x + P_mundo.y*U.y + P_mundo.z*U.z;
	p_vista.y = P_mundo.x *V.x + P_mundo.y*V.y + P_mundo.z*V.z;
	p_vista.z = P_mundo.x *N.x + P_mundo.y*N.y + P_mundo.z*N.z;

	return p_vista;
}

Point_2D calcular_ponto_tela(Point P_vista) {
	Point_2D p_tela;

	float x = (d / hx)*(P_vista.x / P_vista.z);
	float y = (d / hy)*(P_vista.y / P_vista.z);

	p_tela.x = static_cast<int>(((x + 1) * (SCREEN_WIDTH / 2)));

	p_tela.y = static_cast<int>(((1 - y) * (SCREEN_HEIGHT / 2)));

	return p_tela;
}

void coordenadasToVista() { //calcula as coordenadas de vista e de tela de uma vez
	Point_2D p;
	points_tela.push_back(p);
	for (int i = 1; i < points.size(); i++) {
		points[i] = convertToVista(points[i]);
		points_tela.push_back(calcular_ponto_tela(points[i]));
	}

	Pl = convertToVista(Pl);
}

Point normal_triangulo(Point p1, Point p2, Point p3) { //calcula a normal do triangulo
	Point v1; v1.setxyz(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	Point v2; v2.setxyz(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
	Point retorna = produtoVetorial(v1, v2);

	return retorna;
}

void inicializaNormaisPoints() { //inicializa as normais dos pontos com (0,0,0)
	Point p;
	p.setxyz(0, 0, 0);
	normaisPoints.push_back(p);
	for (int i = 1; i <= points.size(); i++) {
		normaisPoints.push_back(p);
	}
}

void normalizaNormaisPoints() {
	for (int i = 1; i < normaisPoints.size(); i++) {
		normaisPoints[i] = normaliza(normaisPoints[i]);
	}
}

void calcular_normal_triangulos() { //calcular normais dos triangulos e pontos
	inicializaNormaisPoints(); //inicializa normais dos pontos)
	Point p;
	normaisTriangulos.push_back(p);
	for (int i = 1; i < triangulos.size(); i++) {
		Triangulo t = triangulos[i];//calcula normal dos triangulos
		Point x = normal_triangulo(points[t.v1], points[t.v2], points[t.v3]);

		normaisTriangulos.push_back(x);
		normaisTriangulos[i] = normaliza(normaisTriangulos[i]);

		//soma a normal do triangulo à cada normal dos seus vertices
		//v1
		normaisPoints[t.v1].x += normaisTriangulos[i].x;
		normaisPoints[t.v1].y += normaisTriangulos[i].y;
		normaisPoints[t.v1].z += normaisTriangulos[i].z;

		//v2
		normaisPoints[t.v2].x += normaisTriangulos[i].x;
		normaisPoints[t.v2].y += normaisTriangulos[i].y;
		normaisPoints[t.v2].z += normaisTriangulos[i].z;

		//v3
		normaisPoints[t.v3].x += normaisTriangulos[i].x;
		normaisPoints[t.v3].y += normaisTriangulos[i].y;
		normaisPoints[t.v3].z += normaisTriangulos[i].z;
	}

	normalizaNormaisPoints();
}

void inicializa_z_buffer() { //inicializa x_buffer com +infinito
	for (int i = 0; i < SCREEN_WIDTH + 1; i++) {
		for (int j = 0; j < SCREEN_HEIGHT + 1; j++) {
			z_buffer[i][j] = 99999999;
		}
	}
}

//Ponto 2d + os 3 vertices do triangulo em 2D
Point calcula_alfa_beta_gama(Point p, Point_2D v1, Point_2D v2, Point_2D v3) {
	float alfa, beta, gama;

	Linha l1 = Linha(v1.x, v2.x, v3.x, p.x);
	Linha l2 = Linha(v1.y, v2.y, v3.y, p.y);
	Linha l3 = Linha(1, 1, 1, 1);
	Linha esc = Escalona(l1, l2, l3).esc();
	alfa = esc.a;
	beta = esc.b;
	gama = esc.c;

	Point retorno; retorno.setxyz(alfa, beta, gama);

	return retorno;
}

Point calcula_alfa_beta(Point p, Point_2D v1, Point_2D v2) {
	Point retorno;

	if (eq(abs(v2.x - v1.x), 0)) {
		retorno.setxyz(1, 0, 0);
		return retorno;
	}

	float beta = (p.x - v1.x) / (v2.x - v1.x);
	float alfa = 1 - beta;
	retorno.setxyz(alfa, beta, 0);

	return retorno;
}

Point vetor_cor(Point p, Point cor) {
	Point retorno;

	retorno.x = p.x *cor.x;
	retorno.y = p.y *cor.y;
	retorno.z = p.z *cor.z;

	return retorno;
}

Point vetor_escalar(float k, Point v) {
	v.x *= k;
	v.y *= k;
	v.z *= k;

	return  v;
}

Point somaVetor(Point p1, Point p2) {
	p1.setxyz(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
	return p1;
}

Point phong(Point L, Point V, Point R, Point normal, float kdAUx, float ksAux) {
	Point  ia, Ib, Ic;

	ia = vetor_escalar(ka, Ia);

	float kdNL = kdAUx*produtoInterno(normal, L);
	Point vetor = vetor_cor(Od, Il);
	Ib = vetor_escalar(kdNL, vetor);

	float ksRV = ksAux*pow(produtoInterno(R, V), n);
	Ic = vetor_escalar(ksRV, Il);

	Point sum = somaVetor(Ib, ia);
	sum = somaVetor(sum, Ic);

	return sum;
}

void drawLine(Triangulo t, Point_2D p1, Point_2D p2, Point_2D p3) {
	int a = t.v1, b = t.v2, c = t.v3;
	/*ordenar os pontos internamente em relacao ao valor de X*/
	if (points_tela[a].x > points_tela[b].x)
		swap(a, b);
	if (points_tela[a].x > points_tela[c].x)
		swap(a, c);
	if (points_tela[b].x > points_tela[c].x)
		swap(b, c);

	Triangulo t1;
	t1.setvs(a, b, c);

	int x1 = points_tela[a].x;
	int sline = points_tela[a].y;

	int x2 = points_tela[b].x;

	while (x1 <= x2) {
		calculoGeral(x1, sline, t1, true);
		x1++;
	}
}

void calculoGeral(int xscan, int scanlineY, Triangulo t, bool line) {
	Point aux, pixel, p_3D, normal, v, l, r;
	float alfa, beta, gama;

	pixel.setxyz(xscan, scanlineY, 0);

	if (line == false) {
		aux = calcula_alfa_beta_gama(pixel, points_tela[t.v1], points_tela[t.v2], points_tela[t.v3]);
		alfa = aux.x;   beta = aux.y;   gama = aux.z;
		
		
		/*aproximacao do ponto 3D*/
		p_3D.setxyz(alfa * points[t.v1].x + beta * points[t.v2].x + gama * points[t.v3].x,
			alfa * points[t.v1].y + beta * points[t.v2].y + gama * points[t.v3].y,
			alfa * points[t.v1].z + beta * points[t.v2].z + gama * points[t.v3].z);
	}
	else {
		aux = calcula_alfa_beta(pixel, points_tela[t.v1], points_tela[t.v3]);

		alfa = aux.x;
		beta = aux.y;
		/*aproximacao do ponto 3D*/
		p_3D.setxyz(alfa * points[t.v1].x + beta  * points[t.v3].x,
			alfa * points[t.v1].y + beta  * points[t.v3].y,
			alfa * points[t.v1].z + beta  * points[t.v3].z);
	}

	int x = static_cast<int>(pixel.x + 0.5);
	int y = static_cast<int>(pixel.y + 0.5);
	float z = p_3D.z;

	/*consulta ao z-buffer*/
	if ((z < z_buffer[x][y]) && (x < (SCREEN_WIDTH + 1)) && y < (SCREEN_HEIGHT + 1)) {
		z_buffer[x][y] = z;
		if (line == false) {
			normal.setxyz(alfa * normaisPoints[t.v1].x + beta * normaisPoints[t.v2].x + gama * normaisPoints[t.v3].x,
				alfa * normaisPoints[t.v1].y + beta * normaisPoints[t.v2].y + gama * normaisPoints[t.v3].y,
				alfa * normaisPoints[t.v1].z + beta * normaisPoints[t.v2].z + gama * normaisPoints[t.v3].z);
		}
		else {
			normal.setxyz(alfa * normaisPoints[t.v1].x + beta  * normaisPoints[t.v3].x,
				alfa * normaisPoints[t.v1].y + beta  * normaisPoints[t.v3].y,
				alfa * normaisPoints[t.v1].z + beta  * normaisPoints[t.v3].z);
		}

		v.setxyz(0 - p_3D.x, 0 - p_3D.y, 0 - p_3D.z);
		l.setxyz(Pl.x - p_3D.x, Pl.y - p_3D.y, Pl.z - p_3D.z);

		normal = normaliza(normal);
		v = normaliza(v);
		l = normaliza(l);

		float kdAux = kd, ksAux = ks;
		if (produtoInterno(v, normal) < 0) normal.setxyz(0 - normal.x, 0 - normal.y, 0 - normal.z);

		if (produtoInterno(normal, l) < 0) {
			kdAux = 0;
			ksAux = 0;
		}

		/*Calcular o vetor R = (2(N . L)N) - L. E normalizar R.*/
		float a = produtoInterno(normal, l);
		Point b = produtoVetorEscalar(2.0, normal);
		b = produtoVetorEscalar(a, b);
		r = subtracaoVetor(b, l);
		r = normaliza(r);

		if ((alfa < 0.1) || (beta < 0.1) || (gama < 0.1)) {
			Od.setxyz(0.9, 0.1, 0.1);
		}
		else if (((alfa < 0.2) && (alfa > 0.1)) || ((beta < 0.2) && (beta > 0.1)) || ((gama < 0.2) && (gama > 0.1))) {
			Od.setxyz(0.1, 0.9, 0.1);
		}
		else {
			Od.setxyz(0.1, 0.1, 0.9);
		}


		if (produtoInterno(r, v) < 0) ksAux = 0;
		Point cor = phong(l, v, r, normal, kdAux, ksAux);

		float red = cor.x, green = cor.y, blue = cor.z;

		if (red > 255) red = 255.0;
		if (green > 255) green = 255.0;
		if (blue > 255) blue = 255.0;

		red /= 255.0; green /= 255.0;  blue /= 255.0;

		glColor3f(red, green, blue);
		drawDot(pixel.x, pixel.y);
	}
}

Point calculoGeral2(int xscan, int scanlineY, Triangulo t, bool line) {
	Point aux, pixel, p_3D, normal, v, l, r;
	float alfa, beta, gama;

	pixel.setxyz(xscan, scanlineY, 0);

	if (line == false) {
		aux = calcula_alfa_beta_gama(pixel, points_tela[t.v1], points_tela[t.v2], points_tela[t.v3]);
		alfa = aux.x;   beta = aux.y;   gama = aux.z;


		/*aproximacao do ponto 3D*/
		p_3D.setxyz(alfa * points[t.v1].x + beta * points[t.v2].x + gama * points[t.v3].x,
			alfa * points[t.v1].y + beta * points[t.v2].y + gama * points[t.v3].y,
			alfa * points[t.v1].z + beta * points[t.v2].z + gama * points[t.v3].z);
	}
	else {
		aux = calcula_alfa_beta(pixel, points_tela[t.v1], points_tela[t.v3]);

		alfa = aux.x;
		beta = aux.y;

		p_3D.setxyz(alfa * points[t.v1].x + beta  * points[t.v3].x,
			alfa * points[t.v1].y + beta  * points[t.v3].y,
			alfa * points[t.v1].z + beta  * points[t.v3].z);
	}

	int x = static_cast<int>(pixel.x + 0.5);
	int y = static_cast<int>(pixel.y + 0.5);
	float z = p_3D.z;
	Point cor;

	if ((x < (SCREEN_WIDTH + 1)) && y < (SCREEN_HEIGHT + 1)) {
		if (line == false) {
			normal.setxyz(alfa * normaisPoints[t.v1].x + beta * normaisPoints[t.v2].x + gama * normaisPoints[t.v3].x,
				alfa * normaisPoints[t.v1].y + beta * normaisPoints[t.v2].y + gama * normaisPoints[t.v3].y,
				alfa * normaisPoints[t.v1].z + beta * normaisPoints[t.v2].z + gama * normaisPoints[t.v3].z);
		}
		else {
			normal.setxyz(alfa * normaisPoints[t.v1].x + beta  * normaisPoints[t.v3].x,
				alfa * normaisPoints[t.v1].y + beta  * normaisPoints[t.v3].y,
				alfa * normaisPoints[t.v1].z + beta  * normaisPoints[t.v3].z);
		}

		v.setxyz(0 - p_3D.x, 0 - p_3D.y, 0 - p_3D.z);
		l.setxyz(Pl.x - p_3D.x, Pl.y - p_3D.y, Pl.z - p_3D.z);

		normal = normaliza(normal);
		v = normaliza(v);
		l = normaliza(l);

		float kdAux = kd, ksAux = ks;
		if (produtoInterno(v, normal) < 0) normal.setxyz(0 - normal.x, 0 - normal.y, 0 - normal.z);

		if (produtoInterno(normal, l) < 0) {
			kdAux = 0;
			ksAux = 0;
		}

		/*Calcular o vetor R = (2(N . L)N) - L. E normalizar R.*/
		float a = produtoInterno(normal, l);
		Point b = produtoVetorEscalar(2.0, normal);
		b = produtoVetorEscalar(a, b);
		r = subtracaoVetor(b, l);
		r = normaliza(r);

		if (produtoInterno(r, v) < 0) ksAux = 0;

		cor = phong(l, v, r, normal, kdAux, ksAux);

		float red = cor.x, green = cor.y, blue = cor.z;

		if (red > 255) red = 255.0;
		if (green > 255) green = 255.0;
		if (blue > 255) blue = 255.0;

		red /= 255.0; green /= 255.0;  blue /= 255.0;

		cor.setxyz(red, green, blue);
		return cor;
	}
	cor.setxyz(0, 0, 0);
	return cor;
}

void fillBottomFlatTriangle(Triangulo t, Point_2D v1, Point_2D v2, Point_2D v3)
{
	float invslope1 = (v2.x - v1.x) / (v2.y - v1.y);
	float invslope2 = (v3.x - v1.x) / (v3.y - v1.y);

	if (invslope1 > invslope2)   swap(invslope1, invslope2);
	float curx1 = v1.x;
	float curx2 = v1.x;

	Point cor_v1, cor_v2, cor_v3, bar, ponto, cor_ponto;
	//cores dos vertices para gouraud
	if (gr) {
		cor_v1 = calculoGeral2(v1.x, v1.y, t, false);
		cor_v2 = calculoGeral2(v2.x, v2.y, t, false);
		cor_v3 = calculoGeral2(v3.x, v3.y, t, false);
	}
	Point p_3D;
	for (int scanlineY = v1.y; scanlineY <= v2.y; scanlineY++)
	{
		for (int xscan = curx1; xscan <= curx2; xscan++) {
			if (gr) {
				ponto.setxyz(xscan, scanlineY, 1);
				bar = calcula_alfa_beta_gama(ponto, v1, v2, v3);

				p_3D.setxyz(bar.x * points[t.v1].x + bar.y * points[t.v2].x + bar.z * points[t.v3].x,
					bar.x * points[t.v1].y + bar.y * points[t.v2].y + bar.z * points[t.v3].y,
					bar.x * points[t.v1].z + bar.y * points[t.v2].z + bar.z * points[t.v3].z);

				float z = p_3D.z;

				if ((z < z_buffer[xscan][scanlineY]) && (xscan < (SCREEN_WIDTH + 1)) && scanlineY < (SCREEN_HEIGHT + 1)) {
					z_buffer[xscan][scanlineY] = z;

					//cor do pixel por combinacao baricentrica das cores dos vertices do triangulo
					cor_ponto.setxyz(bar.x * cor_v1.x + bar.y * cor_v2.x + bar.z * cor_v3.x,
						bar.x * cor_v1.y + bar.y * cor_v2.y + bar.z * cor_v3.y,
						bar.x * cor_v1.z + bar.y * cor_v2.z + bar.z * cor_v3.z);

					//pintar o pixel
					float red = cor_ponto.x, green = cor_ponto.y, blue = cor_ponto.z;

					glColor3f(red, green, blue);
					drawDot(xscan, scanlineY);
				}
			}
			else {
				calculoGeral(xscan, scanlineY, t, false);
			}
		}
		curx1 += invslope1;
		curx2 += invslope2;
	}
}

void fillTopFlatTriangle(Triangulo t, Point_2D v1, Point_2D v2, Point_2D v3)
{
	float invslope1 = (v3.x - v1.x) / (v3.y - v1.y);
	float invslope2 = (v3.x - v2.x) / (v3.y - v2.y);

	if (invslope2 > invslope1)   swap(invslope1, invslope2);
	float curx1 = v3.x;
	float curx2 = v3.x;

	Point cor_v1, cor_v2, cor_v3, bar, ponto, cor_ponto, Ia, Ib, Ip;
	//cores dos vertices para gouraud
	if (gr) {
		cor_v1 = calculoGeral2(v1.x, v1.y, t, false);
		cor_v2 = calculoGeral2(v2.x, v2.y, t, false);
		cor_v3 = calculoGeral2(v3.x, v3.y, t, false);
	}
	Point p_3D;
	for (int scanlineY = v3.y; scanlineY >= v1.y; scanlineY--)
	{
		for (int xscan = curx1; xscan <= curx2; xscan++) {
			if (gr) {
				ponto.setxyz(xscan, scanlineY, 1);
				bar = calcula_alfa_beta_gama(ponto, v1, v2, v3);

				p_3D.setxyz(bar.x * points[t.v1].x + bar.y * points[t.v2].x + bar.z * points[t.v3].x,
					bar.x * points[t.v1].y + bar.y * points[t.v2].y + bar.z * points[t.v3].y,
					bar.x * points[t.v1].z + bar.y * points[t.v2].z + bar.z * points[t.v3].z);

				float z = p_3D.z;

				if ((z < z_buffer[xscan][scanlineY]) && (xscan < (SCREEN_WIDTH + 1)) && scanlineY < (SCREEN_HEIGHT + 1)) {
					z_buffer[xscan][scanlineY] = z;

					//cor do pixel por combinacao baricentrica das cores dos vertices do triangulo
					cor_ponto.setxyz(bar.x * cor_v1.x + bar.y * cor_v2.x + bar.z * cor_v3.x,
						bar.x * cor_v1.y + bar.y * cor_v2.y + bar.z * cor_v3.y,
						bar.x * cor_v1.z + bar.y * cor_v2.z + bar.z * cor_v3.z);

					//pintar o pixel
					float red = cor_ponto.x, green = cor_ponto.y, blue = cor_ponto.z;

					glColor3f(red, green, blue);
					drawDot(xscan, scanlineY);
				}
			}
			else {
				calculoGeral(xscan, scanlineY, t, false);
			}
		}
		curx1 -= invslope1;
		curx2 -= invslope2;
	}
}


void drawTriangle(Triangulo t, Point_2D v1, Point_2D v2, Point_2D v3) {
	Point_2D a, b, c;

	int va = t.v1, vb = t.v2, vc = t.v3;
	/*ordenar os pontos internamente em relacao ao valor de X*/
	if (points_tela[va].y > points_tela[vb].y)
		swap(va, vb);
	if (points_tela[va].y > points_tela[vc].y)
		swap(va, vc);
	if (points_tela[vb].y > points_tela[vc].y)
		swap(vb, vc);

	Triangulo t1;
	t1.setvs(va, vb, vc);

	a = points_tela[va];
	b = points_tela[vb];
	c = points_tela[vc];

	t = t1;

	/* verifica o caso base para bottomTriangle */
	if (b.y == c.y) fillBottomFlatTriangle(t, a, b, c);

	/* verifica o caso base para topTriangle */
	else if (a.y == b.y) fillTopFlatTriangle(t, a, b, c);

	else {
		/*caso geral - dividir o triangulo em dois e aplicar ao bottom e top triangle*/
		Point_2D d;
		//printf("fillTopFlatTriangle && fillBottomFlatTriangle\n");
		d.setxyz((a.x + ((float)(b.y - a.y) / (float)(c.y - a.y)) * (c.x - a.x)), b.y);/*ponto de intermedio para divisao em dois triangulos*/
		fillBottomFlatTriangle(t, a, b, d);
		fillTopFlatTriangle(t, b, d, c);
	}
}

void reshape(int l, int a) {
	glutReshapeWindow(SCREEN_WIDTH, SCREEN_HEIGHT);
}

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(2.0);
	glBegin(GL_POINTS);

	start();

	glEnd();
	glFlush();
	glutSwapBuffers();
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Projeto 2");
	glClearColor(1.0, 1.0, 1.0, 0);

	glutDisplayFunc(myDisplay);
	glutReshapeFunc(reshape);

	myInit();
	glutMainLoop();

	return 0;
}