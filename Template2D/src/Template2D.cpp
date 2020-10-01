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

// Read files
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

class Triangle {
public:
	int  v1, v2, v3;
	void setvs(int a, int b, int c) { v1 = a; v2 = b; v3 = c; }
};

bool eq(float a, float b) {
	if (abs(a - b) < (10e-9))
		return true;
	return false;
}

class Line {
public:

	float a, b, c, d;
	Line() :a(0), b(0), c(0), d(0) {}
	Line(float a, float b, float c, float d) :a(a), b(b), c(c), d(d) {}
	Line(const Line &l) :a(l.a), b(l.b), c(l.c), d(l.d) {}
	Line operator /(float x) const {
		if (eq(x, 0))
			x = 1.0f;
		return Line(a / x, b / x, c / x, d / x);
	}
	Line operator %(const Line &x) const {
		return Line(a - (x.a*a), b - (x.b*a), c - (x.c*a), d - (x.d*a));
	}
	Line operator ^(const Line &x) const {
		return Line(a, b - (x.b*b), c - (x.c*b), d - (x.d*b));
	}
	Line operator +(const Line &x) const {
		return Line(a, b, c - (x.c*c), d - (x.d*c));
	}
	string to_string() {
		char res[150];
		sprintf(res, "(%.03f, %.03f, %.03f, %.03f)", a, b, c, d);
		return res;
	}
};

class RowReduction {
public:

	Line l1, l2, l3;
	RowReduction() :l1(0, 0, 0, 0), l2(0, 0, 0, 0), l3(0, 0, 0, 0) {}
	RowReduction(const Line &l1, const Line &l2, const Line &l3) :l1(l1), l2(l2), l3(l3) {}
	RowReduction(const RowReduction &l) :l1(l.l1), l2(l.l2), l3(l.l3) {}

	Line reduction() {
		l1 = l1 / l1.a;
		l2 = l2%l1;
		l3 = l3%l1;
		l2 = l2 / l2.b;
		l1 = l1^l2;
		l3 = l3^l2;
		l3 = l3 / l3.c;
		l1 = l1 + l3;
		l2 = l2 + l3;
		return Line(l1.d, l2.d, l3.d, 0);
	}
};

Point U, N, C, V, Pl, Ia, Il, Od;
vector<Point> points;
vector<Triangle> triangles;
vector<Point> normalTriangles;
vector<Point> normalPoints;
vector<Point_2D> screenPoints;

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

void calculateTriangleNormalVector();
void calculateScreenAndViewCoordinates();
void FindU();
void drawTriangle(Triangle t, Point_2D p1, Point_2D p2, Point_2D p3);
void initializeZBuffer();
void drawLine(Triangle t, Point_2D p1, Point_2D p2, Point_2D p3);
void ordinaryCalculation(int xscan, int scanlineY, Triangle t, bool line);

void myInit() {
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glOrtho(0.0, float(SCREEN_WIDTH), float(SCREEN_HEIGHT), 0.0, -5.0, 5.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//camera
	string myVector, myNumber;
	ifstream camera("calice2.cfg"); // this file should be in \bin folder of the project :)
	if (camera.is_open())
	{
		getline(camera, myVector);
		C.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(camera, myVector);
		N.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(camera, myVector);
		V.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(camera, myNumber);
		d = atof(split(myNumber)[0].c_str());
		hx = atof(split(myNumber)[1].c_str());
		hy = atof(split(myNumber)[2].c_str());

		camera.close();
	}

	//illumination
	ifstream illumination("illumination.txt"); // this file should be in \bin folder of the project :)
	if (illumination.is_open())
	{
		getline(illumination, myVector);
		Pl.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(illumination, myNumber);
		ka = atof(myNumber.c_str());

		getline(illumination, myVector);
		Ia.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(illumination, myNumber);
		kd = atof(split(myNumber)[0].c_str());

		getline(illumination, myVector);
		Od.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(illumination, myNumber);
		ks = atof(split(myNumber)[0].c_str());

		getline(illumination, myVector);
		Il.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));

		getline(illumination, myNumber);
		n = atof(split(myNumber)[0].c_str());

		illumination.close();
	}

	//objects
	ifstream objects("calice2.byu"); // this file should be in \bin folder of the project :)
	int pts, t;
	if (objects.is_open())
	{
		getline(objects, myVector);
		pts = atoi(split(myVector)[0].c_str());
		t = atoi(split(myVector)[1].c_str());

		Point myPoint;
		myPoint.x = 0.0;
		myPoint.y = 0.0;
		myPoint.z = 0.0;
		points.push_back(myPoint);
		for (int i = 1; i <= pts; i++) {
			getline(objects, myVector);
			myPoint.setxyz(atof(split(myVector)[0].c_str()), atof(split(myVector)[1].c_str()), atof(split(myVector)[2].c_str()));
			points.push_back(myPoint);
		}

		Triangle myTriangle;
		myTriangle.v1 = 0;
		myTriangle.v2 = 0;
		myTriangle.v3 = 0;
		triangles.push_back(myTriangle);
		for (int i = 1; i <= t; i++) {
			getline(objects, myVector);
			myTriangle.setvs(atoi(split(myVector)[0].c_str()), atoi(split(myVector)[1].c_str()), atoi(split(myVector)[2].c_str()));
			triangles.push_back(myTriangle);
		}
		objects.close();
	}
}

void start() {
	FindU();
	calculateScreenAndViewCoordinates();
	calculateTriangleNormalVector();
	initializeZBuffer();
	int r = 0, t = 0;
	for (int i = 1; i < triangles.size(); i++) {
		drawTriangle(triangles[i], screenPoints[triangles[i].v1], screenPoints[triangles[i].v2], screenPoints[triangles[i].v3]);
		t++;
	}
}

void drawDot(int x, int y) {
	glVertex2f(x, y);
}

Point subtractVectors(Point p1, Point p2) {
	p1.setxyz(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
	return p1;
}

Point escalarVectorProduct(float k, Point p1) {
	p1.setxyz(p1.x * k, p1.y * k, p1.z * k);
	return p1;
}

float innerProduct(Point p1, Point p2) {
	float value;
	value = p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
	return value;
}

Point projection(Point v, Point n) {
	float innerProductVN = innerProduct(v, n);
	float innerProductN = innerProduct(n, n);

	float value = innerProductVN / innerProductN;

	n.x *= value;
	n.y *= value;
	n.z *= value;

	return n;
}

Point normalize(Point p) {
	float innerP = innerProduct(p, p);
	float norm = sqrt(innerP);
	p.x /= norm;
	p.y /= norm;
	p.z /= norm;

	return p;
}

Point vectorProduct(Point p1, Point p2) {
	float i = p1.y*p2.z - p1.z*p2.y;
	float j = p1.z*p2.x - p1.x*p2.z;
	float k = p1.x*p2.y - p1.y*p2.x;
	Point p;
	p.setxyz(i, j, k);

	return p;
}

void FindU() {
	N = normalize(N);
	Point projec = projection(V, N);

	V.x = V.x - projec.x;
	V.y = V.y - projec.y;
	V.z = V.z - projec.z;

	V = normalize(V);
	U = vectorProduct(N, V);
}

Point convertToView(Point P_world) {
	P_world.x -= C.x;
	P_world.y -= C.y;
	P_world.z -= C.z;

	Point p_view;
	p_view.x = P_world.x *U.x + P_world.y*U.y + P_world.z*U.z;
	p_view.y = P_world.x *V.x + P_world.y*V.y + P_world.z*V.z;
	p_view.z = P_world.x *N.x + P_world.y*N.y + P_world.z*N.z;

	return p_view;
}

Point_2D calculateScreenPoint(Point P_view) {
	Point_2D p_screen;

	float x = (d / hx)*(P_view.x / P_view.z);
	float y = (d / hy)*(P_view.y / P_view.z);

	p_screen.x = static_cast<int>(((x + 1) * (SCREEN_WIDTH / 2)));

	p_screen.y = static_cast<int>(((1 - y) * (SCREEN_HEIGHT / 2)));

	return p_screen;
}

void calculateScreenAndViewCoordinates() { // calculate screen and view coordinates in one shot
	Point_2D p;
	screenPoints.push_back(p);
	for (int i = 1; i < points.size(); i++) {
		points[i] = convertToView(points[i]);
		screenPoints.push_back(calculateScreenPoint(points[i]));
	}

	Pl = convertToView(Pl);
}

Point triangleNormal(Point p1, Point p2, Point p3) { // triangle normal calculation
	Point v1; v1.setxyz(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	Point v2; v2.setxyz(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
	Point product = vectorProduct(v1, v2);

	return product;
}

void initializeNormalPoints() { // normal point initialization in (0,0,0)
	Point p;
	p.setxyz(0, 0, 0);
	normalPoints.push_back(p);
	for (int i = 1; i <= points.size(); i++) {
		normalPoints.push_back(p);
	}
}

void normalizeNormalPoints() {
	for (int i = 1; i < normalPoints.size(); i++) {
		normalPoints[i] = normalize(normalPoints[i]);
	}
}

void calculateTriangleNormalVector() { // points and triangles normal calculation
	initializeNormalPoints(); // normal points initialization
	Point p;
	normalTriangles.push_back(p);
	for (int i = 1; i < triangles.size(); i++) {
		Triangle t = triangles[i];// normal triangles calculation
		Point x = triangleNormal(points[t.v1], points[t.v2], points[t.v3]);

		normalTriangles.push_back(x);
		normalTriangles[i] = normalize(normalTriangles[i]);

		// sum triangle normal and each one of its vertex
		//v1
		normalPoints[t.v1].x += normalTriangles[i].x;
		normalPoints[t.v1].y += normalTriangles[i].y;
		normalPoints[t.v1].z += normalTriangles[i].z;

		//v2
		normalPoints[t.v2].x += normalTriangles[i].x;
		normalPoints[t.v2].y += normalTriangles[i].y;
		normalPoints[t.v2].z += normalTriangles[i].z;

		//v3
		normalPoints[t.v3].x += normalTriangles[i].x;
		normalPoints[t.v3].y += normalTriangles[i].y;
		normalPoints[t.v3].z += normalTriangles[i].z;
	}

	normalizeNormalPoints();
}

void initializeZBuffer() { // +infinity x_buffer initialization
	for (int i = 0; i < SCREEN_WIDTH + 1; i++) {
		for (int j = 0; j < SCREEN_HEIGHT + 1; j++) {
			z_buffer[i][j] = 99999999;
		}
	}
}

// 2D Point and all 3 triangle vertex in 2D
Point calculateAlphaBetaGamma(Point p, Point_2D v1, Point_2D v2, Point_2D v3) {
	float alpha, beta, gamma;

	Line l1 = Line(v1.x, v2.x, v3.x, p.x);
	Line l2 = Line(v1.y, v2.y, v3.y, p.y);
	Line l3 = Line(1, 1, 1, 1);
	Line esc = RowReduction(l1, l2, l3).reduction();
	alpha = esc.a;
	beta = esc.b;
	gamma = esc.c;

	Point returnPoint;
	returnPoint.setxyz(alpha, beta, gamma);

	return retorno;
}

Point calculateAlphaBeta(Point p, Point_2D v1, Point_2D v2) {
	Point returnPoint;

	if (eq(abs(v2.x - v1.x), 0)) {
		returnPoint.setxyz(1, 0, 0);
		return returnPoint;
	}

	float beta = (p.x - v1.x) / (v2.x - v1.x);
	float alpha = 1 - beta;
	returnPoint.setxyz(alpha, beta, 0);

	return returnPoint;
}

Point colorVector(Point p, Point color) {
	Point returnPoint;

	returnPoint.x = p.x *color.x;
	returnPoint.y = p.y *color.y;
	returnPoint.z = p.z *color.z;

	return returnPoint;
}

Point escalarVector(float k, Point v) {
	v.x *= k;
	v.y *= k;
	v.z *= k;

	return  v;
}

Point sumVector(Point p1, Point p2) {
	p1.setxyz(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
	return p1;
}

Point phong(Point L, Point V, Point R, Point normal, float kdAUx, float ksAux) {
	Point  ia, Ib, Ic;

	ia = escalarVector(ka, Ia);

	float kdNL = kdAUx*innerProduct(normal, L);
	Point myVector = colorVector(Od, Il);
	Ib = escalarVector(kdNL, myVector);

	float ksRV = ksAux*pow(innerProduct(R, V), n);
	Ic = escalarVector(ksRV, Il);

	Point sum = sumVector(Ib, ia);
	sum = sumVector(sum, Ic);

	return sum;
}

void drawLine(Triangle t, Point_2D p1, Point_2D p2, Point_2D p3) {
	int a = t.v1, b = t.v2, c = t.v3;
	/*order the points internally in relation to the value of X*/
	if (screenPoints[a].x > screenPoints[b].x)
		swap(a, b);
	if (screenPoints[a].x > screenPoints[c].x)
		swap(a, c);
	if (screenPoints[b].x > screenPoints[c].x)
		swap(b, c);

	Triangle t1;
	t1.setvs(a, b, c);

	int x1 = screenPoints[a].x;
	int sline = screenPoints[a].y;

	int x2 = screenPoints[b].x;

	while (x1 <= x2) {
		ordinaryCalculation(x1, sline, t1, true);
		x1++;
	}
}

void ordinaryCalculation(int xscan, int scanlineY, Triangle t, bool line) {
	Point aux, pixel, p_3D, normal, v, l, r;
	float alpha, beta, gamma;

	pixel.setxyz(xscan, scanlineY, 0);

	if (line == false) {
		aux = calculateAlphaBetaGamma(pixel, screenPoints[t.v1], screenPoints[t.v2], screenPoints[t.v3]);
		alpha = aux.x;   beta = aux.y;   gamma = aux.z;
		
		
		/*3D point approximation*/
		p_3D.setxyz(alpha * points[t.v1].x + beta * points[t.v2].x + gamma * points[t.v3].x,
			alpha * points[t.v1].y + beta * points[t.v2].y + gamma * points[t.v3].y,
			alpha * points[t.v1].z + beta * points[t.v2].z + gamma * points[t.v3].z);
	}
	else {
		aux = calculateAlphaBeta(pixel, screenPoints[t.v1], screenPoints[t.v3]);

		alpha = aux.x;
		beta = aux.y;
		/*3D point approximation*/
		p_3D.setxyz(alpha * points[t.v1].x + beta  * points[t.v3].x,
			alpha * points[t.v1].y + beta  * points[t.v3].y,
			alpha * points[t.v1].z + beta  * points[t.v3].z);
	}

	int x = static_cast<int>(pixel.x + 0.5);
	int y = static_cast<int>(pixel.y + 0.5);
	float z = p_3D.z;

	/*check z-buffer*/
	if ((z < z_buffer[x][y]) && (x < (SCREEN_WIDTH + 1)) && y < (SCREEN_HEIGHT + 1)) {
		z_buffer[x][y] = z;
		if (line == false) {
			normal.setxyz(alpha * normalPoints[t.v1].x + beta * normalPoints[t.v2].x + gamma * normalPoints[t.v3].x,
				alpha * normalPoints[t.v1].y + beta * normalPoints[t.v2].y + gamma * normalPoints[t.v3].y,
				alpha * normalPoints[t.v1].z + beta * normalPoints[t.v2].z + gamma * normalPoints[t.v3].z);
		}
		else {
			normal.setxyz(alpha * normalPoints[t.v1].x + beta  * normalPoints[t.v3].x,
				alpha * normalPoints[t.v1].y + beta  * normalPoints[t.v3].y,
				alpha * normalPoints[t.v1].z + beta  * normalPoints[t.v3].z);
		}

		v.setxyz(0 - p_3D.x, 0 - p_3D.y, 0 - p_3D.z);
		l.setxyz(Pl.x - p_3D.x, Pl.y - p_3D.y, Pl.z - p_3D.z);

		normal = normalize(normal);
		v = normalize(v);
		l = normalize(l);

		float kdAux = kd, ksAux = ks;
		if (innerProduct(v, normal) < 0) normal.setxyz(0 - normal.x, 0 - normal.y, 0 - normal.z);

		if (innerProduct(normal, l) < 0) {
			kdAux = 0;
			ksAux = 0;
		}

		/*Calculate the vector R = (2(N . L)N) - L. And normalize R.*/
		float a = innerProduct(normal, l);
		Point b = escalarVectorProduct(2.0, normal);
		b = escalarVectorProduct(a, b);
		r = subtractVectors(b, l);
		r = normalize(r);

		if ((alpha < 0.1) || (beta < 0.1) || (gamma < 0.1)) {
			Od.setxyz(0.9, 0.1, 0.1);
		}
		else if (((alpha < 0.2) && (alpha > 0.1)) || ((beta < 0.2) && (beta > 0.1)) || ((gamma < 0.2) && (gamma > 0.1))) {
			Od.setxyz(0.1, 0.9, 0.1);
		}
		else {
			Od.setxyz(0.1, 0.1, 0.9);
		}


		if (innerProduct(r, v) < 0) ksAux = 0;
		Point color = phong(l, v, r, normal, kdAux, ksAux);

		float red = color.x, green = color.y, blue = color.z;

		if (red > 255) red = 255.0;
		if (green > 255) green = 255.0;
		if (blue > 255) blue = 255.0;

		red /= 255.0; green /= 255.0;  blue /= 255.0;

		glColor3f(red, green, blue);
		drawDot(pixel.x, pixel.y);
	}
}

Point ordinaryCalculation2(int xscan, int scanlineY, Triangle t, bool line) {
	Point aux, pixel, p_3D, normal, v, l, r;
	float alpha, beta, gamma;

	pixel.setxyz(xscan, scanlineY, 0);

	if (line == false) {
		aux = calculateAlphaBetaGamma(pixel, screenPoints[t.v1], screenPoints[t.v2], screenPoints[t.v3]);
		alpha = aux.x;   beta = aux.y;   gamma = aux.z;


		/*3D point approximation*/
		p_3D.setxyz(alpha * points[t.v1].x + beta * points[t.v2].x + gamma * points[t.v3].x,
			alpha * points[t.v1].y + beta * points[t.v2].y + gamma * points[t.v3].y,
			alpha * points[t.v1].z + beta * points[t.v2].z + gamma * points[t.v3].z);
	}
	else {
		aux = calculateAlphaBeta(pixel, screenPoints[t.v1], screenPoints[t.v3]);

		alpha = aux.x;
		beta = aux.y;

		p_3D.setxyz(alpha * points[t.v1].x + beta  * points[t.v3].x,
			alpha * points[t.v1].y + beta  * points[t.v3].y,
			alpha * points[t.v1].z + beta  * points[t.v3].z);
	}

	int x = static_cast<int>(pixel.x + 0.5);
	int y = static_cast<int>(pixel.y + 0.5);
	float z = p_3D.z;
	Point color;

	if ((x < (SCREEN_WIDTH + 1)) && y < (SCREEN_HEIGHT + 1)) {
		if (line == false) {
			normal.setxyz(alpha * normalPoints[t.v1].x + beta * normalPoints[t.v2].x + gamma * normalPoints[t.v3].x,
				alpha * normalPoints[t.v1].y + beta * normalPoints[t.v2].y + gamma * normalPoints[t.v3].y,
				alpha * normalPoints[t.v1].z + beta * normalPoints[t.v2].z + gamma * normalPoints[t.v3].z);
		}
		else {
			normal.setxyz(alpha * normalPoints[t.v1].x + beta  * normalPoints[t.v3].x,
				alpha * normalPoints[t.v1].y + beta  * normalPoints[t.v3].y,
				alpha * normalPoints[t.v1].z + beta  * normalPoints[t.v3].z);
		}

		v.setxyz(0 - p_3D.x, 0 - p_3D.y, 0 - p_3D.z);
		l.setxyz(Pl.x - p_3D.x, Pl.y - p_3D.y, Pl.z - p_3D.z);

		normal = normalize(normal);
		v = normalize(v);
		l = normalize(l);

		float kdAux = kd, ksAux = ks;
		if (innerProduct(v, normal) < 0) normal.setxyz(0 - normal.x, 0 - normal.y, 0 - normal.z);

		if (innerProduct(normal, l) < 0) {
			kdAux = 0;
			ksAux = 0;
		}

		/*Calculate the vector R = (2(N . L)N) - L. And normalize R.*/
		float a = innerProduct(normal, l);
		Point b = escalarVectorProduct(2.0, normal);
		b = escalarVectorProduct(a, b);
		r = subtractVectors(b, l);
		r = normalize(r);

		if (innerProduct(r, v) < 0) ksAux = 0;

		color = phong(l, v, r, normal, kdAux, ksAux);

		float red = color.x, green = color.y, blue = color.z;

		if (red > 255) red = 255.0;
		if (green > 255) green = 255.0;
		if (blue > 255) blue = 255.0;

		red /= 255.0; green /= 255.0;  blue /= 255.0;

		color.setxyz(red, green, blue);
		return color;
	}
	color.setxyz(0, 0, 0);
	return color;
}

void fillBottomFlatTriangle(Triangle t, Point_2D v1, Point_2D v2, Point_2D v3)
{
	float invslope1 = (v2.x - v1.x) / (v2.y - v1.y);
	float invslope2 = (v3.x - v1.x) / (v3.y - v1.y);

	if (invslope1 > invslope2)   swap(invslope1, invslope2);
	float curx1 = v1.x;
	float curx2 = v1.x;

	Point color_v1, color_v2, color_v3, bar, point, color_point;
	// vertex colors in gouraud
	if (gr) {
		color_v1 = ordinaryCalculation2(v1.x, v1.y, t, false);
		color_v2 = ordinaryCalculation2(v2.x, v2.y, t, false);
		color_v3 = ordinaryCalculation2(v3.x, v3.y, t, false);
	}
	Point p_3D;
	for (int scanlineY = v1.y; scanlineY <= v2.y; scanlineY++)
	{
		for (int xscan = curx1; xscan <= curx2; xscan++) {
			if (gr) {
				point.setxyz(xscan, scanlineY, 1);
				bar = calculateAlphaBetaGamma(point, v1, v2, v3);

				p_3D.setxyz(bar.x * points[t.v1].x + bar.y * points[t.v2].x + bar.z * points[t.v3].x,
					bar.x * points[t.v1].y + bar.y * points[t.v2].y + bar.z * points[t.v3].y,
					bar.x * points[t.v1].z + bar.y * points[t.v2].z + bar.z * points[t.v3].z);

				float z = p_3D.z;

				if ((z < z_buffer[xscan][scanlineY]) && (xscan < (SCREEN_WIDTH + 1)) && scanlineY < (SCREEN_HEIGHT + 1)) {
					z_buffer[xscan][scanlineY] = z;

					// pixel's color by baricentric combination of triangle vertex colors
					color_point.setxyz(bar.x * color_v1.x + bar.y * color_v2.x + bar.z * color_v3.x,
						bar.x * color_v1.y + bar.y * color_v2.y + bar.z * color_v3.y,
						bar.x * color_v1.z + bar.y * color_v2.z + bar.z * color_v3.z);

					// coloring the pixel
					float red = color_point.x, green = color_point.y, blue = color_point.z;

					glColor3f(red, green, blue);
					drawDot(xscan, scanlineY);
				}
			}
			else {
				ordinaryCalculation(xscan, scanlineY, t, false);
			}
		}
		curx1 += invslope1;
		curx2 += invslope2;
	}
}

void fillTopFlatTriangle(Triangle t, Point_2D v1, Point_2D v2, Point_2D v3)
{
	float invslope1 = (v3.x - v1.x) / (v3.y - v1.y);
	float invslope2 = (v3.x - v2.x) / (v3.y - v2.y);

	if (invslope2 > invslope1) swap(invslope1, invslope2);
	float curx1 = v3.x;
	float curx2 = v3.x;

	Point color_v1, color_v2, color_v3, bar, point, color_point, Ia, Ib, Ip;
	// vertex color in gouraud
	if (gr) {
		color_v1 = ordinaryCalculation2(v1.x, v1.y, t, false);
		color_v2 = ordinaryCalculation2(v2.x, v2.y, t, false);
		color_v3 = ordinaryCalculation2(v3.x, v3.y, t, false);
	}
	Point p_3D;
	for (int scanlineY = v3.y; scanlineY >= v1.y; scanlineY--)
	{
		for (int xscan = curx1; xscan <= curx2; xscan++) {
			if (gr) {
				point.setxyz(xscan, scanlineY, 1);
				bar = calculateAlphaBetaGamma(point, v1, v2, v3);

				p_3D.setxyz(bar.x * points[t.v1].x + bar.y * points[t.v2].x + bar.z * points[t.v3].x,
					bar.x * points[t.v1].y + bar.y * points[t.v2].y + bar.z * points[t.v3].y,
					bar.x * points[t.v1].z + bar.y * points[t.v2].z + bar.z * points[t.v3].z);

				float z = p_3D.z;

				if ((z < z_buffer[xscan][scanlineY]) && (xscan < (SCREEN_WIDTH + 1)) && scanlineY < (SCREEN_HEIGHT + 1)) {
					z_buffer[xscan][scanlineY] = z;

					// pixel's color by baricentric combination of triangle vertex colors
					color_point.setxyz(bar.x * color_v1.x + bar.y * color_v2.x + bar.z * color_v3.x,
						bar.x * color_v1.y + bar.y * color_v2.y + bar.z * color_v3.y,
						bar.x * color_v1.z + bar.y * color_v2.z + bar.z * color_v3.z);

					//coloring the pixel
					float red = color_point.x, green = color_point.y, blue = color_point.z;

					glColor3f(red, green, blue);
					drawDot(xscan, scanlineY);
				}
			}
			else {
				ordinaryCalculation(xscan, scanlineY, t, false);
			}
		}
		curx1 -= invslope1;
		curx2 -= invslope2;
	}
}


void drawTriangle(Triangle t, Point_2D v1, Point_2D v2, Point_2D v3) {
	Point_2D a, b, c;

	int va = t.v1, vb = t.v2, vc = t.v3;
	/*order the points internally in relation to the value of X*/
	if (screenPoints[va].y > screenPoints[vb].y)
		swap(va, vb);
	if (screenPoints[va].y > screenPoints[vc].y)
		swap(va, vc);
	if (screenPoints[vb].y > screenPoints[vc].y)
		swap(vb, vc);

	Triangle t1;
	t1.setvs(va, vb, vc);

	a = screenPoints[va];
	b = screenPoints[vb];
	c = screenPoints[vc];

	t = t1;

	/* checks the base case for bottomTriangle */
	if (b.y == c.y) fillBottomFlatTriangle(t, a, b, c);

	/* checks the base case for topTriangle */
	else if (a.y == b.y) fillTopFlatTriangle(t, a, b, c);

	else {
		/*general case - divide the triangle in two and apply to the bottom and top triangle*/
		Point_2D d;
		//printf("fillTopFlatTriangle && fillBottomFlatTriangle\n");
		d.setxyz((a.x + ((float)(b.y - a.y) / (float)(c.y - a.y)) * (c.x - a.x)), b.y);/*intermediate point for dividing into two triangles*/
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