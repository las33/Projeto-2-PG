/*
-----------------------------------------------------------------------------
OpenGL Tutorial
VOXAR Labs
Computer Science Center - CIn
Federal University of Pernambuco - UFPE
http://www.cin.ufpe.br/~voxarlabs

@author Tullio Lucena - tjsl

-----------------------------------------------------------------------------
*/

#ifndef _OPENGL_TUTORIAL_H_
#define _OPENGL_TUTORIAL_H_

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <gl/glut.h>

using namespace std;

// Constantes do programa
#define FPS			30	// Quantidade de atualizacoes por segundo do programa
#define ESC			27	// Valor da tecla ESC
#define IDLE		-2	// Nada a fazer
#define MODIFIED	-1	// A tela foi modificada

struct Quadrado{
	GLfloat lado, x, y, r, g, b;
	Quadrado(GLfloat lado=0.0, GLfloat x=0.0, GLfloat y=0.0, GLfloat r=0.0, GLfloat g=0.0, GLfloat b=0.0):lado(lado),x(x),y(y),r(r),g(g),b(b){}
};

/**
* Funcao para inicializacao do programa.
*/
void myinit();

/**
* Esta funcao vai ser chamada quando a tela for redimensionada.
* @param w Largura atual da tela.
* @param h Altura atual da tela.
*/
void myreshape (GLsizei w, GLsizei h);

/**
* Esta é a funcao responsavel por pintar os objetos na tela.
*/
void mydisplay();

/**
* Esta funcao vai ser chamada toda vez que o mouse se mover com uma tecla apertada.
* Vai atualizar a posicao do quadrado que esta sendo movido.
* @param x Coordenada x atual do mouse.
* @param y Coordenada y atual do mouse.
*/
void handleMotion(int x, int y);

/**
* Esta funcao vai ser chamada toda vez que uma tecla do mouse seja clicada ou levantada.
* Com o botao direito do mouse, um quadrado pode ser criado na posicao atual do mouse, com uma cor e tamanho aleatorio.
* Com o botao esquerdo do mouse, um quadrado pode ser movimentado.
* @param btn Tecla levantada. Ex: GLUT_LEFT_BUTTON.
* @param state Estado da tecla. Ex: GLUT_DOWN.
* @param x Coordenada x atual do mouse.
* @param y Coordenada y atual do mouse.
*/
void handleMouse(int btn, int state, int x, int y);

/**
* Esta funcao vai ser chamada toda vez que uma tecla comum seja levantada.
* @param key Tecla levantada.
* @param x Coordenada x atual do mouse.
* @param y Coordenada y atual do mouse.
*/
void hadleKeyboard(unsigned char key, int x, int y);

/**
* Esta funcao vai ser chamada toda vez que uma tecla especial seja levantada.
* Botao F5 o programa volta ao estado inicial.
* @param key Tecla levantada. Ex: GLUT_KEY_F5.
* @param x Coordenada x atual do mouse.
* @param y Coordenada y atual do mouse.
*/
void hadleSpecialKeyboard(int key, int x, int y);

/**
* Representa uma maquina de estados do programa. Vai ser executada FPS vezes por segundo.
* @param id Identificador da maquina de estados.
*/
void loop(int id);

#endif //_OPENGL_TUTORIAL_H_
