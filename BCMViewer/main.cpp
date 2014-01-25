/*
 * BCMViewer - BCM mesh viewer
 *
 * Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <ugl.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include "GridBCM.h"

template<typename T> inline T MAX_(T A, T B)
{
	return ( A < B ? B : A);
}

template<typename T> inline T MIN_(T A, T B)
{
	return ( A > B ? B : A);
}

static GridBCM *gridBcm = NULL;


static vec2<float> windowSize(800.0, 600.0);

//static Projection *projection;
//static Control    *control;

static ViewController *control;

static int pmx       = 0;
static int pmy       = 0;
static int mouseMode = 0; // 0 none, 1 : rotate, 2 : trans, 3 : zoom

static bool OnShift = false;
static bool OnCtrl  = false;

static int slicePos[3] = {0};
static GridBCM::AXIS mode = GridBCM::AXIS_X;

void Motion( int x, int y ){
	if(mouseMode == 1){
		float rx = static_cast<float>(y - pmy) / windowSize[0] * 10.0;
		float ry = static_cast<float>(x - pmx) / windowSize[1] * 10.0;
		control->SetRotateX(rx);
		control->SetRotateY(ry);
		pmx = x;
		pmy = y;	
		return;
	}
	
	if(mouseMode == 2){
		float tx = static_cast<float>(x - pmx)  / windowSize[0];
		float ty = static_cast<float>(pmy - y)  / windowSize[1];
		control->SetTranslate(vec3<float>(tx, ty, 0.0));
		pmx = x;
		pmy = y;
		return;
	}

	if(mouseMode == 3){
		float tz = static_cast<float>(y - pmy) / windowSize[1];
		control->SetTranslate(vec3<float>(0.0, 0.0, tz));
		pmx = x;
		pmy = y;
		return;
	}
	
}

void Mouse(int button, int state, int x, int y)
{
	pmx = x;
	pmy = y;

	if( state == GLUT_UP ){
		mouseMode = 0;
		return;
	}
	switch(glutGetModifiers()){
		case GLUT_ACTIVE_SHIFT:
			mouseMode = 2;
			break;
		case GLUT_ACTIVE_CTRL:
			mouseMode = 3;
			break;
		default:
			mouseMode = 1;
			break;
	}
}


void Display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	control->ViewUpdate();

	//glutSolidTeapot(1.0);
	
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);

	gridBcm->Render();

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	glutSwapBuffers();
}



void Keyboard(unsigned char key, int x, int y)
{
	static const char* axisStr[5] = {"", "AXIS_X", "AXIS_Y", "", "AXIS_Z" };

	switch(key){
		case 27: // ESC Key
			exit(0);
			break;

		case 'x':
			mode = GridBCM::AXIS_X;
			gridBcm->SetActiveAxis(mode);
			gridBcm->AddSlicePlane(mode);
			break;

		case 'y':
			mode = GridBCM::AXIS_Y;
			gridBcm->SetActiveAxis(mode);
			gridBcm->AddSlicePlane(mode);
			break;

		case 'z':
			mode = GridBCM::AXIS_Z;
			gridBcm->SetActiveAxis(mode);
			gridBcm->AddSlicePlane(mode);
			break;

		case 'd':
			gridBcm->AddSlicePlane(mode);
			break;

		case 'D':
			gridBcm->DeleteSlicePlane(mode);
			break;

		case 'n':
		{
			printf("mode : %s\n", axisStr[mode]);
			int i = 0;
			if(mode == GridBCM::AXIS_X) i = 0;
			if(mode == GridBCM::AXIS_Y) i = 1;
			if(mode == GridBCM::AXIS_Z) i = 2;
			size_t dn = 1;
			slicePos[i] = gridBcm->SetSlicePosition(slicePos[i] + dn);
			break;
		}

		case 'N':
		{
			printf("mode : %s\n", axisStr[mode]);
			int i = 0;
			if(mode == GridBCM::AXIS_X) i = 0;
			if(mode == GridBCM::AXIS_Y) i = 1;
			if(mode == GridBCM::AXIS_Z) i = 2;
			size_t dn = 10;
			slicePos[i] = gridBcm->SetSlicePosition(slicePos[i] + dn);
			break;
		}

		case 'p':
		{
			printf("mode : %s\n", axisStr[mode]);
			int i = 0;
			if(mode == GridBCM::AXIS_X) i = 0;
			if(mode == GridBCM::AXIS_Y) i = 1;
			if(mode == GridBCM::AXIS_Z) i = 2;
			size_t dp = 1;
			size_t pos = slicePos[i] < dp ? 0 : slicePos[i] - dp;
			slicePos[i] = gridBcm->SetSlicePosition(pos);
			break;
		}

		case 'P':
		{
			printf("mode : %s\n", axisStr[mode]);
			int i = 0;
			if(mode == GridBCM::AXIS_X) i = 0;
			if(mode == GridBCM::AXIS_Y) i = 1;
			if(mode == GridBCM::AXIS_Z) i = 2;
			size_t dp = 10;
			size_t pos = slicePos[i] < dp ? 0 : slicePos[i] - dp;
			slicePos[i] = gridBcm->SetSlicePosition(pos);
			break;
		}

		case 'g':
		{
			gridBcm->ShowGrid(!gridBcm->IsShowGrid());
		}

		default :
			break;
	}
}

void SpecialKey(unsigned char key, int x, int y)
{

}


void Resize(int w, int h)
{
	glViewport(0, 0, w, h);
	windowSize = vec2<float>(w, h);
	control->ResizeWindow(windowSize);
}

void Idle()
{
	glutPostRedisplay();
}


void Init()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	static GLfloat lightpos[]  = { 3.0f, 4.0f, 5.0f, 1.0f };
	static GLfloat lightdiff[] = { 0.6f, 0.6f, 0.6f, 0.6f };

	glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightdiff);
}


int main(int argc, char *argv[])
{
	std::string filename;
	if(argc != 2){
		fprintf(stderr, "Useage : %s <filename> \n", argv[0]);
		return -1;
	}

	filename = std::string(argv[1]);

	glutInitWindowSize(windowSize[0], windowSize[1]);
	glutInit(&argc, argv);
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );

	glutCreateWindow(argv[0]);

	glutDisplayFunc(Display);
	glutReshapeFunc(Resize);
	glutIdleFunc(Idle);
	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);

	glutKeyboardFunc(Keyboard);
	
	int n = 0;
	glGetIntegerv( GL_MAX_TEXTURE_SIZE, &n); printf("GL_MAX_TEXTURE_SIZE, %d\n", n);


	gridBcm = new GridBCM(filename.c_str());

	for(int i = 0; i < 3; i++){
		GridBCM::AXIS axis = static_cast<GridBCM::AXIS>((0x1) << i);
		slicePos[i] = gridBcm->GetCellCount(axis) / 2;
	}
	
	mode = GridBCM::AXIS_Z;
	gridBcm->SetActiveAxis(mode);
	gridBcm->AddSlicePlane(mode);
	gridBcm->DeleteSlicePlane(GridBCM::AXIS_X);
	gridBcm->DeleteSlicePlane(GridBCM::AXIS_Y);
	
	const Vec3r org = gridBcm->GetGlobalOrigin();
	const Vec3r rgn = gridBcm->GetGlobalRegion();
	
	BBox box( vec3<float>(org.x, org.y, org.z), vec3<float>(rgn.x, rgn.y, rgn.z) );

	control = new ViewController(windowSize, box);
	
	Init();

	glutMainLoop();

	return 0;
}

