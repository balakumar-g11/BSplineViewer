// This is a template for assignment 3 of DM6122: 3D Modeling and Reconstruction.
// NTU
// August 2009
//
// Open a new project "Win32 Console Application" and add sample.c to Source Files 

#include <stdio.h>
#include <iostream>
//#include <GL/glut.h>
#include "glut.h"
#include <math.h>
#include "BSpline.h"
#include "DataTypes.h"

using namespace std;
using namespace BSplineViewer;

// B-spline curve
BSpline*  bSpline;

// The parameters are used to define a visible window in this application's World Coordinate System. 
double	winLLx = 0.0;
double  winLLy = 0.0;
double	winLen = 100.0;

// 
int    displayCP = 1;    // flag for whether the control polygon is displayed or not
int    adaptivePlot = 1;    // flag for whether the plot is adaptive or uniform
int    samplingPnt = 0;    // flag for whether the sampling points are displayed or not
int    tessNum = 10;   // number of sampling points
double tessEps = 2.;     // approximation error for tessellation


//============================================================
static void Init(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);         // set display-window color to white
	glMatrixMode(GL_PROJECTION);	           // set projection parameters
	gluOrtho2D(winLLx, winLLx + winLen, winLLy, winLLy + winLen);   // set an orthogonal projection
}

//============================================================
static void drawCurve(void)
{
	glClear(GL_COLOR_BUFFER_BIT);	// clear display window
	glColor3f(1.0, 0.0, 0.0);   // set line segment color to red


	// Draw the control polygon
	if (displayCP != 0) {
		bSpline->drawControlPolygon();
	}

	// Draw the curve
	if (adaptivePlot) {  // plot adaptively
		bSpline->adaptiveRender();
	}

	else {  // plot uniformly
		bSpline->uniformRender();
	}


	glFlush();		    // process all openGL routines as quickly as possible	         
	glutSwapBuffers();  // swap buffers to display the current frame
}

//============================================================
static void idle(void)
{
	drawCurve();
}


//============================================================
static void hotkey(unsigned char k, int x, int y)
{
	// Here we are processing keyboard events.
	switch (k)
	{
	case 27:
		delete bSpline;
		exit(0);
		break;

		// Toggle plotting the control polygon
	case 'C':
	case 'c':
		bSpline->toggleControlPoints();
		break;

		// Toggle sampling points
	case 'P':
	case 'p':
		bSpline->toggleSamplingPoints();
		break;

		// Toggle adaptive/uniform plotting
	case 'A':
	case 'a':
		adaptivePlot = !adaptivePlot;
		break;

		// Increase tessellation
	case '+':
	case '=':
		if (adaptivePlot) {
			bSpline->setApproximationLimit(bSpline->getApproximationLimit() * 0.7);
			if (bSpline->getApproximationLimit() < 0.5)
				bSpline->setApproximationLimit(0.01);
		}
		else {
			bSpline->setTessellationNum(bSpline->getTessellationNum() + 1);
			if (bSpline->getTessellationNum() > 100)
				bSpline->setTessellationNum(100);
		}
		break;

		// Decrease tessellation
	case '-':
	case '_':
		if (adaptivePlot) {
			bSpline->setApproximationLimit(bSpline->getApproximationLimit() * 1.4);
			if (bSpline->getApproximationLimit() > 50)  
				bSpline->setApproximationLimit(100);
		}
		else {
			bSpline->setTessellationNum(bSpline->getTessellationNum() - 1);
			if (bSpline->getTessellationNum() < 2)
				bSpline->setTessellationNum(1);
		}
		break;
	}
}

//============================================================
void main(int argc, char *argv[])
{
	// load the curve from a file
	string fileName;

	bSpline = new BSpline();

	cout << endl << "Please enter a filename : ";
	cin >> fileName;

	if (!bSpline->loadFile(fileName)) return;

	bSpline->chooseWindow(winLen, winLLx, winLLy);

	// help information
	printf("\n\nB-spline curve plotting\n");
	printf("NTU, September 2006\n");
	printf("\n");
	printf(" ESC      - Quit program\n");
	printf("\n");
	printf(" A/a : Toggle adaptive/uniform plotting (Default adaptive)\n");
	printf(" C/c : Toggle plotting the control polygon (Default On)\n");
	printf(" P/p : Toggle sampling points (Default Off)\n");
	printf(" +   : Increase tessellation\n");
	printf(" -   : Decrease tessellation\n");
	printf("\n");


	// set up graphics window
	glutInit(&argc, argv);                         // initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);  // set display mode
	glutInitWindowSize(650, 650);                 // set display window width and height
	glutInitWindowPosition(100, 100);             // set top-left display window position
	glutCreateWindow("Balakumar Marimuthu G1401999A:        use +, -, c, a, p, and Esc keys.");

	Init();                        // execute initialization procedure
	glutIdleFunc(idle);            // enables us to make interaction.
	glutDisplayFunc(drawCurve);    // send graphics to display window
	glutKeyboardFunc(hotkey);

	glutMainLoop();                // display everything and wait
}
