#include "Bezier.h"
#include <algorithm>
#include <iostream>
#include "glut.h"

using namespace BSplineViewer;
using namespace std;

Bezier::Bezier(double tessEps, vector<Point2d>& vecPoints) 
{
	dTessEps = tessEps;
	vecBezierPoints = vecPoints;
}

Bezier::~Bezier()
{
}

Point2d Bezier::deCasteljau(double tVal, int index, int lerpLevel)
{
	if (lerpLevel == 0)
		return vecBezierPoints[index];

	Point2d pointAtTVal = (deCasteljau(tVal, index, lerpLevel - 1)*(1 - tVal)) + 
							(deCasteljau(tVal, index + 1, lerpLevel - 1)*tVal);

	maplerpSteps[lerpLevel].push_back(pointAtTVal);

	return pointAtTVal;
}

void Bezier::midSubdivide(vector<Point2d>& vecBez, vector<Point2d>& vecLeftBez, vector<Point2d>& vecRightBez)
{
	double tVal = 0.5f;
	
	maplerpSteps.clear();

	deCasteljau(tVal, 0, vecBez.size()-1);
	
	vecLeftBez.push_back(vecBez.front());
	vecRightBez.push_back(vecBez.back());

	for (int i = vecBez.size()-1; i > 0; i--)
	{
		vecLeftBez.push_back(maplerpSteps[i].front());
		vecRightBez.push_back(maplerpSteps[i].back());
	}
}

void Bezier::adaptiveRender(vector<Point2d>& vecBez)
{
	double height = maxDistance(vecBez);
	if (height < dTessEps)
	{
		drawLine(vecBez.front(), vecBez.back());
		return;
	}
	else
	{
		vector<Point2d> leftBez;
		vector<Point2d> rightBez;
		midSubdivide(vecBez, leftBez, rightBez);
		adaptiveRender(leftBez);
		adaptiveRender(rightBez);
	}
}

void Bezier::formVector(Point2d startPoint, Point2d endPoint, Point2d& lineVector)
{
	lineVector.x = endPoint.x - startPoint.x;
	lineVector.y = endPoint.y - startPoint.y;
}

double Bezier::crossProduct2D(Point2d vector1, Point2d vector2)
{
	return fabs((vector1.x * vector2.y) - (vector1.y * vector2.x));
}

double Bezier::lengthOfVector2D(Point2d vect)
{
	return sqrt(vect.x * vect.x + vect.y * vect.y);
}

double Bezier::maxDistance(vector<Point2d>& vecBez)
{
	double maxDist = -1.0f;

	if (vecBez.size() <= 0)
		return maxDist;

	vector<double> vecDistance;
	double baseLineLen;

	vector<Point2d>::iterator vecCntPtsIter = vecBez.begin();

	Point2d baseLine;
	Point2d startPoint = vecBez.front();
	Point2d endPoint = vecBez.back();

	formVector(startPoint, endPoint, baseLine);
	baseLineLen = lengthOfVector2D(baseLine);

	vecCntPtsIter++;

	for (; vecCntPtsIter != vecBez.end(); vecCntPtsIter++)
	{
		Point2d ctrlLine;

		formVector(startPoint, (*vecCntPtsIter), ctrlLine);

		double crossProd = crossProduct2D(baseLine, ctrlLine);

		double distance = crossProd / baseLineLen;

		vecDistance.push_back(distance);
	}

	if (vecDistance.size() > 0)
	{
		sort(vecDistance.begin(), vecDistance.end());
		maxDist = vecDistance.back();
	}

	return maxDist;
}

void Bezier::lerp(Point2d& point1, Point2d& point2, Point2d& resultPoint, double tVal)
{
	resultPoint.x = (((1 - tVal) * point1.x) + (tVal * point2.x));
	resultPoint.y = (((1 - tVal) * point1.y) + (tVal * point2.y));
}

void Bezier::drawLine(Point2d startPoint, Point2d endPoint)
{
	glColor3f(0.0, 1.0, 0.0);
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glVertex2f(startPoint.x, startPoint.y);
	glVertex2f(endPoint.x, endPoint.y);
	glEnd();
}
