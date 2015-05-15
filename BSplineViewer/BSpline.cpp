#include "BSpline.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include "Bezier.h"
#include "glut.h"

using namespace BSplineViewer;

BSpline::BSpline()
{
	tessNum = 10;
	tessEps = 2.;
	samplingPoint = false;
	controlPoint = true;
}

BSpline::~BSpline()
{
	vecKnots.clear();
	vecCntPts.clear();
}

void BSpline::extractBezier(vector<Point2d>& bez, int ind)
{
	int     i, j;
	int     k;
	double  knots[50];
	Point2d cnt[30];

	k = iDegree;

	// copy one segment
	for (i = ind - k, j = 0; i <= ind; i++) {
		cnt[j].x = vecCntPts[i].x;
		cnt[j].y = vecCntPts[i].y;
		j++;
	}
	for (i = ind - k, j = 0; i <= ind + k + 1; i++) {
		knots[j] = vecKnots[i];
		j++;
	}

	// insert knots to make the left end be Bezier end
	while (1) {
		for (i = k - 1; i>0; i--) {
			if (knots[i] < knots[k]) {
				j = i;
				break;
			}
			j = 0;
		}

		if (j == 0) break;

		// update control points
		for (i = 0; i<j; i++) {
			cnt[i].x = ((knots[k + 1 + i] - knots[k]) / (knots[k + i + 1] - knots[i + 1]))*cnt[i].x
				+ ((knots[k] - knots[i + 1]) / (knots[k + i + 1] - knots[i + 1]))*cnt[i + 1].x;
			cnt[i].y = ((knots[k + 1 + i] - knots[k]) / (knots[k + i + 1] - knots[i + 1]))*cnt[i].y
				+ ((knots[k] - knots[i + 1]) / (knots[k + i + 1] - knots[i + 1]))*cnt[i + 1].y;
		}
		// update knots
		for (i = 0; i<j; i++)
			knots[i] = knots[i + 1];
		knots[j] = knots[k];
	}

	// insert knots to make the right end be Bezier end
	while (1) {
		for (i = k + 2; i< k + k + 1; i++) {
			if (knots[i] > knots[k + 1]) {
				j = i;
				break;
			}
			j = 0;
		}

		if (j == 0) break;

		// update control points
		for (i = k; i >= j - k; i--) {
			cnt[i].x = ((knots[k + i] - knots[k + 1]) / (knots[k + i] - knots[i]))*cnt[i - 1].x
				+ ((knots[k + 1] - knots[i]) / (knots[k + i] - knots[i]))*cnt[i].x;
			cnt[i].y = ((knots[k + i] - knots[k + 1]) / (knots[k + i] - knots[i]))*cnt[i - 1].y
				+ ((knots[k + 1] - knots[i]) / (knots[k + i] - knots[i]))*cnt[i].y;
		}
		// update knots
		for (i = k + k + 1; i>j; i--)
			knots[i] = knots[i - 1];
		knots[j] = knots[k + 1];
	}

	// return the Bezier control points
	for (i = 0; i <= iDegree; i++) {
		Point2d pnt;
		pnt.x = cnt[i].x;
		pnt.y = cnt[i].y;
		bez.push_back(pnt);
	}
}

void BSpline::formVector(Point2d startPoint, Point2d endPoint, Point2d& lineVector)
{
	lineVector.x = endPoint.x - startPoint.x;
	lineVector.y = endPoint.y - startPoint.y;
}

double BSpline::crossProduct2D(Point2d vector1, Point2d vector2)
{
	return fabs((vector1.x * vector2.y) - (vector1.y * vector2.x));
}

double BSpline::lengthOfVector2D(Point2d vect)
{
	return sqrt(vect.x * vect.x + vect.y * vect.y);
}

double BSpline::maxDistance(vector<Point2d>& vecBez)
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

void BSpline::lerp(Point2d& point1, Point2d& point2, Point2d& resultPoint, double tVal)
{
	resultPoint.x = (((1 - tVal) * point1.x) + (tVal * point2.x));
	resultPoint.y = (((1 - tVal) * point1.y) + (tVal * point2.y));
}

void BSpline::deCasteljau(vector<Point2d>& vecBez, vector<Point2d>& vecLerpPoints, double tVal)
{
	vector<vector<Point2d>> vecLerpSteps;
	vecLerpSteps.push_back(vecBez);

	for (int i = 0; i < iDegree; i++)
	{
		vector<Point2d> lerpStep;
		vector<Point2d>::iterator lerpPtsIter = vecLerpSteps[i].begin();
		for (; lerpPtsIter != vecLerpSteps[i].end(); lerpPtsIter++)
		{
			Point2d point1 = *(lerpPtsIter);
			lerpPtsIter++;

			if (lerpPtsIter == vecLerpSteps[i].end())
				break;

			Point2d point2 = *(lerpPtsIter);

			Point2d resultPoint;

			lerp(point1, point2, resultPoint, tVal);

			lerpStep.push_back(resultPoint);

			lerpPtsIter--;
		}
		vecLerpSteps.push_back(lerpStep);
	}

	vector<vector<Point2d>>::iterator lerpStepIter = vecLerpSteps.begin();
	for (; lerpStepIter != vecLerpSteps.end(); lerpStepIter++)
	{
		Point2d point;
		point.x = (*lerpStepIter).front().x;
		point.y = (*lerpStepIter).front().y;
		vecLerpPoints.push_back(point);
	}

	lerpStepIter = vecLerpSteps.end();
	lerpStepIter--;
	lerpStepIter--;
	for (; lerpStepIter != vecLerpSteps.begin(); lerpStepIter--)
	{
		Point2d point;
		point.x = (*lerpStepIter).back().x;
		point.y = (*lerpStepIter).back().y;
		vecLerpPoints.push_back(point);
	}

	Point2d point;
	point.x = vecBez.back().x;
	point.y = vecBez.back().y;
	vecLerpPoints.push_back(point);
}


void BSpline::midSubDivide(vector<Point2d>& vecBez, vector<Point2d>& vecLeftBez, vector<Point2d>& vecRightBez)
{
	vector<Point2d> subDividedPoints;
	double tVal = 0.5f;

	deCasteljau(vecBez, subDividedPoints, tVal);

	int midPoint = subDividedPoints.size() / 2;

	for (int i = 0; i <= midPoint; i++)
	{
		vecLeftBez.push_back(subDividedPoints[i]);
	}

	for (size_t i = midPoint; i < subDividedPoints.size(); i++)
	{
		vecRightBez.push_back(subDividedPoints[i]);
	}
}

void BSpline::drawLine(Point2d startPoint, Point2d endPoint)
{
	glColor3f(0.0, 1.0, 0.0);
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glVertex2f(startPoint.x, startPoint.y);
	glVertex2f(endPoint.x, endPoint.y);
	glEnd();

	if (samplingPoint)
	{
		glBegin(GL_POINTS);
		glPointSize(4.0f);
		glColor3f(0.0, 1.0, 1.0);
		glVertex2f(startPoint.x, startPoint.y);
		glVertex2f(endPoint.x, endPoint.y);
		glEnd();
	}
}

void BSpline::drawControlPolygon()
{
	if (!controlPoint)
		return;

	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(3.0);
	glBegin(GL_LINE_STRIP);      // display the control polygon
	for (int i = 0; i< iCntNum; i++)
		glVertex2f(vecCntPts[i].x, vecCntPts[i].y);
	glEnd();

	glPointSize(6.0);            // display the control points
	glBegin(GL_POINTS);
	for (int i = 0; i<iCntNum; i++)
		glVertex2f(vecCntPts[i].x, vecCntPts[i].y);
	glEnd();
}

void BSpline::plotBezier(vector<Point2d>& vecBez)
{
	double height = maxDistance(vecBez);
	if (height < tessEps)
	{
		drawLine(vecBez.front(), vecBez.back());
		return;
	}
	else
	{
		vector<Point2d> leftBez;
		vector<Point2d> rightBez;
		midSubDivide(vecBez, leftBez, rightBez);
		plotBezier(leftBez);
		plotBezier(rightBez);
	}
}

void BSpline::adaptiveRender()
{
	vector<Point2d> bez;
	int i;

	for (i = iDegree; i< iCntNum; i++) {
		if (fabs(vecKnots[i] - vecKnots[i + 1]) < 0.00001) continue;  // no segment, skip over
		bez.clear();
		extractBezier(bez, i); 
		plotBezier(bez);
	}
}

bool BSpline::loadFile(string strFileName)
{
	bool bResult = false;
	
	FILE *fp;
	int  i;

	if ((fp = fopen(strFileName.c_str(), "r")) == NULL) return bResult;  // fail to open the file

	fscanf(fp, "%d%d", &(iDegree), &(iCntNum));

	for (i = 0; i <= iCntNum + iDegree; i++)
	{
		double knot;
		fscanf(fp, "%lf", &(knot));
		vecKnots.push_back(knot);
	}
		
	for (i = 0; i < iCntNum; i++)
	{
		Point2d cntPt;
		fscanf(fp, "%lf%lf", &cntPt.x, &cntPt.y);
		vecCntPts.push_back(cntPt);
	}
	fclose(fp);
	
	bResult = true;
	return bResult;
}

//============================================================
void BSpline::chooseWindow(double& winLen, double& winLLx, double& winLLy)
{
	int    i;
	double left, right, bottom, top;


	left = right = vecCntPts.front().x;
	for (i = 1; i< iCntNum; i++) {
		if (left > vecCntPts[i].x)  left = vecCntPts[i].x;
		if (right < vecCntPts[i].x) right = vecCntPts[i].x;
	}

	bottom = top = vecCntPts.front().y;
	for (i = 1; i< iCntNum; i++) {
		if (bottom > vecCntPts[i].y)  bottom = vecCntPts[i].y;
		if (top < vecCntPts[i].y) top = vecCntPts[i].y;
	}

	winLen = top - bottom;
	if (winLen < right - left) winLen = right - left;

	winLen += 100;
	winLLy = bottom - 50;
	winLLx = left - 50;
}

void BSpline::drawLines(vector<Point2d>& vecLineStrip)
{
	vector<Point2d>::iterator pointIter = vecLineStrip.begin();
	glBegin(GL_LINE_STRIP);
	glLineWidth(2.0f);
	glColor3f(0.0f, 0.0f, 1.0f);
	for (; pointIter != vecLineStrip.end(); pointIter++)
	{
		glVertex2f((*pointIter).x, (*pointIter).y);
	}
	glEnd();
	if (samplingPoint)
	{
		glBegin(GL_POINTS);
		glPointSize(4.0f);
		glColor3f(0.0f, 0.1f, 1.0f);
		pointIter = vecLineStrip.begin();
		for (; pointIter != vecLineStrip.end(); pointIter++)
		{
			glVertex2f((*pointIter).x, (*pointIter).y);
		}
		glEnd();
	}
}

void BSpline::uniformRender()
{
	vector<Point2d> vecInterpolation;

	double tInterval = (vecKnots[iCntNum] - vecKnots[iDegree]) / tessNum;
	
	for (int i = 0; i <= tessNum; i++)
	{
		double pointAtCurve = vecKnots[iDegree] + (i * tInterval);
		
		int curveSeg = findKnot(pointAtCurve);

		if (curveSeg < 0)
			continue;

		vecInterpolation.push_back(deBoor(pointAtCurve, curveSeg, iDegree));
	}
	
	drawLines(vecInterpolation);
}

int BSpline::findKnot(double tVal)
{
	int curveSeg = -1;

	for (int i = iDegree; i <= iCntNum+1; i++)
	{
		if (tVal >= vecKnots[i] && tVal <= vecKnots[i+1])
		{
			curveSeg = i;
			break;
		}
	}

	return curveSeg;
}

Point2d BSpline::deBoor(double tVal, int point, int lerpLevel)
{
	if (lerpLevel == 0){
		return vecCntPts[point];
	}

	Point2d point1 = deBoor(tVal, point - 1, lerpLevel - 1);
	Point2d point2 = deBoor(tVal, point, lerpLevel - 1);

	float blend = ((tVal - vecKnots[point]) / (vecKnots[point + iDegree + 1 - lerpLevel] - vecKnots[point]));

	return ((point1 * (1 - blend)) + (point2 * blend));
}
