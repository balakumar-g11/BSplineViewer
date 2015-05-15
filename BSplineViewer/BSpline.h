#ifndef _BSPLINES_H
#define _BSPLINES_H

#include <vector>
#include <string>
#include "DataTypes.h"

using namespace std;

namespace BSplineViewer
{
	class BSpline
	{
	public:
		BSpline();

		~BSpline();

		bool loadFile(string);

		void adaptiveRender();

		void uniformRender();

		void drawControlPolygon();

		void setTessellationNum(int num)
		{
			tessNum = num;
		}

		void setApproximationLimit(double num)
		{
			tessEps = num;
		}

		int getTessellationNum()
		{
			return tessNum;
		}	

		double getApproximationLimit()
		{
			return tessEps;
		}

		void toggleSamplingPoints()
		{
			samplingPoint = !samplingPoint;
		}

		void toggleControlPoints()
		{
			controlPoint = !controlPoint;
		}

		void chooseWindow(double& winLen, double& winLLx, double& winLLy);

	private:
		void extractBezier(vector<Point2d>& vecBez, int ind);

		void plotBezier(vector<Point2d>& vecBez);

		double maxDistance(vector<Point2d>& vecBez);

		void formVector(Point2d startPoint, Point2d endPoint, Point2d& lineVector);

		double crossProduct2D(Point2d vector1, Point2d vector2);

		double lengthOfVector2D(Point2d vect);

		void lerp(Point2d& point1, Point2d& point2, Point2d& resultPoint, double tVal);

		void deCasteljau(vector<Point2d>& vecBez, vector<Point2d>& vecLerpPoints, double tVal);

		void midSubDivide(vector<Point2d>& vecBez, vector<Point2d>& vecLeftBez, vector<Point2d>& vecRightBez);

		void drawLine(Point2d startPoint, Point2d endPoint); 
				
		Point2d deBoor(double tVal, int point, int lerpLevel);

		void analyseBSpline();

		int findKnot(double tVal);

		void drawLines(vector<Point2d>& vecLineStrip);

		int iDegree;   // degree of the B-spline curve

		int iCntNum;   // number of the deBoor points of the B-spline curve

		vector<double> vecKnots;   // knot vector of the B-spline curve
		
		vector<Point2d> vecCntPts;     // control points of the B-spline curve

		int tessNum;   // number of sampling points

		double tessEps;     // approximation error for tessellation

		bool samplingPoint;

		bool controlPoint;
	};
}

#endif
