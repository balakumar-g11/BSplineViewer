#ifndef _BEZIER_H
#define _BEZIER_H

#include <vector>
#include <map>
#include <string>
#include "DataTypes.h"

using namespace std;

namespace BSplineViewer
{
	class Bezier
	{
	public:
		Bezier(double tessEps, vector<Point2d>& vecPoints);

		~Bezier();

		void adaptiveRender(vector<Point2d>& vecBez);

	private:
		void midSubdivide(vector<Point2d>& vecBez, vector<Point2d>& vecLeftBez, vector<Point2d>& vecRightBez);

		void formVector(Point2d startPoint, Point2d endPoint, Point2d& lineVector);

		double crossProduct2D(Point2d vector1, Point2d vector2);

		double lengthOfVector2D(Point2d vect);

		double maxDistance(vector<Point2d>& vecBez);

		Point2d deCasteljau(double, int, int);

		void lerp(Point2d& point1, Point2d& point2, Point2d& resultPoint, double tVal);

		void drawLine(Point2d startPoint, Point2d endPoint);

		map<int, vector<Point2d>> maplerpSteps;

		double dTessEps;

		vector<Point2d> vecBezierPoints;
	};
}

#endif
