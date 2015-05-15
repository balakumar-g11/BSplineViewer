#ifndef _DATA_TYPES_H
#define _DATA_TYPES_H

struct Point2d{
	double x, y;      // x, y coordinates of a 2D point

public:
	Point2d operator*(double mul) const
	{
		Point2d temp = (*this);
		temp.x = x*mul;
		temp.y = y*mul;
		return temp;
	}

	Point2d operator+(Point2d add) const
	{
		Point2d temp = (*this);
		temp.x += add.x;
		temp.y += add.y;
		return temp;
	}
};


/*typedef struct{
	int      degree;   // degree of the B-spline curve
	int      cntNum;   // number of the deBoor points of the B-spline curve
	double   *knots;   // knot vector of the B-spline curve
	Point2d  *cnt;     // control points of the B-spline curve
} Bspline;*/

#endif

