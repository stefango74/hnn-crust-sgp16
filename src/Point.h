/*
	Copyright 2015, 2016 Stefan Ohrhallinger
	This file is part of TestReconstruct2D.

    TestReconstruct2D is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TestReconstruct2D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TestReconstruct2D.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef POINT_H_
#define POINT_H_

#include <array>
#include <math.h>

using namespace std;

#ifdef OLD
class Point
{
	double coord[2];

public:
	Point()
	{
		coord[0] = 0.0;
		coord[1] = 0.0;
	}

	Point(double x, double y)
	{
		coord[0] = x;
		coord[1] = y;
	}

	double& operator[](int const& index)
	{
		return coord[index];
	}

	bool operator==(const Point& p) const
	{
		return ((coord[0] == p.coord[0]) && (coord[1] == p.coord[1]));
	}

	/*
	 * return Euclidean distance between points
	 */
	double distance(Point p)
	{
		double d[2] = { p[0] - coord[0], p[1] - coord[1] };

		return sqrt(d[0]*d[0] + d[1]*d[1]);
	}
};
#endif

struct Point
{
	array<float, 2> data;

	Point() {}

	Point(const float& x, const float& y)
	{
		data[0] = x;
		data[1] = y;
	}

	float & operator[](int ind) {return data[ind];}

	Point operator-(const Point &rhs) const
	{
		return Point(this->x()-rhs.x(), this->y()-rhs.y());
	}

	Point operator+(const Point &rhs) const
	{
		return Point(this->x()+rhs.x(), this->y()+rhs.y());
	}

	//dot product
	float operator*(const Point &rhs) const
	{
		return (this->x()*rhs.x()) + (this->y()*rhs.y());
	}

	//scale
	Point operator*(const float &s) const
	{
		return Point(this->x()*s, this->y()*s);
	}

	// dot product
	float dot(Point p)
	{
		return x()*p.x() + y()*p.y();
	}

	float squared_length()
	{
		return x()*x() + y()*y();
	}

	float squared_distance(const Point &p0)
	{
		auto d = p0 - *this;
		return d.squared_length();
	}

	float distance(const Point &p0)
	{
		return sqrt(squared_distance(p0));
	}

	void normalize()
	{
		float l = sqrt(squared_length());
		data[0] /= l;
		data[1] /= l;
	}

	const float &x() const  {return data[0];}
	const float &y() const  {return data[1];}
};


#endif /* POINT_H_ */
