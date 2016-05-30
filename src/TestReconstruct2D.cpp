/*
	Copyright 2015, 2016 Stefan Ohrhallinger.
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

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <png++/png.hpp>
#include <GL/glew.h>
#include <GL/glut.h>

#include "Reconstruct2D.h"

using namespace std;

/*
 * read a pts-file (x, y)
 *
 * connect where conforming to sampling condition
 * classify points as conforming, non-conforming
 * display points and connecting edges
 *
 * perturb with noise
 * determine feature for non-conforming points
 * classify points as conforming, noisy, outliers
 * determine connecting edges between selected points (how to select deterministically?)
 * display projected points, spherical fits, connecting edges
 *
 * output SVG (write primitives both to openGL and SVG)
 */

#define INIT_WIDTH 1000
#define INIT_HEIGHT 1000
int viewWidth = INIT_WIDTH;
int viewHeight = INIT_HEIGHT;

template<typename T> inline T CUB(T t) { return (SQR(t)*(t)); };

enum class PointState: char { EPSHALF, EPSONE, EPSLARGE, EPSHALF_NM, EPSONE_NM, EPSLARGE_NM };

class TestReconstruct2D
{
private:
	int viewWidth, viewHeight;
	float scale;
	float offset[2], lfsMin, lfsMax;
	bool isClosed;
	vector<Point> points, projPoints, curvePoints, extremaPoints, mPoints;
	vector<float> lfsCurvePoints, kCurvePoints;
	vector<PointState> pointStates;
	vector<Point> normals;
	vector<PointsEnum> pClasses;
	map<pair<int, int>, EdgeEnum> edgeMap;
	vector<Circle> circles;
	vector<pair<int, int> > arcs;

	void drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap);
	void drawPoints(vector<Point> &points);
	void drawPointsWithStates(vector<Point> &points, vector<PointState> &pointStates);
	void drawCurvePoints(vector<Point> &points, vector<float> &lfsCurvePoints, float minLfs, float maxLfs);
	void drawCircle(Circle &circle);
	void drawCircles(vector<Circle> &circles);
	void drawArc(Circle &c, pair<int, int> arc);
	void drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arc);
	void drawLabels(vector<Point> &points);
	void drawNormals(vector<Point> &points);
	Point transform(Point);
	float scaleOfPoints(float *translate, vector<Point> &points);

public:
	int sampleCount;
	TestReconstruct2D();
	virtual ~TestReconstruct2D();
	void reshape(int width, int height);
	void draw(string fileName, bool isDrawLines);
	void loadPointSet(string name);
	void loadCurvePointsFromPNG(char *name);
	void generateCurveFromBezierString(string str);
	void generatePointSetCosinesOnCircleWithVaryingAmplitude(float minAmp, float maxAmp);
	void generatePointSetCosinesOnCircleWithVaryingFrequency(float minFreq, float maxFreq);
	void sampleCurveByEps(float maxEps, float error, bool isClosed);
	void sampleCurveByReach(float maxReach, float error, bool isClosed);
	void generatePointSetBySamplingCurveEpsAndManifold(float error);
	void generateCurveCosinesOnCircleWithVaryingFrequency(int count, float minFreq, float maxFreq, float amplitude);
	void reconstruct(bool isClosed);
	void computeScale();
	void determinePointStates(float maxEps);
	float scaleOfCurvePoints();
};

/*
 * determine scale of curve points for unit square
 */
float TestReconstruct2D::scaleOfCurvePoints()
{
	float offset[2];

	return scaleOfPoints(offset, curvePoints);
}

/*
 * determine scale of points for unit square
 */
float TestReconstruct2D::scaleOfPoints(float *offset, vector<Point> &points)
{
	int i;
	float min[2] = { numeric_limits<float>::max(), numeric_limits<float>::max() };
	float max[2] = { -numeric_limits<float>::max(), -numeric_limits<float>::max() };

	for (auto p:points)
	{
		for (i = 0; i < 2; i++)
		{
			if (p[i] < min[i])
				min[i] = p[i];

			if (p[i] > max[i])
				max[i] = p[i];
		}
	}

	float dim[2] = { max[0] - min[0], max[1] - min[1] };

	for (i = 0; i < 2; i++)
		offset[i] = min[i];

	i = (dim[0] > dim[1]) ? 0 : 1;

	offset[1 - i] -= (dim[i] - dim[1 - i])/2;

	return dim[i];
}

void TestReconstruct2D::reconstruct(bool isClosed)
{
	Reconstruct2D *instance = new Reconstruct2D(points, isClosed);
	instance->reconstruct();
	edgeMap = instance->getEdgeMap();
	pClasses = instance->getPointClassification();
	circles = instance->getCircles();
	arcs = instance->getArcs();
	normals = instance->getNormals();
	projPoints = instance->getProjectedPoints();
	sampleCount = points.size();

	delete instance;
}

void TestReconstruct2D::computeScale()
{
	scale = scaleOfPoints(offset, points);
}

/*
 * load point set from file with name
 */
void TestReconstruct2D::loadPointSet(string name)
{
	int i;

	points.clear();
	ifstream file(name);

	if (file)
	{
		while (file)
		{
			float x, y;
			file >> x >> y;

			if (file)
				points.push_back(Point(x, y));
		}
	}
	else
	{
		cerr << "ERROR: input file " << name << " could not be read." << endl;
		exit(2);
	}

	if (points.size() < 3)
	{
		cerr << "ERROR: input file " << name << " contains less than 3 points." << endl;
		exit(3);
	}

	for (i = 0; i < (int)points.size(); i++)
	{
		pClasses.push_back(CONFORMING);
		pointStates.push_back(PointState::EPSLARGE_NM);
	}

	computeScale();
}

/*
 * load point set from PNG file with name (white pixels as points)
 */
void TestReconstruct2D::loadCurvePointsFromPNG(char *name)
{
	string filename(name);
	png::image< png::gray_pixel_1 > image(filename, png::require_color_space< png::gray_pixel_1>());

	for (png::uint_32 y = 0; y < image.get_height(); ++y)
	{
	    for (png::uint_32 x = 0; x < image.get_width(); ++x)
	    {
	    	if (image.get_pixel(x, y) == 1)
	    		curvePoints.push_back(Point(x, y));
	    }
	}
}

/*
 * converts a string two comma-separated floats into a Point
 */
Point str2point(string str)
{

	int ofs = str.find(",");

	return Point(stof(str.substr(0, ofs)), stof(str.substr(ofs + 1, str.length() - ofs - 1)));
}

/*
 * parses SVG curveto string and extracts bezier control points
 */
void parseSVGCurveToString(string str, vector<Point> &bezierVec)
{
//	"m 100,309.50506 c 282.85714,91.42857 397.14286,-74.28571 480,-11.42857"
	int i;
	Point p[2];
	assert((str[0] == 'M') || (str[0] == 'm'));
	string::size_type ofs = str.find(' ', 2);
	string pStr = str.substr(2, ofs - 1);
	p[0] = str2point(pStr);
	char cc = str[ofs + 1];
	bool isRelative = (cc == 'c');
	assert(isRelative || (cc == 'C'));
	bool isEnd = false;
	ofs += 2;

	do
	{
		bezierVec.push_back(p[0]);

		for (i = 0; i < 3; i++)
		{
			int prevOfs = ofs;
			ofs = str.find(' ', ofs + 1);

			if (ofs == string::npos)
			{
				isEnd = true;
				ofs = str.length();
			}

			pStr.clear();
			pStr = str.substr(prevOfs, ofs - prevOfs);
			p[1] = str2point(pStr);

			if (isRelative)
				p[1] = p[1] + p[0];

			bezierVec.push_back(p[1]);
		}

		if (!isEnd)
			p[0] = p[1];

	} while (!isEnd);
}

/*
 * generate curve point set from SVG curveto string
 */
void TestReconstruct2D::generateCurveFromBezierString(string str)
{
	int i, j;
	const int SAMPLE_COUNT = 300;
	string bezierStr(str);
	vector<Point> b;
	parseSVGCurveToString(bezierStr, b);

	// iterate all cubic bezier curves
	for (j = 0; j < (int)b.size()/4; j++)
	{
		// sample the cubic bezier curve by evaluating with parameter t [0..1]
		for (i = 0; i < SAMPLE_COUNT; i++)
		{
			float t = (float)i/SAMPLE_COUNT;
			Point p = b[j*4]*CUB(1 - t) + b[j*4 + 1]*3*t*SQR(1 - t) + b[j*4 + 2]*3*SQR(t)*(1 - t) + b[j*4 + 3]*CUB(t);
			curvePoints.push_back(p);
		}
	}
}

/*
 * generate point set: cosines on circle, with varying amplitude
 */
void TestReconstruct2D::generatePointSetCosinesOnCircleWithVaryingAmplitude(float minAmp, float maxAmp)
{
	int i;
	const int SAMPLE_COUNT = 100;
	const float frequency = 10;

	points.clear();

	// sample a regular-frequency, changing-amplitude sine function on a circle
	for (i = 0; i < SAMPLE_COUNT; i++)
	{
		float deg = i*PI*2/SAMPLE_COUNT;
		float cx = cos(deg);
		float cy = sin(deg);
		float currAmp = minAmp + (float)i/SAMPLE_COUNT*(maxAmp - minAmp);
		float f = cos(frequency*i*PI*2/SAMPLE_COUNT)*currAmp;
		float dx = f*cx;
		float dy = f*cy;
		Point p(cx + dx, cy + dy);
		points.push_back(p);
	}

	for (i = 0; i < (int)points.size(); i++)
		pClasses.push_back(CONFORMING);

	computeScale();
}

/*
 * return radius for circle through point p with normalized normal n and point q
 * q can be mirrored on the line through n, therefore the radius is the circumradius of the triangle pqq'
 */
float radiusForCircleThrough2PointsandNormal(Point p, Point n, Point q)
{
	float a = p.distance(q);
	Point n2(-n[1], n[0]);
	Point pq = p - q;
	float dist = abs(pq*n2);
	float b = 2*dist;

	if (b == 0)
		return 0.5*sqrt(pq.squared_length());	// distance pq = diameter

	float e = (2*a + b)*SQR(b)*(2*a - b);

	if (e <= 0.0)
		return numeric_limits<float>::max();	// triangle points are collinear, infinite radius

	float d = sqrt(e);
	return abs(SQR(a)*b/d);	// circumradius of triangle adapted to isosceles version
}

/*
 * compute LFS values for curve points
 */
void computeLFSForCurve(vector<Point> &curvePoints, vector<float> &lfsCurvePoints, vector<Point> &mPoints, bool isClosed)
{
	int i, j;

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		Point prevP, nextP, currP = curvePoints[i];

		if (i > 0)
			prevP = curvePoints[i - 1];
		else
		{
			if (isClosed)
				prevP = curvePoints[curvePoints.size() - 1];
			else
				prevP = currP;
		}

		if (i < (int)curvePoints.size() - 1)
			nextP = curvePoints[i + 1];
		else
		{
			if (isClosed)
				nextP = curvePoints[0];
			else
				nextP = currP;
		}

		Point normal = prevP - nextP;
		normal.normalize();
		swap(normal[0], normal[1]);
		normal[0] = -normal[0];
		float minR = numeric_limits<float>::max();

		for (j = 0; j < (int)curvePoints.size(); j++)
			if (i != j)
			{
				// at point p, determine radius r for maximum empty circle with center through normal n (one neighbor q on its boundary)
				Point curr2P = curvePoints[j];
				float radius = radiusForCircleThrough2PointsandNormal(currP, normal, curr2P);

				if (radius < minR)
				{
					minR = radius;
					float direction = (normal*(curr2P - currP) < 0) ? -1.0 : 1.0;
					mPoints[i] = currP + normal*(radius*direction);
				}
			}
	}

	// compute lfs from nearest medial axis points
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points;

	ann_points = annAllocPts(mPoints.size(), 2);

	for(i = 0; i < (int)mPoints.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = mPoints[i][0];
		p[1] = mPoints[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, mPoints.size(), 2);
	ANNpointArray search_point = annAllocPts(1, 2);
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray distances = new ANNdist[1];

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		// get nearest neighbor in medial axis
		search_point[0][0] = curvePoints[i][0];
		search_point[0][1] = curvePoints[i][1];
		kdTree->annkSearch(search_point[0], 1, nnIdx, distances);
		lfsCurvePoints[i] = sqrt(distances[0]);
	}

	delete nnIdx;
	delete distances;
	annDeallocPts(ann_points);
}

/*
 * compute curvature from radius of circle fitting 3 points
 */
float computeCurvature3Points(Point p0, Point p1, Point p2)
{
	// compute slopes of lines
	float m1 = (p1[1] - p0[1])/(p1[0] - p0[0]);
	float m2 = (p2[1] - p1[1])/(p2[0] - p1[0]);

	// compute circle center
	float x = (m1*m2*(p0[1] - p2[1]) + m2*(p0[0] + p1[0]) - m1*(p1[0] + p2[0]))/(2*(m2 - m1));
	float y = -1/m1*(x - (p0[0] + p1[0])/2) + (p0[1] + p1[1])/2;
	Point cc(x, y);

	float radius = cc.distance(p0);

	if (radius >= numeric_limits<float>::max())
		return 0.0;
	else
		return 1/radius;
}

/*
 * compute curvature values for curve points
 */
void computeKForCurve(vector<Point> &curvePoints, vector<float> &kCurvePoints, vector<Point> &extremaPoints)
{
	int i;

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		Point currP = curvePoints[i];
		Point prevP = curvePoints[(i + curvePoints.size() - 1) % curvePoints.size()];
		Point nextP = curvePoints[(i + 1) % curvePoints.size()];
		float k = computeCurvature3Points(prevP, currP, nextP);

		kCurvePoints.push_back(k);
	}

	// consider tolerance to determine real extrema?
	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		if ((kCurvePoints[i] > kCurvePoints[(i + curvePoints.size() - 1) % curvePoints.size()]) &&
			(kCurvePoints[i] > kCurvePoints[(i + 1) % curvePoints.size()]))
			extremaPoints.push_back(curvePoints[i]);
	}
}

const int K_MAX = 15;

/*
 * generate curve: cosines on circle, with varying frequency
 */
void TestReconstruct2D::generateCurveCosinesOnCircleWithVaryingFrequency(int count, float minFreq, float maxFreq, float amplitude)
{
	int i;

	// sample a frequency-varying, constant amplitude cosine function on a circle
	for (i = 0; i < count; i++)
	{
		float deg = i*PI*2/count;
		float cx = cos(deg);
		float cy = sin(deg);
		float currFreq = minFreq + (float)i/count*(maxFreq - minFreq);
//		float f = cos(currFreq*PI*2/2*i/count)*amplitude;	// divide by 2 for extra factor from differentiation
		float f = cos(currFreq*PI*2*i/count)*amplitude;
		float dx = f*cx;
		float dy = f*cy;
		Point p(cx + dx, cy + dy);
		curvePoints.push_back(p);
	}
}

/*
 * return distance of p0 from line p1-p2
 */
float distancePFromLine(Point p0, Point p1, Point p2)
{
	Point normal(p2 - p1);
	swap(normal[0], normal[1]);
	normal[0] = -normal[0];
	normal.normalize();
	return abs(normal*(p0 - p1));
}

/*
 * determine pointState and pClasses
 */
void TestReconstruct2D::determinePointStates(float maxEps)
{
	int i, kMax = (K_MAX > points.size()) ? points.size() : K_MAX;

	// construct kd-tree
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points;

	ann_points = annAllocPts(points.size(), 2);
	ANNidxArray nnIdx = new ANNidx[kMax];
	ANNdistArray distances = new ANNdist[kMax];

	for (i = 0; i < (int)points.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = points[i][0];
		p[1] = points[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, points.size(), 2);

	for (i = 0; i < (int)points.size(); i++)
	{
		pClasses.push_back(CONFORMING);
		int orderedN[2] = { (int)((i + points.size() - 1) % points.size()), (int)((i + 1) % points.size()) };

		// determine whether point is uniquely connected
		Point currP = points[i];
		Point prevP = points[orderedN[0]];
		Point nextP = points[orderedN[1]];
		float nDist = prevP.distance(nextP);

		if ((currP.distance(prevP) > nDist) && (currP.distance(nextP) > nDist))
			cout << "failed assumption: eps<" << maxEps << " makes samples uniquely connected" << endl;

		// determine whether point is bijective/manifold
		sort(orderedN, &orderedN[2]);

		kdTree->annkSearch(ann_points[i], kMax, nnIdx, distances);
		int nn[2] { nnIdx[1], -1 };

		// add neighbors until feature condition fulfilled
		int k = 3;
		bool featureFitted = false;

		do
		{
			// check feature size
			nn[1] = nnIdx[k - 1];
			float nnDist = points[nn[0]].distance(points[nn[1]]);
			featureFitted = ((nnDist > points[i].distance(points[nn[0]])) && (nnDist > points[i].distance(points[nn[1]])));
			k++;
		} while (!featureFitted && (k < kMax));

		// check if feature neighbors are actual neighbors = manifold condition
		sort(nn, &nn[2]);
		bool isBijective = ((nn[0] == orderedN[0]) && (nn[1] == orderedN[1]));

		if (isBijective)
			if (maxEps <= 0.5)
				pointStates.push_back(PointState::EPSHALF);
			else
				pointStates.push_back(PointState::EPSONE);
		else
			if (maxEps <= 0.5)
				pointStates.push_back(PointState::EPSHALF_NM);
			else
				pointStates.push_back(PointState::EPSONE_NM);
	}

	delete nnIdx;
	delete distances;
	annDeallocPts(ann_points);
}

/*
 * generate point set by sampling with condition < epsMax, max error to original curve and optionally manifold condition
 */
void TestReconstruct2D::sampleCurveByEps(float maxEps, float error, bool isClosed)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	vector<int> curveIndices;
	i = 0;

	do
	{
		int prevI = i;
		Point prevP = curvePoints[prevI];
		curveIndices.push_back(prevI);
		points.push_back(prevP);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];

			j = prevI + 1;
			isConforming = true;

			// test eps and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for epsilon condition (a sample within dist/lfs < maxEps)
					float lfs = lfsCurvePoints[j];

					float dist = currP.distance(prevP);

					if (dist/lfs < maxEps)
						isConforming = true;
					else
					{
						dist = currP.distance(p);
						isConforming = (dist/lfs < maxEps);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

	determinePointStates(maxEps);
	computeScale();
}

/*
 * generate point set by sampling with condition < maxReach, max error to original curve
 */
void TestReconstruct2D::sampleCurveByReach(float maxReach, float error, bool isClosed)
{
	int i, j;
	lfsCurvePoints.resize(curvePoints.size());
	mPoints.resize(curvePoints.size());
	computeLFSForCurve(curvePoints, lfsCurvePoints, mPoints, isClosed);
	computeKForCurve(curvePoints, kCurvePoints, extremaPoints);
	vector<int> curveIndices;
	i = 0;

	do
	{
		int prevI = i;
		Point prevP = curvePoints[prevI];
		curveIndices.push_back(prevI);
		points.push_back(prevP);
		i++;

		// test candidate sample if it conforms to the sampling condition
		bool isConforming = true;

		while ((i < (int)curvePoints.size()) && isConforming)
		{
			Point p = curvePoints[i];
			float reach = numeric_limits<float>::max();

			for (j = prevI; j <= i; j++)
			{
				float lfs = lfsCurvePoints[j];

				if (lfs < reach)
					reach = lfs;
			}

			j = prevI + 1;
			isConforming = true;

			// test reach and error conditions for each curve point between samples
			while ((j < i) && isConforming)
			{
				Point currP = curvePoints[j];

				// test error from curve (of chord prevP-p)
				if (error > 0.0)
					isConforming = (distancePFromLine(currP, prevP, p) < error);

				if (isConforming)
				{
					// test for reach condition (a sample within dist/reach < maxReach)
					// use prev/next points so that discrete interval of curve points does not impact
					float dist = curvePoints[j + 1].distance(prevP);

					if (dist/reach < maxReach)
						isConforming = true;
					else
					{
						dist = curvePoints[j - 1].distance(p);
						isConforming = (dist/reach < maxReach);
					}
				}

				j++;
			}

			if (isConforming)
				i++;
		}
	} while (i < (int)curvePoints.size());

	if (!isClosed)
		points.push_back(curvePoints[curvePoints.size() - 1]);

	determinePointStates(maxReach);
	computeScale();
}

TestReconstruct2D::TestReconstruct2D()
{
	viewWidth = 1;
	viewHeight = 1;
	scale = 1.0;
	offset[0] = 0.0;
	offset[1] = 0.0;
	lfsMin = 0.0;
	lfsMax = 0.0;
}

TestReconstruct2D::~TestReconstruct2D()
{
}

void TestReconstruct2D::reshape(int width, int height)
{
	this->viewWidth = width;
	this->viewHeight = height;
}

/*
 * transform point into unit square
 */
Point TestReconstruct2D::transform(Point p)
{
	int i;
	Point tp;

	for (i = 0; i < 2; i++)
		tp[i] = 0.05 + (p[i] - offset[i])/scale*0.9;	// add 5% border

	return tp;
}

void TestReconstruct2D::drawPoints(vector<Point> &points)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawPointsWithStates(vector<Point> &points, vector<PointState> &pointStates)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		glColor3f(0.0, 0.0, 0.0);	// black

		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawCurvePoints(vector<Point> &points, vector<float> &lfsCurvePoints, float lfsMin, float lfsMax)
{
	int i;

	glBegin(GL_POINTS);

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawEdges(vector<Point> &points, map<pair<int, int>, EdgeEnum> &edgeMap)
{
	int i, j;

	glBegin(GL_LINES);

	for (j = 1; j < 2; j++)
	{
		for (auto edgeItem:edgeMap)
			for (i = 0; i < 2; i++)
			{
				int index = ((i == 0) ? edgeItem.first.first : edgeItem.first.second);
				Point tp = transform(points[index]);
				glVertex2f(tp[0], tp[1]);
			}
	}

	glEnd();
}

void TestReconstruct2D::drawCircle(Circle &c)
{
	int i;
	glBegin(GL_LINE_LOOP);

	for (i = 0; i < 360; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a + cos(degInRad)*c.r, c.b + sin(degInRad)*c.r);
		Point tp = transform(p);
		glVertex2f(tp[0], tp[1]);
	}

	glEnd();
}

void TestReconstruct2D::drawArc(Circle &c, pair<int, int> arc)
{
	int i;
	glBegin(GL_LINES);
	Point oldTP;

	if (c.r == 0.0)
		return;

	if (arc.second < arc.first)
		arc.second += 360;

	for (i = arc.first; i < arc.second; i++)
	{
		float degInRad = i*3.14159/180;
		Point p(c.a - cos(degInRad)*c.r, c.b - sin(degInRad)*c.r);
		Point tp = transform(p);

		if (i != arc.first)
		{
			glVertex2f(oldTP[0], oldTP[1]);
			glVertex2f(tp[0], tp[1]);
		}

		oldTP = tp;
	}

	glEnd();
}

void TestReconstruct2D::drawCircles(vector<Circle> &circles)
{
	int i;

	for (i = 0; i < (int)circles.size(); i++)
			drawCircle(circles[i]);
}

void TestReconstruct2D::drawArcs(vector<Circle> &circles, vector<pair<int, int> > &arcs)
{
	int i;

	for (i = 0; i < (int)circles.size(); i++)
			drawArc(circles[i], arcs[i]);
}

void TestReconstruct2D::drawLabels(vector<Point> &points)
{
	int i, j;

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tp = transform(points[i]);
		glRasterPos2f(tp[0] + 0.01, tp[1] - 0.00);
		char str[5];
		sprintf(str, "%d", i);

		for (j = 0; j < (int)strlen(str); j++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, str[j]);
	}
}

void drawNormal(Point tp, Point n)
{
	glVertex2f(tp[0], tp[1]);
	Point endP = tp + n*0.05;
	glVertex2f(endP[0], endP[1]);

	// draw arrow hooks
	Point hookVec;
	float cs = cos(0.75*PI);
	float sn = sin(0.75*PI);
	hookVec[0] = 0.01*(cs*n[0] - sn*n[1]);
	hookVec[1] = 0.01*(sn*n[0] + cs*n[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] + hookVec[0], endP[1] + hookVec[1]);
	glVertex2f(endP[0], endP[1]);
	glVertex2f(endP[0] - hookVec[1], endP[1] + hookVec[0]);
}

void TestReconstruct2D::drawNormals(vector<Point> &points)
{
	int i;

	glBegin(GL_LINES);

	for (i = 0; i < (int)points.size(); i++)
	{
		// draw arrow line
		Point tp = transform(points[i]);
		drawNormal(tp, normals[i]);
	}

	glEnd();
}

/*
 * cubic solver: input the 4 coefficients c[4], output the 3 solutions sol[3] and return number of real roots
 */
int cubicSolve(double (&c)[4], double (&sol)[3])
{
    if (c[0] == 0)
    	assert("use quadratic solver");

    if (c[3] == 0)
        assert("One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken.");

    c[1] /= c[0];
    c[2] /= c[0];
    c[3] /= c[0];

    double q = (3.0*c[2] - SQR(c[1]))/9.0;
    double r = -(27.0*c[3]) + c[1]*(9.0*c[2] - 2.0*SQR(c[1]));
    r /= 54.0;
    double disc = q*q*q + r*r;
    double term1 = c[1]/3.0;

    if (disc > 0)
    {
    	// one root real, two are complex
        double s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
        double t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
        sol[0] = -term1 + s + t;

        return 1;
    }

    // the remaining options are all real
    if (disc == 0)
    {
    	// all roots real, at least two are equal.
        double r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
        sol[0] = -term1 + 2.0*r13;
        sol[1] = -(r13 + term1);

        if (sol[0] == sol[1])
        	return 1;
        else
        	return 2;
    }

    // only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    double dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    double r13 = 2.0*sqrt(q);
    sol[0] = -term1 + r13*cos(dum1/3.0);
    sol[1] = -term1 + r13*cos((dum1 + 2.0*PI)/3.0);
    sol[2] = -term1 + r13*cos((dum1 + 4.0*PI)/3.0);

    return 3;
}

const float A = -1.0, B = 1.0;	// weight function coefficients
const float W = 1.0;	// weight coefficient (larger -> more interpolating)

/*
 * compute implicit function value F(x) = sum_i w_i(x)f_i(x) for the input point (x, y)
 */
float computeImplicit(float x, float y, vector<Point> &points, vector<Point> &normals)
{
	int i;

	// w_i(x) = 1/(1+Wd²)(was: Ad²+B), d=|x, p_i|=sqrt((x-p_i(x))²+(y-p_i(y))²) - only positive weights are considered
	// f_i(x) = (x-p_i)n_i = n_i(y)y + n_i(x)(x-p_i(x))-p_i(y)

	// compute implicit function value
	float value = 0.0;

	for (i = 0; i < (int)points.size(); i++)
	{
		Point p = points[i];
		Point n = normals[i];
		float d = sqrt(SQR(x - p[0]) + SQR(y - p[1]));
		float w = 1.0/(1.0 + W*SQR(d));	// degressive curve with maximum 1 at d=0
		float f = (x - p[0])*n[0] + (y - p[1])*n[1];
		value += w*f;
	}

	return value;
}

/*
 * sign function
 */
template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

#ifdef OBSOLETE
/*
 * evaluate zero crossing of implicit function
 */
float solveImplicit(float x, vector<Point> &points, vector<Point> &normals)
{
	const float THRESHOLD = 0.01;
	float f[2], y[2] = { -5.0, 5.0 }, newY = -100.0, oldY;
	f[0] = computeImplicit(x, y[0], points, normals);
	f[1] = computeImplicit(x, y[1], points, normals);

	do
	{
		oldY = newY;
		newY = y[0] - f[0]/(f[1] - f[0])*(y[1] - y[0]);
		float newF = computeImplicit(x, newY, points, normals);

		if (sgn(newF) == sgn(f[0]))
		{
			y[0] = newY;
			f[0] = newF;
		}
		else
		{
			y[1] = newY;
			f[1] = newF;
		}
	} while (abs(oldY - newY) > THRESHOLD);

	return newY;
}

/*
 * intersect explicit cubic function for quadratic weighting of points on a baseline: F(x) = sum_i w_i(x)f_i(x) = 0
 */
float solveExplicit(float x, vector<Point> &points, vector<Point> &normals)
{
	int i;

	// NEW (buggy -> would require O(N^2) complexity for N points, by multiplying with all but one denominator)
	// F(x) = sum_i w_i(x)f_i(x)
	// w_i(x) = 1/(1+Wd²), d=|x, p_i|=sqrt((x-p_x)²+(y-p_y)²)
	// f_i(x) = (x-p_i)n_i = n_y*y + n_x*(x-p_x)-p_y

	// OLD (buggy):
	// w_i(x) = Ad²+B, d=|x, p_i|=sqrt((x-p_i(x))²+(y-p_i(y))²) ->0
	// w_i(x) = A((x-p_i(x))²+(y-p_i(y))²)+B = Ay² - 2Ap_i(y)y + A((x-p_i(x))²+p_i(y)²)+B
	// f_i(x) = (x-p_i)n_i = n_i(y)y + n_i(x)(x-p_i(x))-p_i(y)

	// compute coefficients for explicit cubic function y=x0x³+x1x²+x2x+x3: compute w_i(x)f_i(x) for all points and sum it up
	double c[4] = { 0.0, 0.0, 0.0, 0.0 };

	for (i = 0; i < (int)points.size(); i++)
	{
		Point p = points[i];
		Point n = normals[i];
		float w[3] = { W, -2*p[1]*W, W*(SQR(x) - 2*p[0]*x + SQR(p[0]) + SQR(p[1])) + 1 };
		float f[2] = { n[1], (x - p[0])*n[0] - p[1]*n[1] };

//		float w[3] = { A, -2*A*p[1], A*(SQR(x - p[0]) + SQR(p[1])) + B };	// w_i(x)
//		float f[2] = { n[1], n[0]*(x - p[0]) - p[1] };		// f_i(x)0

		// sum up w_i(x)f_i(x)
		c[0] += w[0]*f[0];
		c[1] += w[0]*f[1] + w[1]*f[0];
		c[2] += w[1]*f[1] + w[2]*f[0];
		c[3] += w[2]*f[1];
	}

	// solve cubic
	double s[3];
	int solutionCount = cubicSolve(c, s);
	assert(solutionCount == 1);

	return s[0];
}

/*
 * test explicit cubic function for quadratic weighting of points on a baseline
 */
void drawTest()
{
	int i;
/*
	// test cubic root solver
	float c[4] = { 1, -2, -3, 4 }, s[3];
	int count = cubicSolve(c, s);
*/

	vector<Point> points, normals;
	points.push_back(Point(1.0, 1.0));
/*
	points.push_back(Point(2.0, 3.0));
	points.push_back(Point(4.0, 2.0));
	points.push_back(Point(5.0, 1.0));
	points.push_back(Point(7.0, 2.0));
*/
	normals.push_back(Point(-0.2, 1.0));
	normals.push_back(Point(0.0, 1.0));
	normals.push_back(Point(0.4, 1.0));
	normals.push_back(Point(-0.1, 1.0));
	normals.push_back(Point(-0.3, 1.0));

	for (i = 0; i < (int)points.size(); i++)
		normals[i].normalize();

	// test
	computeImplicit(2.0, 2.0, points, normals);

	// test
// 	float testy1 = solveExplicit(2.0, points, normals);
// 	float testy2 = solveImplicit(2.0, points, normals);

	glColor3f(0.0, 0.0, 0.0);	// black
	glBegin(GL_POINTS);

	for (auto p:points)
		glVertex2f(p[0]*0.05, (5.0 - p[1])*0.05);

	glEnd();

	glBegin(GL_LINES);

	for (i = 0; i < (int)points.size(); i++)
	{
		Point tn = normals[i];
		tn[1] = -tn[1];
		Point tp(points[i][0]*0.05, (5.0 - points[i][1])*0.05);
		drawNormal(tp, tn);
	}

	glEnd();

	glColor3f(1.0, 0.0, 0.0);	// red
	glBegin(GL_LINE_STRIP);

	for (i = 0; i < 100; i++)
	{
		float x = i*0.1;
		float y = solveExplicit(x, points, normals);
//		float y = solveImplicit(x, points, normals);
		glVertex2f(x*0.05, (5.0 - y)*0.05);
	}

	glEnd();

	glBegin(GL_LINES);
	glVertex2f(0.0, 0.25);
	glVertex2f(0.5, 0.25);
	glEnd();

	// second baseline: 22.5 degrees rotated about its origin, translated by (0.0, -2)
	glColor3f(0.0, 0.0, 1.0);	// blue

	float sinT = sin(-PI/8), cosT = cos(-PI/8), ofsX = 0.0, ofsY = 0.0;

	glBegin(GL_LINES);
	glVertex2f(ofsX, (5.0 - ofsY)*0.05);
	glVertex2f((10.0*cos(PI/8) + ofsX)*0.05, (5.0 - (10.0*sin(PI/8) + ofsY))*0.05);
	glEnd();

	for (i = 0; i < (int)points.size(); i++)
	{
		// transform point
		Point tp(cosT*points[i][0] - sinT*points[i][1] - ofsX, sinT*points[i][0] + cosT*points[i][1] - ofsY);
		points[i] = tp;

		// transform normal
		Point tn(cosT*normals[i][0] - sinT*normals[i][1], sinT*normals[i][0] + cosT*normals[i][1]);
		normals[i] = tn;
	}

	glBegin(GL_LINE_STRIP);

	for (i = 0; i < 100; i++)
	{
		float x = i*0.1;
		float y = solveImplicit(x, points, normals);
		float tx = cos(PI/8)*x - sin(PI/8)*y;
		float ty = sin(PI/8)*x + cos(PI/8)*y;
		glVertex2f(tx*0.05, (5.0 - ty)*0.05);
	}

	glEnd();
}
#endif

/*
 * write screen to PNG file
 */
void writePNG(int width, int height, string fileName)
{
	int x, y, npixels = width*height;
	GLfloat* pixels = new GLfloat[npixels*3];
	glReadPixels(0.0, 0.0, width, height, GL_RGB, GL_FLOAT, pixels);
	png::image<png::rgb_pixel> image(width, height);
	double R, G, B;

	for (y = height - 1; y >= 0; y--)
		for (x = 0; x < width; x++)
		{
			R = *pixels++;
			G = *pixels++;
			B = *pixels++;
			image[y][x] = png::rgb_pixel(255*R, 255*G, 255*B);     // set pixel to position (x, y)
		}

	image.write(fileName);
}

void TestReconstruct2D::draw(string fileName, bool isDrawLines)
{
	// clear screen to white
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	// draw arcs
	glColor3f(0.0, 1.0, 0.0);	// green
	glLineWidth(1.0);
	drawArcs(circles, arcs);

	// draw discs at curve points
	glColor3f(1.0, 1.0, 0.8);	// very light yellow
	float noiseBandwidth = 0.5;
	glPointSize(noiseBandwidth/scale*0.9*viewWidth);
	glEnable(GL_POINT_SMOOTH);
//	drawPoints(curvePoints);

	// draw curve points
	glColor3f(0.0, 0.0, 0.0);	// black
	glPointSize(1.0);
	glDisable(GL_POINT_SMOOTH);
	drawCurvePoints(curvePoints, lfsCurvePoints, lfsMin, lfsMax);

	// draw medial axis points
	glColor3f(0.0, 0.0, 0.0);	// black
	glPointSize(1.0);
	glDisable(GL_POINT_SMOOTH);
//	drawPoints(mPoints);

	// draw curvature extrema
	glColor3f(0.5, 0.0, 0.5);	// lilac
	glPointSize(8.0);
	glEnable(GL_POINT_SMOOTH);
//	drawPoints(extremaPoints);

	if (isDrawLines)
	{
		// draw lines
		glColor3f(1.0, 0.0, 0.0);	// red
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
		drawEdges(projPoints, edgeMap);
	}

	// draw points
	glPointSize(16.0);
	glEnable(GL_POINT_SMOOTH);

	if (pointStates.size() > 0)
		drawPointsWithStates(points, pointStates);

	// draw averaged points
	glColor3f(0.0, 0.0, 1.0);	// blue
//	drawPoints(projPoints);

	// draw normals
	glColor3f(0.0, 0.0, 0.0);	// black
	glLineWidth(1.0);
//	drawNormals(projPoints);

	// draw labels
	glColor3f(0.0, 0.0, 0.0);	// black
//	drawLabels(points);

//	drawTest();

	writePNG(viewWidth, viewHeight, fileName);
}

static TestReconstruct2D *instance = NULL;

///////////////////////////// GL FUNCTIONS /////////////////////////////////

// GLUT idle function
void idle()
{
}

// GLUT display function
void display(string fileName, bool isDrawLines)
{
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);
	glViewport(0, 0, viewWidth, viewHeight);
	instance->draw(fileName, isDrawLines);

	// restore the stored viewport dimensions
	glViewport(vp[0], vp[1], vp[2], vp[3]);
	glutSwapBuffers();
}

// GLUT reshape function
void reshape(int newWidth, int newHeight)
{
	if (newWidth == 0)
		newWidth = 1;

	if (newHeight == 0)
		newHeight = 1;

	viewWidth = newWidth;
	viewHeight = newHeight;

	glViewport(0, 0, viewWidth, viewHeight);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, 1, 1, 0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// update application
	instance->reshape(viewWidth, viewHeight);
}

void mouseMotionFunc(int x, int y)
{
}

void mouseButtonFunc(int b, int s, int x, int y)
{
	glutPostRedisplay();
}

void processNormalKeys(unsigned char key, int x, int y)
{
	if (key == 27)
		exit(0);
	else
		cout << "key " << (int)key << " pressed, not handled" << endl;
}

void processSpecialKeys(int key, int x, int y)
{
}

string teaser = "m 459.41947,59.854719 c 73.85195,-20.931012 115.59266,1.472244 128.48931,38.923482 10.04049,29.157109 -10.47593,80.018759 -65.49878,95.369409 -50.87306,14.19292 -117.60902,-12.28741 -167.84646,6.656 -48.81755,18.40801 -90.57522,54.8466 -117.17769,104.55079 -41.15894,76.90154 -121.64283,100.91335 -145.461964,64.14468 -26.21363,-40.46495 26.458014,-62.9026 49.497464,-73.23606 42.35902,-18.9985 82.29881,-56.17798 75.76144,-82.32743 -5.05077,-20.20305 7.60004,-39.6756 -103.16379,-90.78029 -42.34202,-19.53594 -47.179322,-83.122834 52.829,-61.484748 71.07716,15.378456 73.53591,64.130198 110.9465,69.361498 53.08911,7.42371 107.77302,-50.246318 181.62497,-71.177331";
string cat = "m -359.20637,154.76853 c 1.31737,20.90141 -22.5042,19.16747 -41.04941,22.90222 -36.84896,7.42087 -7.86444,80.19238 -9.25282,111.85808 -2.12028,48.35848 -28.61898,61.69306 -27.94568,74.52182 0.6733,12.82876 8.9512,41.36747 10.55726,53.40731 5.60901,42.04793 -3.10508,42.22903 -1.86305,59.61745 1.24203,17.38843 18.32402,31.88623 24.21959,55.27036 6.72462,26.67244 1.55865,69.6502 -21.11451,64.58558 -18.25714,-4.07819 0.684,-14.08244 -0.80291,-34.97167 -1.52671,-21.44846 -53.77362,-66.03567 -60.05658,-72.46396 -9.61626,-9.83868 -59.14815,116.13232 -97.80076,98.87918 -16.4618,-7.34797 14.25215,-19.90371 21.08331,-29.35669 9.27355,-12.83277 8.40577,-79.45873 4.05867,-82.56381 -4.34711,-3.10508 -5.45141,8.96877 -50.96738,8.95145 -30.14494,-0.0115 -51.50013,-14.54059 -55.84724,-12.67754 -4.3471,1.86304 18.76438,123.84914 -20.4935,124.82405 -12.96318,0.32192 -14.90436,-2.48406 -16.14639,-11.17827 -1.24203,-8.69421 11.39897,-7.14505 12.52685,-17.17533 1.58967,-14.13693 -24.32614,-82.1871 -31.15731,-82.1871 -6.83116,0 -19.40215,19.76594 -31.67177,27.94568 -14.90437,9.93624 -39.12396,16.76741 -48.43919,26.70366 -9.31523,9.93624 -6.24026,38.90074 -24.84061,39.12395 -22.00738,0.26409 -10.55725,-37.88192 4.34711,-57.1334 14.90437,-19.25147 64.58558,-54.64933 62.72254,-58.99644 -1.86305,-4.34711 -8.69422,-2.48406 -22.97757,-14.90436 -14.28335,-12.42031 -45.33114,-68.90534 -50.92324,-68.31168 -4.49133,0.4768 -5.58914,8.69422 -16.1464,8.0732 -10.55726,-0.62101 -16.76741,-4.96812 -18.63045,-8.0732 -1.86305,-3.10507 3.10507,-12.4203 3.10507,-19.25147 0,-6.83116 -0.83634,-15.34496 0.62102,-24.8406 1.78298,-11.61726 12.4203,-15.52538 12.4203,-20.49351 0,-2.48406 -6.49401,-3.45136 -11.61739,-13.0771 -4.37675,-8.22299 -8.75093,-19.0134 -6.42326,-21.51785 4.23313,-4.55462 19.97645,11.836 22.38776,9.75435 2.66394,-2.29975 -13.19657,-15.37013 -4.96812,-18.63046 1.44336,-0.5719 11.38224,2.15149 20.82608,7.09525 9.49993,4.97313 16.43483,15.2613 20.16092,16.50333 7.45218,2.48406 16.7674,0 28.56669,6.21015 11.79929,6.21016 19.54274,21.8625 50.30224,32.91381 51.85476,18.63045 140.57649,-10.41292 203.69298,-9.93625 63.11649,0.47667 70.64811,10.92575 79.48995,-18.00944 12.74644,-41.71311 -20.14572,-100.99802 -7.11961,-130.59508 10.04418,-22.82171 37.54935,-43.91019 61.76895,-43.91019 24.21959,0 34.77685,9.93624 35.39786,21.11452";
string parallel = "m 212.84016,42.209659 c 89.05973,0 2190.85664,-0.0947 2292.84504,-0.0947 101.9884,0 168.9495,89.872011 168.9495,176.871401 0,86.99938 -75.0613,168.91644 -168.948,168.91644 -44.4735,0 -2237.03561,0.0316 -2290.06862,0.0316 -106.06602,0 -171.725931,-90.66119 -171.725931,-177.78685 0,-87.12566 76.691511,-167.937865 168.948011,-167.937861";
string bunny = "m -755.54157,439.19348 c -2.81588,-35.3868 -73.42744,-49.1442 -84.52047,-72.1131 -16.34716,-33.84797 2.26058,-62.71611 17.45083,-81.27764 14.10537,-17.23587 13.61005,-19.65993 13.66234,-70.79573 0.0523,-51.13581 4.00051,-61.97973 16.1464,-62.10152 18.06052,-0.1811 11.86373,29.49677 12.70874,59.17833 1.04073,36.55644 -6.06677,54.03577 2.78541,55.27036 6.26796,0.87418 6.94735,-22.5034 11.8305,-59.79935 3.49259,-26.67529 5.60268,-54.70426 21.11452,-52.16528 15.83216,2.59141 15.67466,26.3087 8.40577,56.15545 -8.6868,35.66883 -11.40314,65.14933 10.84569,78.60485 46.36972,28.0432 87.88088,-40.45104 156.49582,-9.93625 51.81346,23.04275 60.58667,55.5695 62.10153,73.90081 4.46432,54.02268 -7.29574,55.14578 -1.24203,73.9008 6.05371,18.75502 19.00256,11.9741 19.25148,36.63989 0.40003,39.6392 -52.42394,41.64734 -156.16325,40.1841 -126.77474,-1.78816 -149.40364,-0.43075 -149.37621,-23.41669 0.0333,-27.93684 40.95673,-11.39249 38.50293,-42.22903";

void createFigure1a()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(teaser);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(1.0/3.0, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure1a.png", true);
}

void createFigure1b()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(teaser);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(0.47, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure1b.png", true);
}

void createFigure1c()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(teaser);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByReach(0.9, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure1c.png", true);
}

void createFigure3a()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("crocodile-conf.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure3a.png", false);
}

void createFigure3c()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("crocodile-conf.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure3c.png", true);
}

void createFigure4a()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("ventilator.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure4a.png", false);
}

void createFigure4c()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("ventilator.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure4c.png", true);
}

void createFigure9aleft()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(parallel);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(1.0/3.0, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9aleft.png", true);
}

void createFigure9aright()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(parallel);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByReach(0.9, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9aright.png", true);
}

void createFigure9bleft()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(bunny);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(1.0/3.0, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9bleft.png", true);
}

void createFigure9bright()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(bunny);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByReach(0.9, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9bright.png", true);
}

void createFigure9cleft()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(cat);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(1.0/3.0, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9cleft.png", true);
}

void createFigure9cright()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(cat);
	instance->scaleOfCurvePoints();
	instance->sampleCurveByReach(0.9, 0.0, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure9cright.png", true);
}

void createFigure10a()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(bunny);
	float maxDim = instance->scaleOfCurvePoints();
	instance->sampleCurveByEps(1.0/3.0, 0.01*maxDim, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure10a.png", true);
}

void createFigure10b()
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(bunny);
	float maxDim = instance->scaleOfCurvePoints();
	instance->sampleCurveByReach(0.9, 0.01*maxDim, true);
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure10b.png", true);
}

void createFigure11b()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("ventilator.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure11b.png", true);
}

void createFigure11d()
{
	instance = new TestReconstruct2D();
	instance->loadPointSet("mushroom.p2d");
	instance->scaleOfCurvePoints();
	instance->reconstruct(true);
	reshape(viewWidth, viewHeight);
	display("figure11d.png", true);
}

/*
 * minimum angle between 3 samples for given rho in degrees
 */
float rhoAngle(float rho)
{
	return (PI - 4*asin(rho/2))/PI*180.0; // Lemma 8
}

/*
 * minimum angle between 3 samples for given eps in degrees
 */
float epsAngle(float eps)
{
	return (PI - 4*asin(eps/(2 - 2*eps)))/PI*180.0;	// Lemma 9
}

void createTable1()
{
	int i;
	float bounds[] = { 0.5, 0.2, 1.0/3.0, 0.4, 0.48, 0.47, 0.9 };

	cout << "Table 1:" << endl;
	cout << "Bound\t\t| min alpha | circle | par" << endl;
	cout << "==========================================" << endl;

	for (i = 0; i < 7; i++)
	{
		float bound = bounds[i];
		float angle;

		if (i == 0)
			angle = 150.0;
		else
			angle = rhoAngle(bound);

		float circleSamples = 360.0/(180.0 - angle);
		float parallelSamples = 1.0/bound;

		cout << ((i == 6) ? "rho < " : "eps < ") << bound << "\t| >" << angle << "° | " << circleSamples << " | " << parallelSamples << endl;
	}
}

int countSamples(string fileName, bool isRho, float bound, float error)
{
	instance = new TestReconstruct2D();
	instance->generateCurveFromBezierString(fileName);
	float maxDim = instance->scaleOfCurvePoints();

	if (isRho)
		instance->sampleCurveByReach(bound, error*maxDim, true);
	else
		instance->sampleCurveByEps(bound, error*maxDim, true);

	instance->reconstruct(true);

	return instance->sampleCount;
}

void createTable2()
{
	int i;
	cout << "Table 2:" << endl;
	cout << "Model\t\t| rho < 0.9\t| eps < 0.47\t| e < 0.33" << endl;
	cout << "==========================================================" << endl;

	string models[4] = { parallel, teaser, bunny, cat };
	string modelNames[4] = { "Parallel", "Teaser  ", "Bunny   ", "Cat     " };

	for (i = 0; i < 4; i++)
	{
		float countRho09 = countSamples(models[i], true, 0.9, 0.0);
		float countEps047 = countSamples(models[i], false, 0.47, 0.0);
		float countEps033 = countSamples(models[i], false, 1.0/3.0, 0.0);
		cout << modelNames[i] << "\t| " << countRho09 << "\t\t| " <<
				countEps047 << "(" << (countEps047 - countRho09)/countRho09*100 << "%)\t| " <<
				countEps033 << "(" << (countEps033 - countRho09)/countRho09*100 << "%)\t" << endl;
	}
}

void createTable3()
{
	int i;
	cout << "Table 3:" << endl;
	cout << "Hausdorff dist\t| rho < 0.9\t| eps < 0.47\t| e < 0.33" << endl;
	cout << "==========================================================" << endl;

	float bounds[] = { 0.0, 0.01, 0.003, 0.001, 0.0003 };

	for (i = 0; i < 5; i++)
	{
		float bound = bounds[i];
		float countRho09 = countSamples(bunny, true, 0.9, bound);
		float countEps047 = countSamples(bunny, false, 0.47, bound);
		float countEps033 = countSamples(bunny, false, 1.0/3.0, bound);

		if (i == 0)
			cout << "inf";
		else
			cout << bound << "%";

		cout << "\t\t| " << countRho09 << "\t\t| " <<
				countEps047 << "(" << (countEps047 - countRho09)/countRho09*100 << "%)\t| " <<
				countEps033 << "(" << (countEps033 - countRho09)/countRho09*100 << "%)\t" << endl;
	}
}

void printUsage(char *argv[])
{
	cout << "Usage: " << argv[0] << " fig1a|fig1b|fig1c|fig3a|fig3c|fig4a|fig4c|" <<
		"fig9aleft|fig9aright|fig9bleft|fig9bright|fig9cleft|fig9cright" <<
		"fig10a|fig10b|fig11b|fig11d|tab1|tab2|tab3" << endl;
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitWindowSize(viewWidth, viewHeight);
	glutCreateWindow("Reconstruct from 2D Points - SGP Reproducibility Stamp");

	if (argc <= 1)
		printUsage(argv);
	else
	{
		string a = argv[1];

		if (!a.compare("fig1a"))
			createFigure1a();
		else
		if (!a.compare("fig1b"))
			createFigure1b();
		else
		if (!a.compare("fig1c"))
			createFigure1c();
		else
		if (!a.compare("fig3a"))
			createFigure3a();
		else
		if (!a.compare("fig3c"))
			createFigure3c();
		else
		if (!a.compare("fig4a"))
			createFigure4a();
		else
		if (!a.compare("fig4c"))
			createFigure4c();
		else
		if (!a.compare("fig9aleft"))
			createFigure9aleft();
		else
		if (!a.compare("fig9aright"))
			createFigure9aright();
		else
		if (!a.compare("fig9bleft"))
			createFigure9bleft();
		else
		if (!a.compare("fig9bright"))
			createFigure9bright();
		else
		if (!a.compare("fig9cleft"))
			createFigure9cleft();
		else
		if (!a.compare("fig9cright"))
			createFigure9cright();
		else
		if (!a.compare("fig10a"))
			createFigure10a();
		else
		if (!a.compare("fig10b"))
			createFigure10b();
		else
		if (!a.compare("fig11b"))
			createFigure11b();
		else
		if (!a.compare("fig11d"))
			createFigure11d();
		else
		if (!a.compare("tab1"))
			createTable1();
		else
		if (!a.compare("tab2"))
			createTable2();
		else
		if (!a.compare("tab3"))
			createTable3();
		else
		{
			cout << "No valid parameter chosen" << endl;
			printUsage(argv);
		}
	}

	return 0;
}
