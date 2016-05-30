/*
	Copyright 2015, 2016 Stefan Ohrhallinger, except where otherwise noted.
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

#include <list>
#include <map>
#include "Reconstruct2D.h"

using namespace std;

const int K_MAX = 15;

Reconstruct2D::Reconstruct2D(const vector<Point> &p_points, bool p_isClosed)
{
	int i;

	points = p_points;
	isClosed = p_isClosed;

	for (i = 0; i < (int)points.size(); i++)
		pClasses.push_back(CONFORMING);
}

/*
 * determine kNN
 */
void Reconstruct2D::determineNeighbors(ANNkd_tree *kdTree, ANNpointArray ann_points, int kMax,
		vector<pair<int, int> > &manifoldNeighbors)
{
	int i;

	// determine feature size
	for (i = 0; i < (int)points.size(); i++)
	{
		Point currP = points[i];
		ANNidxArray nnIdx = new ANNidx[kMax];
		ANNdistArray distances = new ANNdist[kMax];
		kdTree->annkSearch(ann_points[i], kMax, nnIdx, distances);
		int nn1 = nnIdx[1], nn2 = -1;
		manifoldNeighbors.push_back(make_pair(nn1, -1));
		bool featureFitted = false;

		// for open curves, do not connect their end points
		if (isClosed || ((i > 0) && (i < (int)points.size() - 1)))
		{
			// add neighbors until feature condition fulfilled
			int k = 2;

			do
			{
				// check feature size
				nn2 = nnIdx[k];
				float nnDist = points[nn1].distance(points[nn2]);
				featureFitted = ((nnDist > currP.distance(points[nn1])) && (nnDist > currP.distance(points[nn2])));
				k++;
			} while (!featureFitted && (k < kMax));

			if (featureFitted)
				manifoldNeighbors[i].second = nn2;
		}

		if (!featureFitted)
			pClasses[i] = NONCONFORMING;

		delete nnIdx;
		delete distances;
	}
}

/*
 * reconstruct points assuming no noise
 */
void Reconstruct2D::reconstruct()
{
	int i, j, kMax = (K_MAX > points.size()) ? points.size() : K_MAX;
	vector<pair<int, int> > manifoldNeighbors;
	vector<pair<vector<int>, vector<int> > > neighbors(points.size());
	ANNkd_tree *kdTree = NULL;
	ANNpointArray ann_points;

	ann_points = annAllocPts(points.size(), 2);

	for(i = 0; i < (int)points.size(); i++)
	{
		auto p = ann_points[i];
		p[0] = points[i][0];
		p[1] = points[i][1];
	}

	kdTree = new ANNkd_tree(ann_points, points.size(), 2);

	determineNeighbors(kdTree, ann_points, kMax, manifoldNeighbors);

	// check local consistency
	for (i = 0; i < (int)points.size(); i++)
	{
		int nn[2] = { manifoldNeighbors[i].first, manifoldNeighbors[i].second };

		for (j = 0; j < 2; j++)
		{
			if (nn[j] != -1)
			{
				pair<int, int> edge(nn[j], -1);
				bool isConsistent = true;

				if (i == manifoldNeighbors[nn[j]].first)
					edge.second = manifoldNeighbors[nn[j]].first;
				else
				if (i == manifoldNeighbors[nn[j]].second)
					edge.second = manifoldNeighbors[nn[j]].second;
				else
					isConsistent = false;

				// only insert consistent edges
				if (isConsistent)
				{
					if (edge.first > edge.second)
						swap(edge.first, edge.second);

					visEdgeMap[edge] = BIJECTIVE;
				}
			}
		}
	}

	// normals and connected components not necessary for visual reconstruction:

	// orient direction between manifold neighbors: follow and flip until all points handled
	bool isOriented[points.size()];

	for (i = 0; i < (int)points.size(); i++)
		isOriented[i] = false;

	// DEBUG
	int ccCount = 0;

	for (i = 0; i < (int)points.size(); i++)
		if (!isOriented[i] && (pClasses[i] == CONFORMING))
		{
			isOriented[i] = true;

			// DEBUG
			ccCount++;
			int ccSize = 1;

			// follow both incident edges until leaf vertex or reaching back to it and orient consistently
			for (j = 0; j < 2; j++)
			{
				int prevI = -1, nextI = i;
				bool isFinished = false;

				do
				{
					prevI = nextI;
					nextI = (j == 0) ? manifoldNeighbors[nextI].first : manifoldNeighbors[nextI].second;

					if ((nextI != -1) && (pClasses[nextI] == CONFORMING) && (!isOriented[nextI]))
					{
						int n = (j == 0) ? manifoldNeighbors[nextI].first : manifoldNeighbors[nextI].second;

						if (n == prevI)
						{
							swap(manifoldNeighbors[nextI].first, manifoldNeighbors[nextI].second);
							swap(neighbors[nextI].first, neighbors[nextI].second);
						}

						isOriented[nextI] = true;

						// DEBUG
						ccSize++;
					}
					else
						isFinished = true;
				} while (!isFinished);
			}
		}

	// compute normals
	normals.resize(points.size());

	for (i = 0; i < (int)points.size(); i++)
		if (pClasses[i] == CONFORMING)
		{
			Point tangent = points[manifoldNeighbors[i].first] - points[manifoldNeighbors[i].second];
			normals[i] = Point(-tangent[1], tangent[0]);
			normals[i].normalize();
		}
		else
			normals[i] = Point(0.0, 0.0);

	projPoints = points;
}

map<pair<int, int>, EdgeEnum> Reconstruct2D::getEdgeMap()
{
	return visEdgeMap;
}

vector<Point> Reconstruct2D::getProjectedPoints()
{
	return projPoints;
}

vector<Point> Reconstruct2D::getNormals()
{
	return normals;
}

vector<Circle> Reconstruct2D::getCircles()
{
	return circles;
}

vector<pair<int, int> > Reconstruct2D::getArcs()
{
	return arcs;
}

vector<PointsEnum> Reconstruct2D::getPointClassification()
{
	return pClasses;
}
