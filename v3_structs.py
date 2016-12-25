from __future__ import print_function
import os,sys
from sympy.geometry import *
import heapq
import sympy
import numpy as np
import matplotlib.pyplot as plt

'''	Supporting structures needed for Voronoi Diagram of Points
'''

'''	The class Site is neatly abstracted as follows:
	In case of Voronoi of Line Segments, Sites can either be points or Line Segments
	In case of Voronoi of Points, Sites are points.
'''
class Site:
	''' pt2=None => Site is a point. '''
	def __init__(self, pt1, pt2=None, isUpperEndpoint=True, isInterior=False):
		if pt2 is None:
			self.points = sorted([pt1, pt1], lambda p: (-p.y, p.x))
			self.isPoint= True
		else:
			self.points = sorted([pt1, pt2], lambda p: (-p.y, p.x))
			self.isPoint = (pt1==pt2)
		self.isUpperEndPoint=isUpperEndpoint
		if not isInterior:
			self.isPoint=True
		self.isInterior=isInterior
	def point(self):
		return self.points[0] if self.isUpperEndpoint else self.points[1]
	def line(self):
		return Line(self.points[0], self.points[1])
	def addToPlot(self):
		if self.isPoint:
			pt = self.points[0 if self.isUpperEndpoint else 1]
			plt.scatter( [pt.x], [pt.y] )
		else:
			plt.plot([self.points[0].x, self.points[1].x],[self.points[0].y, self.points[1].y])

''' An Arc is the set of points that are equidistant to a site and the sweepline.
	The set of points on an Arc are also closest to the site that induces the arc and nothing else.
	When two Arcs intersect, that intersection point is equidistant from the sites of the two Arcs, and also equidistant from the sweepline.
	Arcs are eliminated by circleEvents which can happen whenever.

	Notice that a site can induce multiple Arcs. Each arc would have a different circleEvent
'''
class Arc:
	def __init__(self, site):
		self.site=site
		self.circleEvent=None
		self.isBreakpoint=False

''' A breakpoint is the intersection of two arcs.
	Notice that the order of the two arcs is important.
	Hence, we represent a Breakpoint as an ordered pair of two Sites

	As per Fortune's algorithm, we can identify a few types of Breakpoints:
		1. Between two "site endpoints" - results in locus being a line segment.
			=> Here, "site endpoint" can also be a normal point
		2. Between two "site interiors" - results in locus being a line segment.
		3. Between "site endpoint" and "site interior" - locus is a parabola.
			=> Here, "site endpoint" can also be a normal point
		4. Closest to a "site endpoint" - locus is a line segment perpendicular the the "site interior"
			=> NOTE Here, the "site endpoint" CANNOT be normal point
		5. Intersection of "site interior" and sweepline
'''
class Breakpoint:
	def __init__(self, leftSite, rightSite):
		self.sites =[leftSite, rightSite]
		self.voronoiEdge=None
		self.isBreakpoint=True

	''' The exact position of a Breakpoint depends on the position of the sweepline.
	'''
	def eval(self, sweeplinePosition):
		if self.sites[0].isPoint and self.sites[1].isPoint:
			x0 = self.sites[0].point().x
			x1 = self.sites[1].point().x
			y0 = self.sites[0].point().y
			y1 = self.sites[1].point().y
			ys = sweeplinePosition
			# Some quadractic equations.
			a2c = (y1 - y0)
			a1c = ((y1-ys)*(-2*x0)) - ((y0-ys)*(-2*x1))
			a0c = ((y1-ys)*(x0*x0)) + (((y0*y0) - (ys*ys))*(y1-ys)) - ((y0-ys)*(x1*x1)) - (((y1*y1)-(ys*ys))*(y0-ys))
			# A little bit of geometry trickery.
			if a2c.evalf()==0:
				aOption= ((x0+x1)/2.0) 
			else:
				aOption= ((-a1c + (sympy.sqrt(a1c*a1c - 4*a2c*a0c)))/(2*a2c))
			return Triangle(Point(x0,y0), Point(x1,y1), Point(aOption,ys)).circumcenter
		else:
			raise NotImplementedError # TODO ... handle other cases
	def locus(self):
		if self.sites[0].isPoint and self.sites[1].isPoint: # The locus becomes the perpendicular bisector of the two points.
			return Segment(self.sites[0].point(), self.sites[1].point()).perpendicular_bisector()
		else:
			raise NotImplementedError # TODO ... handle other cases 
	def type(self):
		if self.sites[0].isPoint and self.sites[1].isPoint:
			return 1
		else:
			return 0
	def convergence(self, other):
		if self.type()==1 and other.type()==1:
			return self.locus().intersection(other.locus())
		else:
			raise NotImplementedError


''' In fortune's algorithm, two kinds of events happen
	SiteEvents, wherein a new Arc is created upon inserting a new Site
	CircleEvent, wherein an Arc is deleted because two breakpoints converged...i.e., their loci intersect
'''
class Event:
	def __init__(self, eventSite):
		self.site=eventSite
		self.isSiteEvent=True
		# fields needed for circleEvents
		self.isFalseAlarm=False
		self.arc=None
		self.breakpoints=None
		self.center=None
	def site(self):
		return self.site
	def circle(self):
		return Circle(self.center, abs(self.center.y - self.site.point().y) )
''' A quick function to make a circle event. Also ensures that the targetArc's circle event is set correctly.

	@params atPoint the sympy.geometry.Point at which the circleEvent should happen. only care about the y-coordinate, really.
	@params withCenter the center of the circle event. keeping it for miscellaneous reasons.
	@params forBreakpoints the list of Breakpoints (pair of Breakpoints) which converged ... in order from left to right.
	@params targetArc the arc which will be removed by the circleEvent
		=> imposed ASSERT: the targetArc is between forBreakpoints[0] and forBreakpoints[1]
	@returns event The corresponding circle event.
'''
def makeCircleEvent(atPoint, withCenter, forBreakpoints, targetArc):
	assert (len(forBreakpoints)==2),"Circle Events must be caused by convergence of two breakpoints."
	# assert ((forBreakpoints[0].rightArc==targetArc) and (forBreakpoints[1].leftArc==targetArc)), "The targetArc isn't matching up with the breakpoints provided."
	assert ((targetArc.circleEvent is None) or (targetArc.circleEvent.isFalseAlarm)), "An Arc can't have two circle events."
	event = Event(Site(atPoint))
	event.arc=targetArc
	targetArc.circleEvent=event
	event.center=withCenter
	event.breakpoints=forBreakpoints
	event.isSiteEvent=False
	return event


''' A VoronoiEdge is effectively the path traced by a Breakpoint
	in the final VoronoiDiagram
'''
class VoronoiEdge:
	def __init__(self, breakpoint, sweeplinePositionStart):
		self.breakpoint=breakpoint
		assert (breakpoint.voronoiEdge is None),"A breakpoint can't have two Voronoi Edges mate."
		self.start=sweeplinePositionStart
		self.end=None
	def close(self, sweeplinePositionEnd):
		self.end=sweeplinePositionEnd
	def addToPlot(self):
		if breakpoint.type()==1:
			if self.end is not None:
				start=breakpoint.eval(self.start)
				end=breakpoint.eval(self.end)
				plt.plot([start.x, end.x],[start.y, end.y])
			else:
				raise NotImplementedError
		else:
			raise NotImplementedError

'''	A VoronoiDiagram is a collection of VoronoiEdges and sites
	Ideally, it'd be a DCEL. But for now, this will do.
'''
class VoronoiDiagram:
	def __init__(self):
		self.sites = []
		self.edges = []
	def plot(self):
		for s in self.sites:
			s.addToPlot()
		for e in edges:
			e.addToPlot()
		plt.show()
	def addSite(s): # @params s is a Site
		self.sites.append(s)
	def addEdge(e): # @params e is a VoronoiEdge
		self.edges.append(e)


''' TODO create a class Beachline.
	it should be a balanced binary search tree which supports the following:
		stores Breakpoints and Arcs.
		insertion
		deletion
		replacement by subtree
		finding the neighbouring arcs/breakpoints

	It's a TODO because for now, we live with a plain old list.
'''
class Beachline:
	def __init__(self, breakpoint=None, arc=None):
		assert (breakpoint is None or arc is None),"A Beachline node cannot both be a breakpoint and an arc"
		self.par = None
		self.left= None
		self.right=None
		self.size=1
		self.breakpoint=breakpoint
		self.arc=arc

	def isArc(self):
		return ((self.left is None) and (self.right is None))

	def isBreakpoint(self):
		return (not ((self.left is None) and (self.right is None)))

	def search(self, linesweepPosition, point):
		if self.isArc():
			return self
		else:
			loc = self.breakpoint.eval(linesweepPosition)
			if loc.x<point.x:
				return self.right.search(linesweepPosition, point)
			else:
				return self.left.search(linesweepPosition, point)

	# Need a function that lets you get the leftward arc and the rightward arc