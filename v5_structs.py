from __future__ import print_function
import os,sys
from sympy.geometry import *
import heapq
import sympy
import numpy as np
import matplotlib.pyplot as plt

''' Data structures for Fortune's Algorithm
	We're using the algorithm on the wikipedia page.
'''

'''	The class Site is neatly abstracted as follows:
	In case of Voronoi of Line Segments, Sites can either be points or Line Segments
	In case of Voronoi of Points, Sites are points.
'''
class Site:
	''' pt2=None => Site is a point. '''
	def __init__(self, pt1, pt2=None, isUpperEndpoint=True, isInterior=False):
		if pt2 is None:
			self.points = sorted([pt1, pt1], key=lambda p: (-p.y, p.x))
			self.isPoint= True
		else:
			self.points = sorted([pt1, pt2], key=lambda p: (-p.y, p.x))
			self.isPoint = (pt1==pt2)
		self.isUpperEndpoint=isUpperEndpoint
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
	def eval(self, ys):
		if self.isPoint:
			return self.point()
		else:
			raise NotImplementedError
	def __str__(self):
		if self.isPoint:
			return str(self.point().evalf())
		else:
			return str(self.line().evalf())
	def addToPlot(self):
		if self.isPoint:
			plt.scatter([self.points[0].evalf().x],[self.points[0].evalf().y])
		else:
			plt.scatter([self.points[0].evalf().x, self.points[1].evalf().x],[self.points[0].evalf().y, self.points[1].evalf().y])

''' We could make it O(logn) search etc, but this will do for now.
'''
class BeachlineElement:
	def __init__(self, site1, site2=None):
		if site2 is None: # Its a region of a single site.
			self.sites = [site1]
			self.isArc = True
			self.isBreakpoint = False
		else:
			self.sites = [site1, site2]
			self.isBreakpoint = True
			self.isArc = False
		self.circleEvent=None
		self.voronoiEdge=None
	def type(self):
		if self.sites[0].isPoint and self.sites[1].isPoint:
			return 1
	def eval(self, ys):
		if self.isArc:
			raise NotImplementedError
		else:
			if self.sites[0].isPoint and self.sites[1].isPoint:
				x0 = self.sites[0].point().x
				x1 = self.sites[1].point().x
				y0 = self.sites[0].point().y
				y1 = self.sites[1].point().y
				# Some quadractic equations.
				a2c = (y1 - y0)
				a1c = ((y1-ys)*(-2*x0)) - ((y0-ys)*(-2*x1))
				a0c = ((y1-ys)*(x0*x0)) + (((y0*y0) - (ys*ys))*(y1-ys)) - ((y0-ys)*(x1*x1)) - (((y1*y1)-(ys*ys))*(y0-ys))
				# A little bit of geometry trickery.
				if a2c.evalf()==0:
					aOption= ((x0+x1)/2.0) 
				else:
					aOption= ((-a1c + (sympy.sqrt(a1c*a1c - 4*a2c*a0c)))/(2*a2c))
				try:
					ans = Triangle(Point(x0,y0), Point(x1,y1), Point(aOption,ys)).circumcenter.evalf()
				except:
					ans = Triangle(Point(x0,y0), Point(x1,y1), Point(aOption,ys)).midpoint.evalf()
				return ans
			else:
				raise NotImplementedError # TODO ... handle other cases, where sites might be Segments etc.
	def __str__(self):
		if self.isArc:
			return '*R(' + str(self.sites[0]) + ')'
		else:
			return 'C( ['+ str(self.sites[0]) +'], ['+ str(self.sites[1]) +'] )'
class Event:
	def __init__(self, eventSite):
		self.site = eventSite
		self.isSiteEvent = True
		self.isFalseAlarm=False
		self.center=None
		self.beachlineElements=None # Beachline elements
	def circle(self):
		return Circle(self.center, abs(self.center.y - self.site.point().y) )
	def sweeplinePosition(self):
		return self.site.point().y

	def __str__(self):
		ans = 'Event[ @(' + str(self.site) + ') ]'
		if self.isSiteEvent:
			ans += ' is siteEvent'
		else:
			ans += ' is circleEvent, falseAlarm=' + str(self.isFalseAlarm) + ' '
			ans += str(self.center.evalf())
		return ans
''' A quick function to make a circle event. Also ensures that the targetArc's circle event is set correctly.

	@params atPoint the sympy.geometry.Point at which the circleEvent should happen. only care about the y-coordinate, really.
	@params withCenter the center of the circle event. keeping it for miscellaneous reasons.
	@params forBreakpoints the list of Breakpoints (pair of Breakpoints) which converged ... in order from left to right.
	@params targetArc the arc which will be removed by the circleEvent
		=> imposed ASSERT: the targetArc is between forBreakpoints[0] and forBreakpoints[1]
	@returns event The corresponding circle event.
'''
def makeCircleEvent(atPoint, withCenter, beachlineElements):
	# assert (len(forBreakpoints)==2),"Circle Events must be caused by convergence of two breakpoints."
	# assert ((forBreakpoints[0].rightArc==targetArc) and (forBreakpoints[1].leftArc==targetArc)), "The targetArc isn't matching up with the breakpoints provided."
	# assert ((targetArc.circleEvent is None) or (targetArc.circleEvent.isFalseAlarm)), "An Arc can't have two circle events."
	event = Event(Site(atPoint))
	if beachlineElements[1].circleEvent is not None:
		beachlineElements[1].circleEvent.isFalseAlarm = True
	beachlineElements[1].circleEvent=event
	
	event.center=withCenter
	event.beachlineElements=beachlineElements
	event.isSiteEvent=False
	return event

''' class VoronoiEdge provides a way to store edges traced out by BeachlineElements/Breakpoints
'''
class VoronoiEdge:
	def __init__(self, beachlineElement, ysStart):
		assert(beachlineElement.isBreakpoint), "Can't make a VoronoiEdge out of a non-Breakpoint"
		self.breakpoint=beachlineElement
		self.ysStart = ysStart
		self.ysEnd = None
	def close(self, ysEnd):
		self.ysEnd = ysEnd

	def addToPlot(self, boundary):
		if self.breakpoint.type()==1:
			start = self.breakpoint.eval(self.ysStart).evalf()
			if self.ysEnd is None:
				end = self.breakpoint.eval(self.ysStart - boundary).evalf()
			else:
				end = self.breakpoint.eval(self.ysEnd).evalf()
			plt.plot([start.x, end.x],[start.y, end.y])
		else:
			raise NotImplementedError

	def __str__(self):
		start = self.breakpoint.eval(self.ysStart).evalf()
		ans = 'Edge( ' + str(start) + ' , '
		if self.ysEnd is not None:
			end = self.breakpoint.eval(self.ysEnd).evalf()
			ans += str(end)
		else:
			end = self.breakpoint.eval(self.ysStart-50).evalf()
			ans += str(end)
		ans += ' )'
		return ans

class VoronoiDiagram:
	def __init__(self, lowerBoundary=20):
		self.edges = []
		self.sites = []
		self.lowerBoundary = lowerBoundary

	# @params e is an instance of VoronoiEdge.
	def addEdge(self, e):
		self.edges.append(e)
	def addSite(self, s):
		self.sites.append(s)
	def plot(self):
		for e in self.edges:
			e.addToPlot(self.lowerBoundary)
		for s in self.sites:
			s.addToPlot()
		plt.show()