from __future__ import print_function
import os,sys
from sympy.geometry import *
import heapq
import sympy
import numpy as np
import matplotlib.pyplot as plt
'''
	Implementation of Fortune's Algorithm for constructing Voronoi Diagram of Points.
	We do a line sweep from top to bottom.

'''

class Site:
	def __init__(self, pt1, pt2=None):
		if pt2 is None:
			self.point=pt1
			self.points=[pt1,pt1]
			self.isPoint=True
		else:
			self.point=None
			self.points=[pt1,pt2]
			self.points.sort(key=lambda p:(-p.y,p.x))
			self.isPoint=False
			self.isUpperEndpoint=None
			self.isLowerEndpoint=None
			self.isInterior=None

'''
	Need to generalize this to use
	Sites instead of just points.

	Breakpoints are defined by the arc to its left and the arc to its right.
	Each arc is defined by a single site.
	Hence, we define breakpoints using the sites that define the arcs to the left and right.
'''
class Breakpoint:
	def __init__(self, leftSite, rightSite):
		self.sites = [leftSite, rightSite]
		# self.points.sort(key=lambda p: (-p.y, p.x)) # Top to bottom and left to right. :)
		self.isBreakpoint=True

	def eval(self, sweeplinePosition):
		# Calculate the point that is equidistant from self.sites and Line(Point(0,sweeplinePosition), slope=0)
		if self.sites[0].isPoint and self.sites[1].isPoint:

			x0 = self.sites[0].point.x
			x1 = self.sites[1].point.x
			y0 = self.sites[0].point.y
			y1 = self.sites[1].point.y
			ys = sweeplinePosition

			# Say the answer is (a,b)
			# We construct a quadractic equation for a.
			# Then we get b by using the value of a... heck, let's just feed the stuff into the circumcentre 
			a2c = (y1 - y0)
			a1c = ((y1-ys)*(-2*x0)) - ((y0-ys)*(-2*x1))
			a0c = ((y1-ys)*(x0*x0)) + (((y0*y0) - (ys*ys))*(y1-ys)) - ((y0-ys)*(x1*x1)) - (((y1*y1)-(ys*ys))*(y0-ys))
			
			''' The following commented code would've returned all intersections ... but we want the one as defined by ordering of sites. '''
			# if a2c.evalf()==0.0:
			# 	aOptions = [ ((x0+x1)/ 2) ]
			# else:
			# 	aOptions = [ ((-a1c + p*(sympy.sqrt(a1c*a1c - 4*a2c*a0c)))/(2*a2c)) for p in [1.0,-1.0] ]
			# centers = [ Triangle(Point(x0,y0), Point(x1,y1), Point(a,ys)).circumcenter for a in aOptions ]
			# return centers
			if a2c.evalf()==0:
				aOption= ((x0+x1)/2.0) 
			else:
				aOption= ((-a1c + (sympy.sqrt(a1c*a1c - 4*a2c*a0c)))/(2*a2c)) # FIXME see if this is okay.
			return Triangle(Point(x0,y0), Point(x1,y1), Point(aOption,ys)).circumcenter
		else: # All the other cases of being an endpoint, being an interior etc... 5 cases -_-
			raise NotImplementedError
	def line(self):
		if self.leftSite.isPoint and self.rightSite.isPoint:
			return Segment(self.leftSite.point, self.rightSite.point).perpendicular_bisector()
def convergence(brk1, brk2):
	# returns the point at which two breakpoints will converge. This will act as the 
	l1 = brk1.line()
	l2 = brk2.line()
	return l1.intersection(l2)

class Event:
	def __init__(self, loc):
		self.point=loc # The location at which the event is.
		self.isSiteEvent = False
		self.isFalseAlarm = False
		self.arc = None # Used for circles.
		self.sites = None
		self.breakpoints = None
		self.center = None
class Arc:
	def __init__(self, cause):
		self.point=cause
		self.isBreakpoint=False


class BeachlineNode: # Balanced Binary Search 
	def __init__(self):
		# structural
		self.par = None
		self.left= None
		self.right=None
		# No need to maintain size, really.
		# data
		self.isLeaf = False
		self.isInterior=False
		self.breakpoint = None
		self.point = None
		self.circleEvent = None
		self.edgeInVoronoiDiagram=None
	def balance(self):
		pass # TODO Not a balanced tree. Ideally, it should be. But it isn't.

	def inorderBreakpoints(self):
		if self.isLeaf:
			return []
		else:
			return self.left.inorderBreakpoints() + [self.breakpoint] + self.right.inorderBreakpoints()


	def search(self, linesweepLocation, point):
		if self.isLeaf:
			return self
		else:
			loc = self.breakpoint.eval(linesweepLocation)
			if loc.y<point.y: # Go right.
				return self.right.search(linesweepLocation, point)
			else: # Go left.
				return self.left.search(linesweepLocation, point)

class VoronoiEdge:
	def __init__(self, brk):
		self.breakpoint=brk
		self.sweeplinePositionStart=None
		self.sweeplinePositionEnd=None
	def trace(self, stepSize=-0.01):
		pts = []
		for ys in np.arange(self.sweeplinePositionStart, self.sweeplinePositionEnd, stepSize):
			pts.append(self.breakpoint.eval(ys).evalf())
		return pts
	def draw(self):
		start = self.breakpoint.eval(self.sweeplinePositionStart).evalf()
		end = self.breakpoint.eval(self.sweeplinePositionEnd).evalf()
		plt.plot([start.x, end.x],[start.y, end.y])
class VoronoiDiagram:
	def __init__(self):
		self.edges=[]
	def addEdge(self, edge):
		self.edges.append(edge)
	def plot(self):
		# Create a matplotlib figure
		for e in self.edges:
			# trace e from e.sweeplinePositionStart to e.sweeplinePositionEnd.
			# add the trace to the matplotlib figure
			e.draw()
		# display the plot or (save&close) it
		plt.show()

class VoronoiBuilder:
	def build(self, S):
		self.Q = []

		for s in S:
			e = Event(s)
			e.isSiteEvent = True
			heapq.heappush(self.Q, (-s.y, e))

		self.beachline = None
		self.vd = VoronoiDiagram()

		while( len(self.Q)!=0 ):
			event = heapq.heappop(self.Q)
			if event.isSiteEvent:
				# Find the breakpoint and the correct arc.
				self.handleSiteEvent(event)
			else:
				if not event.isFalseAlarm:
					self.handleCircleEvent(event)
				else: # It's a false alarm.
					pass
		return self.vd

	def handleSiteEvent(self, event):
		if self.beachline is None: # The size is 0.
			self.beachline = BeachlineNode()
			self.beachline.isLeaf = True
			self.beachline.point = event.point
			return
		else:
			arc = self.beachline.search(event.point)
			if arc.circleEvent is not None:
				arc.circleEvent.isFalseAlarm = True
			# Following notation of the book.
			pj = arc.point
			pi = event.point

			# Modify arc to pj-<pj,pi>-(pi-<pi,pj>-pj)
				# call them: pjl-<pjpi>-(pi_node-<pipj>-pjr)
				# arc is #0, and <pjpi> IS arc
			pjl = BeachlineNode() #1
			pipj= BeachlineNode() #2
			pi_node  = BeachlineNode() #3
			pjr = BeachlineNode() #4

				#0
			arc.left = pj
			arc.right= pipj
			arc.isLeaf=False
			arc.isInterior=True
			arc.breakpoint=Breakpoint(Site(pj), Site(pi))

				#1
			pjl.par = arc
			pjl.isInterior=False
			pjl.isLeaf=True
			pjl.point=pj

				#2
			pipj.par=arc
			pipj.isLeaf=False
			pipj.isInterior=True
			pipj.breakpoint=Breakpoint(Site(pi), Site(pj))
			pipj.left=pi_node
			pipj.right=pjr
				#3
			pi_node.par = pipj
			pi_node.isLeaf=True
			pi_node.isInterior=False
			pi_node.point=pi

				#4
			pjr.par=pipj
			pjr.isLeaf=True
			pjr.isInterior=False
			pjr.point=pj

			# Add edges to output
			ij_edge=VoronoiEdge(pipj.breakpoint)
			ij_edge.sweeplinePositionStart=event.point.y
			self.vd.addEdge(ij_edge)
			pipj.edgeInVoronoiDiagram=ij_edge
			ji_edge=VoronoiEdge(pjpi.breakpoint)
			ji_edge.sweeplinePositionStart=event.point.y
			self.vd.addEdge(ji_edge)
			pjpi.edgeInVoronoiDiagram=ji_edge


			# Add circle events
			# TODO This is presently O(n). Make it O(logn) by implementing a predecessor, successor function
			#		in class BeachlineNode.
			inorder = self.beachline.inorderNodes()
			
			leftArcBreakPoints = None
			rightArcBreakPoints= None

			for i in range(0, len(inorder)):
				if inorder[i].breakpoint==arc.breakpoint: # right arc.
					rightArcBreakPoints = [ inorder[i-1].breakpoint , arc.breakpoint] if (i-1)>0 else []
				if inorder[i].breakpoint==pipj.breakpoint:
					leftArcBreakPoints = [ pipj.breakpoint, inorder[i+1].breakpoint ] if (i+1)<len(inorder) else []
					
			if len(leftArcBreakPoints)>1:
				# Circle event which causes disappearence of arc to the left of leftArcBreakPoints[1]
				circleCenter = convergence(leftArcBreakPoints[0], leftArcBreakPoints[1]) 
				radius = circleCenter.distance(leftArcBreakPoints[1].sites[0].point)
				newEventPoint = Point(circleCenter.x, circleCenter.y - radius)
				if newEventPoint.y < event.point.y: # create a new circle event.
					circleEvent = Event(newEventPoint)
					heapq.heappush(self.Q, (-newEventPoint.y,circleEvent))
					leftArcBreakPoints[1].left.circleEvent = circleEvent
					circleEvent.arc = leftArcBreakPoints[1].left
					circleEvent.breakpoints = leftArcBreakPoints
					circleEvent.sites = [ pipj.breakpoint.sites[0], pipj.breakpoint.sites[1], leftArcBreakPoints[1].sites[1] ]
					circleEvent.center=circleCenter
			if len(rightArcBreakPoints)>1:
				# Circle event which causes disappearence of arc to the right of rightArcBreakPoints[0]
				circleCenter = convergence(rightArcBreakPoints[0], rightArcBreakPoints[1])
				radius = circleCenter.distance(rightArcBreakPoints[1].sites[0].point)
				newEventPoint = Point(circleCenter.x, circleCenter.y - radius)
				newEventPoint = Point(circleCenter.x, circleCenter.y - radius)
				if newEventPoint.y < event.point.y: # create a new circle event.
					circleEvent = Event(newEventPoint)
					heapq.heappush(self.Q, (-newEventPoint.y,circleEvent))
					rightArcBreakPoints[0].right.circleEvent = circleEvent
					circleEvent.arc = rightArcBreakPoints[0]
					circleEvent.breakpoints = rightArcBreakPoints
					circleEvent.sites = [ rightArcBreakPoints[0].sites[0], rightArcBreakPoints[0].sites[1], arc.breakpoint.sites[1] ]


	def handleCircleEvent(self, event):
		
		

if __name__ == '__main__':
	pt1 = Site(Point(4,3))
	pt2 = Site(Point(5,3.5))
	ys = 0.0
	b1 = Breakpoint(pt1,pt2)
	b2 = Breakpoint(pt2,pt1)
	e = VoronoiEdge(b1)
	e.sweeplinePositionStart=2.0
	e.sweeplinePositionEnd=0.0
	e2 = VoronoiEdge(b2)
	e2.sweeplinePositionStart=2.0
	e2.sweeplinePositionEnd=0.0
	vd = VoronoiDiagram()
	vd.addEdge(e)
	vd.addEdge(e2)
	vd.plot()
	
