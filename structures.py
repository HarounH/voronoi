'''
	This file contains basic geometric structures needed to implement 
	Fortune's Algorithm to construct Voronoi Diagram of Line Segment.


'''

# Depends heavily on the Geometry module provided by Sympy.
# We use the Point, Segment structures.

from sympy.geometry import *
import sympy
import math # For sqrt

'''
	A site is the input
	It is either a point or a line segment.
		If it is a point, it contians two points which are equal
		If it is a line segment, then it contains two points.
'''
class Site:
	def __init__(self, p1, p2):
		self.points = [p1,p2]
		self.points.sort(key=lambda p:-p.y) # Upper endpoint above.
		
		self.isPoint = p1==p2
		self.isSegment= not self.isPoint


		self.isUpperEndpoint = None
		self.isLowerEndpoint = None
		self.isInterior = None

	def isEndpoint(self):
		return self.isUpperEndpoint or self.isLowerEndpoint
	'''
		This function is called only if the Site is a point, or the endpoint of a line segment.
	'''
	def point(self):
		if self.isUpperEndpoint:
			return self.points[0]
		elif self.isLowerEndpoint:
			return self.points[1]
		else:
			# FIXME throw an exception and exit.
			pass
	def interior(self):
		if self.isInterior:
			return Segment(self.points[0], self.points[1])
		else:
			# FIXME throw an exception and exit.
			pass
	def line(self):
		if self.isInterior:
			return Line(self.points[0], self.points[1])
		else:
			# FIXME throw an exception and exit.
			pass
'''
	In Fortune's, we maintain a beachline. The beachline's interior nodes contain breakpoints, which are intersections
	between two arcs.
	Depending on what induces an arc, the breakpoints can be of different types.

	Each breakpoint is defined as the set of points equidistant from two different sites (except for type 4,5)
	Then, given a sweepline position, the breakpoint is defined as the point equidistant from the
	sites and the sweepline as well.



	We have 5 types of breakpoints
		1. Between two site endpoints - results in locus being a line segment.
			=> Here, site endpoint can also be a normal point
		2. Between two site interiors - results in locus being a line segment.
		3. Between site endpoint and site interior - locus is a parabola.
			=> Here, site endpoint can also be a normal point
		4. Closest to a site endpoint - locus is a line segment perpendicular the the site interior
			=> NOTE Here, the site endpoint CANNOT be normal point
		5. Intersection of site interior and sweepline
'''
class Breakpoint:
	def __init__(self, site1, site2=None):
		# Need to store the sites and also like... sort them, in some sense.
		if site2 is None:
			self.sites = [site1]
			self.type = 4 if site1.isEndpoint() else 5
		else:
			self.sites = [ site1, site2 ]
			if site1.isEndpoint():
				if site2.isEndpoint():
					self.type = 1
				else:
					self.type = 3
			else:
				if site2.isEndpoint():
					self.type = 3
				else:
					self.type = 2
	def eval(self,sweeplinePosition):
		if self.type==1: # the breakpoint follows the perpendicular bisector of the two sites.
			points = [self.sites[0].point(), self.sites[1].point()]
			x0 = points[0].x
			x1 = points[1].x
			y0 = points[0].y
			y1 = points[1].y
			ys = sweeplinePosition

			# Say the answer is (a,b)
			# We construct a quadractic equation for a.
			# Then we get b by using the value of a... heck, let's just feed the stuff into the circumcentre 
			a2c = y0 - y1
			a1c = 2*(0 -x1*y0 + x1*ys + x0*y1 - x0*ys)
			a0c = ((x1*x1 + y1*y1 - ys*ys)*(y0-ys)) - ((x0*x0 + y0*y0 - ys*ys)*(y1-ys))

			aOptions = [ ((-a1c + p*(math.sqrt(a1c*a1c - 4*a2c*a0c)))/(2*a2c)) for p in [1.0,-1.0] ]
			if aOptions[0] <= max(x0,x1) and aOptions[0]>=min(x0,x1):
				a = aOptions[0]
			else:
				a = aOptions[1]

			return Triangle(points[0], points[1], Point(a,ys)).circumcenter

		elif self.type==2: # Be smart. This is the incircle of three points
			l0 = self.sites[0].line()
			l1 = self.sites[1].line()
			l2 = Line(Point(0,sweeplinePosition), slope=0)
			p0 = l0.intersection(l1)[0]
			p1 = l1.intersection(l2)[0]
			p2 = l2.intersection(l0)[0]
			return Triangle(p0,p1,p2).incircle.center.evalf() # FIXME ideally, we'd just pass this around as a symbolic expression.

		elif self.type==3: # between site endpoint and site interior
			temp = (self.sites[0],self.sites[1]) if self.sites[0].isEndpoint() else (self.sites[1], self.sites[0])
			point,line = temp[0].point(), temp[1].line()

			coeffs = line.coefficients

			s = sympy.sqrt(coeffs[0]*coeffs[0] + coeffs[1]*coeffs[1])

			a = (coeffs[0]/s)
			b = (coeffs[1]/s)
			c = (coeffs[2]/s)

			x0 = point.x
			y0 = point.y
			ys = sweeplinePosition

			# TODO
		elif self.type==4: # It is a site endpoint connected to a site interior.
			point = self.sites[0].point()
			line = Line( self.sites[0].points[0], self.sites[0].points[1] )
			perp1 = line.perpendicular(point)

			flat = Line(Point(0, sweeplinePosition), slope=0)
			lf = line.intersection(flat)
			d = lf.distance(point)
			drop = Point( (lf.x + (1.0 if point.x>lf.x else -1.0)*d) , sweeplinePosition)

			perp2 = flat.perpendicular(drop)

			return perp2.intersection(perp1)
		
		elif self.type==5:
			interior = self.sites[0].interior()
			line=Line(interior)
			flat = Line(Point(0, sweeplinePosition), slope=0)
			return line.intersection(flat)
		else:
			pass
			# FIXME raise an exception



'''
	Fortune's algorithm has two kinds of Events - Site Events and Circle/Vertex Events.
'''
# class Event:

if __name__ == '__main__':
	p1 = Point(2,3)
	p2 = Point(4,4)
	p3 = Point(5,5)
	p4 = Point(6,6)
	p5 = Point(5.5, 4.5)
	p6 = Point(7, 3.5)
	
	site1 = Site(p1,p2)
	site1.isInterior = True
	site2 = Site(p3,p4)
	site2.isLowerEndpoint = True
	site3 = Site(p5,p6)
	site3.isUpperEndpoint = True

	b1 = Breakpoint(site1, site2)
	# b2 = Breakpoint(site2, site3)

	print(b1.eval(0.0).evalf())
	# print(b2.eval(0.0).evalf())