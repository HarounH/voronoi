from v5_structs import *
import pdb

class VorMaker:
	''' @params S is a set of points.
	'''
	def build(self, S):
		self.Q = [] # Priority Queue
		self.vd = VoronoiDiagram() # FIXME pass it a bounding box.
		
		for s in S:
			self.vd.addSite(Site(s))
			heapq.heappush(self.Q, ( -s.y, Event(Site(s))))
			
		p1 = heapq.heappop(self.Q)[1]
		
		# ASSERT: self.T always has the form of ...Arc, Breakpoint, Arc...
		# T[0] is always an Arc
		# T[-1] is always an Arc
		self.T = [ BeachlineElement(p1.site) ]
		while len(self.Q)!=0:

			e =  heapq.heappop(self.Q)[1] # This is an event.
			# pdb.set_trace()
			if e.isSiteEvent:
				self.handleSiteEvent(e)
			else:
				self.handleCircleEvent(e)
		return self.vd
	''' @params e is an Event, for which isSite is True.
		It represents a single site.
	'''
	def handleSiteEvent(self, e):
		ys = e.sweeplinePosition()
		# print('handling site event e=' + str(e) + ' and ys=' + str(ys))
		p = e.site
		leftT = None
		rightT= None
		Rq = None
		# Find Crq, Rq, Cqs which contains p

		# Four cases...
		# case 0 : there are no breakpoints to compare against.
		if len(self.T)==1:
			Rq = self.T[0]
			leftT = []
			rightT= []
		# case 1 : p is entirely to left of self.T
		elif (p.eval(ys).x < self.T[1].eval(ys).x):
			leftT = []
			Rq = self.T[0]
			rightT = self.T[1:]
		# case 2 : p is entirely to right of self.T
		elif (p.eval(ys).x > self.T[-2].eval(ys).x):
			leftT = self.T[:-1]
			Rq = self.T[-1]
			rightT = []
		# case 3 : some arc lies above p.
		else:
			for i in range(0, len(self.T)-2):
				if i%2==0: # It's an arc.
					continue
				else:
					if (p.eval(ys).x < self.T[i+2].eval(ys).x) and (p.eval(ys).x > self.T[i].eval(ys).x):
						leftT = self.T[:(i+1)]
						Rq = self.T[i+1]
						rightT = self.T[(i+2):]
						break
			
		# We have the left list, and the right list and the item to be replaced.
		if Rq.circleEvent is not None:
			Rq.circleEvent.isFalseAlarm = True

		Crq = leftT[-1] if len(leftT)>0 else None
		Cqs = rightT[0] if len(rightT)>0 else None
		
		q = Rq.sites[0]
		
		subtree = [ BeachlineElement(q), BeachlineElement(q,p), BeachlineElement(p), BeachlineElement(p,q), BeachlineElement(q) ]
		

		Cqp = subtree[1]
		Cpq = subtree[3]
		# Step 4
		qpEdge = VoronoiEdge(Cqp, ys)
		Cqp.voronoiEdge=qpEdge
		pqEdge = VoronoiEdge(Cpq, ys)
		Cpq.voronoiEdge=pqEdge
		self.vd.addEdge(qpEdge)
		self.vd.addEdge(pqEdge)

		# Step 5
		if Crq is not None and Cqp is not None and subtree[0] is not None:
			self.checkAndAddCircleEvent([Crq, subtree[0], Cqp], ys)

		if Cpq is not None and Cqs is not None and subtree[4] is not None:
			self.checkAndAddCircleEvent([Cpq, subtree[4], Cqs], ys)

		self.T = leftT + subtree + rightT

	'''	@params beachlineElements is [Crq, Rq, Cqp]
		@params ys is y-coordinate of sweepline when function is called.
		it checks if Crq and Cqp intersect. If they do, then it adds a circleEvent.
		Notice - if Crq and Cqp intersect, then a circle passes through r,q,p.
	'''
	def checkAndAddCircleEvent(self, beachlineElements, ys):
		Crq = beachlineElements[0]
		Rq = beachlineElements[1]
		Cqp = beachlineElements[2]
		assert (Crq.isBreakpoint and Rq.isArc and Cqp.isBreakpoint),"Invalid input to checkAndAddCircleEvent - BeachlineElements aren't of correct type."
		r = Crq.sites[0]
		q = Crq.sites[1]
		qPrime = Cqp.sites[0]
		qPrime2= Rq.sites[0]
		assert (q==qPrime and q==qPrime2),"Invalid input to checkAndAddCircleEvent - middle arc doesn't match"
		p = Cqp.sites[1]
		if (r.isPoint and q.isPoint and p.isPoint):
			if (Triangle(r.point(),q.point(),p.point()).area<0): # rqp is a left hand turn => circle is down under.
				circ = Circle(r.point(), q.point(), p.point())
				c = circ.center
				radius = circ.radius
				eventPoint = Point(c.x, c.y-radius)
				if eventPoint.y <= ys:
					e = makeCircleEvent(eventPoint, c, beachlineElements)
					heapq.heappush(self.Q, (-eventPoint.y, e))
				else:
					assert False,"We have a problem, right turn triangle."
		else: # There's a lot of combinations.
			raise NotImplementedError


	def handleCircleEvent(self, e):
		if e.isFalseAlarm:
			return
		ys = e.sweeplinePosition()
		arcIdx = self.T.index(e.beachlineElements[1])
		
		# e contains [Cqr Rr Crs]. arcIdx is the index of Rq
		Cqr = e.beachlineElements[0]
		Rr = e.beachlineElements[1]
		Crs = e.beachlineElements[2]
		assert ((Cqr==self.T[arcIdx-1]) and (Crs==self.T[arcIdx+1])), "something is not right"
		leftT = self.T[:(arcIdx-1)]
		rightT= self.T[(arcIdx+2):]
		Cuq = leftT[-2] if len(leftT)>1 else None
		Csv = rightT[1] if len(rightT)>1 else None
		Rq = leftT[-1] if len(leftT)>0 else None
		Rs = rightT[0] if len(rightT)>0 else None
		q = Cqr.sites[0]
		s = Crs.sites[1]
		Cqs = BeachlineElement(q,s)
		self.T = leftT + [ Cqs ] + rightT
		if ((Rq is not None) and (Rq.circleEvent is not None)):
			Rq.circleEvent .isFalseAlarm=True
		if ((Rs is not None) and (Rs.circleEvent is not None)):
			Rs.circleEvent .isFalseAlarm=True

		if Cuq is not None:
			self.checkAndAddCircleEvent([Cuq, Rq, Cqs], ys)
		if Csv is not None:
			self.checkAndAddCircleEvent([Cqs, Rs, Csv], ys)

		Cqr.voronoiEdge.close(ys)
		Crs.voronoiEdge.close(ys)
		qsEdge = VoronoiEdge(Cqs, ys)
		self.vd.addEdge(qsEdge)
		Cqs.voronoiEdge = qsEdge


	def handleCircleEvent_depre(self, e):
		# print('handling circle event e=' + str(e))		
		if e.isFalseAlarm:
			return
		ys = e.sweeplinePosition()
		arcIdx = None
		try:
			arcIdx = self.T.index(e.beachlineElements[1])
			lidx = arcIdx-1
			ridx = arcIdx+1
		except ValueError:
			print("Called for handleCircleEvent for an arc which is no longer present. ")
			raise ValueError

		leftT=self.T[:lidx]
		rightT=self.T[(ridx+1):]

		Cqr = e.beachlineElements[0]
		Crs = e.beachlineElements[2]
		
		assert (Cqr==self.T[lidx]),"left is not ok"
		assert (Crs==self.T[ridx]),"right is not ok"
		Cuq = leftT[-2] if len(leftT)>2 else None
		Csv = rightT[1] if len(rightT)>2 else None

		Cqs = BeachlineElement( Cqr.sites[0], Crs.sites[1] ) # New beachline element.
		self.T = leftT + [Cqs] + rightT

		# Goes from [Cuq Rq Cqr Rr Crs Rs Csv] to [Cuq Rq Cqs Rs Csv]

		# Step 1: delete some circles.
		
		# FIXME buggy bug bug. :P
		Rq = leftT[-1] if len(leftT)>0 else None
		if Rq is not None:
			if Rq.circleEvent is not None:
				Rq.circleEvent.isFalseAlarm = True
		Rs = rightT[0] if len(rightT)>0 else None
		if Rs is not None:
			if Rs.circleEvent is not None:
				Rs.circleEvent.isFalseAlarm = True

		# Step 2: Close off the edges, 
		Cqr.voronoiEdge.close(ys)
		Crs.voronoiEdge.close(ys)
		qsEdge = VoronoiEdge(Cqs, ys)
		Cqs.voronoiEdge = qsEdge

		# Step 3:
		if Cuq is not None:
			# Create a new circle event between leftT[-2] and Cqs
			self.checkAndAddCircleEvent([Cuq, Rq, Cqs], ys)
		if Csv is not None:
			self.checkAndAddCircleEvent([Cqs, Rs, Csv], ys)


def main():
	S = [ Point(1,-2), Point(5.0, 10.0), Point(1.5,1), Point(0,0) ]
	vm = VorMaker()
	vd = vm.build(S)
	for e in vd.edges:
		print e
	vd.plot()
	return vd

if __name__ == '__main__':
	main()
