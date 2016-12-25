'''
	Voronoi Diagram of Points.

	Inefficiencies:
		Beachline is simply a list of of breakpoints instead of 
		a binary search tree
		Arcs are maintained in the same fashion.

'''
from v3_structs import *


''' A class that builds the Voronoi Diagram
'''
class VoronoiMaker:
	
	''' @params S a list of sympy.geometry.Point ... 
	'''
	def build(self, S):
		self.vd = VoronoiDiagram()
		self.events=[]
		for s in S:
			heapq.heappush(self.events, (-s.y,Event(Site(s))))

		self.beachline=None

		while len(self.events)!=0:
			event=heapq.heappop(self.events)
			if event.isSiteEvent:
				self.handleSiteEvent(event)
			else:
				self.handleCircleEvent(event)
		return self.vd

	def constructSubTree(self, si, sj): # Dirty function for Step3 of handleSiteEvent()
		larc = Arc(sj)
		rarc = Arc(sj)
		marc = Arc(si)
		lbrk = Breakpoint(sj, si)
		rbrk = Breakpoint(si, sj)
		subtree = Beachline(breakpoint=lbrk)
		subtree.left = Beachline(arc=larc)
		rightTree = Beachline(breakpoint=rbrk)
		subtree.right = rightTree
		rightTree.left = Beachline(arc=marc)
		rightTree.right= Beachline(arc=rarc)
		subtree.left.par = subtree
		subtree.right.par = subtree
		rightTree.left.par = rightTree
		rightTree.right.par = rightTree
		return subtree, lbrk, rbrk

	def handleSiteEvent(self, event):
		eventPoint = event.site.point()
		if self.beachline is None: # Step 1
			self.beachline = Beachline(arc=Arc(event.site))
		else:
			# Step 2
			arc = self.beachline.search(eventPoint.y, event.site.point())
			assert (arc is not None),"Unknown error here. arc was None"
			if arc.circleEvent is not None:
				arc.circleEvent.isFalseAlarm = True # We aren't deleting it. Instead, we'll simply ignore it.

			# Step 3
			arcPar = arc.par # Parent of arc, if changes need to be made.
			sj = arc.site
			si = event.site
			newSubTree, jiBrk, ijBrl = self.constructSubTree(si, sj) # Three new leaf nodes, and 2 new interior nodes.
			if arcPar is None:
				self.beachline = newSubTree
			else:
				if arcPar.left==arc:
					arcPar.left=newSubTree
				else:
					arcPar.right=newSubTree
				newSubTree.par=arcPar
				arc.par=None # Its gone now.


			# Step 4: Add lbrk and rbrk to our Voronoi Diagrams.
			ledge = VoronoiEdge(lbrk, eventPoint.y)
			lbrk.voronoiEdge = ledge
			redge = VoronoiEdge(rbrk, eventPoint.y)
			rbrk.voronoiEdge = redge
			self.vd.addEdge(ledge)
			self.vd.addEdge(redge)

			# Step 5: Check the consecutive thingies.
			# TODO
			

	def handleCircleEvent(self, event):
		if event.isFalseAlarm:
			return
		arc = event.arc
		linesweepPosition=event.site.point()
