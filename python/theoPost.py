import math
import sys

#if ( sys.version_info < (3,0) ):
#	raise "Must use python version 3 or later. See /broad/software/free/Linux/redhat_5_x86_64/pkgs/python_3.1.2/bin/python3.1"

class controls:
	def __init__(self,l_x,l_y,c_x,c_y,r_x,r_y):
		self.x_a = l_x - 2*c_x + r_x
		self.y_a = l_y - 2*c_y + r_y
		self.x_b = -2*l_x + 2*c_x
		self.y_b = -2*l_y + 2*l_x
		self.x_c = l_x
		self.y_c = l_y
		
class IntegrationCollection:
	def __init__(self,start,stop,err,num_ints):
		self.start = start
		self.stop = stop
		self.err = err
		self.num_ints = num_ints

lfmap = dict()
def logfact(a):
	global lfmap
	if ( a < 2 ):
		return 0.0
	if ( a not in lfmap ):
		lfmap[a] = math.log10(a) + logfact(a-1)
	return lfmap[a]

def logchoose(a,b):
	return logfact(a)-logfact(b)-logfact(a-b)

def logbinomial(success,trials,prob):
	return logchoose(trials,success) + success*math.log10(prob) + (trials-success)*math.log10(1-prob)

def quad(a,b,c,x):
	return a*x*x + b*x + c

def qformula(a,b,c,equivVal):
	return (-b + math.sqrt(b*b-4*a*c))/(2*a)

def cbezierf(cts,pt):
	t = qformula(cts.x_a,cts.x_b,cts.x_c,pt)
	y = quad(cts.y_a,cts.y_b,cts.y_c,t)
	return y

bez_cts = controls(-7,10.99919,-2.849154,0.1444735,-0.0043648054,-1.559080) # based on previous gradient descent


def simpson(f,ic):
	class DeprecationError(Exception):
		def __init__(self,val):
			self.value = val
		def __str__(self):
			return repr(self.value)
	raise DeprecationError("Simpson is deprecated. Do not use it.")

def simpAux(f,a,b,eps,s,fa,fb,fc,cap):
	if ( s == 0 ):
		return []
	c = ( a + b )/2
	h = b-a
	d = (a + c)/2
	e = (c + b)/2
	fd = f(d)
	fe = f(e)
	s_l = (h/12)*(fa + 4*fd + fc)
	s_r = (h/12)*(fc + 4*fe + fb)
	s_2 = s_l + s_r
	if ( cap <= 0 or abs(s_2 - s) <= 15*eps ):
		try:
			return [math.log10(s_2 + (s_2 - s)/15.0)]
		except OverflowError:
			print(s_2)
			print(s_2-s)
			return [-350]
	return simpAux(f,a,c,eps/2,s_l,fa,fc,fd,cap-1) + simpAux(f,c,b,eps/2,s_r,fc,fb,fe,cap-1)

def adaptiveSimpson(f,start,stop,error,cap):
	mid = (start + stop)/2
	size = stop - start
	fa = f(start)
	fb = f(mid)
	fc = f(stop)
	s = (size/6)*(fa + 4*fc + fb)
	h = simpAux(f,start,stop,error,s,fa,fb,fc,int(cap))
	h.sort()
	#print("first: "+str(h[0]))
	#print("last: "+str(h[len(h)-1]))
	return sum(map(lambda x: 10**x,h))

def neutral(x):
	return -1.0*math.log10(x)

def twoState(x):
	if ( x < 0.04 ):
		return -1.5*math.log10(x)
	else:
		return -1.0*math.log10(x)

def bezier(x):
	return cbezierf(bez_cts,math.log10(x))

norm_cache = (None,None)
def resampleProbability(logshape,ic,ac,ns,ac_new,ns_new):
	global norm_cache
	logpost = lambda x: logshape(x) + logbinomial(ac,2*ns,x)
	if ( norm_cache[1] == None or norm_cache[0] != (ac,ns,logshape) ):
		print("Caching posterior norm")
		norm_cache = ((ac,ns,logshape),math.log10(adaptiveSimpson( lambda v: math.pow(10,logpost(v)), ic.start,ic.stop,ic.err,ic.num_ints)))
	logpost_normed = lambda v: logpost(v) - norm_cache[1]
	newshape = lambda y: math.pow(10,logpost_normed(y) + logbinomial(ac_new, 2*ns_new, y))
	return adaptiveSimpson(newshape,ic.start,ic.stop,ic.err,ic.num_ints)

sim_ic = IntegrationCollection(4e-7,0.999,1e-200,16)
sys.setrecursionlimit(int(2e6))
neutral_post = map( lambda v: resampleProbability(neutral,sim_ic,1,900,v,900), range(0,21) ) 
twostate_post = list(map( lambda v: resampleProbability(twoState,sim_ic,1,900,v,900), range(0,21) ))
g = open("n_ts.txt",'w')
idx = 0
for e in neutral_post:
	g.write(str(idx))
	g.write("\t"+str(e)+"\t"+str(twostate_post[idx])+"\n")
	idx += 1
