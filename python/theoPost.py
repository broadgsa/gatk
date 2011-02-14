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
		return -1.5*math.log10(x) + 0.5*math.log10(0.04)
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

def getPost(logshape,ic,ac,ns):
	global norm_cache
	logpost = lambda x: logshape(x) + logbinomial(ac,2*ns,x)
	if ( norm_cache[1] == None or norm_cache[0] != (ac,ns,logshape) ):
		print("Caching posterior norm")
		norm_cache = ((ac,ns,logshape),math.log10(adaptiveSimpson(lambda v: math.pow(10,logpost(v)),ic.start,ic.stop,ic.err,ic.num_ints)))
	return lambda v: logpost(v) - norm_cache[1]

sim_ic = IntegrationCollection(5e-8,0.999,1e-2000,22)
sys.setrecursionlimit(int(2e6))
#neutral_post = map( lambda v: resampleProbability(neutral,sim_ic,1,900,v,900), range(0,21) ) 
#twostate_post = list(map( lambda v: resampleProbability(twoState,sim_ic,1,900,v,900), range(0,21) ))
#g = open("n_ts.txt",'w')
#idx = 0
#for e in neutral_post:
#	g.write(str(idx))
#	g.write("\t"+str(e)+"\t"+str(twostate_post[idx])+"\n")
#	idx += 1

DO_1 = False
if ( DO_1 ):
	eomiautism_ac_1 = 317763
	eomiautism_ac_2 = 78844
	eomiautism_ac_3p = 239526 # all of these go on chip by default

	new_set = 10000-917-998

	num_unseen_sites = 125*new_set

	unseen_unseen = resampleProbability(twoState,sim_ic,0,917+998,0,new_set)
	unseen_1 = resampleProbability(twoState,sim_ic,0,917+998,1,new_set)/(1-unseen_unseen)
	unseen_2 = resampleProbability(twoState,sim_ic,0,917+998,2,new_set)/(1-unseen_unseen)
	ac1_unseen = resampleProbability(twoState,sim_ic,1,917+998,0,new_set)
	ac1_ac1 = resampleProbability(twoState,sim_ic,1,917+998,1,new_set)
	ac2_unseen = resampleProbability(twoState,sim_ic,2,917+998,0,new_set)

	total = 636133 + num_unseen_sites
	ac1 = unseen_1*num_unseen_sites + ac1_unseen*eomiautism_ac_1
	ac2 = unseen_2*num_unseen_sites + ac1_unseen*eomiautism_ac_1 + ac2_unseen*eomiautism_ac_2

	print("\t".join(map(lambda u: str(u), [unseen_unseen,unseen_1,unseen_2])))
	print("\t".join(map(lambda u: str(u), [total,ac1,ac2])))

	ea_ns = 343877
	ea_ns_ac1 = 204223
	ea_ns_ac2 = 42280

	ns_new_ac1 = ea_ns_ac1*ac1_unseen + unseen_1*num_unseen_sites*(1.7/(1+1.7))
	ns_new_ac2 = ea_ns_ac2*ac2_unseen + unseen_2*num_unseen_sites*(1.4/(1+1.4))
	ns_new_total = ea_ns + ns_new_ac1 + ns_new_ac2 + (num_unseen_sites*(1-unseen_1-unseen_2))*(0.6/(1+0.6))

	print("\t".join(map(lambda u: str(u), [ns_new_total,ns_new_ac1,ns_new_ac2])))
	print(1-resampleProbability(twoState,sim_ic,2,1000,0,10000)-resampleProbability(twoState,sim_ic,2,1000,1,10000))
	print(1-resampleProbability(twoState,sim_ic,1,100,0,2000)-resampleProbability(twoState,sim_ic,1,100,1,2000))
	print(1-resampleProbability(twoState,sim_ic,2,100,0,2000)-resampleProbability(twoState,sim_ic,2,100,1,2000))
	print(1-resampleProbability(twoState,sim_ic,20,1000,0,2000)-resampleProbability(twoState,sim_ic,20,1000,1,2000))

def emitPosterior(ac):
	post_ac = getPost(twoState,sim_ic,ac,10000)
	o = open("post_%d.txt" % ac,'w')
	pt = sim_ic.start
	while ( pt < 0.2 ):
		o.write("%e\t%e\n" % (pt,post_ac(pt)))
		pt = 1.015*pt
	o.close()

emitPosterior(2)
emitPosterior(3)
emitPosterior(10)
emitPosterior(25)

#o = open("test2s.txt",'w')
#pt = sim_ic.start
#while ( pt<0.4 ):
#	o.write("%e\t%e\n" % (pt,twoState(pt)))
#	pt = 1.015*pt
#o.close()
