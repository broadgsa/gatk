#!/usr/bin/env python

import time
class timer:
    def __init__(self):
        self.start()
    def start(self):
        self.start_time = time.time()
    def stop(self):
        return time.time() - self.start_time

class time_func:
    def __init__(self, fn):
        self._fn = fn
    def __call__(self, *args, **keywords):
        tmr = timer()
        result = self._fn(*args, **keywords)
        elapsed = tmr.stop()
        print "%.5fs elapsed." % elapsed
        return result

# DiskMemoize caches the results of function calls on disk for later, speedy retrieval
import cPickle, os, sys, math
class DiskMemoize:
    """Memoize(fn) - an instance which acts like function fn given at instantiation
    but memoizes the returned result of fn according to the function's arguments.
    Memoized results are stored on disk 
    """
    # Supposedly will only work on functions with non-mutable arguments
    def __init__(self, fn, fn_name, global_deps):
        """fn - name of function to call when object is called as a function
          fn_name - string corresponding to fn that is used to name memoization disk file
          global_deps - global variables that affect the outcome of the function call provided
            as a list of variable names as strings"""
        self.fn = fn
        self.fn_name = fn_name
        self.global_deps = global_deps
        #self.memo = {}
    def __call__(self, global_deps = [], volatile_keywords = {}, skip_cache=False, verbose=True, *args, **keywords):
        """Calls the function assigned at initialization - self.fn - attempting to 
        use a memoized result if it has been stored on disk.

        The function arguments are (in order of most likely use):
          args - non-keyword arguments are DISABLED and will cause an exception to be raised
          keywords - keyword arguments provided in standard function call style
          volatile_keywords - hash of keywords (var_name, value) pairs to be passed to the 
            function call that do NOT affect the cached result in a meaningful way (e.g. verbose flag)
          skip_cache - skips checking the cache for a memoized result, only memoize a new result 
            to disk (default: True)"""
               
        #print "BONUS:",BONUS
        #print globals()["BONUS"]
        
        # Raise an exception if any non-keyword arguments were provided
        if len(args) > 0:
            raise Exception("DiskMemoize.__call__ encountered non-keyword arguments (args); it can only be called with keyword arguments (keywords)")
        
        # Pickling filename specifying all dependencies
        pkl_filename = self.fn_name+"__"
        pkl_filename += "__".join([str(arg)+"."+str(val) for arg,val in keywords.items()] + \
                                  ["global_"+arg_name+"."+str(globals()[arg_name]) for arg_name in global_deps])
        pkl_filename += ".pkl"

        if os.path.exists(pkl_filename) and not skip_cache:
            # Just read the memoized result
            pkl_file = open(pkl_filename, "rb")
            print "DiskMemoize: Loading result from",pkl_filename, ;sys.stdout.flush()
            tmr = timer()
            result = cPickle.load(pkl_file)
            calc_time = cPickle.load(pkl_file)
            elapsed = tmr.stop()
            speedup = calc_time/elapsed if elapsed != 0 else 9999.9999
            print ": %.5fs loading, %.5fs original calc. time, %.2f X speedup" % (elapsed, calc_time, speedup)
            return result
        else:
            # Calculate the result and memoize it to disk
            pkl_file = open(pkl_filename, "wb")

            print "DiskMemoize: Calculating result for:",
            #print self.fn_name+" ("+(", ".join(map(str,args))) + ")"
            keywords.update(volatile_keywords)
            tmr = timer()
            print self.fn_name+" ("+(", ".join(["=".join(map(str, it)) for it in keywords.items()])) + ")",
            sys.stdout.flush()
            result = self.fn(*args,**keywords)
            calc_time = tmr.stop()
            print ": %.5fs calculating" % calc_time

            tmr = timer()
            print "DiskMemoize: Writing result to",pkl_filename, ;sys.stdout.flush()
            cPickle.dump(result, pkl_file, protocol=2)
            cPickle.dump(calc_time, pkl_file, protocol=2)
            print ": %.5fs writing" % tmr.stop()

            return result

def primes(first_prime, last_prime, verbose=False):
    """Generates lists of prime numbers from first_prime to last_prime.
    Used to test DiskMemoize"""

    #if prime < 2: return 
    primes = range(first_prime, last_prime)
    if verbose:
        print "I'm being verbose"
        print "Prime range is",first_prime,"to",last_prime
    last_div = int(math.sqrt(last_prime)+1)
    for div in range(2, last_div): # (len(primes)/2.0)+1):
        next_primes = []
        for x in primes:
            if x % div != 0 or x == div:
                #primes[x] = 0
                next_primes.append(x)
                #primes.pop(i)
                #pop_list.append(i)
        primes = next_primes
        if "BONUS" in globals():
            primes += [0] * BONUS
    return primes * 10000

BONUS = 0  

if __name__ == "__main__":

    print_checks = False

    print "Running DiskMemoize test suite."
    primes = DiskMemoize(primes, "primes", global_deps = ['BONUS'])

    prime_list = primes(first_prime = 20, last_prime = 100, skip_cache = True, volatile_keywords = {'verbose' : True})
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 20, last_prime = 100)
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    print "skip_check = True"
    prime_list = primes(first_prime = 20, last_prime = 100, skip_cache = True)
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 20, last_prime = 100)
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 20, last_prime = 100, skip_cache = True)
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    BONUS = 5
    print "BONUS:",BONUS
    prime_list = primes(first_prime = 20, last_prime = 100, skip_cache = True)
    if print_checks: print "primes:", prime_list, "\nLength:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 20, last_prime = 10000, skip_cache = True)
    if print_checks: print "Length:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 20, last_prime = 10000)
    if print_checks: print "Length:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 1, last_prime = 30000, skip_cache = True)
    if print_checks: print "Length:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"

    prime_list = primes(first_prime = 1, last_prime = 30000)
    if print_checks: print "Length:", len(prime_list), "\nHash:", hash(str(prime_list)), "\n"
 




#INTJ????
#INTP????
#extreme P

# cookies with oatmeal, coconut, chocolate chip


















































# For functions taking mutable arguments, use the cPickle module, as
# in class MemoizeMutable:

#class MemoizeMutable:
#    """Memoize(fn) - an instance which acts like fn but memoizes its arguments
#       Will work on functions with mutable arguments (slower than Memoize)
#    """
#    def __init__(self, fn):
#        self.fn = fn
#        self.memo = {}
#    def __call__(self, *args):
#        import cPickle
#        str = cPickle.dumps(args)
#        if not self.memo.has_key(str):
#            self.memo[str] = self.fn(*args)
#        return self.memo[str]
