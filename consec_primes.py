def consec_primes():
	LIMIT = 3628800

	primes = [False, False] + [True] * (LIMIT+1)
	ans = [2]
	for i in xrange(3, LIMIT/2 +1, 2):
		if primes[i]:
			ans.append(i)
			for j in xrange(i+i, LIMIT/2, i):
				primes[j] = False
	
	for k in xrange(1, len(ans)):
		if ans[k]-ans[k-1] > 9:
			print "%d-%d" % (ans[k-1]+1, ans[k]+1)

if __name__ == "__main__":
	consec_primes()
