#!/usr/bin/python

# Find all the positions of stop codons plus 0 and len(s)+1
def findall(s, value='*', start=0):
	i = start - 1
	pos = []
	pos.append(0)
	try:
		while 1:
			i = s.index(value, i+1)
			pos.append(i+1)
	except ValueError:
		pass
	pos.append(len(s)+1)
	
	return pos


# Return the biggest distance and the index of start and end pos
def max_dis(L):	
	m = 0
	start = 0
	end = 0
	for idx, obj in enumerate(L):
		if idx == 0 :
			pass
		else :
			if obj - L[idx-1] > m :
				m = obj - L[idx-1]
				start = L[idx-1]
				end = obj
	return [m, start, end]