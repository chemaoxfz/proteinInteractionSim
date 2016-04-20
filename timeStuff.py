# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 22:29:28 2015

@author: xfz
"""

import time

start=time.time()
for i in xrange(100000000):
    x=1+1e-8
    x=x*x

end=time.time()
print(end-start)