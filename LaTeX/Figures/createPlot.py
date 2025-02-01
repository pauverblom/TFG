#!/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(5,5))
afx=np.arange(0.0,1.0,0.01)
afsin=np.sin(afx)
plt.plot(afx,afsin)
plt.savefig('sin.pdf',bbox_inches='tight')

plt.clf()
plt.figure(figsize=(5,5))
afx=np.arange(0.0,1.0,0.01)
afsin=np.sin(afx)
plt.plot(afx,afsin)
plt.savefig('sin.png',bbox_inches='tight')

plt.clf()
plt.figure(figsize=(8,8))
afx=np.arange(0.0,1.0,0.01)
afsin=np.sin(afx)
plt.plot(afx,afsin)
plt.title("Sin function")
plt.xlabel(r'$x$')
plt.ylabel(r'$\sin(x)$')
plt.savefig('sinBig.pdf',bbox_inches='tight')

plt.clf()
plt.figure(figsize=(4,4))
afx=np.arange(0.0,1.0,0.01)
afsin=np.sin(afx)
plt.plot(afx,afsin)
plt.title("Sin function")
plt.xlabel(r'$x$')
plt.ylabel(r'$\sin(x)$')
plt.savefig('sinSmall.pdf',bbox_inches='tight')
