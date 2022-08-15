# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 14:02:13 2022

@author: syang
"""


import subprocess
import os

try:
    os.mkdir('Plots')
except:
    pass

subprocess.call('python dip_domain_view.py', shell = True)
subprocess.call('python nondip_domain_view.py', shell = True)
subprocess.call('python nondip_lithology_view.py', shell = True)
subprocess.call('python dip_lithology_view.py', shell = True)
