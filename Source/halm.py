# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 13:40:19 2022

@author: syang
"""

import subprocess
import os

try:
    os.mkdir('Output')
except:
    pass

subprocess.call('python transform_horizon_top.py', shell = True)
subprocess.call('python transform_horizon_bot.py', shell = True)
subprocess.call('python map_logs.py', shell = True)
subprocess.call('python Interpolation_nondip.py', shell = True)
subprocess.call('python map_back_lithology.py', shell = True)
subprocess.call('python generate_grid_lithology.py', shell = True)
