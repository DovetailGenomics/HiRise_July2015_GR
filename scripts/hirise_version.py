#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys
import argparse
import subprocess
import shlex

parser = argparse.ArgumentParser()

parser.add_argument('-d','--debug',default=False,action="store_true")


args = parser.parse_args()

import os

script_path=os.path.realpath(__file__)
scripts_dir = os.path.dirname(script_path)
#print((script_path,os.path.dirname(script_path)))

cmd = "git describe --dirty"
output,error = subprocess.Popen(shlex.split(cmd),stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.PIPE, cwd=scripts_dir).communicate()
print( scripts_dir )
print( output.decode('utf-8').strip() )
