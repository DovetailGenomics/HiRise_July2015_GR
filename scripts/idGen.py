#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from builtins import range

import random

def id(length=5):
    id=""
    for i in range(length):
        id += random.choice( "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz" )
    return(id)

