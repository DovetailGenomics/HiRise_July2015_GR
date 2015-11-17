#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#

from builtins import object

class BamTags(object):

    @staticmethod
    def mate_mapq(line):
        try:
            return line.opt("xm")
        except KeyError:
            return line.opt("MQ")

    @staticmethod
    def mate_start(line):
        return line.opt("xs")

    @staticmethod
    def orientation(line):
        return line.opt("xt")

    @staticmethod
    def junction(line):
        return line.opt("xj")
