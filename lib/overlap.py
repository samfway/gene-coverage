#!/usr/bin/env python

__author__ = "The QIIME Development Team"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Sam Way", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"

import bisect

class overlap_detector:
    """ Detect overlap of genes/sequence features using binary searches 
        Each feature has a start and stop point, on a particular strand.

        The purpose of this class is to provide a means of efficiently 
        searching through those features to find overlap.  To do this,
        the start and end points of the feature are treated as "marks"
        along the sequence.  These marks are sorted and then fed into 
        a binary search to narrow the range of possible overlapping 
        features.  

        Intended use:  pass, as feature_locations, the features found 
        in contigs.  Then use find_max_overlap to determine the maximum
        fraction of coverage of a gene/feature, covered by a single contig.

     """ 
    def __init__(self, feature_locations):
        self.mark_pairs = []
        for f in feature_locations:
            self.mark_pairs.append( (f.start, f.end, f.strand) )
            self.mark_pairs.append( (f.end, f.start, f.strand) )
        self.mark_pairs.sort()
        self.marks = [mp[0] for mp in self.mark_pairs ]
            
    def find_max_overlap(self, start, end, strand):
        """ Find maximum (fraction of) overlap by a single contig """ 
        range_start = bisect.bisect_left(self.marks, start)
        range_end = bisect.bisect_right(self.marks, end)
        max_overlap = 0

        for i in xrange(range_start, range_end):
            if self.mark_pairs[i][2] != strand: continue # Not on same strand
            if self.mark_pairs[i][0] < self.mark_pairs[i][1]:  
                # Found start point
                overlap = min(self.mark_pairs[i][1], end) - self.mark_pairs[i][0] + 1
            else:                                  
                # Found end point 
                overlap = self.mark_pairs[i][0] - max(self.mark_pairs[i][1], start) + 1

            if overlap > max_overlap:
                max_overlap = overlap 

        return float(max_overlap)/(end-start+1)
    
