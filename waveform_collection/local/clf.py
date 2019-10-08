# -*- coding: utf-8 -*-

from .common import check_file_exists, check_file_extension
from os.path import splitext
from obspy import UTCDateTime
from datetime import datetime
import os
from glob import glob


def find_latest_clf_file(array):
    clf_dir = os.path.join(os.sep,'cdrecv')
    latest_clf = sorted(glob(os.path.join(clf_dir,'{}*.clf'.format(array))))[-1]
    
    return latest_clf


class clf():
    '''
    clf = clf(clf_file)
    
    The clf object creates a list containing dictionaries describing each data frame
    from the specified clf file. For details about frame contents, consult the CDTools Software User Guide.
    
    '''
    
    def __init__(self, clf_file):
        self.clf_file = clf_file
        self._check_file()
        
        self.frames = list()
        self._read_frames()
    
    def _check_file(self):
        assert check_file_exists(self.clf_file)
        assert check_file_extension(self.clf_file, '.clf')
    
    def _parse_frame(self,line):
        f = dict()
        columns = line.split()
        f['frameCreator'], f['frameDestination'], f['frameNumber'] = columns[0].split('/')
        dataFileExtension = columns[1]
        f['dataFile'] = '.'.join([splitext(self.clf_file)[0],dataFileExtension])
        f['frameOffsetBytes'] = int(columns[2])
        f['frameSize'] = int(columns[3])
        f['frameDirection'], f['frameType'] = list(columns[4])
        f['frameTime'] = UTCDateTime(float(columns[5]))
        if f['frameType'] == 'd':
            f['frameNominalTime'] = UTCDateTime(float(columns[6]))
        else:
            f['frameNominalTime'] = columns[6]
        f['frameParsingCode'] = columns[7]
        f['frameSpanTime'] = float(columns[8])
        f['frameChannelMask'] = columns[9].strip('|')
        if f['frameType'] == 'h':
            f['frameChannelMaskIndex'] = int(columns[10])
        else:
            f['frameTimeLengthMillisec'] = int(columns[10])
        
        self.frames.append(f)
    
    def _read_frames(self):
        with open(self.clf_file, 'r') as fp:
            lines = [line.strip() for line in fp.readlines()]
        for line in lines:
            self._parse_frame(line)
    
    def calc_timeliness(self):
        '''
        calculate timeliness of each data frame by subtracting the frameTime (when we received
        the frame) from the frameNominalTime (time of the first data point in the frame), minus
        the frameSpanTime (to remove transmission time from IDC sender to WATC receiver)
        
        '''
        
        for frame in self.frames:
            if 'd' in frame['frameType']:
                frame['timeliness'] = frame['frameTime'] - frame['frameNominalTime'] - frame['frameSpanTime']
#        self.timeliness = {frame['frameNumber']: frame['frameTime'] - frame['frameNominalTime'] - frame['frameSpanTime'] 
#                                for frame in self.frames if 'd' in frame['frameType']}
