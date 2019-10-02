# -*- coding: utf-8 -*-

'''
This section, smart24read, is used to read uncompressed CD-1.1 formatted data 
files written by a Geotech SMART-24 digitizer. 

This code is based on the Matlab version written by Dave Whitoff, 
Duncan Marriott, Jay Helmericks, David Fee 
    University of Alaska Geophysical Institute
    last updated on 14 October, 2010
Andrew Winkelman
    University of Alaska Geophysical Institute
    last updated 2 October, 2019

Explaination of function:
    The Geotech SMART-24 digitizer, in the usual WATC use case and configuration,
    writes 15-minute files that roughly follow the CD-1.1 format. There is a single
    file header, equivalent to the frame header specified by CD-1.1. Following the 
    frame header is a single channel subframe header, Followed channel subframes
    each containing a minute of data.
    
    Fields follow the IEEE standard for numerical representation. This implies 
    forward, or big-endian, byte ordering of multibyte values. Additionally, all
    variable-length fields that are not a multiple of four bytes long are padded
    to the next multiple of four with with trailing null (0) bytes.
    
    Frame header:
        Field                   Format      Number      Description
        frame_type              int32       1           numeric identifier of this frame type (5 for data frame)
        trailer_offset          int32       1           byte offset from first byte of the frame to the beginning of the trailer
        frame_creator           char        8           assigned identifier of the creator of the frame
        frame_destination       char        8           identifier of the destination of the frame
        sequence_number         int64       1           sequence number assigned by the frame creator
        series_number           int32       1           series number assigned by the frame creator
    
    Channel subframe header:
        Field                   Format      Number      Description
        num_channels            int32       1           number of channels in this frame
        frame_time_length       int32       1           time in milliseconds this frame encompasses
        nominal_time            char        20          nominal UTC start time of all channels in frame; yyyyddd_hh:mm:ss.ttt
        channel_string_count    int32       1           unpadded length in bytes of the channel string; must be ten times the number of channels field
        channel_string          char        N           channel string listing of the channel subframes to follow, 10 bytes per subframe

    Channel subframe of data frame:
        Field                   Format      Number      Description
        channel_length          int32       1           length, in bytes and divisible by four, of this channel subframe, not counting this integer
        authentication_offset   int32       1           byte offset from the first byte of the frame to the authentication key identifier
        authentication          bool        1           channel_description byte 1: authentication (0 = off; 1 = on)
        transformation          uint8       1           channel_description byte 2: transformation (compression; see spec for details)
        sensor_type             uint8       1           channel_description byte 3: sensor type; (0 = seismic; 1 = hydroacoustic; 2 = infrasonic; 3 = weather; 4 = other; 5 = velocity; 6 = acceleration)
        calibration_provided    bool        1           channel_description byte 4: option flag (0 = unused; 1 = calib and calper provided in bytes 17-24)
        site_name               char        5           channel_description bytes 5-9: site name; left justified, padded with ASCII null bytes as required
        channel_name            char        3           channel_description bytes 10-12: channel name; left justified, padded with ASCII null bytes as required
        location_name           char        2           channel_description bytes 13-14: location name; left justified, padded with ASCII null bytes as required
        data_format             char        2           channel_description bytes 15-16: uncompressed data format (CSS 3.0 data type) ASCII characters, set before signature if frame is signed
        calibration_factor      int32       1           channel_description bytes 17-20: CSS 3.0 calibration factor when byte 4 = 1 (IEEE float)
        calibration_period      int32       1           channel_description bytes 21-24: CSS 3.0 calibration period when byte 4 = 1 (IEEE float)
        timestamp               char        20          UTC start time for first sample of this channel
        duration                int32       1           time in milliseconds spanned by this channel data (always 60,000 for geotech)
        samples                 int32       1           number of samples in this channel subframe
        channel_status_size     int32       1           unpadded length, in bytes, of next field - needs to be 48 bytes for geotech cd11; see additional comments below regarding format
        channel_status_data     bytes       48          status for channel, padded to be divisible by four bytes; see additional comments below regarding format (48 bytes for geotech)
        channel_data_size       int32       1           unpadded length in bytes of next field - needs to be 4 bytes times the number of samples
        channel_data            data        N           data for channel, padded to be divisible by four bytes, of type specified in data_format
        subframe count          int32       1           subframe count as assigned by digitizer; zero for digitizers that do not support this count
        auth_key_identifier     int32       1           pointer to the certificate with the public key to be used for verifying the authentication value field
        authentication_size     int32       1           unpadded length in bytes of the next field
        authentication_value    data        N           DSS signature, padded as necessary to be divisible by four bytes, of type specified in data; see additional comments regarding format
         
        Channel status bytes:
            byte 1:     format of channel status field (1 = this format)
            byte 2:     data status byte:
                bit 1:      future use
                bit 2:      zeroed data
                bit 3:      clipped
                bit 4:      calibration underway
                bit 5-8:    future use
            byte 3:     channel security byte:
                bit 1:      equipment housing open
                bit 2:      digitizing equipment open (tied to bit 4)
                bit 3:      vault door open
                bit 4:      authentication seal broken (tied to bit 2)
                bit 5-8:    future use
            byte 4:     miscellaneous status byte:
                bit 1:      clock unlocked
                bit 2:      GPS receiver off
                bit 3:      GSP receiver unlocked
                bit 4:      digitizer analog input shorted
                bit 5:      digitizer calibration loop back
                bit 6-8:    future use
            byte 5-8:   future use
            byte 9-29:  time of last GPS syncronization [20-byte ASCII chars]
            byte 30-34: clock differential in microseconds [4-byte integer / int32]
            byte 33-36: GPS reported latitude (float)
            byte 37-40: GPS reported longitude (float)
            byte 41-44: GPS reported altitude (float)
            byte 45-48: LSB bitweight volts/count (float)
    
        CD-1.1 allows the data type of channel_data and authentication_value to
        be specified as one of the following formats.
            Data type   Size        Description
            s4          4 bytes     IEEE integer (default)
            s3          3 bytes     IEEE integer, packed
            s2          2 bytes     IEEE integer, short
            i4          4 bytes     4-byte integer
            i2          2 bytes     2-byte integer
            CD          N/A         Encapsulated CD-1 data
        Note that the Geotech SMART-24 only uses the s4 "IEEE integer" data type.
    
        CD-1.1 and the Geotech SMART-24 allows the use of compression, presumably Canadian compression
        as specified by CD-1.1. WATC habit is to use uncompressed data at the digitizer level.
'''


import os
import struct
import obspy
import numpy as np


class Smart24:
    '''
    read an uncompressed cd11_file produced from a Geotech Smart24 digitizer
    
    inputs:
        cd11_file   [string]    path to cd11 file to read; must have .cd11 file extension
    output:
        smart24     [class]     smart24 class containing the following objects
                                (See CD-1.1 specification and Geotech SMART24 User's manual for details)
            cd11_file           [string]    path to input cd11_file
            file_header         [dict]      dictionary containing file header info
            subframe_header     [dict]      dictionary containing subframe header info
            channel_subframe    [list]      list of dictionaries containing read channel subframes
            st                  [Stream]    obspy.Stream object holding the waveform data
            
    '''
    
    def __init__(self, cd11_file, verbose=False):
        
        # parse inputs
        self.cd11_file = cd11_file
        self._verbose = verbose
        
        # set some constants
        self._bytes_size = {'int32': 4,
                      'float': 4,
                      'int64': 8,
                      'char':  1,
                      'uint8': 1,
                      'bool':  1,
                      'bytes': 1}
        
        self._bytes_code = {'int32': 'i',
                      'float': 'f',
                      'int64': 'q',
                      'char':  'c',
                      'uint8': 'b',
                      'bool':  '?'}
        
        # verify file exists and read binary data from file
        self._test_file(extension='.cd11')
        self._binarydata = self._read_binary_cd11_file()
        if self._verbose:
            print(self._binarydata[:150])
        
        # begin reading file header and subframe header
        self._position = 0
        if self._verbose:
            print('Reading file header...')
        self._parse_file_header()
        if self._verbose:
            print('reading subframe header...')
        self._parse_subframe_header()
        
        # read channel subframes until end of binary data is reached
        self.channel_subframe = list()
        # if there's an odd number of channels, skip two bytes
        if self.subframe_header['number_of_channels'] % 2 > 0:
            self._position += 2
        while self._position < len(self._binarydata):
            if self._verbose:
                print('reading channel subframe...')
            self._parse_channel_subframe()
        
        # parse read channel subframes into useful information
        self._build_stream()
    
    
    def _test_file(self, extension):
        '''Test if cd11_file exists and has '.cd11' file extension'''
        
        f = os.path.abspath(self.cd11_file)
        assert os.path.isfile(self.cd11_file), \
            "Error, file {} doesn't appear to exist".format(f)
        assert os.path.splitext(self.cd11_file)[-1] == extension, \
            "Error, file {} doesn't appear to have the '{}' file extension".format(f, extension)
    
    
    def _read_binary_cd11_file(self):
        '''Read cd11_file bytes'''
        
        with open(self.cd11_file, 'rb') as fb:
            return fb.read()
    
    
    def _pad_to_four(self,n):
        '''
        if n is divisible by 4, returns n. if not, adds one until the resulting number is divisible by 4 and returns
        '''
        
        while n % 4:
            n += 1
        return n
    
    
    def _decode_bytes(self, b, n, bytes_format):
        '''decode bytes in b as n values of format b_format'''
        
        if bytes_format in ['int32', 'int64','bool','uint8','float']:
            format_code = str(n) + self._bytes_code[bytes_format]
            if n == 1:
                result = struct.unpack(format_code, b)[0]
            else:
                result = struct.unpack(format_code, b)
        elif bytes_format in ['char']:
            format_code = str(n) + self._bytes_code[bytes_format]
            result = b''.join(struct.unpack(format_code, b)).decode('ISO-8859-1').replace('\x00', '')
        elif bytes_format in ['bytes']:
            result = b
        return result
    
    
    def _parse_next_bytes_as(self, n, bytes_format, pad=False):
        '''parse bytestring from self._binarydata as amount n of bytes_format starting
        from the current self._position, where bytes_format is one of int32, int64,
        float, char, uint8, and bool'''
        

        
        nextbytes = self._binarydata[self._position:self._position+n*self._bytes_size[bytes_format]]
        
        result = self._decode_bytes(nextbytes, n, bytes_format)
        
        if self._verbose:
            print('\tDebug: reading bytes {} (starting at position {}]) as {} {}: {}'.format(nextbytes, self._position, n, bytes_format, result))
        
        if pad:
            self._position += self._pad_to_four(len(nextbytes))
        else:
            self._position += len(nextbytes)
        return result
    
    
    def _parse_file_header(self):
        '''read the first 36 bytes 'file' frame header'''
        
        self.file_header = {'frame_type':        self._parse_next_bytes_as(1, 'int32'),
                            'trailer_offset':    self._parse_next_bytes_as(1, 'int32'),
                            'frame_creator':     self._parse_next_bytes_as(8, 'char'),
                            'frame_destination': self._parse_next_bytes_as(8, 'char'),
                            'sequence_number':   self._parse_next_bytes_as(1, 'int64'),
                            'series_number':     self._parse_next_bytes_as(1, 'int32')}
        
        assert self.file_header['frame_type'] == 5, "Error, frame format not data frame (5)!"


    def _parse_subframe_header(self):
        '''read payload/subframe header'''
        
        self.subframe_header = {'number_of_channels':   self._parse_next_bytes_as(1, 'int32'),
                                'frame_time_length':    self._parse_next_bytes_as(1, 'int32'),
                                'nominal_time':         self._parse_next_bytes_as(20, 'char'),
                                'channel_string_count': self._parse_next_bytes_as(1, 'int32'),
                                'channel_string': []}
        # the subframe_header channel_string length is determined by the channel_string_count
        # so do this after the subframe header exists
        channel_string = self._parse_next_bytes_as(self.subframe_header['channel_string_count'], 'char', pad=True)
        for i in range(self.subframe_header['number_of_channels']):
            self.subframe_header['channel_string'].append(channel_string[10*i:(10*i)+10])
        
        assert self.subframe_header['channel_string_count'] == 10 * self.subframe_header['number_of_channels'],\
            "Error, channel_string-count not 10 * number of channels!"
        
        # let's interpret the nominal_time as obspy.UTCDateTime object
        self.subframe_header['nominal_time'] = obspy.UTCDateTime.strptime(self.subframe_header['nominal_time']+'000','%Y%j %H:%M:%S.%f')
        
        # let's interpret the channel strings in SEED-compliant form
        for i, channel_string in enumerate(self.subframe_header['channel_string']):
            net = ''
            site = channel_string[:5]
            chan = channel_string[5:8]
            loc = channel_string[8:]
            self.subframe_header['channel_string'][i] = '.'.join([net, site, loc, chan])


    def _parse_channel_subframe(self):
        '''parse channel subframe which contains the subframe waveform data and metadata'''
        
        
        channel_subframe = {'channel_length': self._parse_next_bytes_as(1, 'int32'),
                            'authentication_offset': self._parse_next_bytes_as(1, 'int32'),
                            'channel_description': {'authentication': self._parse_next_bytes_as(1, 'uint8'),
                                                    'transformation': self._parse_next_bytes_as(1, 'uint8'),
                                                    'sensor_type': self._parse_next_bytes_as(1, 'uint8'),
                                                    'option_flag': self._parse_next_bytes_as(1, 'uint8'),
                                                    'site_name': self._parse_next_bytes_as(5, 'char'), 
                                                    'channel_name': self._parse_next_bytes_as(3, 'char'),
                                                    'location_name': self._parse_next_bytes_as(2, 'char'),
                                                    'uncompressed_data_format': self._parse_next_bytes_as(2, 'char'),
                                                    'calibration_factor': self._parse_next_bytes_as(1, 'float'),
                                                    'calibration period': self._parse_next_bytes_as(1, 'float')},
                            'time_stamp': obspy.UTCDateTime.strptime(self._parse_next_bytes_as(20, 'char')+'000','%Y%j %H:%M:%S.%f'),
                            'subframe_time_length': self._parse_next_bytes_as(1, 'int32'),
                            'samples': self._parse_next_bytes_as(1, 'int32'),
                            'channel_status_size': self._parse_next_bytes_as(1, 'int32')}
        # now that we know how big the channel_status_data is, let's read that in
        channel_subframe['channel_status_data'] = self._parse_next_bytes_as(channel_subframe['channel_status_size'],'bytes', pad=True)
        channel_subframe['channel_status_data'] = self._decode_channel_status_data(channel_subframe['channel_status_data'])
        channel_subframe['channel_data_size'] = self._parse_next_bytes_as(1, 'int32')
        # now that we know how big the channel data is, lets read that in
        channel_subframe['channel_data'] = self._parse_next_bytes_as(channel_subframe['channel_data_size'], 'bytes', pad=True)
        channel_subframe['subframe_count'] = self._parse_next_bytes_as(1, 'int32')
        channel_subframe['auth_key_id'] = self._parse_next_bytes_as(1, 'int32')
        channel_subframe['auth_size'] = self._parse_next_bytes_as(1, 'int32')
        channel_subframe['auth_value'] = self._parse_next_bytes_as(channel_subframe['auth_size'], 'bytes', pad=True)
        
        assert channel_subframe['subframe_time_length'] == 60000, 'Error: channel_subframe subframe_time_length not 60000!'
        assert channel_subframe['channel_status_size'] == 48, 'Error:  channel_subframe channel_status_size not 48!'
        
        self.channel_subframe.append(channel_subframe)
    
    
    def _decode_channel_status_data(self,channel_status_data):
        '''
        channel status data is up to 32 bytes with content represented as:
            channel status  [8-byte ASCII]  
                byte 1:     format of channel status field (2 = Geotech SMART24)
                byte 2:     data status byte
                                bit 1: 1 = dead sensor
                                bit 2: 1 = zeroed data
                                bit 3: 1 = clipped
                                bit 4: 1 = calibration underway
                                bits 5-8: for future use
                byte 3:     channel security byte
                                bit 1: 1 = equipment housing open
                                bit 2: 1 = digitizing equipment open
                                bit 3: 1 = vault door opened
                                bit 4: 1 = authentication seal broken
                                bit 5: 1 = equipment moved
                                bits 6-8: for future use
                byte 4:     miscellaneous status byte
                                bit 1: 1 = clock differential too large
                                bit 2: 1 = GPS receiver off
                                bit 3: 1 = GPS receiver unlocked
                                bit 4: 1 = digitizer analog input shorted
                                bit 5: 1 = digitizer calibration loop back
                                bits 6-8: for future use
                byte 5:     voltage indicator byte
                                bit 1: 1 = main power failure
                                bit 2: 1 = backup power unstable
                                bits 3-8: for future use
                bytes 6-8:  undefined
                20-byte ASCII:  time of last GPS synchronization
                IEEE int:       clock differential in microseconds
        '''
        
        assert len(channel_status_data) == 48, 'Error, channel_status_data not 48 bytes long'
        
        channel_status = dict()
        cs_pos = 0
        
        # byte 1: format of channel status data (2 = Geotech SMART24)
        channel_status_format = self._decode_bytes(channel_status_data[cs_pos:cs_pos+1], 1, 'uint8')
        cs_pos += self._bytes_size['uint8']
        assert channel_status_format == 2, "Error, channel_status_data['channel_status_format'] != 2 (Geotech SMART24 format)!"
        
        # byte 2: data status
        data_status = channel_status_data[cs_pos:cs_pos+1]
        channel_status['dead_sensor']          = self._is_bit_set(data_status, 1)
        channel_status['zeroed_data']          = self._is_bit_set(data_status, 2)
        channel_status['clipped']              = self._is_bit_set(data_status, 3)
        channel_status['calibration_underway'] = self._is_bit_set(data_status, 4)
        cs_pos += 1
        
        # byte 3: channel security
        channel_security = channel_status_data[cs_pos:cs_pos+1]
        channel_status['equipment_housing_open']     = self._is_bit_set(channel_security, 1)
        channel_status['digitizing_equipment_open']  = self._is_bit_set(channel_security, 2)
        channel_status['vault_door_open']            = self._is_bit_set(channel_security, 3)
        channel_status['authentication_seal_broken'] = self._is_bit_set(channel_security, 4)
        channel_status['equipment_moved']            = self._is_bit_set(channel_security, 5)
        cs_pos += 1
        
        # byte 4: miscellaneous
        miscellaneous = channel_status_data[cs_pos:cs_pos+1]
        channel_status['clock_differential_too_large']    = self._is_bit_set(miscellaneous, 1)
        channel_status['GPS_receiver_off']                = self._is_bit_set(miscellaneous, 2)
        channel_status['GPS_receiver_unlocked']           = self._is_bit_set(miscellaneous, 3)
        channel_status['digitizer_analog_input_shorted']  = self._is_bit_set(miscellaneous, 4)
        channel_status['digitizer_calibration_loop_back'] = self._is_bit_set(miscellaneous, 5)
        cs_pos += 1
        
        # byte 5: voltage indicator
        voltage_indicator = channel_status_data[cs_pos:cs_pos+1]
        channel_status['main_power_failure']    = self._is_bit_set(voltage_indicator, 1)
        channel_status['backup_power_unstable'] = self._is_bit_set(voltage_indicator, 2)
        cs_pos += 1
        
        # bytes 6, 7, 8 undefined
        cs_pos += 3
        
        # bytes 9-28: time of last GPS sync
        channel_status['last_GPS_sync'] = channel_status_data[cs_pos:cs_pos+20].decode()
        cs_pos += 20
        
        # bytes 29-32: clock differential in microseconds
        channel_status['clock_differential'] = self._decode_bytes(channel_status_data[cs_pos:cs_pos+self._bytes_size['int32']], 1, 'int32')
        cs_pos += self._bytes_size['int32']
        
        # bytes 33-36: lattitude (float)
        channel_status['latitude'] = self._decode_bytes(channel_status_data[cs_pos:cs_pos+self._bytes_size['float']], 1, 'float')
        cs_pos += self._bytes_size['float']
        # bytes 37:40: longitude (float)
        channel_status['longitude'] = self._decode_bytes(channel_status_data[cs_pos:cs_pos+self._bytes_size['float']], 1, 'float')
        cs_pos += self._bytes_size['float']
        # bytes 41-44: altitude (float)
        channel_status['altitude'] = self._decode_bytes(channel_status_data[cs_pos:cs_pos+self._bytes_size['float']], 1, 'float')
        cs_pos += self._bytes_size['float']
        # bytes 45-48: LSB volts/count (float)
        channel_status['LSB_bitweight'] = self._decode_bytes(channel_status_data[cs_pos:cs_pos+self._bytes_size['float']], 1, 'float')
        cs_pos += self._bytes_size['float']
        
        assert cs_pos == len(channel_status_data), 'Error, wrong amount of bytes read'
        
        return channel_status
    
    
    def _is_bit_set(self,byte,bitposition):
        '''
        checks if bit in bitposition is set in given byte, returns bool
        '''
        
        bitmask = list('00000000')
        bitmask[-bitposition] = '1'  # negative bitposition because the bits we're checking count from right, Big Endian
        bitmask = ''.join(bitmask)
        
    
        return bool(int(bitmask,2) & int.from_bytes(byte,byteorder='big'))
    
    
    def _build_stream(self):
        '''builds obspy.Stream object containing waveform data'''
        
        DTYPE = {# Big-endian integers
                 's4': '>i',
                 's2': '>h',
                 # Little-endian integers
                 'i4': '<i',
                 'i2': '<h',
                 # ASCII integers
                 'c0': ('S12', np.int),
                 'c#': ('S12', np.int),
                 # Big-endian floating point
                 't4': '>f',
                 't8': '>d',
                 # Little-endian floating point
                 'f4': '<f',
                 'f8': '<d',
                 # ASCII floating point
                 'a0': ('S15', np.float32),
                 'a#': ('S15', np.float32),
                 'b0': ('S24', np.float64),
                 'b#': ('S24', np.float64)}
        
        traces = []
        for csf in self.channel_subframe:
            # trace data
            dtype = DTYPE[csf['channel_description']['uncompressed_data_format']]
            if isinstance(dtype, tuple):
                read_fmt = np.dtype(dtype[0])
                fmt = dtype[1]
            else:
                read_fmt = np.dtype(dtype)
                fmt = read_fmt
            data = obspy.core.compatibility.from_buffer(csf['channel_data'],dtype=read_fmt)
            data = np.require(data, dtype=fmt)
            
            # trace header
            header = {'station': csf['channel_description']['site_name'],
                      'location': csf['channel_description']['location_name'],
                      'channel': csf['channel_description']['channel_name'],
                      'starttime': csf['time_stamp'],
                      'sampling_rate': csf['samples'] / csf['subframe_time_length'] * 1000,
                      'calib': csf['channel_description']['calibration_factor'],
                      'calper': csf['channel_description']['calibration period']}
            tr = obspy.Trace(data, header=header)
            traces.append(tr)
            
        self.st = obspy.Stream(traces)
