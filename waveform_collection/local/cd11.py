# -*- coding: utf-8 -*-

from .clf import clf

import numpy as np


def read_cd11(clf_file,amount='all'):
    '''
    read clf file, find CD-1.1 frames in accompanying CD-1.1 binary data file specified in clf_file, 
    read CD-1.1 frames from binary data file
    
    input:
        clf_file    [string]    path to clf file to read
        amount      [str|int]   how many frames to format output data for
                                'all' means do all the frames in CLF
                                N, where type(N) == int means do first or last N frames
                                whether N > 0 or N < 0, respectively
                                note that if amount is not 'all', then only data frames included
    output:
        frames      [list]      list of frame objects, see cd11_frame class below
    
    '''
    
    cd11_frames = []
    
    clf_frames  = clf(clf_file).frames
    
    if amount == 'all':
        pass
    elif type(amount) == int and amount < 0:
        clf_frames = [f for f in clf_frames if 'd' in f['frameType']][amount:]
    elif type(amount) == int and amount > 0:
        clf_frames = [f for f in clf_frames if 'd' in f['frameType']][:amount]
    
    with open(clf_frames[0]['dataFile'],'rb') as bf:
        for clf_frame_dict in clf_frames:
            # move cursor to offset
            bf.seek(clf_frame_dict['frameOffsetBytes'])
            # read bytes
            cd11_frame_bytes = bf.read(clf_frame_dict['frameSize'])
            # read framebytes into frame object
            # note reading the waveform data within is currently not supported
            cd11_frames.append(cd11_frame(cd11_frame_bytes,skipdata=True))
    

    
    return cd11_frames


class cd11_frame:
    '''
    every cd-1.1 formatted frame consists of a header frame (36 bytes), a frame payload (variable bytes),
    and a frame trailer (variable bytes). This class reads and decodes these portions.
    
    inputs:
        frame       [bytestring]    bytes representing a complete cd-1.1 frame
        skipdata    [bool]          flag to skip reading the channel subframe data into memory
    '''
    
    # sizes of IEEE representations, note these are Big Endian
    _sizes = {'IEEE int': 4,
              'IEEE short': 2,
              'IEEE long': 8,
              'IEEE float': 4}
    
    
    def __init__(self, frame, skipdata=False):
        self._read_frame_header(frame)
        self._read_frame_trailer(frame)
        
        if self.frameheader['frame_type'] == 5:
            self._read_data_frame(frame, skipdata)
    
    
    def __repr__(self):
        s = '{}'.format(self.frameheader)
        return s
    
    
    def _read_frame_header(self,frame):
        '''
        extracts the header from a cd-1.1 formatted frame and returns it as a dictionary in self.frameheader
        
        frame header fields are:
            frame type:         numeric identifier of this frame type
            trailer offset:     byte offset from first byte of the frame to the beginning of the trailer
            frame creator:      assigned identifier of the creator of the frame
            frame destination:  identifier of the destination of the frame
            sequence number:    sequence number assigned by the frame creator
            series: series      number assinged by the frame creator
        
        input:
            frame   [bytes]     cd-1.1 frame to read (note IEEE representations are Big Endian)
                                    frame type          [IEEE integer]
                                    trailer offset      [IEEE integer]
                                    frame creator       [8-byte ASCII]
                                    frame destination   [8-byte ASCII]
                                    sequence number     [IEEE long]
                                    series              [IEEE integer]
        output:
            header  [dict]      header with keys:
                                    frame_type          [int]
                                    trailer_offset      [int]
                                    frame_creator       [str]
                                    frame_destination   [str]
                                    sequence_number     [int]
                                    series              [int]
        '''
        
        sizes = {'frame_type': 4,
                 'trailer_offset': 4,
                 'frame_creator': 8,
                 'frame_destination': 8,
                 'sequence_number': 8,
                 'series': 4}
        
        position = 0
        self.frameheader = dict()
        self.frameheader['frame_type'] = int.from_bytes(frame[position:position+sizes['frame_type']], byteorder='big')
        position += sizes['frame_type']
        self.frameheader['trailer_offset'] = int.from_bytes(frame[position:position+sizes['trailer_offset']], byteorder='big')
        position += sizes['trailer_offset']
        self.frameheader['frame_creator'] = frame[position:position+sizes['frame_creator']].split(b'\x00')[0].decode()
        position += sizes['frame_creator']
        self.frameheader['frame_destination'] = frame[position:position+sizes['frame_destination']].split(b'\x00')[0].decode()
        position += sizes['frame_destination']
        self.frameheader['sequence_number'] = int.from_bytes(frame[position:position+sizes['sequence_number']], byteorder='big')
        position += sizes['sequence_number']
        self.frameheader['series'] = int.from_bytes(frame[position:position+sizes['series']], byteorder='big')
    
    
    def _read_frame_trailer(self,frame):
        '''
        extracts the trailer from a cd-1.1 formatted frame and returns it as a dictionary in self.frametrailer
        
        frame trailer fields are:
            authentication key identifier:  identifier of the certificate with the public key required to verify the authentication value field; if non-zero, then authentication is used to verify communications
            authentication size:            unpadded length in bytes of next field; zero when no authentication; that is, authentication key is zero
            authentication value:           Digital Signature Standard (DSS) signature, padded to be divisible by four. DSA is 40 bytes. This field may have a length of zero.
            comm verification:              error detection for transmission verification
        
        input:
            frame   [bytes]     cd1.1 frame to read (note IEEE representations are Big Endian)
                                    auth_key_id     [IEEE int]
                                    auth_size       [IEEE int]
                                    auth_value      [N data bytes]
                                    comm_ver        [IEEE long]
        output:
            trailer     [dict]  trailer with keys:
                                    auth_key_id     [int]
                                    auth_size       [int]
                                    auth_value      [int]
                                    comm_ver        [int]
        '''
        
        sizes = {'auth_key_id': 4,
                 'auth_size': 4,
                 'comm_ver': 8}
        
        try:
            position = self.frameheader['trailer_offset']
        except (KeyError, AttributeError) as E:
            #print('frame header not read yet or trailer offset not defined; {}'.format(E))
            raise SystemExit('frame header not read yet or trailer offset not defined; {}'.format(E))
        
        self.frametrailer = dict()
        self.frametrailer['auth_key_id'] = int.from_bytes(frame[position:position+sizes['auth_key_id']], byteorder='big')
        position += sizes['auth_key_id']
        self.frametrailer['auth_size'] = int.from_bytes(frame[position:position+sizes['auth_size']], byteorder='big')
        position += sizes['auth_size']
        self.frametrailer['auth_value'] = int.from_bytes(frame[position:position+self.frametrailer['auth_size']], byteorder='big')
        position += self.frametrailer['auth_size']
        self.frametrailer['comm_ver'] = int.from_bytes(frame[position:position+sizes['comm_ver']], byteorder='big')
    
    
    def _read_data_frame(self,frame,skipdata):
        '''extracts the frame payload for data frames (frame type 5)
        
        data frame constituents:
        channel subframe header:
            number of channels:     number of channels (channel subframes) in this frame
            frame time length:      time in milliseconds this frame encompasses
            nominal time:           nominal UTC start time of all channels in frame
            channel string count:   unpadded length in bytes of the channel string; must be 10 times the number of channels field
            channel string:         channl string listing the channel subframes to follow, 10 bytes per subframe, null-padded to a multiple of 4 bytes. 
                                    each 10-byte string is formatted as follows:
                                        5 bytes for the site name, left justified and padded with ASCII null (if necessary)
                                        3 bytes for the channel name, left justified and padded with ASCII null (if necessary)
                                        2 bytes for the location name, left justified and padded with ASCII null (if necessary)
                                    note that IDC software does not currently read, use, or store the location field.
        channel subframe:
            channel length:         length in bytes and divisible by 4 of this channel subframe, not counting this field (IEEE int)
            authentication offset:  byte offset from the first byte of the frame to the authentication key identifier
            channel description:    flags and identifiers for this channel
                                        byte 1:      authentication (0 = off; 1 = on)
                                        byte 2:      transformation (0 = no transformation; 
                                                                     1 = Canadian compression applied before signature;
                                                                     2 = Canadian compression applied after signature;
                                                                     3 = Steim compression applied before signature;
                                                                     4 = Steim compression applied after signature)
                                        byte 3:      sensor type (0 = seismic; 1 = hydroacoustic; 2 = infrasonic; 3 = weather; >3 = other)
                                        byte 4:      option flag (0 = unused; 1 = calib and calper provided in bytes 17-24)
                                        bytes 5-9:   site name; left-justified, padded with ASCII null bytes as required
                                        bytes 10-12: channel name; left justified, padded wth ASCII null bytes as required
                                        bytes 13-14: location name; left justified, padded with ASCII null bytes as required
                                        bytes 15-16: uncompressed data format (CSS 3.0 data type) ASCII characters, set before signature if frame is signed
                                        bytes 17-20: CSS 3.0 calibration factor when byte 4 = 1 (IEEE float)
                                        bytes 21-24: CSS 3.0 calibration period when byte 4 = 1 (IEEE float)
            time stamp:             UTC start time for first sample of this channel
            subframe time length:   time in milliseconds spanned by this channel data
            samples:                number of samples in this channel subframe
            channel status size:    unpadded length in bytes of next field
            channel status data:    status data for channel, padded to be divisible by four
            channel data size:      unpadded length in bytes of next field
            channel data:           data for channel, padded to be divisible by four
            subframe count:         subframe count as assigned by digitizer; zero for digitizers that do not support this count
            authentication key id:  pointer to the certificate with the public key to be used for verifying the authentication value field
            authentication size:    unpadded length in bytes of next field
            authentication value:   DSS signature over the following fields: 
                                        channel description, timestamp, subframe time length, samples, 
                                        channel status size, channel status data, channel data size, 
                                        channel data, and subframe count; 
                                    this frame is padded as necessary to be divisible by four
        
        input:
            frame   [bytes]     cd1.1 frame to read (note IEEE representations are Big Endian)
                                    number of channels:     [IEEE int]
                                    frame time length:      [IEEE int]
                                    nominal time:           [20-byte ASCII]
                                    channel string count:   [IEEE int]
                                    channel string:         [N-byte ASCII]
        
        output:
            channel_subframe_header     [dict]  channel subframe header with keys:
                                                    n_channels          [int]
                                                    frame_time_length   [int]
                                                    nominal_time        [str]
                                                    channel_str_count   [int]
                                                    channel_str         [list] of strings of 10 characters
            channel_subframes           [list]  list of channel subframe dictionaries with keys:
                                                    channel_length:         [IEEE int]
                                                    auth_offset:            [IEEE int]
                                                    channel_description:    [24 data bytes]
                                                    time_stamp:             [20-byte ASCII]
                                                    subframe_time_length:   [IEEE int]
                                                    samples:                [IEEE int]
                                                    channel_status_size:    [IEEE int]
                                                    channel_status_data:    [N data bytes]
                                                    channel_data_size:      [IEEE int]
                                                    channel_data:           [N data bytes]
                                                    subframe_count:         [IEEE int]
                                                    auth_key_id:            [IEEE int]
                                                    auth_size:              [IEEE int]
                                                    auth_value:             [N data bytes]
        '''
        
        try: self.frameheader
        except AttributeError: raise SystemExit('frame header not read yet, do that first')
        
        try: assert self.frameheader['frame_type'] == 5
        except AssertionError: raise SystemExit('frame type not data frame')
        
        position = 36  # length of frame header
        self.channel_subframe_header = {}
        self.channel_subframes = []
        
        # read channel subframe header
        self.channel_subframe_header['n_channels'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
        position += self._sizes['IEEE int']
        self.channel_subframe_header['frame_time_length'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
        position += self._sizes['IEEE int']
        self.channel_subframe_header['nominal_time'] = frame[position:position+20].decode()
        position += 20
        self.channel_subframe_header['channel_str_count'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
        position += self._sizes['IEEE int']
        channel_str = frame[position:position+self.channel_subframe_header['channel_str_count']]
        self.channel_subframe_header['channel_str'] = []
        while channel_str:
            self.channel_subframe_header['channel_str'].append(self._clean_channel_str(channel_str[:10]))
            channel_str = channel_str[10:]
        position += self._pad_to_four(self.channel_subframe_header['channel_str_count'])
        
        # read channel subframes
        for subframe_index in range(self.channel_subframe_header['n_channels']):
            channel_subframe = {}
            channel_subframe['channel_length'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['auth_offset'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['channel_description'] = self._read_channel_description(frame[position:position+24])
            position += 24
            channel_subframe['time_stamp'] = frame[position:position+20].decode()
            position += 20
            channel_subframe['subframe_time_length'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['samples'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['channel_status_size'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            # need to break apart channel_status_data
            channel_status_data = frame[position:position+channel_subframe['channel_status_size']]
            channel_subframe['channel_status'] = self._decode_channel_status_data(channel_status_data)
            position += self._pad_to_four(channel_subframe['channel_status_size'])
            channel_subframe['channel_data_size'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            if not skipdata:
                channel_subframe['channel_data'] = frame[position:position+channel_subframe['channel_data_size']]
            else:
                channel_subframe['channel_data'] = b''
            position += self._pad_to_four(channel_subframe['channel_data_size'])
            channel_subframe['subframe_count'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['auth_key_id'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['auth_size'] = int.from_bytes(frame[position:position+self._sizes['IEEE int']], byteorder='big')
            position += self._sizes['IEEE int']
            channel_subframe['auth_value'] = int.from_bytes(frame[position:position+channel_subframe['auth_size']], byteorder='big')
            position += channel_subframe['auth_size']
            
            self.channel_subframes.append(channel_subframe)
    
    
    def _pad_to_four(self,n):
        '''
        if n is divisible by 4, returns n. if not, adds one until the resulting number is divisible by 4 and returns
        '''
        
        while n % 4:
            n += 1
        return n
    
    
    def _clean_channel_str(self,s):
        '''
        cleans the ASCII null characters from the channel strings and separates station, channel, location with forward slashes
        '''
        
        station = s[:5].strip(b'\x00').decode()
        channel = s[5:8].strip(b'\x00').decode()
        location = s[8:].strip(b'\x00').decode()
        
        return '/'.join([station,channel,location])
    
    
    def _read_channel_description(self,channel_description):
        '''
        reads channel description (24 bytes) and returns a dict with keys:
            auth            (1 byte)    [bool]  False for off; True for on
            transformation  (1 byte)    [int]   0 = no transformation;
                                                1 = Canadian compression applied before signature
                                                2 = Canadian compression applied after signature
                                                3 = Steim compression applied before signature
                                                4 = Steim compression applied after signature
            sensor_type     (1 byte)    [int]   0 = seismic
                                                1 = hydroacoustic
                                                2 = infrasound
                                                3 = weather
                                                >3 = other
            option_flag     (1_byte)    [int]   0 = unused
                                                1 = calib and calper provided in bytes 17-24
            site_name       (5 bytes)   [str]   site name, left justified, padded with ASCII null bytes as required
            channel_name    (3 bytes)   [str]   channel name, left justified, padded with ASCII null bytes as required
            location_name   (2 bytes)   [str]   location name, left justified, padded with ASCII null bytes as required
            data_format     (2 bytes)   [str]   uncompressed data format (CSS 3.0 data type), set before signature frame is signed
            calib           (4 bytes)   [float] CSS 3.0 calibration factor when byte 4 is 1, IEEE float
            calper          (4 bytes)   [float] CSS 3.0 calibration period when byte 4 is 1, IEEE float
        '''
        
#        print(channel_description,'\n')
        
        dt = np.dtype(np.single)
        dt = dt.newbyteorder('B')
        
        cd_pos = 0
        cd = {}
        cd['auth'] = bool(int.from_bytes(channel_description[cd_pos:cd_pos+1],byteorder='big'))
#        print('auth',cd_pos,cd_pos+1,channel_description[cd_pos:cd_pos+1],cd['auth'])
        cd_pos += 1
        cd['transformation'] = int.from_bytes(channel_description[cd_pos:cd_pos+1],byteorder='big')
#        print('transformation',cd_pos,cd_pos+1,channel_description[cd_pos:cd_pos+1],cd['transformation'])
        cd_pos += 1
        cd['sensor_type'] = int.from_bytes(channel_description[cd_pos:cd_pos+1],byteorder='big')
#        print('sensor_type',cd_pos,cd_pos+1,channel_description[cd_pos:cd_pos+1],cd['sensor_type'])
        cd_pos += 1
        cd['option_flag'] = int.from_bytes(channel_description[cd_pos:cd_pos+1],byteorder='big')
#        print('option_flag',cd_pos,cd_pos+1,channel_description[cd_pos:cd_pos+1],cd['option_flag'])
        cd_pos += 1
        cd['site_name'] = channel_description[cd_pos:cd_pos+5].strip(b'\x00').decode()
#        print('site_name',cd_pos,cd_pos+5,channel_description[cd_pos:cd_pos+5],cd['site_name'])
        cd_pos += 5
        cd['channel_name'] = channel_description[cd_pos:cd_pos+3].strip(b'\x00').decode()
#        print('channel_name',cd_pos,cd_pos+3,channel_description[cd_pos:cd_pos+3],cd['channel_name'])
        cd_pos += 3
        cd['location_name'] = channel_description[cd_pos:cd_pos+2].strip(b'\x00').decode()
#        print('location_name',cd_pos,cd_pos+2,channel_description[cd_pos:cd_pos+2],cd['location_name'])
        cd_pos += 2
        cd['data_format'] = channel_description[cd_pos:cd_pos+2].strip(b'\x00').decode()
#        print('data_format',cd_pos,cd_pos+2,channel_description[cd_pos:cd_pos+2],cd['data_format'])
        cd_pos += 2
        if cd['option_flag'] == 1:
            cd['calib'] = np.frombuffer(channel_description[cd_pos:cd_pos+self._sizes['IEEE float']], dtype=dt).item()
#            print('calib',cd_pos,cd_pos+self._sizes['IEEE float'],channel_description[cd_pos:cd_pos+self._sizes['IEEE float']],cd['calib'])
            cd_pos += self._sizes['IEEE float']
            cd['calper'] = np.frombuffer(channel_description[cd_pos:cd_pos+self._sizes['IEEE float']], dtype=dt).item()
#            print('calper',cd_pos,cd_pos+self._sizes['IEEE float'],channel_description[cd_pos:cd_pos+self._sizes['IEEE float']],cd['calper'])
            cd_pos += self._sizes['IEEE float']
        else:
            cd['calib'] = 1.0
            cd['calper'] = 1.0
            cd_pos += 2*self._sizes['IEEE float']
        
#        print('\n')
        
        return cd
    
    
    def _decode_channel_status_data(self,channel_status_data):
        '''
        channel status data is up to 32 bytes with content represented as:
            channel status  [8-byte ASCII]  
                byte 1:     format of channel status field (1 = this format)
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
        
        channel_status = dict()
        cs_pos = 0
        
        channel_status_format = int.from_bytes(channel_status_data[cs_pos:cs_pos+1], byteorder='big')
        cs_pos += 1
        try: assert channel_status_format == 1
        except AssertionError:
            raise SystemExit('channel status format unknown')
        
        data_status = channel_status_data[cs_pos:cs_pos+1]
        channel_status['dead_sensor']          = self._is_bit_set(data_status, 1)
        channel_status['zeroed_data']          = self._is_bit_set(data_status, 2)
        channel_status['clipped']              = self._is_bit_set(data_status, 3)
        channel_status['calibration_underway'] = self._is_bit_set(data_status, 4)
        cs_pos += 1
        channel_security = channel_status_data[cs_pos:cs_pos+1]
        channel_status['equipment_housing_open']     = self._is_bit_set(channel_security, 1)
        channel_status['digitizing_equipment_open']  = self._is_bit_set(channel_security, 2)
        channel_status['vault_door_open']            = self._is_bit_set(channel_security, 3)
        channel_status['authentication_seal_broken'] = self._is_bit_set(channel_security, 4)
        channel_status['equipment_moved']            = self._is_bit_set(channel_security, 5)
        cs_pos += 1
        miscellaneous = channel_status_data[cs_pos:cs_pos+1]
        channel_status['clock_differential_too_large']    = self._is_bit_set(miscellaneous, 1)
        channel_status['GPS_receiver_off']                = self._is_bit_set(miscellaneous, 2)
        channel_status['GPS_receiver_unlocked']           = self._is_bit_set(miscellaneous, 3)
        channel_status['digitizer_analog_input_shorted']  = self._is_bit_set(miscellaneous, 4)
        channel_status['digitizer_calibration_loop_back'] = self._is_bit_set(miscellaneous, 5)
        cs_pos += 1
        voltage_indicator = channel_status_data[cs_pos:cs_pos+1]
        channel_status['main_power_failure']    = self._is_bit_set(voltage_indicator, 1)
        channel_status['backup_power_unstable'] = self._is_bit_set(voltage_indicator, 2)
        cs_pos += 1
        #skip three undefined bytes
        cs_pos += 3
        channel_status['last_GPS_sync'] = channel_status_data[cs_pos:cs_pos+20].decode()
        cs_pos += 20
        channel_status['clock_differential'] = int.from_bytes(channel_status_data[cs_pos:cs_pos+self._sizes['IEEE int']], byteorder='big')
        cs_pos += self._sizes['IEEE int']
        
        return channel_status
    
    
    def _is_bit_set(self,byte,bitposition):
        '''
        checks if bit in bitposition is set in given byte, returns bool
        '''
        
        bitmask = list('00000000')
        bitmask[-bitposition] = '1'  # negative bitposition because the bits we're checking count from right, Big Endian
        bitmask = ''.join(bitmask)
        
    
        return bool(int(bitmask,2) & int.from_bytes(byte,byteorder='big'))