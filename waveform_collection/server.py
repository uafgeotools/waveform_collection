from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.earthworm import Client as EW_Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.geodetics import gps2dist_azimuth
from obspy import Stream
import numpy as np
import os
import fnmatch
import warnings
from . import CollectionWarning
from .local.common import load_json_file


# Get location of AVO JSON files
json_dir = os.path.join(os.path.dirname(__file__), '..', 'avo_json')

# Load AVO infrasound station calibration values (units are Pa/ct)
AVO_INFRA_CALIBS = load_json_file(os.path.join(json_dir, 'avo_infra_calibs.json'))

# Load AVO station coordinates (elevation units are meters)
AVO_COORDS = load_json_file(os.path.join(json_dir, 'avo_coords.json'))

# Default infrasound channels - covering all the bases here!
INFRASOUND_CHANNELS = 'BDF,HDF,CDF,DDF'

# Define some conversion factors
KM2M = 1000     # [m/km]
SEC2MIN = 1/60  # [min/s]


def gather_waveforms(source, network, station, location, channel, starttime,
                     endtime, time_buffer=0, merge=True, remove_response=False,
                     return_failed_stations=False, watc_url=None,
                     watc_username=None, watc_password=None):
    """
    Gather seismic/infrasound waveforms from IRIS or WATC FDSN, or AVO Winston,
    and output a :class:`~obspy.core.stream.Stream` with station/element
    coordinates attached. Optionally remove the sensitivity.

    **NOTE**

    Usual RTM usage is to specify a starttime/endtime that brackets the
    estimated source origin time. Then time_buffer is used to download enough
    extra data to account for the time required for an infrasound signal to
    propagate to the farthest station.

    Args:
        source (str): Which source to gather waveforms from. Options are:

            * `'IRIS'` – IRIS FDSN
            * `'WATC'` – WATC FDSN
            * `'AVO'` – AVO Winston

        network (str): SEED network code [wildcards (``*``, ``?``) accepted]
        station (str): SEED station code [wildcards (``*``, ``?``) accepted]
        location (str): SEED location code [wildcards (``*``, ``?``) accepted]
        channel (str): SEED channel code [wildcards (``*``, ``?``) accepted]
        starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time for
            data request
        endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): End time for
            data request
        time_buffer (int or float): Extra amount of data to download after
            `endtime` [s]
        merge (bool): Toggle merging of :class:`~obspy.core.trace.Trace` objects
            with identical IDs
        remove_response (bool): Toggle response removal via
            :meth:`~obspy.core.trace.Trace.remove_sensitivity` or a simple
            scalar multiplication
        return_failed_stations (bool): If `True`, returns a list of station
            codes that were requested but not downloaded. This disables the
            standard failed station warning message
        watc_url (str): URL for WATC FDSN server
        watc_username (str): Username for WATC FDSN server
        watc_password (str): Password for WATC FDSN server

    Returns:
        :class:`~obspy.core.stream.Stream` containing gathered waveforms. If
        `return_failed_stations` is `True`, additionally returns a list
        containing station codes that were requested but not downloaded
    """

    print('--------------')
    print('GATHERING DATA')
    print('--------------')

    # IRIS FDSN
    if source == 'IRIS':

        client = FDSN_Client('IRIS')
        print('Reading data from IRIS FDSN...')
        try:
            st_out = client.get_waveforms(network, station, location, channel,
                                          starttime, endtime + time_buffer,
                                          attach_response=True)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # WATC FDSN
    elif source == 'WATC':

        print('Connecting to WATC FDSN...')
        client = FDSN_Client(base_url=watc_url, user=watc_username,
                             password=watc_password)
        print('Successfully connected. Reading data from WATC FDSN...')
        try:
            st_out = client.get_waveforms(network, station, location, channel,
                                          starttime, endtime + time_buffer,
                                          attach_response=True)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # AVO Winston
    elif source == 'AVO':

        client = EW_Client('pubavo1.wr.usgs.gov', port=16023)  # 16023 is long-term
        print('Reading data from AVO Winston...')
        st_out = Stream()  # Make empty Stream object to populate

        # Brute-force "dynamic grid search" over network/station/channel/location codes
        for nw in _restricted_matching('network', network, client):
            for sta in _restricted_matching('station', station, client, network=nw):
                for cha in _restricted_matching('channel', channel, client, network=nw, station=sta):
                    for loc in _restricted_matching('location', location, client, network=nw, station=sta, channel=cha):
                        try:
                            st_out += client.get_waveforms(nw, sta, loc, cha, starttime, endtime + time_buffer)
                        except KeyError:
                            pass

    else:
        raise ValueError('Unrecognized source. Valid options are \'IRIS\', '
                         '\'WATC\', or \'AVO\'.')

    if merge:
        st_out.merge()  # Merge Traces with the same ID
    st_out.sort()

    # Check that all requested stations are present in Stream
    requested_stations = station.split(',')
    downloaded_stations = [tr.stats.station for tr in st_out]
    failed_stations = []
    for sta in requested_stations:
        # The below check works with wildcards, but obviously cannot detect if
        # ALL stations corresponding to a given wildcard (e.g., O??K) were
        # downloaded. Thus, if careful station selection is desired, specify
        # each station explicitly and the below check will then be effective.
        if not fnmatch.filter(downloaded_stations, sta):
            if not return_failed_stations:
                # If we're not returning the failed stations, then show this
                # warning message to alert the user
                warnings.warn(f'Station {sta} not downloaded from {source} '
                              'server for this time period.', CollectionWarning)
            failed_stations.append(sta)

    # If the Stream is empty, then we can stop here
    if st_out.count() == 0:
        print('No data downloaded.')
        if return_failed_stations:
            return st_out, failed_stations
        else:
            return st_out

    # Otherwise, show what the Stream contains
    print(st_out.__str__(extended=True))  # This syntax prints the WHOLE Stream

    # Add zeros to ensure all Traces have same length
    st_out.trim(starttime, endtime + time_buffer, pad=True, fill_value=0)

    print('Assigning coordinates...')

    # Use IRIS inventory info for AVO data source
    if source == 'AVO':
        client = FDSN_Client('IRIS')

    try:
        inv = client.get_stations(network=network, station=station,
                                  location=location, channel=channel,
                                  starttime=starttime,
                                  endtime=endtime + time_buffer,
                                  level='channel')
    except FDSNNoDataException:
        inv = []

    for tr in st_out:
        for nw in inv:
            for sta in nw:
                for cha in sta:
                    # Being very thorough to check if everything matches!
                    if (tr.stats.network == nw.code and
                            tr.stats.station == sta.code and
                            tr.stats.location == cha.location_code and
                            tr.stats.channel == cha.code):

                        tr.stats.longitude = cha.longitude
                        tr.stats.latitude = cha.latitude
                        tr.stats.elevation = cha.elevation

    # Check if any Trace did NOT get coordinates assigned, and try to use JSON
    # coordinates if available
    for tr in st_out:
        try:
            tr.stats.longitude, tr.stats.latitude, tr.stats.elevation
        except AttributeError:
            try:
                tr.stats.latitude, tr.stats.longitude,\
                    tr.stats.elevation = AVO_COORDS[tr.id]
                warnings.warn(f'Using coordinates from JSON file for {tr.id}.',
                              CollectionWarning)
            except KeyError:
                print(f'No coordinates available for {tr.id}. Stopping.')
                raise

    # Remove sensitivity
    if remove_response:

        print('Removing sensitivity...')

        for tr in st_out:
            try:
                # Just removing sensitivity for now. remove_response() can lead
                # to errors. This should be sufficient for now. Plus some
                # IRIS-AVO responses are wonky.
                tr.remove_sensitivity()
            except ValueError:  # No response information found
                # This is only set up for infrasound calibration values
                try:
                    calib = AVO_INFRA_CALIBS[tr.id]
                    tr.data = tr.data * calib
                    warnings.warn('Using calibration value from JSON file for '
                                  f'{tr.id}.', CollectionWarning)
                except KeyError:
                    print(f'No calibration value available for {tr.id}. '
                          'Stopping.')
                    raise

    print('Done')

    # Return the Stream with coordinates attached (and responses removed if
    # specified)
    if return_failed_stations:
        return st_out, failed_stations
    else:
        return st_out


def gather_waveforms_bulk(lon_0, lat_0, max_radius, starttime, endtime,
                          channel, network='*', station='*', location='*',
                          time_buffer=0, merge=True, remove_response=False,
                          watc_url=None, watc_username=None, watc_password=None):
    """
    Bulk gather infrasound waveforms within a specified maximum radius of a
    specified location. Waveforms are gathered from IRIS (and optionally WATC)
    FDSN, and AVO Winston. Outputs a :class:`~obspy.core.stream.Stream` with
    station/element coordinates attached. Optionally removes the sensitivity.
    (Output :class:`~obspy.core.stream.Stream` has the same properties as output
    :class:`~obspy.core.stream.Stream` from :func:`gather_waveforms`.)

    **NOTE 1**

    WATC database will NOT be used for station search NOR data download unless
    BOTH `watc_username` and `watc_password` are set.

    **NOTE 2**

    Usual RTM usage is to specify a starttime/endtime that brackets the
    estimated source origin time. Then time_buffer is used to download enough
    extra data to account for the time required for an infrasound signal to
    propagate to the farthest station.

    Args:
        lon_0 (int or float): Longitude of search center [deg.]
        lat_0 (int or float): Latitude of search center [deg.]
        max_radius (int or float): Maximum radius to search for stations within
            [km]
        starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time for
            data request
        endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): End time for
            data request
        channel (str): SEED channel code [wildcards (``*``, ``?``) accepted]
            (REQUIRED PARAMETER!)
        network (str): SEED network code [wildcards (``*``, ``?``) accepted]
        station (str): SEED station code [wildcards (``*``, ``?``) accepted]
        location (str): SEED location code [wildcards (``*``, ``?``) accepted]
        time_buffer (int or float): Extra amount of data to download after
            `endtime` [s]
        merge (bool): Toggle merging of :class:`~obspy.core.trace.Trace` objects
            with identical IDs
        remove_response (bool): Toggle response removal via
            :meth:`~obspy.core.trace.Trace.remove_sensitivity` or a simple
            scalar multiplication
        watc_url (str): URL for WATC FDSN server
        watc_username (str): Username for WATC FDSN server
        watc_password (str): Password for WATC FDSN server

    Returns:
        :class:`~obspy.core.stream.Stream` containing bulk gathered waveforms
    """

    print('-------------------')
    print('BULK GATHERING DATA')
    print('-------------------')

    print('Creating station list...')

    inventories = []  # Create empty list of inventories

    # Grab IRIS inventory
    iris_client = FDSN_Client('IRIS')
    try:
        iris_inv = iris_client.get_stations(starttime=starttime,
                                            endtime=endtime + time_buffer,
                                            network=network, station=station,
                                            location=location, channel=channel,
                                            level='channel')

        inventories.append(iris_inv)  # Add IRIS inventory to list

    except FDSNNoDataException:
        print('No stations found on IRIS FDSN.')

    # If the user supplied both a WATC password and WATC username, then search
    # through WATC database
    if watc_username and watc_password:

        # Grab WATC inventory
        print('Connecting to WATC FDSN...')
        watc_client = FDSN_Client(base_url=watc_url, user=watc_username,
                                  password=watc_password)
        print('Successfully connected.')
        try:
            watc_inv = watc_client.get_stations(starttime=starttime,
                                                endtime=endtime + time_buffer,
                                                network=network, station=station,
                                                location=location, channel=channel,
                                                level='channel')

            inventories.append(watc_inv)  # Add WATC inventory to list

        except FDSNNoDataException:
            print('No stations found on WATC FDSN.')

    requested_station_list = []  # Initialize list of stations to request

    # Big loop through all channels in all inventories!
    for inv in inventories:
        for nw in inv:
            for stn in nw:
                for cha in stn:
                    dist, _, _ = gps2dist_azimuth(lat_0, lon_0, cha.latitude,
                                                  cha.longitude)  # [m]
                    if dist <= max_radius * KM2M:
                        requested_station_list.append(stn.code)

    # Loop through each entry in AVO station coordinates JSON file
    for tr_id, coord in AVO_COORDS.items():

        nw, sta, loc, cha = tr_id.split('.')  # Extract codes from Trace.id

        # Only add station to requested stations list if it satisfies the
        # user-supplied query restrictions
        if (_matching([nw], network) and
                _matching([sta], station) and
                _matching([loc], location) and
                _matching([cha], channel)):

            dist, _, _ = gps2dist_azimuth(lat_0, lon_0, *coord[0:2])  # [m]
            if dist <= max_radius * KM2M:
                requested_station_list.append(sta)

    if not requested_station_list:
        raise ValueError('Station list is empty. Expand the station search '
                         'and try again.')

    # Put into the correct format for ObsPy (e.g., 'HOM,O22K,DLL')
    requested_stations = ','.join(np.unique(requested_station_list))

    print('Done')

    if time_buffer != 0:
        print(f'Using time buffer of {time_buffer:.1f} s '
              f'(~{time_buffer * SEC2MIN:.0f} min)')

    print('Making calls to gather_waveforms()...')

    st_out = Stream()  # Initialize empty Stream to populate

    # Gather waveforms from IRIS
    iris_st, iris_failed = gather_waveforms(source='IRIS', network=network,
                                            station=requested_stations,
                                            location=location, channel=channel,
                                            starttime=starttime,
                                            endtime=endtime,
                                            time_buffer=time_buffer,
                                            remove_response=remove_response,
                                            return_failed_stations=True)
    st_out += iris_st

    # If IRIS couldn't grab all stations in requested station list, try WATC
    # (if the user set username and password)
    if iris_failed:

        if watc_username and watc_password:
            # Gather waveforms from WATC
            watc_st, watc_failed = gather_waveforms(source='WATC', network=network,
                                                    station=','.join(iris_failed),
                                                    location=location,
                                                    channel=channel,
                                                    starttime=starttime,
                                                    endtime=endtime,
                                                    time_buffer=time_buffer,
                                                    remove_response=remove_response,
                                                    return_failed_stations=True,
                                                    watc_url=watc_url,
                                                    watc_username=watc_username,
                                                    watc_password=watc_password)
        else:
            # Return an empty Stream and same failed stations
            watc_st = Stream()
            watc_failed = iris_failed

        st_out += watc_st

        # If WATC couldn't grab all stations missed by IRIS, try AVO
        if watc_failed:

            # Gather waveforms from AVO
            avo_st, remaining_failed = gather_waveforms(source='AVO',
                                                        network=network,
                                                        station=','.join(watc_failed),
                                                        location=location,
                                                        channel=channel,
                                                        starttime=starttime,
                                                        endtime=endtime,
                                                        time_buffer=time_buffer,
                                                        remove_response=remove_response,
                                                        return_failed_stations=True)

            if remaining_failed:
                print('--------------')
                for sta in remaining_failed:
                    warnings.warn(f'Station {sta} found in radius search but '
                                  'no data found.', CollectionWarning)

            st_out += avo_st

    if merge:
        st_out.merge()  # Merge Traces with the same ID
    st_out.sort()

    print('--------------')
    print('Finishing gathering waveforms from station list. Check warnings '
          'for any missed stations.')

    return st_out


def _restricted_matching(code_type, requested_codes, avo_client,
                         **restriction_kwargs):
    """
    Find all SEED network/station/location/channel codes on AVO Winston that
    match a user-supplied query string. Optionally constrain the search to a
    particular network/station/location/channel using keyword arguments passed
    on to ``avo_client.get_availability()``.

    Args:
        code_type (str): One of `'network'`, `'station'`, `'location'`, or
            `'channel'`
        requested_codes (str): Comma-separated SEED code string (wildcards
            accepted)
        avo_client (:class:`~obspy.clients.earthworm.client.Client`): AVO
            Winston client instance
        **restriction_kwargs: Query restrictions to be passed on to
            ``avo_client.get_availability()``

    Returns:
        list: A list of SEED codes for `code_type`, subject to the query
        restrictions given in `**restriction_kwargs` AND matching the patterns
        in `requested_codes`
    """

    # Get availability on AVO Winston subject to optional network/station/
    # location/channel restrictions
    inv = np.array(avo_client.get_availability(**restriction_kwargs))

    if inv.size == 0:
        # Nothing available; create empty code container
        all_codes = np.empty((4, 0))
    else:
        # Discard the time info (5th and 6th entries) and convert to strings
        all_codes = inv[:,0:4].T.astype(str)

    # Gather unique codes present on AVO Winston given restrictions
    nw_unique, sta_unique, \
    loc_unique, cha_unique = [np.unique(code).tolist() for code in all_codes]

    if code_type == 'network':
        restricted_matching_codes = _matching(nw_unique, requested_codes)
    elif code_type == 'station':
        restricted_matching_codes = _matching(sta_unique, requested_codes)
    elif code_type == 'location':
        restricted_matching_codes = _matching(loc_unique, requested_codes)
    elif code_type == 'channel':
        restricted_matching_codes = _matching(cha_unique, requested_codes)
    else:
        raise ValueError(f'Code type \'{code_type}\' not recognized!')

    return restricted_matching_codes


def _matching(unique_code_list, requested_codes):
    """
    Takes a comma-separated SEED code string (e.g., ``'BD?,HDF'``) and returns
    the subset of an input list of unique codes (e.g.,
    ``['BDF', 'EHZ', 'DDF']``) that matches the patterns in the comma-separated
    SEED code string.

    Args:
        unique_code_list (list): List of unique code strings
        requested_codes (str): Comma-separated SEED code string (wildcards
            accepted)

    Returns:
        list: Subset of `unique_code_list` that matches the patterns in
        `requested_codes`
    """

    matching_codes = []

    # Split requested codes/patterns into a list of strings
    for pattern in requested_codes.split(','):

        # Return the subset of unique codes matched by the pattern
        matching_codes += fnmatch.filter(unique_code_list, pattern)

    return np.unique(matching_codes).tolist()
