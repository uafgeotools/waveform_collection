from obspy.clients.fdsn import Client as FDSN_Client
from obspy.clients.earthworm import Client as EW_Client
from obspy.clients.fdsn.header import FDSNNoDataException
from obspy.geodetics import gps2dist_azimuth
from obspy import Inventory, Stream
from multiprocess import Pool
import numpy as np
import os
import fnmatch
import warnings
from . import CollectionWarning

# Default infrasound channels - covering all the bases here!
INFRASOUND_CHANNELS = 'BDF,HDF,CDF,DDF'

# Define some conversion factors
KM2M = 1000     # [m/km]
SEC2MIN = 1/60  # [min/s]


def gather_waveforms(source, network, station, location, channel, starttime,
                     endtime, time_buffer=0, merge_fill_value=0,
                     trim_fill_value=0, remove_response=False,
                     return_failed_stations=False, watc_url=None,
                     watc_username=None, watc_password=None, verbose=True, parallel=False, cores=None):
    """
    Gather seismic/infrasound waveforms from any ObsPy-supported FDSN
    (see https://docs.obspy.org/packages/obspy.clients.fdsn.html), WATC FDSN, 
    or AVO Winston, and output a :class:`~obspy.core.stream.Stream` with 
    station/element coordinates attached. Optionally remove the response.

    **NOTE**

    Usual RTM usage is to specify a starttime/endtime that brackets the
    estimated source origin time. Then time_buffer is used to download enough
    extra data to account for the time required for an infrasound signal to
    propagate to the farthest station.

    Args:
        source (str): Which source to gather waveforms from. Options are:

            * Any ObsPy-supported FDSN (e.g., `'IRIS'`, `'NCEDC'`). For full 
              list, see https://docs.obspy.org/packages/obspy.clients.fdsn.html
            * `'WATC'` – WATC FDSN (requires `watc_url`, `watc_username`, and
              `watc_password`)
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
        merge_fill_value (bool, int, float, str, or None): Controls merging of
            :class:`~obspy.core.trace.Trace` objects with identical IDs. If
            `False`, no merging is performed. Otherwise, a merge is performed
            with the ``fill_value`` provided to this parameter. For details,
            see the docstring of :meth:`obspy.core.stream.Stream.trim`
        trim_fill_value (bool, int, float, or None): Controls trimming of the
            output :class:`~obspy.core.stream.Stream`, useful if precisely
            uniform start and end times are desired. If `False`, no trimming is
            performed. Otherwise, a trim is performed with the ``fill_value``
            provided to this parameter. For details, see the docstring of
            :meth:`obspy.core.stream.Stream.merge`
        remove_response (bool or str): Response removal via full frequency deconvolution
            (`'full'`) or single frequency sensitivity (`'sens'`) / a simple
            scalar multiplication. Default is `False` to return stream in
            counts.
        return_failed_stations (bool): If `True`, returns a list of station
            codes that were requested but not downloaded. This disables the
            standard failed station warning message
        watc_url (str): URL for WATC FDSN server (``http://10.30.6.3:8080``, or
            ``http://10.30.5.10:8080`` if using VPN)
        watc_username (str): Username for WATC FDSN server
        watc_password (str): Password for WATC FDSN server
        verbose (bool): If `False`, all print statements will be blocked. 
            Default is `True`.
        parallel (bool): If `True`, use increments and parallel processing to download data
        cores (int): Number of cores to use for parallel processing

    Returns:
        :class:`~obspy.core.stream.Stream` containing gathered waveforms. If
        `return_failed_stations` is `True`, additionally returns a list
        containing station codes that were requested but not downloaded
    """
    # log() does nothing if `verbose=False`
    log = print if verbose else lambda *args, **kwargs: None
    
    # Check for issues with fill value args
    if merge_fill_value is True or trim_fill_value is True:
        raise ValueError('Cannot provide True to fill value parameters.')

    log('--------------')
    if parallel:
        log('GATHERING DATA IN PARALLEL')
    else:
        log('GATHERING DATA')
    log('--------------')

    # WATC FDSN
    if source == 'WATC':

        # Check that all three arguments required for the WATC server are present
        if watc_url is None or watc_username is None or watc_password is None:
            raise ValueError('WATC source requires watc_url, watc_username, and '
                             'watc_password.')

        log('Connecting to WATC FDSN...')
        client = FDSN_Client(base_url=watc_url, user=watc_username,
                             password=watc_password)
        log('Successfully connected. Reading data from WATC FDSN...')
        try:
            if parallel:
                st_out = parallel_gather(client, network, station, location,
                                         channel, starttime, endtime + time_buffer, cores)
            else:
                st_out = client.get_waveforms(network, station, location, channel,
                                              starttime, endtime + time_buffer,
                                              attach_response=True)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # AVO Winston
    elif source == 'AVO':

        client = EW_Client('pubavo1.wr.usgs.gov', port=16023)  # 16023 is long-term
        log('Reading data from AVO Winston...')
        st_out = Stream()  # Make empty Stream object to populate

        # Brute-force "dynamic grid search" over network/station/channel/location codes
        for nw in _restricted_matching('network', network, client):
            for sta in _restricted_matching('station', station, client, network=nw):
                for cha in _restricted_matching('channel', channel, client, network=nw, station=sta):
                    for loc in _restricted_matching('location', location, client, network=nw, station=sta, channel=cha):
                        try:
                            if parallel:
                                st_out = parallel_gather(client, nw, sta, loc,
                                                         cha, starttime, endtime + time_buffer, cores)
                            else:
                                st_out += client.get_waveforms(nw, sta, loc, cha, starttime, endtime + time_buffer)
                        except KeyError:
                            pass

    # ObsPy-supported FDSN
    else:
        client = FDSN_Client(source)
        log('Reading data from %s FDSN...' % source)
        try:
            if parallel:
                st_out = parallel_gather(client, network, station, location,
                                         channel, starttime, endtime + time_buffer, cores)
            else:
                st_out = client.get_waveforms(network, station, location, channel,
                                              starttime, endtime + time_buffer,
                                              attach_response=True)
        except FDSNNoDataException:
            st_out = Stream()  # Just create an empty Stream object

    # Merge, if specified
    if merge_fill_value is not False:
        _safe_merge(st_out, merge_fill_value)  # Merge Traces with same ID
        warnings.warn(f'Merging with "fill_value={merge_fill_value}"',
                      CollectionWarning)

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
        log('No data downloaded.')
        if return_failed_stations:
            return st_out, failed_stations
        else:
            return st_out

    # Otherwise, show what the Stream contains
    log(st_out.__str__(extended=True))  # This syntax prints the WHOLE Stream

    # Trim, if specified
    if trim_fill_value is not False:
        st_out.trim(starttime, endtime + time_buffer, pad=True,
                    fill_value=trim_fill_value)
        warnings.warn(f'Trimming with "fill_value={trim_fill_value}"',
                      CollectionWarning)

    log('Assigning coordinates...')

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
        inv = Inventory()  # Make an empty inv
        warnings.warn('Creating empty inventory.', CollectionWarning)

    for tr in st_out:
        try:
            coords = inv.get_coordinates(tr.id)
            tr.stats.longitude = coords['longitude']
            tr.stats.latitude = coords['latitude']
            tr.stats.elevation = coords['elevation']
        except Exception as e:
            if str(e) == 'No matching channel metadata found.':
                warnings.warn(f'No metadata for {tr.id} found in inventory.',
                              CollectionWarning)
            else:
                raise

    # Check if any Trace did NOT get coordinates assigned
    for tr in st_out:
        try:
            tr.stats.longitude, tr.stats.latitude, tr.stats.elevation
        except AttributeError:
            log(f'No coordinates available for {tr.id}. Stopping.')
            raise

    # Remove response
    if remove_response:
    
        # Backwards compatibility
        if remove_response is True:
          remove_response = 'sens'
        
        log(f'Removing response, method={remove_response}')

        for tr in st_out:
            try:
                # Remove full frequency response
                if remove_response == 'full':
                    # pre-filter for response removal, these values should
                    # work for most cases
                    pre_filt = [0.0005, 0.001, (tr.stats.sampling_rate / 2) - 2,
                                tr.stats.sampling_rate/2]
                    tr.remove_response(pre_filt=pre_filt, output='VEL',
                                       water_level=None)

                # Just apply single freq sensitivity
                elif remove_response == 'sens':
                    tr.remove_sensitivity()

            except ValueError:  # No response information found
                log(f'No calibration value available for {tr.id}. '
                        'Stopping.')
                raise

    log('Done')

    # Return the Stream with coordinates attached (and responses removed if
    # specified)
    if return_failed_stations:
        return st_out, failed_stations
    else:
        return st_out


def gather_waveforms_bulk(lon_0, lat_0, max_radius, starttime, endtime,
                          channel, network='*', station='*', location='*',
                          time_buffer=0, merge_fill_value=0, trim_fill_value=0,
                          remove_response=False, watc_url=None,
                          watc_username=None, watc_password=None, iris_only=True,
                          verbose=True, parallel=False, cores=None):
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
        merge_fill_value (bool, int, float, str, or None): Controls merging of
            :class:`~obspy.core.trace.Trace` objects with identical IDs. If
            `False`, no merging is performed. Otherwise, a merge is performed
            with the ``fill_value`` provided to this parameter. For details,
            see the docstring of :meth:`obspy.core.stream.Stream.trim`
        trim_fill_value (bool, int, float, or None): Controls trimming of the
            output :class:`~obspy.core.stream.Stream`, useful if precisely
            uniform start and end times are desired. If `False`, no trimming is
            performed. Otherwise, a trim is performed with the ``fill_value``
            provided to this parameter. For details, see the docstring of
            :meth:`obspy.core.stream.Stream.merge`
        remove_response (bool): Toggle response removal via
            :meth:`~obspy.core.trace.Trace.remove_sensitivity` or a simple
            scalar multiplication
        watc_url (str): URL for WATC FDSN server
        watc_username (str): Username for WATC FDSN server
        watc_password (str): Password for WATC FDSN server
        iris_only (bool): If `True`, only the IRIS FDSN source is used
        verbose (bool): If `False`, all print statements will be blocked. 
            Default is `True`.

    Returns:
        :class:`~obspy.core.stream.Stream` containing bulk gathered waveforms
    """
    # log() does nothing if `verbose=False`
    log = print if verbose else lambda *args, **kwargs: None
    
    # Check for issues with fill value args
    if merge_fill_value is True or trim_fill_value is True:
        raise ValueError('Cannot provide True to fill value parameters.')

    log('-------------------')
    log('BULK GATHERING DATA')
    log('-------------------')

    log('Creating station list...')

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
        log('No stations found on IRIS FDSN.')

    # If the user supplied both a WATC password and WATC username, then search
    # through WATC database
    if watc_username and watc_password:

        # Grab WATC inventory
        log('Connecting to WATC FDSN...')
        watc_client = FDSN_Client(base_url=watc_url, user=watc_username,
                                  password=watc_password)
        log('Successfully connected.')
        try:
            watc_inv = watc_client.get_stations(starttime=starttime,
                                                endtime=endtime + time_buffer,
                                                network=network, station=station,
                                                location=location, channel=channel,
                                                level='channel')

            inventories.append(watc_inv)  # Add WATC inventory to list

        except FDSNNoDataException:
            log('No stations found on WATC FDSN.')

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

    if not requested_station_list:
        raise ValueError('Station list is empty. Expand the station search '
                         'and try again.')

    # Put into the correct format for ObsPy (e.g., 'HOM,O22K,DLL')
    requested_stations = ','.join(np.unique(requested_station_list))

    log('Done')

    if time_buffer != 0:
        log(f'Using time buffer of {time_buffer:.1f} s '
              f'(~{time_buffer * SEC2MIN:.0f} min)')

    log('Making calls to gather_waveforms()...')

    st_out = Stream()  # Initialize empty Stream to populate

    # Gather waveforms from IRIS
    iris_st, iris_failed = gather_waveforms(source='IRIS', network=network,
                                            station=requested_stations,
                                            location=location, channel=channel,
                                            starttime=starttime,
                                            endtime=endtime,
                                            time_buffer=time_buffer,
                                            merge_fill_value=False,
                                            trim_fill_value=False,
                                            remove_response=remove_response,
                                            return_failed_stations=True,
                                            verbose=verbose, parallel=parallel, cores=cores)
    st_out += iris_st

    # If IRIS couldn't grab all stations in requested station list, try WATC
    # (if the user set username and password)
    if iris_failed:

        if iris_only:
            log('--------------')
            for sta in iris_failed:
                warnings.warn(f'Station {sta} found in radius search but '
                              'no data found. Try specifying `iris_only=False` to look '
                              'in additional data archives.', CollectionWarning)
        else:
            if watc_username and watc_password:
                # Gather waveforms from WATC
                watc_st, watc_failed = gather_waveforms(source='WATC', network=network,
                                                        station=','.join(iris_failed),
                                                        location=location,
                                                        channel=channel,
                                                        starttime=starttime,
                                                        endtime=endtime,
                                                        time_buffer=time_buffer,
                                                        merge_fill_value=False,
                                                        trim_fill_value=False,
                                                        remove_response=remove_response,
                                                        return_failed_stations=True,
                                                        watc_url=watc_url,
                                                        watc_username=watc_username,
                                                        watc_password=watc_password,
                                                        verbose=verbose, parallel=parallel, cores=cores)
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
                                                            merge_fill_value=False,
                                                            trim_fill_value=False,
                                                            remove_response=remove_response,
                                                            return_failed_stations=True,
                                                            verbose=verbose, parallel=parallel, cores=cores)

                if remaining_failed:
                    log('--------------')
                    for sta in remaining_failed:
                        warnings.warn(f'Station {sta} found in radius search but '
                                      'no data found.', CollectionWarning)

                st_out += avo_st

    # Merge, if specified
    if merge_fill_value is not False:
        _safe_merge(st_out, merge_fill_value)  # Merge Traces with same ID
        warnings.warn(f'Merging with "fill_value={merge_fill_value}"',
                      CollectionWarning)

    # Trim, if specified
    if trim_fill_value is not False:
        st_out.trim(starttime, endtime + time_buffer, pad=True,
                    fill_value=trim_fill_value)
        warnings.warn(f'Trimming with "fill_value={trim_fill_value}"',
                      CollectionWarning)

    st_out.sort()

    log('--------------')
    log('Finishing gathering waveforms from station list. Check warnings '
          'for any missed stations.')

    return st_out


def _safe_merge(st, fill_value):
    """
    Merge Traces with same ID, modifying data types and rounding off non-integer sampling rates if necessary.
    Modified from code by Aaron Wech.

    Args:
        st (:class:`~obspy.core.stream.Stream`): Input Stream (modified in-place!)
        fill_value (int, float, str, or None): Passed on to
            :meth:`obspy.core.stream.Stream.merge`
    """

    try:
        st.merge(fill_value=fill_value)
    except Exception:  # ObsPy raises an Exception if data types are not all identical
        for tr in st:
            if tr.data.dtype != np.dtype(np.int32):
                tr.data = tr.data.astype(np.int32, copy=False)
    try:
        st.merge(fill_value=fill_value)
    except Exception:  # ObsPy also raises an Exception if traces with the same ids have different sampling rates
        for tr in st:
            if tr.stats.sampling_rate != np.round(tr.stats.sampling_rate):
                warnings.warn('Rounding off %s sampling rate from %f Hz to %.1f Hz for merge compatibility.' % (
                              tr.id, tr.stats.sampling_rate, np.round(tr.stats.sampling_rate)), CollectionWarning)
                tr.stats.sampling_rate = np.round(tr.stats.sampling_rate)
        st.merge(fill_value=fill_value)  # Try merging with rounded sampling rates


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

def get_increment_waveforms(args):
    """
    Function required for parallel processing

    Args:
        args (tuple): A tuple containing the following elements:
            client (obspy.clients.fdsn.Client): The FDSN client to use for data retrieval.
            network (str): SEED network code [wildcards (``*``, ``?``) accepted].
            station (str): SEED station code [wildcards (``*``, ``?``) accepted].
            location (str): SEED location code [wildcards (``*``, ``?``) accepted].
            channel (str): SEED channel code [wildcards (``*``, ``?``) accepted].
            starttime (obspy.core.utcdatetime.UTCDateTime): Start time for data request.
            endtime (obspy.core.utcdatetime.UTCDateTime): End time for data request.

    Returns:
        obspy.core.stream.Stream: Stream of gathered waveforms for the specified time increment.
    """
    client, network, station, location, channel, starttime, endtime = args
    try:
        st_inc = client.get_waveforms(network=network, station=station, location=location, channel=channel,
                                  starttime=starttime, endtime=endtime, attach_response=True)
        return st_inc
    except FDSNNoDataException:
        return Stream()  # Just create an empty Stream object


def parallel_gather(client, network, station, location, channel, starttime,
                     endtime, cores, increment_s=12*60*60):
    """
    Gather waveforms incrementally and parallelize by time, station, or location using multiple cores.

    Args:
        client (obspy.clients.fdsn.Client): The FDSN client to use for data retrieval.
        network (str): SEED network code [wildcards (``*``, ``?``) accepted]
        station (str): SEED station code [wildcards (``*``, ``?``) accepted]
        location (str): SEED location code [wildcards (``*``, ``?``) accepted]
        channel (str): SEED channel code [wildcards (``*``, ``?``) accepted]
        starttime (:class:`~obspy.core.utcdatetime.UTCDateTime`): Start time for
            data request
        endtime (:class:`~obspy.core.utcdatetime.UTCDateTime`): End time for
            data request
        cores (int): Number of cores to use for parallel processing.
        increment_s (int, optional): Time increment in seconds for each parallel task. Default is 12 hours (43200 seconds).

    Returns:
        obspy.core.stream.Stream: Combined stream of gathered waveforms.
    """
    if cores is None:  # Check if cores is specified
        warnings.warn("Number of cores not specified. Using available cores - 2.", CollectionWarning)
        cores = os.cpu_count() -2

    if cores > os.cpu_count():  # Check if requested cores exceeds available cores
        warnings.warn(f"Number of cores requested ({cores}) exceeds number of available cores ({os.cpu_count()}). "
                      f"Using {os.cpu_count()} cores instead.", CollectionWarning)
        cores = os.cpu_count()

    if len(station.split(',')) == 1:  # For a single station, parallelize by either time windows or locations.
        inv = client.get_stations(network=network, station=station,
                                  starttime=starttime, endtime=endtime, level="channel")
        for net in inv:
            for sta in net:
                locs = [channel.location_code for channel in sta]  # Get location codes

        if len(locs) > 2:  # For more than 2 locations (meant for arrays), parallelize by locations instead of time
            with Pool(processes=cores) as pool:
                results = pool.map(get_increment_waveforms,
                                   [(client, network, station, locations, channel, starttime, endtime) for locations in locs])

        else:  # For a single location, parallelize by time windows
            time_ranges = []
            increment_start = starttime
            while increment_start < endtime:
                increment_end = min(increment_start + increment_s, endtime)  # selects end time of the increment
                if increment_end > increment_start:
                    time_ranges.append((increment_start, increment_end))  # appends the start and end time of the increment
                increment_start = increment_end  # updates new increment start time

            with Pool(processes=cores) as pool:  # Use multiprocess to gather data in parallel
                results = pool.map(get_increment_waveforms,
                                   [(client, network, station, location, channel, times[0], times[1]) for times in time_ranges])

    else:  # For multiple stations, parallelize by stations
        station = station.split(',')  # Split the station string into a list of strings

        with Pool(processes=cores) as pool:
            results = pool.map(get_increment_waveforms,
                               [(client, network, stations, location, channel, starttime, endtime) for stations in station])

    combined_stream = Stream()  # Initialize empty Stream object
    for st_inc in results:  # Combine all the streams into one
        combined_stream += st_inc

    return combined_stream
