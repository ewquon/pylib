#!/usr/bin/env python
"""
Utilities to get information about a running jupyter notebook server

Based on https://stackoverflow.com/questions/23886989/match-a-running-ipython-notebook-to-a-process
"""
import sys,os
import json
import posixpath
import subprocess
from urllib.request import urlopen

import pandas as pd
import psutil


def notebooks_table(host='localhost',port=8888,token=None):
    """Show table with info about running jupyter notebooks.

    Args:
        host: host of the jupyter server.
        port: port of the jupyter server.

    Returns:
        DataFrame with rows corresponding to running notebooks and following columns:
            * index: notebook kernel id.
            * path: path to notebook file.
            * pid: pid of the notebook process.
            * memory: notebook memory consumption in percentage.
    """
    notebooks = get_running_notebooks(host, port, token)
    df = pd.DataFrame(notebooks)
    df = df.set_index('kernel_id')
    df['pid'] = df.apply(lambda row: get_process_id(row.name), axis=1)
    # same notebook can be run in multiple processes
    df = expand_column(df, 'pid')
    df['memory'] = df.pid.apply(memory_usage_psutil)
    return df.sort_values('memory', ascending=False)


def get_running_notebooks(host='localhost',port=8888,token=None):
    sessions_url = posixpath.join('http://{:s}:{:d}'.format(host, int(port)), 'api', 'sessions')
    if token is not None:
        sessions_url += '?token='+token
    response = urlopen(sessions_url).read()
    res = json.loads(response)
    notebooks = [{'kernel_id': notebook['kernel']['id'],
                  'path': notebook['notebook']['path']} for notebook in res]
    return notebooks


def get_process_id(name):
    """Return process ids found by (partial) name or regex.

    Source: https://stackoverflow.com/a/44712205/304209.
    >>> get_process_id('kthreadd')
    [2]
    >>> get_process_id('watchdog')
    [10, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61]  # ymmv
    >>> get_process_id('non-existent process')
    []
    """
    child = subprocess.Popen(['pgrep', '-f', name], stdout=subprocess.PIPE, shell=False)
    response = child.communicate()[0]
    return [int(pid) for pid in response.split()]


def memory_usage_psutil(pid=None):
    """Get memory usage percentage by current process or by process specified by id, like in top.

    Source: https://stackoverflow.com/a/30014612/304209.

    Args:
        pid: pid of the process to analyze. If None, analyze the current process.

    Returns:
        memory usage of the process, in percentage like in top, values in [0, 100].
    """
    if pid is None:
        pid = os.getpid()
    process = psutil.Process(pid)
    return process.memory_percent()


def long_substr(strings):
    """Find longest common substring in a list of strings.

    Source: https://stackoverflow.com/a/2894073/304209.

    Args:
        strings: list of strings.

    Returns:
        longest substring which is found in all of the strings.
    """
    substr = ''
    if len(strings) > 1 and len(strings[0]) > 0:
        for i in range(len(strings[0])):
            for j in range(len(strings[0])-i+1):
                if j > len(substr) and all(strings[0][i:i+j] in x for x in strings):
                    substr = strings[0][i:i+j]
    print(strings)
    print(substr)
    return substr


def expand_column(dataframe, column):
    """Transform iterable column values into multiple rows.

    Source: https://stackoverflow.com/a/27266225/304209.

    Args:
        dataframe: DataFrame to process.
        column: name of the column to expand.

    Returns:
        copy of the DataFrame with the following updates:
            * for rows where column contains only 1 value, keep them as is.
            * for rows where column contains a list of values, transform them
                into multiple rows, each of which contains one value from the list in column.
    """
    tmp_df = dataframe.apply(
        lambda row: pd.Series(row[column]), axis=1).stack().reset_index(level=1, drop=True)
    tmp_df.name = column
    return dataframe.drop(column, axis=1).join(tmp_df)


#==============================================================================
if __name__ == '__main__':
    Nargs = len(sys.argv) - 1
    if Nargs == 1:
        host = 'localhost'
        port = int(sys.argv[1])
        token = None
    elif Nargs == 2:
        host = 'localhost'
        port = int(sys.argv[1])
        token = sys.argv[2]
    elif Nargs == 3:
        host = sys.argv[1]
        port = int(sys.argv[2])
        token = sys.argv[3]
    else:
        print('USAGE:')
        print(' ',sys.argv[0],'[port]')
        print(' ',sys.argv[0],'[port] [token]')
        print(' ',sys.argv[0],'[hostname] [port] [token]')
        sys.exit(1)
    #print(get_running_notebooks(host,port,token))
    print(notebooks_table(host,port,token))


